# ── 1. Imports + constants ─────────────────────────────────────────────────────
from __future__ import annotations

from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

import asyncio
import math as _math
import os

import pandas as pd
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from pydantic import BaseModel
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

# Optional: RDKit for 2D structure rendering + Tanimoto similarity
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs, inchi as rdInchi
    from rdkit.Chem.Draw import rdMolDraw2D as _rdMolDraw2D
    import numpy as np
    _RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    AllChem = None
    DataStructs = None
    rdInchi = None
    _rdMolDraw2D = None
    _RDKIT_AVAILABLE = False

# Full-chain reference values (oracle — never alter)
FULLCHAIN_SI = 13.4463
FULLCHAIN_APPROVAL_RATE = 16.23

DATA_DIR = Path(os.environ.get("SYNERGPCR_DATA", "./data"))


def _norm_phase(v):
    """Normalise clinical_phase values from the CSV at read time."""
    if not v or str(v).strip().lower() in ('nan', 'none', '', 'null'):
        return "Investigational"
    s = str(v).strip()
    abbrev = {
        'Inv.': 'Investigational',
        'inv.': 'Investigational',
        'investigational': 'Investigational',
        'Phase 1': 'Phase 1',
        'Phase 2': 'Phase 2',
        'Phase 3': 'Phase 3',
        'Phase 4': 'Phase 4',
        'Approved': 'Approved',
        'approved': 'Approved',
        'Withdrawn': 'Withdrawn',
    }
    return abbrev.get(s, s)

# Chain-pattern priority ranking (higher = more complete assay coverage)
_PATTERN_RANK: dict[str, int] = {
    "Full_chain": 7,
    "L1+L2+L3+L4": 7,
    "L1+L2+L4": 6,
    "L2+L3+L4": 5,
    "L2+L4": 4,
    "L1+L4": 4,
    "L3+L4": 3,
    "L1+L2+L3": 3,
    "L1+L2": 2,
    "L2+L3": 2,
    "L2_only": 1,
    "L1_only": 1,
    "L3_only": 1,
    "L4_only": 1,
}


def _compute_pair_chain_pattern(l1: str, l2: str,
                                 l3: str, l4: str) -> str:
    """Compute pair-level chain pattern from L1-L4 states.
    Only 'Active' counts; NaN / Inactive / None are excluded.
    """
    active = []
    for label, val in [('L1', l1), ('L2', l2), ('L3', l3), ('L4', l4)]:
        s = str(val or '').strip()
        if s == 'Active':
            active.append(label)
    if not active:
        return 'No_active_layer'
    if active == ['L1', 'L2', 'L3', 'L4']:
        return 'Full_chain'
    if len(active) == 1:
        return active[0] + '_only'
    return '+'.join(active)


def _s(val: object) -> str:
    """Safely coerce a value to str, mapping NaN/None to ''."""
    if val is None:
        return ""
    if isinstance(val, float) and _math.isnan(val):
        return ""
    s = str(val).strip()
    return "" if s.lower() in ("nan", "none") else s


def _sanitize_records(records: list[dict]) -> list[dict]:
    """Replace float NaN values with None so JSON serialization produces null."""
    return [
        {k: (None if (isinstance(v, float) and _math.isnan(v)) else v)
         for k, v in row.items()}
        for row in records
    ]


# ── Fingerprint index (built once at startup) ──────────────────────────────────
FP_INDEX: dict = {
    "ikeys":  [],   # list[str]  — InChIKey strings
    "fps":    [],   # list[ExplicitBitVect]  — RDKit FP objects
    "ready":  False,
}

# ── 2. Global data store ────────────────────────────────────────────────────────
DATA: dict = {}


# ── 3. Lifespan context manager ─────────────────────────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Load all CSV datasets into memory on startup."""
    print("=== SynerGPCR startup: loading data ===")

    # ── Core tables ──
    try:
        DATA["compounds"] = pd.read_csv(
            DATA_DIR / "compound_lookup_final.csv", dtype=str, low_memory=False,
        ).set_index("InChIKey")
        _n_assay = int(
            (DATA["compounds"]["in_assay_universe"].astype(str) == "True").sum()
        )
        DATA["assay_ikeys"] = set(
            DATA["compounds"].index[
                DATA["compounds"]["in_assay_universe"].astype(str) == "True"
            ]
        )
        print(
            f"  compounds loaded: {len(DATA['compounds']):,} total "
            f"({_n_assay:,} with assay data)"
        )
    except Exception as e:
        print(f"  compounds: NOT FOUND ({e})")

    try:
        DATA["targets"] = pd.read_csv(
            DATA_DIR / "target_lookup_filled.csv", dtype=str,
        ).set_index("UniProt_AC")
        print(f"  targets: {len(DATA['targets']):,} rows")
    except Exception as e:
        print(f"  targets: NOT FOUND ({e})")

    try:
        chains = pd.read_csv(
            DATA_DIR / "compound_chain_summary.csv", dtype=str, low_memory=False,
        )
        DATA["chains"] = chains
        DATA["chains_by_ikey"] = {k: v for k, v in chains.groupby("InChIKey")}
        DATA["chains_by_target"] = {k: v for k, v in chains.groupby("UniProt_AC")}
        print(f"  chains: {len(chains):,} rows")
    except Exception as e:
        print(f"  chains: NOT FOUND ({e})")

    try:
        preds = pd.read_csv(
            DATA_DIR / "gpcract_predictions.csv", dtype=str, low_memory=False,
        )
        preds = preds[preds["status"] != "AWAITING_GPCRACT_V2_INFERENCE"].copy()
        DATA["predictions"] = preds
        DATA["pred_by_ikey"] = {k: v for k, v in preds.groupby("InChIKey")}
        n_real_preds = len(preds)
        if n_real_preds == 0:
            print("  WARNING: gpcract_predictions.csv contains ONLY placeholder rows.")
            print("  GPCRact v2 inference has not run yet. AI nodes will show as 'no_data'.")
        else:
            print(f"  predictions: {n_real_preds:,} rows (real inference results)")
    except Exception as e:
        print(f"  predictions: NOT FOUND ({e})")

    # ── Synergy tables ──
    try:
        _syn_raw = pd.read_csv(DATA_DIR / "synergy_stats.csv")
        _syn_raw = _syn_raw.dropna(subset=["bliss_SI"]).copy()
        DATA["synergy"] = _syn_raw.astype(str)
        DATA["synergy_map"] = DATA["synergy"].set_index("chain_pattern").to_dict("index")
        print(f"  synergy_stats: {len(DATA['synergy']):,} rows (bliss_SI NaN rows dropped)")
    except Exception as e:
        print(f"  synergy_stats: NOT FOUND ({e})")

    try:
        DATA["single_layer"] = pd.read_csv(DATA_DIR / "single_layer_table.csv", dtype=str)
        print(f"  single_layer: {len(DATA['single_layer']):,} rows")
    except Exception as e:
        print(f"  single_layer: NOT FOUND ({e})")

    # ── Disease / indication ──
    try:
        DATA["diseases"] = pd.read_csv(DATA_DIR / "drug_indication.csv", dtype=str)
        print(f"  drug_indication: {len(DATA['diseases']):,} rows")
    except Exception as e:
        DATA["diseases"] = None
        print(f"  drug_indication: NOT FOUND ({e})")

    try:
        DATA["target_indication"] = pd.read_csv(DATA_DIR / "target_indication.csv", dtype=str)
        print(f"  target_indication: {len(DATA['target_indication']):,} rows")
    except Exception as e:
        print(f"  target_indication: NOT FOUND ({e})")

    try:
        DATA["disease_cats"] = pd.read_csv(DATA_DIR / "disease_category_lookup.csv", dtype=str)
        print(f"  disease_cats: {len(DATA['disease_cats']):,} rows")
    except Exception as e:
        print(f"  disease_cats: NOT FOUND ({e})")

    # ── AI / repurposing ──
    try:
        DATA["candidates"] = pd.read_csv(DATA_DIR / "candidate_shortlist.csv", dtype=str)
        print(f"  candidates: {len(DATA['candidates']):,} rows")
    except Exception as e:
        print(f"  candidates: NOT FOUND ({e})")

    try:
        DATA["repurposing_crossfam"] = pd.read_csv(DATA_DIR / "repurposing_view_crossfam.csv", dtype=str)
        print(f"  repurposing_crossfam: {len(DATA['repurposing_crossfam']):,} rows")
    except Exception as e:
        print(f"  repurposing_crossfam: NOT FOUND ({e})")

    try:
        DATA["repurposing_full"] = pd.read_csv(DATA_DIR / "repurposing_view_full.csv", dtype=str)
        print(f"  repurposing_full: {len(DATA['repurposing_full']):,} rows")
    except Exception as e:
        print(f"  repurposing_full: NOT FOUND ({e})")

    try:
        DATA["gpcract_3point"] = pd.read_csv(DATA_DIR / "gpcract_3point_si.csv", dtype=str)
        print(f"  gpcract_3point: {len(DATA['gpcract_3point']):,} rows")
    except Exception as e:
        print(f"  gpcract_3point: NOT FOUND ({e})")

    # ── Target enrichment ──
    try:
        DATA["target_enrichment"] = pd.read_csv(
            DATA_DIR / "target_clinical_enrichment.tsv", sep="\t", dtype=str,
        )
        DATA["target_enrichment_map"] = DATA["target_enrichment"].set_index("UniProt").to_dict("index")
        print(f"  target_enrichment: {len(DATA['target_enrichment']):,} rows")
    except Exception as e:
        print(f"  target_enrichment: NOT FOUND ({e})")

    try:
        DATA["target_disease"] = pd.read_csv(
            DATA_DIR / "target_disease_enrichment.tsv", sep="\t", dtype=str,
        )
        print(f"  target_disease: {len(DATA['target_disease']):,} rows")
    except Exception as e:
        print(f"  target_disease: NOT FOUND ({e})")

    # ── UMLS disease categories ──
    try:
        DATA["umls_cats"] = pd.read_csv(DATA_DIR / "UMLS_CUI_to_Category_v2.csv", dtype=str)
        print(f"  umls_cats: {len(DATA['umls_cats']):,} rows")
    except Exception as e:
        print(f"  umls_cats: NOT FOUND ({e})")

    # ── GPCR fallback ──
    try:
        DATA["gpcr_fallback"] = pd.read_csv(DATA_DIR / "gpcr_fallback_lookup.csv", dtype=str)
        print(f"  gpcr_fallback: {len(DATA['gpcr_fallback']):,} rows")
    except Exception as e:
        print(f"  gpcr_fallback: NOT FOUND ({e})")

    # NOTE: assay_stats_precomputed.csv is an empty stub — do NOT load it.

    # ── Pre-compute Morgan fingerprints for Tanimoto search ────────────
    # Run in a thread pool so the event loop is not blocked during startup.
    if _RDKIT_AVAILABLE:
        def _build_fp_index() -> None:
            """CPU-bound FP build — executed in a thread pool executor."""
            try:
                _df_fp = DATA["compounds"]
                _pool = _df_fp[
                    _df_fp["in_assay_universe"].astype(str) == "True"
                ][["smiles"]].copy()

                _ikeys_tmp, _fps_tmp = [], []
                for _ikey, _row in _pool.iterrows():
                    _smi = str(_row.get("smiles", "") or "").strip()
                    if not _smi or _smi.lower() in ("nan", "none", ""):
                        continue
                    _mol = Chem.MolFromSmiles(_smi)
                    if _mol is None:
                        continue
                    _ikeys_tmp.append(_ikey)
                    _fps_tmp.append(
                        AllChem.GetMorganFingerprintAsBitVect(_mol, 2, 2048)
                    )

                FP_INDEX["ikeys"] = _ikeys_tmp
                FP_INDEX["fps"]   = _fps_tmp
                FP_INDEX["ready"] = True
                print(
                    f"[SynerGPCR] Fingerprint cache built: "
                    f"{len(_ikeys_tmp):,} assay-active compounds indexed."
                )
            except Exception as _e:
                print(f"[SynerGPCR] WARNING: FP cache build failed: {_e}")
                FP_INDEX["ready"] = False

        # Schedule FP build as a background task; server is already
        # accepting requests before this completes.
        asyncio.get_running_loop().run_in_executor(None, _build_fp_index)
        print("[SynerGPCR] Fingerprint index build started in background.")

    print(f"  rdkit: {'available' if _RDKIT_AVAILABLE else 'NOT INSTALLED — structure rendering disabled'}")
    print("=== SynerGPCR ready on port 8091 ===")
    yield


# ── 4. FastAPI app ──────────────────────────────────────────────────────────────
app = FastAPI(title="SynerGPCR", lifespan=lifespan)

# ── 5. Static files ─────────────────────────────────────────────────────────────
app.mount("/static", StaticFiles(directory="static"), name="static")

# ── 6. Templates ────────────────────────────────────────────────────────────────
templates = Jinja2Templates(directory="templates")
templates.env.auto_reload = True

from urllib.parse import quote as _url_quote
templates.env.filters["urlencode"] = lambda s: _url_quote(str(s or ""), safe="")


# ── Synergy helper functions ───────────────────────────────────────────────────

def get_synergy_row(pattern: str) -> dict | None:
    """Return synergy stats for a chain_pattern, or None if not found."""
    return DATA.get("synergy_map", {}).get(pattern)



def _best_pattern(chain_records: list[dict]) -> str | None:
    """Return compound-level best chain pattern as union of active layers
    across all GPCR pairs. Matches the oracle definition in synergy_stats.csv
    (validated in webexport_step3_synergy.py / fix_compound_chain_summary_pairlevel.py).

    Each element of chain_records must have keys L1, L2, L3, L4 with
    values 'Active', 'Inactive', or empty / None.
    """
    active: set[str] = set()
    for cr in chain_records:
        for layer in ("L1", "L2", "L3", "L4"):
            if str(cr.get(layer, "") or "").strip() == "Active":
                active.add(layer)
    if not active:
        return "No_active_layer"
    ordered = [L for L in ("L1", "L2", "L3", "L4") if L in active]
    if ordered == ["L1", "L2", "L3", "L4"]:
        return "Full_chain"
    if len(ordered) == 1:
        return ordered[0] + "_only"
    return "+".join(ordered)


# ── 7. Page routes ─────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def home(request: Request) -> HTMLResponse:
    """Render the home page with summary statistics."""
    stats = {
        "n_assay_compounds": int((DATA["compounds"]["in_assay_universe"] == "True").sum()) if "compounds" in DATA else 0,
        "n_approved_db": int((DATA["compounds"]["is_approved"] == "True").sum()) if "compounds" in DATA else 0,
        "n_approved_assay": int(
            ((DATA["compounds"]["in_assay_universe"] == "True") &
             (DATA["compounds"]["is_approved"] == "True")).sum()
        ) if "compounds" in DATA else 0,
        "n_targets": len(DATA.get("targets", [])),
        "n_ai_targets": int(
            DATA["predictions"][DATA["predictions"]["is_high_confidence"] == "True"]["UniProt_AC"].nunique()
        ) if "predictions" in DATA and len(DATA["predictions"]) > 0 else 0,
        "n_disease_cats": len(DATA.get("disease_cats", [])),
        "n_enriched_targets": int(
            (DATA["target_enrichment"]["enough_n_exp"] == "True").sum()
        ) if "target_enrichment" in DATA else 0,
        "n_candidates": len(DATA.get("candidates", [])),
        "n_repurposing_web": DATA["repurposing_crossfam"]["inchikey"].nunique()
                             if "repurposing_crossfam" in DATA else 0,
    }
    synergy_rows = DATA["synergy"].to_dict("records") if "synergy" in DATA else []
    return templates.TemplateResponse(
        request,
        "home.html",
        {"stats": stats, "synergy_rows": synergy_rows},
    )


@app.get("/browse", response_class=HTMLResponse)
async def browse(
    request: Request,
    target: Optional[str] = None,
    disease: Optional[str] = None,
    filter: Optional[str] = None,
    tab: Optional[str] = None,
) -> HTMLResponse:
    """Render the browse page, optionally pre-filtered."""
    return templates.TemplateResponse(
        request,
        "browse.html",
        {
            "filter_target": target,
            "filter_disease": disease,
            "filter_param": filter,
            "tab_param": tab,
        },
    )


@app.get("/about", response_class=HTMLResponse)
async def about(request: Request) -> HTMLResponse:
    """Render the about page."""
    return templates.TemplateResponse(request, "about.html")


@app.get("/help", response_class=HTMLResponse)
async def help_page(request: Request) -> HTMLResponse:
    """Render the Help and Usage Guide page."""
    return templates.TemplateResponse(request, "help.html", {})


@app.get("/compound/{ikey}", response_class=HTMLResponse)
async def compound_view(request: Request, ikey: str) -> HTMLResponse:
    """Render the compound detail page, or 404 if InChIKey is unknown."""
    if ikey not in DATA["compounds"].index:
        raise HTTPException(
            status_code=404,
            detail=f"Compound '{ikey}' not found in SynerGPCR"
        )
    compound = DATA["compounds"].loc[ikey].to_dict()
    compound["clinical_phase"] = _norm_phase(compound.get("clinical_phase", ""))
    return templates.TemplateResponse(
        request,
        "compound.html",
        {"ikey": ikey, "compound": compound},
    )


@app.get("/ai-prediction", response_class=HTMLResponse)
async def ai_prediction(request: Request) -> HTMLResponse:
    """Render the AI prediction page."""
    return templates.TemplateResponse(request, "ai_prediction.html", {})


@app.get("/download", response_class=HTMLResponse)
async def download(request: Request) -> HTMLResponse:
    """Render the download page."""
    return templates.TemplateResponse(request, "download.html", {})


@app.get("/disease/{disease_name:path}", response_class=HTMLResponse)
async def disease_page(request: Request, disease_name: str) -> HTMLResponse:
    """Render disease detail page."""
    return templates.TemplateResponse(
        request,
        "disease.html",
        {"disease_name": disease_name},
    )


# ── 8. API routes ──────────────────────────────────────────────────────────────

@app.get("/api/search")
async def api_search(
    q: str = Query(default="", min_length=0),
    type: str = Query(default="all"),
    limit: int = Query(default=10, ge=1, le=50),
) -> list[dict]:
    """Search compounds, targets, and diseases; returns autocomplete results."""
    q_strip = q.strip()
    if len(q_strip) < 2:
        return []

    results: list[dict] = []
    q_lower = q_strip.lower()
    seen_ikeys: set[str] = set()
    seen_uniprots: set[str] = set()
    seen_diseases: set[str] = set()

    def _compound_subtitle(row, ikey):
        approved = str(row.get("is_approved", "False")).strip() == "True"
        if approved:
            return "Approved drug"
        phase = _norm_phase(row.get("clinical_phase", ""))
        return phase

    def _add_compound(ikey, row):
        if ikey in seen_ikeys:
            return
        seen_ikeys.add(ikey)
        results.append({
            "label": str(row.get("preferred_name", ikey)),
            "value": ikey,
            "category": "compound",
            "subtitle": _compound_subtitle(row, ikey),
        })

    def _add_target(uniprot_ac, row):
        if uniprot_ac in seen_uniprots:
            return
        seen_uniprots.add(uniprot_ac)
        results.append({
            "label": row.get("gene_name", uniprot_ac),
            "value": uniprot_ac,
            "category": "target",
            "subtitle": (row.get("protein_name", "") or "")[:60],
        })

    def _add_disease(dname, row):
        if dname in seen_diseases:
            return
        seen_diseases.add(dname)
        results.append({
            "label": dname,
            "value": dname,
            "category": "disease",
            "subtitle": row.get("category_name", ""),
        })

    # Only expose assay-active compounds (277K) in search results.
    # Compounds without assay data are not part of SynerGPCR's DB.
    if type in ("compound", "all") and "compounds" in DATA:
        cpd_df = DATA["compounds"]
    else:
        cpd_df = None
    tgt_df  = DATA["targets"]     if type in ("target",   "all") else None
    umls_df = DATA.get("umls_cats") if type in ("disease",  "all") else None

    # ── PASS 1: gene_name prefix (targets) ──────────────────────────────
    if tgt_df is not None:
        mask = tgt_df["gene_name"].str.lower().str.startswith(q_lower, na=False)
        for uniprot_ac, row in tgt_df[mask].head(limit).iterrows():
            if len(results) >= limit: break
            _add_target(uniprot_ac, row)

    # ── PASS 2: preferred_name prefix (compounds) ────────────────────────
    if cpd_df is not None:
        mask = cpd_df["preferred_name"].str.lower().str.startswith(q_lower, na=False)
        candidates = cpd_df[mask].copy()
        # Deprioritise entries whose preferred_name IS the InChIKey (no resolved name)
        candidates["_is_ikey"] = candidates.index == candidates["preferred_name"]
        candidates = candidates.sort_values("_is_ikey")
        for ikey, row in candidates.head(limit).iterrows():
            if len(results) >= limit: break
            _add_compound(ikey, row)

    # ── PASS 3: disease_name prefix ──────────────────────────────────────
    if umls_df is not None:
        mask = umls_df["disease_name"].str.lower().str.startswith(q_lower, na=False)
        for _, row in umls_df[mask].head(limit).iterrows():
            if len(results) >= limit: break
            _add_disease(row.get("disease_name", ""), row)

    # ── PASS 4: protein_name contains (targets) ──────────────────────────
    if tgt_df is not None and len(results) < limit:
        mask = tgt_df["protein_name"].str.lower().str.contains(
            q_lower, na=False, regex=False
        )
        for uniprot_ac, row in tgt_df[mask].head(limit).iterrows():
            if len(results) >= limit: break
            _add_target(uniprot_ac, row)

    # ── PASS 5: preferred_name contains (compounds) ──────────────────────
    if cpd_df is not None and len(results) < limit:
        mask = cpd_df["preferred_name"].str.lower().str.contains(
            q_lower, na=False, regex=False
        )
        candidates = cpd_df[mask].copy()
        candidates["_is_ikey"] = candidates.index == candidates["preferred_name"]
        candidates = candidates.sort_values("_is_ikey")
        for ikey, row in candidates.head(limit).iterrows():
            if len(results) >= limit: break
            _add_compound(ikey, row)

    # ── PASS 6: disease_name contains ────────────────────────────────────
    if umls_df is not None and len(results) < limit:
        mask = umls_df["disease_name"].str.lower().str.contains(
            q_lower, na=False, regex=False
        ) & ~umls_df["disease_name"].str.lower().str.startswith(q_lower, na=False)
        for _, row in umls_df[mask].head(limit).iterrows():
            if len(results) >= limit: break
            _add_disease(row.get("disease_name", ""), row)

    # ── PASS 7: InChIKey exact match only ────────────────────────────────
    # Only triggers when user pastes a full InChIKey (27 chars, two hyphens)
    if cpd_df is not None and len(results) < limit:
        q_upper = q_strip.upper()
        if (len(q_upper) == 27 and q_upper.count("-") == 2
                and q_upper in cpd_df.index):
            _add_compound(q_upper, cpd_df.loc[q_upper])

    # ── PASS 8: UniProt AC exact match ────────────────────────────────
    # Triggered when user pastes a UniProt accession (e.g. P14416).
    # UniProt ACs are 6 or 10 chars; format: letter + 5 alphanum
    # (reviewed) or letter + number + 3 alphanum + number + 5 alphanum
    if tgt_df is not None and len(results) < limit:
        q_upper8 = q_strip.upper()
        if q_upper8 in tgt_df.index:
            _add_target(q_upper8, tgt_df.loc[q_upper8])

    return results[:limit]


@app.get("/api/compound/{ikey}")
async def api_compound(ikey: str) -> dict:
    """Return full compound detail: chains, AI predictions, and synergy scoreboard."""
    if ikey not in DATA["compounds"].index:
        raise HTTPException(
            status_code=404,
            detail=f"Compound '{ikey}' not found in SynerGPCR"
        )

    row = DATA["compounds"].loc[ikey]
    targets_df = DATA["targets"]

    result: dict = {
        "name": row.get("preferred_name", ""),
        "smiles": row.get("smiles", ""),
        "is_approved": row.get("is_approved", "False"),
        "clinical_phase": _norm_phase(row.get("clinical_phase", "")),
        "chains": [],
        "ai_predictions": [],
        "synergy": {},
    }

    # ── Chain data ──
    current_pattern: str | None = None
    if ikey in DATA.get("chains_by_ikey", {}):
        chain_df = DATA["chains_by_ikey"][ikey]
        chain_records = []
        for _, cr in chain_df.iterrows():
            uniprot_ac = cr.get("UniProt_AC", "")
            gene_name = (
                targets_df.loc[uniprot_ac, "gene_name"]
                if uniprot_ac in targets_df.index
                else uniprot_ac
            )
            _pat = _compute_pair_chain_pattern(
                cr.get("L1_state", ""),
                cr.get("L2_state", ""),
                cr.get("L3_state", ""),
                cr.get("L4_state", ""),
            )
            _syn_row = DATA.get("synergy_map", {}).get(_pat, {})
            chain_records.append({
                "uniprot_ac": uniprot_ac,
                "gene_name": gene_name,
                "L1": cr.get("L1_state", ""),
                "L2": cr.get("L2_state", ""),
                "L3": cr.get("L3_state", ""),
                "L4": cr.get("L4_state", ""),
                "chain_pattern": _pat,
                "current_bliss_SI": str(_syn_row.get("bliss_SI", "")),
                "current_approval_rate_pct": str(_syn_row.get("approval_rate_pct", "")),
            })
        def _chain_sort_key(rec):
            vals = [rec["L1"], rec["L2"], rec["L3"], rec["L4"]]
            n_active = sum(1 for v in vals if str(v).strip() == "Active")
            n_any = sum(1 for v in vals if str(v).strip() not in ("", "nan", "None"))
            return (-n_active, -n_any, str(rec.get("gene_name", "")))
        chain_records.sort(key=_chain_sort_key)
        result["chains"] = chain_records
        current_pattern = _best_pattern(chain_records)

    # ── AI predictions ──
    predicted_pattern: str | None = None
    if ikey in DATA.get("pred_by_ikey", {}):
        pred_df = DATA["pred_by_ikey"][ikey]
        valid = pred_df[
            (pred_df["is_high_confidence"] == "True")
            & (pred_df["predicted_moa"] != "non-binder")
        ]
        for _, pr in valid.iterrows():
            result["ai_predictions"].append({
                "uniprot_ac": pr.get("UniProt_AC", ""),
                "layer": pr.get("layer", "L2"),
                "predicted_moa": pr.get("predicted_moa", ""),
                "confidence_raw": pr.get("confidence_raw", ""),
                "is_high_confidence": pr.get("is_high_confidence", "False"),
                "upgrades_to_L1L2L4": pr.get("upgrades_to_L1L2L4", "False"),
            })

        upgrades = pred_df[
            (pred_df["is_high_confidence"] == "True")
            & (pred_df["upgrades_to_L1L2L4"] == "True")
        ]
        if len(upgrades) > 0:
            predicted_pattern = "L1+L2+L4"

    # ── Synergy scoreboard ──
    current_syn = get_synergy_row(current_pattern) if current_pattern else None
    predicted_syn = get_synergy_row(predicted_pattern) if predicted_pattern else None

    result["synergy"] = {
        "current_pattern": current_pattern,
        "current_si": current_syn.get("bliss_SI") if current_syn else None,
        "current_ci_lo": current_syn.get("bliss_SI_ci_lo") if current_syn else None,
        "current_ci_hi": current_syn.get("bliss_SI_ci_hi") if current_syn else None,
        "current_approval_rate": current_syn.get("approval_rate_pct") if current_syn else None,
        "current_interpretation": current_syn.get("bliss_interpretation") if current_syn else None,
        "predicted_pattern": predicted_pattern,
        "predicted_si": predicted_syn.get("bliss_SI") if predicted_syn else None,
        "predicted_approval_rate": predicted_syn.get("approval_rate_pct") if predicted_syn else None,
        "fullchain_si": FULLCHAIN_SI,
        "fullchain_approval_rate": FULLCHAIN_APPROVAL_RATE,
    }

    return result


@app.get("/api/target/{uniprot_ac}")
async def api_target(uniprot_ac: str) -> dict:
    """Return target detail with compound counts."""
    if uniprot_ac not in DATA["targets"].index:
        raise HTTPException(status_code=404, detail=f"Target '{uniprot_ac}' not found")

    row = DATA["targets"].loc[uniprot_ac]

    if uniprot_ac in DATA.get("chains_by_target", {}):
        target_chains = DATA["chains_by_target"][uniprot_ac]
        n_compounds = len(target_chains)
        n_approved = int((target_chains["is_approved"] == "True").sum())
    else:
        n_compounds = 0
        n_approved = 0

    return {
        "uniprot_ac": uniprot_ac,
        "gene_name": _s(row.get("gene_name", "")),
        "protein_name": _s(row.get("protein_name", "")),
        "gpcr_class": _s(row.get("gpcr_class_name", "")),
        "gpcr_class_letter": _s(row.get("gpcr_class_letter", "")),
        "gpcr_family": _s(row.get("gpcr_family", "")),
        "organism": _s(row.get("organism", "")),
        "uniprot_entry_name": _s(row.get("uniprot_entry_name", "")),
        "chembl_target_id": _s(row.get("chembl_target_id", "")),
        "n_compounds": n_compounds,
        "n_approved": n_approved,
    }


@app.get("/api/structure/{ikey}")
async def api_structure(ikey: str) -> Response:
    """Return RDKit 2D SVG for a compound SMILES (image/svg+xml)."""
    if ikey not in DATA["compounds"].index:
        raise HTTPException(status_code=404, detail="Compound not found")

    smiles = str(DATA["compounds"].loc[ikey].get("smiles", "") or "").strip()
    if not smiles or smiles.lower() in ("nan", "none", ""):
        raise HTTPException(status_code=404, detail="No SMILES available")

    if not _RDKIT_AVAILABLE:
        raise HTTPException(status_code=503, detail="RDKit not installed")

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        drawer = _rdMolDraw2D.MolDraw2DSVG(360, 280)
        drawer.drawOptions().padding = 0.12
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace("encoding='iso-8859-1'", "encoding='UTF-8'")
        return Response(
            content=svg,
            media_type="image/svg+xml",
            headers={"Cache-Control": "max-age=86400"},
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Structure rendering failed: {exc}")


@app.get("/api/target/{uniprot_ac}/compounds")
async def api_target_compounds(
    uniprot_ac: str,
    limit: int = Query(default=500, ge=1, le=5000),
) -> list[dict]:
    """Return top compounds for a specific GPCR target, approved first then by chain rank."""
    if uniprot_ac not in DATA["targets"].index:
        raise HTTPException(status_code=404, detail=f"Target '{uniprot_ac}' not found")

    if uniprot_ac not in DATA.get("chains_by_target", {}):
        return []

    tc = DATA["chains_by_target"][uniprot_ac].copy()
    tc["_rank"] = tc["chain_pattern"].map(lambda p: _PATTERN_RANK.get(str(p), 0))
    tc["_approved_flag"] = (tc["is_approved"] == "True").astype(int)
    tc = tc.sort_values(["_approved_flag", "_rank"], ascending=False)

    compounds_df = DATA["compounds"]
    results: list[dict] = []
    for _, row in tc.head(limit).iterrows():
        ikey = row.get("InChIKey", "")
        name = ikey
        phase = ""
        if ikey in compounds_df.index:
            cpd = compounds_df.loc[ikey]
            name = str(cpd.get("preferred_name", ikey))
            phase = _norm_phase(cpd.get("clinical_phase", ""))
        _pat2 = _compute_pair_chain_pattern(
            row.get("L1_state", ""),
            row.get("L2_state", ""),
            row.get("L3_state", ""),
            row.get("L4_state", ""),
        )
        _syn2 = DATA.get("synergy_map", {}).get(_pat2, {})
        results.append({
            "ikey": ikey,
            "name": name,
            "is_approved": row.get("is_approved", "False"),
            "clinical_phase": phase,
            "chain_pattern": _pat2,
            "L1": row.get("L1_state", ""),
            "L2": row.get("L2_state", ""),
            "L3": row.get("L3_state", ""),
            "L4": row.get("L4_state", ""),
            "current_bliss_SI": str(_syn2.get("bliss_SI", "")),
            "current_approval_rate_pct": str(_syn2.get("approval_rate_pct", "")),
        })
    return results


@app.get("/api/targets")
async def api_targets() -> list[dict]:
    """Return all GPCR targets with compound counts (for browse page)."""
    df = DATA["targets"].reset_index()
    chains_by_target = DATA.get("chains_by_target", {})
    rows = []
    for _, row in df.iterrows():
        uniprot_ac = row.get("UniProt_AC", "")
        if uniprot_ac in chains_by_target:
            tc = chains_by_target[uniprot_ac]
            n_cpds = len(tc)
            n_appr = int((tc["is_approved"] == "True").sum())
        else:
            n_cpds = 0
            n_appr = 0
        rows.append({
            "uniprot_ac": uniprot_ac,
            "gene_name": _s(row.get("gene_name", "")),
            "protein_name": _s(row.get("protein_name", ""))[:60],
            "gpcr_class_letter": _s(row.get("gpcr_class_letter", "")),
            "gpcr_class_name": _s(row.get("gpcr_class_name", "")),
            "gpcr_family": _s(row.get("gpcr_family", "")),
            "organism": _s(row.get("organism", "")),
            "uniprot_entry_name": _s(row.get("uniprot_entry_name", "")),
            "n_compounds": n_cpds,
            "n_approved": n_appr,
        })
    return rows


@app.get("/api/compounds/approved")
async def api_compounds_approved(
    limit: int = Query(default=100, ge=1, le=5000),
    offset: int = Query(default=0, ge=0),
) -> dict:
    """Return paginated approved compounds with best chain pattern."""
    compounds_df = DATA["compounds"]
    approved = compounds_df[compounds_df["is_approved"] == "True"].copy()
    total = len(approved)
    n_with_assay = int((approved["in_assay_universe"] == "True").sum())

    # Sort: assay-bearing first, then named before InChIKey-only, then A→Z
    approved["_no_assay"] = (approved["in_assay_universe"] != "True").astype(int)
    approved["_name_is_ikey"] = (approved["preferred_name"] == approved.index).astype(int)
    approved["_sort_key"] = approved["preferred_name"].str.lower()
    approved = approved.sort_values(["_no_assay", "_name_is_ikey", "_sort_key"])

    page = approved.iloc[offset : offset + limit]
    results: list[dict] = []
    for ikey, row in page.iterrows():
        in_assay = str(row.get("in_assay_universe", "False")).strip() == "True"
        best_pat = None
        n_targets = None
        if in_assay and ikey in DATA.get("chains_by_ikey", {}):
            cdf = DATA["chains_by_ikey"][ikey]
            n_targets = len(cdf)
            bp = _best_pattern([
                {"L1": r.get("L1_state",""), "L2": r.get("L2_state",""),
                 "L3": r.get("L3_state",""), "L4": r.get("L4_state","")}
                for _, r in cdf.iterrows()
            ])
            best_pat = bp if bp else "—"
        elif in_assay:
            n_targets = 0
            best_pat = "—"
        name = str(row.get("preferred_name", ikey))
        phase = _norm_phase(row.get("clinical_phase", ""))
        results.append({
            "ikey": ikey,
            "name": name,
            "clinical_phase": phase,
            "in_assay_universe": in_assay,
            "chain_pattern": best_pat,
            "n_targets": n_targets,
        })
    return {
        "total": total,
        "n_with_assay": n_with_assay,
        "offset": offset,
        "limit": limit,
        "compounds": results,
    }


@app.get("/api/synergy_stats")
async def api_synergy_stats() -> list[dict]:
    """Return synergy_stats table with NaN cleaned to null."""
    return _sanitize_records(DATA["synergy"].to_dict("records"))


@app.get("/api/stats")
async def api_stats() -> dict:
    """Return global database statistics for the browse page header."""
    n_assay = int((DATA["compounds"]["in_assay_universe"] == "True").sum()) \
              if "compounds" in DATA else 0
    n_approved = int(
        ((DATA["compounds"]["in_assay_universe"] == "True") &
         (DATA["compounds"]["is_approved"] == "True")).sum()
    ) if "compounds" in DATA else 0
    n_approved_db = int((DATA["compounds"]["is_approved"] == "True").sum()) \
                   if "compounds" in DATA else 0
    n_targets = len(DATA.get("targets", []))
    n_ai_targets = int(
        DATA["predictions"][DATA["predictions"]["is_high_confidence"] == "True"]["UniProt_AC"].nunique()
    ) if "predictions" in DATA and len(DATA["predictions"]) > 0 else 0
    n_disease_cats = len(DATA.get("disease_cats", []))
    n_enriched_targets = int(
        (DATA["target_enrichment"]["enough_n_exp"] == "True").sum()
    ) if "target_enrichment" in DATA else 0
    n_candidates = len(DATA.get("candidates", []))
    n_repurposing_web = DATA["repurposing_crossfam"]["inchikey"].nunique() \
                        if "repurposing_crossfam" in DATA else 0
    return {
        "n_assay_compounds": n_assay,
        "n_approved_db": n_approved_db,
        "n_approved_assay": n_approved,
        "n_targets": n_targets,
        "n_ai_targets": n_ai_targets,
        "n_disease_cats": n_disease_cats,
        "n_enriched_targets": n_enriched_targets,
        "n_candidates": n_candidates,
        "n_repurposing_web": n_repurposing_web,
    }


# ── 8b. New API routes ────────────────────────────────────────────────────────

@app.get("/api/single_layer")
async def api_single_layer() -> list[dict]:
    """Return single-layer table with NaN cleaned to null."""
    if "single_layer" not in DATA:
        raise HTTPException(status_code=503, detail="single_layer data not loaded")
    return _sanitize_records(DATA["single_layer"].to_dict("records"))


@app.get("/api/gpcract_3point")
async def api_gpcract_3point() -> list[dict]:
    """Return GPCRact 3-point SI panel data (L2 and L3 rows)."""
    if "gpcract_3point" not in DATA:
        raise HTTPException(status_code=503, detail="gpcract_3point data not loaded")
    return _sanitize_records(DATA["gpcract_3point"].to_dict("records"))


@app.get("/api/candidates")
async def api_candidates(
    layer: Optional[str] = None,
    is_approved: Optional[str] = None,
    limit: int = Query(default=100, ge=1, le=2000),
    offset: int = Query(default=0, ge=0),
) -> dict:
    """Return paginated AI candidate shortlist with optional filters."""
    if "candidates" not in DATA:
        raise HTTPException(status_code=503, detail="candidates data not loaded")
    df = DATA["candidates"].copy()
    if layer:
        df = df[df["ai_filled_layers"] == layer]
    if is_approved == "true":
        df = df[df["is_approved"] == "True"]
    elif is_approved == "false":
        df = df[df["is_approved"] != "True"]
    total = len(df)
    page = df.iloc[offset:offset + limit]
    return {"total": total, "offset": offset, "candidates": _sanitize_records(page.to_dict("records"))}


@app.get("/api/repurposing")
async def api_repurposing(
    repurposing_class: Optional[str] = None,
    limit: int = Query(default=100, ge=1, le=2000),
    offset: int = Query(default=0, ge=0),
) -> dict:
    """Return paginated repurposing candidates with optional class filter."""
    if "repurposing_crossfam" not in DATA:
        raise HTTPException(status_code=503, detail="repurposing data not loaded")
    df = DATA["repurposing_crossfam"].copy()
    if repurposing_class:
        df = df[df["repurposing_class"] == repurposing_class]
    total = len(df)
    page = df.iloc[offset:offset + limit]
    return {"total": total, "offset": offset, "repurposing": _sanitize_records(page.to_dict("records"))}


@app.get("/api/compound/{ikey}/repurposing")
async def api_compound_repurposing(ikey: str) -> list[dict]:
    """Return repurposing rows for a specific compound."""
    if "repurposing_full" not in DATA:
        return []
    df = DATA["repurposing_full"]
    rows = df[df["inchikey"] == ikey]
    if rows.empty:
        return []
    return _sanitize_records(rows.to_dict("records"))


@app.get("/api/target/{uniprot_ac}/enrichment")
async def api_target_enrichment(uniprot_ac: str) -> dict:
    """Return clinical enrichment row for a target."""
    enr_map = DATA.get("target_enrichment_map", {})
    row = enr_map.get(uniprot_ac)
    if row is None:
        return {}
    return _sanitize_records([row])[0]


@app.get("/api/target/{uniprot_ac}/disease")
async def api_target_disease(uniprot_ac: str) -> list[dict]:
    """Return disease enrichment rows for a target."""
    if "target_disease" not in DATA:
        return []
    df = DATA["target_disease"]
    rows = df[df["UniProt"] == uniprot_ac]
    if rows.empty:
        return []
    return _sanitize_records(rows.to_dict("records"))


@app.get("/api/target/{uniprot_ac}/disease_cats")
async def api_target_disease_cats(uniprot_ac: str) -> list[str]:
    """Return disease categories this target has enrichment data for."""
    if "target_disease" not in DATA:
        return []
    df = DATA["target_disease"]
    rows = df[df["UniProt"] == uniprot_ac]
    if rows.empty:
        return []
    return sorted(rows["disease_category"].dropna().unique().tolist())


@app.get("/api/disease_category")
async def api_disease_category(
    category: str = Query(..., min_length=1),
) -> dict:
    """Return summary stats and disease list for a disease category."""
    di = DATA.get("diseases")
    ti = DATA.get("target_indication")
    umls = DATA.get("umls_cats")

    if di is None:
        raise HTTPException(status_code=503, detail="drug_indication data not loaded")

    # Filter by category
    cat_di = di[di["category_name"] == category] if "category_name" in di.columns else di.iloc[0:0]
    cat_ti = ti[ti["category_name"] == category] if (
        ti is not None and "category_name" in ti.columns
    ) else pd.DataFrame()

    n_drugs    = int(cat_di["InChIKey"].nunique())   if not cat_di.empty else 0
    n_targets  = int(cat_ti["UniProt_AC"].nunique()) if not cat_ti.empty else 0

    # Build per-disease summary
    disease_rows: list[dict] = []
    if not cat_di.empty and "disease_name" in cat_di.columns:
        for dname, grp in cat_di.groupby("disease_name"):
            n_d_drugs = int(grp["InChIKey"].nunique())
            # Count targets for this disease from target_indication
            n_d_tgts = 0
            if not cat_ti.empty and "disease_name" in cat_ti.columns:
                t_hits = cat_ti[
                    cat_ti["disease_name"].str.lower() == str(dname).lower()
                ]
                n_d_tgts = int(t_hits["UniProt_AC"].nunique())
            disease_rows.append({
                "disease_name": _s(dname),
                "n_drugs":   n_d_drugs,
                "n_targets": n_d_tgts,
            })

    # Sort by n_drugs descending, take top 30
    disease_rows.sort(key=lambda x: x["n_drugs"], reverse=True)

    # Unique disease count from UMLS lookup if available
    n_diseases = len(disease_rows)

    return {
        "category":   category,
        "n_diseases": n_diseases,
        "n_drugs":    n_drugs,
        "n_targets":  n_targets,
        "diseases":   disease_rows[:30],
    }


@app.get("/api/disease_targets")
async def api_disease_targets(
    category: str = Query(..., min_length=1),
) -> list[dict]:
    """Return all target enrichment rows for a disease category."""
    if "target_disease" not in DATA:
        return []
    df = DATA["target_disease"]
    rows = df[df["disease_category"] == category]
    if rows.empty:
        return []
    cols = ["UniProt", "GPCR_name", "disease_category", "n_co_active_exp",
            "n_approved_cat", "approval_rate_cond", "log2_enrichment",
            "ci_robust_backend"]
    available = [c for c in cols if c in rows.columns]
    return _sanitize_records(rows[available].to_dict("records"))


@app.get("/api/disease/{disease_name:path}")
async def api_disease_detail(disease_name: str) -> dict:
    """Return disease detail: UMLS CUI, category, drugs, targets."""
    ind_df  = DATA.get("diseases")       # drug_indication.csv
    ti_df   = DATA.get("target_indication")
    umls_df = DATA.get("umls_cats")
    cpd_df  = DATA.get("compounds")      # compound_lookup_final (indexed by InChIKey)

    if ind_df is None:
        raise HTTPException(status_code=503, detail="indication data not loaded")

    # ── Resolve display name and UMLS CUI ───────────────────────────────
    disease_display = disease_name.strip()
    umls_cui = ""
    category = ""

    if umls_df is not None:
        dn_lower = disease_display.lower()
        exact    = umls_df[umls_df["disease_name"].str.lower() == dn_lower]
        prefix   = umls_df[umls_df["disease_name"].str.lower().str.startswith(dn_lower, na=False)]
        contains = umls_df[umls_df["disease_name"].str.lower().str.contains(dn_lower, na=False, regex=False)]
        match    = exact if not exact.empty else (
                   prefix if not prefix.empty else contains)
        if not match.empty:
            row0 = match.iloc[0]
            disease_display = row0.get("disease_name", disease_display)
            umls_cui        = str(row0.get("umls_cui", "") or "")
            category        = str(row0.get("category_name", "") or "")

    if umls_df is not None and not category:
        raise HTTPException(status_code=404, detail=f"Disease '{disease_name}' not found")

    # ── Drug indication rows for this disease ───────────────────────────
    dn_lower = disease_display.lower()
    ind_hits = ind_df[
        ind_df["disease_name"].str.lower() == dn_lower
    ]
    if ind_hits.empty:
        ind_hits = ind_df[
            ind_df["disease_name"].str.lower().str.contains(
                dn_lower, na=False, regex=False)
        ]

    # ── Build drug list: SynerGPCR compounds only ───────────────────────
    known_ikeys = set(cpd_df.index) if cpd_df is not None else set()

    drug_rows: list[dict] = []
    seen_ikeys: set[str] = set()

    for _, row in ind_hits.iterrows():
        ikey = str(row.get("InChIKey", "") or "").strip()
        if not ikey or ikey in seen_ikeys:
            continue
        if ikey not in known_ikeys:
            continue   # skip compounds not in SynerGPCR
        seen_ikeys.add(ikey)

        # Use compound_lookup_final for accurate name and phase
        cpd_row = cpd_df.loc[ikey]
        drug_name = str(cpd_row.get("preferred_name", ikey) or ikey)
        is_approved = str(cpd_row.get("is_approved", "False")).strip() == "True"
        phase = _norm_phase(cpd_row.get("clinical_phase", ""))

        if is_approved:
            clinical_status = "Approved"
        elif phase:
            # Normalise: strip any leading "Phase " prefix, then rebuild cleanly
            phase_clean = phase.strip()
            if phase_clean.lower().startswith("phase "):
                phase_clean = phase_clean[6:].strip()
            # If remaining value is a word like "Investigational" → use as-is
            if phase_clean.lower() in ("investigational", "preclinical",
                                        "discovery", "research"):
                clinical_status = phase_clean.capitalize()
            else:
                clinical_status = f"Phase {phase_clean}"
        else:
            clinical_status = "Investigational"

        drug_rows.append({
            "inchikey":        ikey,
            "drug_name":       drug_name,
            "clinical_status": clinical_status,
            "is_approved":     is_approved,
            "in_assay":        str(cpd_row.get("in_assay_universe", "False")).strip() == "True",
        })

    # Sort: Approved first, then alphabetically
    drug_rows.sort(key=lambda x: (0 if x["is_approved"] else 1,
                                   x["drug_name"].lower()))
    n_drugs = len(drug_rows)

    # ── Target list ─────────────────────────────────────────────────────
    target_rows: list[dict] = []
    if ti_df is not None:
        ti_hits = ti_df[
            ti_df["disease_name"].str.lower() == dn_lower
        ]
        if ti_hits.empty:
            ti_hits = ti_df[
                ti_df["disease_name"].str.lower().str.contains(
                    dn_lower, na=False, regex=False)
            ]
        seen_uniprots: set[str] = set()
        tgt_df = DATA.get("targets")
        for _, row in ti_hits.iterrows():
            uniprot = str(row.get("UniProt_AC", "") or "").strip()
            if not uniprot or uniprot in seen_uniprots:
                continue
            seen_uniprots.add(uniprot)
            gene = ""
            if tgt_df is not None and uniprot in tgt_df.index:
                gene = _s(tgt_df.loc[uniprot].get("gene_name", ""))
            target_rows.append({
                "uniprot_ac": uniprot,
                "gene_name":  gene or uniprot,
            })
        n_targets = len(target_rows)
    else:
        n_targets = 0

    return {
        "disease_name":  disease_display,
        "umls_cui":      umls_cui,
        "category":      category,
        "n_drugs":       n_drugs,
        "n_targets":     n_targets,
        "drugs":         drug_rows,
        "targets":       target_rows,
    }


@app.get("/api/target/{uniprot_ac}/fallback")
async def api_target_fallback(uniprot_ac: str) -> dict:
    """Return fallback info for a target not in the SynerGPCR-600 set."""
    # If target IS in our database, no fallback needed
    if "targets" in DATA and uniprot_ac in DATA["targets"].index:
        return {}

    if "gpcr_fallback" not in DATA:
        return {"in_database": False, "same_family_targets": []}

    fb = DATA["gpcr_fallback"]
    match = fb[fb["uniprot_ac"] == uniprot_ac]
    if match.empty:
        return {"in_database": False, "same_family_targets": []}

    query_row = match.iloc[0]
    gene_name = str(query_row.get("gene_name", "") or "").strip()
    family = str(query_row.get("gpcr_family", "") or "").strip()

    # Find same-family targets that ARE in SynerGPCR-600
    same_fam = fb[
        (fb["gpcr_family"] == family)
        & (fb["in_synergpcr_600"] == "True")
    ]
    family_targets = []
    for _, r in same_fam.iterrows():
        family_targets.append({
            "uniprot_ac": r.get("uniprot_ac", ""),
            "gene_name": r.get("gene_name", ""),
            "protein_name": r.get("protein_name", ""),
        })

    result: dict = {
        "query_uniprot": uniprot_ac,
        "query_gene": gene_name if gene_name.lower() not in ("nan", "none", "") else "",
        "query_family": family,
        "in_database": False,
        "same_family_targets": family_targets,
    }
    if gene_name and gene_name.lower() not in ("nan", "none", ""):
        result["gpcrdb_url"] = f"https://gpcrdb.org/protein/{gene_name}/"
    return result


class TanimotoRequest(BaseModel):
    smiles: str
    top_n:  int   = 10
    threshold: float = 0.5


@app.post("/api/tanimoto")
async def api_tanimoto(body: TanimotoRequest) -> dict:
    """SMILES → ECFP4 Tanimoto similarity search.

    Returns an exact_match redirect when the query molecule is
    already in the database, otherwise returns the top-N similar
    compounds above the similarity threshold.
    """
    if not _RDKIT_AVAILABLE:
        raise HTTPException(
            status_code=503,
            detail="RDKit is not installed on this server."
        )
    if not FP_INDEX["ready"]:
        raise HTTPException(
            status_code=503,
            detail="Fingerprint index is not ready yet. "
                   "Please retry in a few seconds."
        )

    # ── Parse query SMILES ─────────────────────────────────────────
    smiles = body.smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid SMILES string. "
                   "Please check your input and try again."
        )

    # ── Exact-match check via InChIKey ─────────────────────────────
    try:
        inchi_str = rdInchi.MolToInchi(mol)
        if inchi_str:
            ikey = rdInchi.InchiToInchiKey(inchi_str)
            _assay_ikeys = set(FP_INDEX["ikeys"])
            if ikey and ikey in _assay_ikeys:
                row = DATA["compounds"].loc[ikey]
                return {
                    "exact_match": True,
                    "ikey": ikey,
                    "name": str(row.get("preferred_name", "") or ""),
                    "similar_compounds": [],
                }
    except Exception:
        pass   # fall through to similarity search

    # Compute query fingerprint in the async context (main thread)
    # before entering the executor — mol is NOT thread-safe.
    query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

    def _run_tanimoto(qfp) -> list[dict]:
        """Thread-safe: receives a pre-computed fingerprint, not mol."""
        sims = DataStructs.BulkTanimotoSimilarity(
            qfp, FP_INDEX["fps"]
        )
        threshold = body.threshold
        top_n     = body.top_n
        cpd_df    = DATA["compounds"]

        hits = []
        for ikey, sim in zip(FP_INDEX["ikeys"], sims):
            if sim >= threshold:
                row = cpd_df.loc[ikey]
                hits.append({
                    "ikey":           ikey,
                    "name":           str(row.get("preferred_name", "") or ""),
                    "similarity":     round(float(sim), 4),
                    "is_approved":    str(row.get("is_approved", "False")),
                    "clinical_phase": _norm_phase(
                                          row.get("clinical_phase", "")
                                      ),
                })

        hits.sort(key=lambda x: -x["similarity"])
        return hits[:top_n]

    hits = await asyncio.get_running_loop().run_in_executor(
        None, _run_tanimoto, query_fp
    )

    return {
        "exact_match":       False,
        "ikey":              None,
        "similar_compounds": hits,
    }


DOWNLOAD_DIR = DATA_DIR / "download"


@app.get("/download/file/{filename}")
async def download_file(filename: str):
    """Serve a curated download file from the download directory."""
    ALLOWED = {
        "compound_lookup.csv",
        "target_lookup.csv",
        "compound_chain_summary.csv",
        "gpcract_predictions.csv",
        "candidate_shortlist.csv",

        "synergy_stats.csv",
        "single_layer_table.csv",
    }
    if filename not in ALLOWED:
        raise HTTPException(status_code=404, detail="File not found")
    fpath = DOWNLOAD_DIR / filename
    if not fpath.exists():
        raise HTTPException(status_code=503,
                            detail="File not yet generated")
    return FileResponse(
        path=str(fpath),
        media_type="text/csv",
        filename=filename,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'}
    )


# ── 9. Exception handlers ──────────────────────────────────────────────────────

@app.exception_handler(404)
async def not_found_handler(request: Request, exc: HTTPException) -> HTMLResponse:
    """Render custom 404 page."""
    detail = getattr(exc, "detail", "Page not found")
    return templates.TemplateResponse(
        request,
        "404.html",
        {"detail": detail},
        status_code=404,
    )
