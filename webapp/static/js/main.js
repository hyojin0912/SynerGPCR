/* SynerGPCR — main.js */
'use strict';

// ── Autocomplete core ──────────────────────────────────────────────────────────

const CATEGORY_CLASS = {
  compound: 'badge-compound',
  target:   'badge-target',
  disease:  'badge-disease',
};

function hideDropdown(dd) {
  if (!dd) return;
  dd.innerHTML = '';
  dd.classList.remove('sg-search-dropdown--open');
}

function showDropdown(dd) {
  if (!dd) return;
  dd.classList.add('sg-search-dropdown--open');
}

function renderDropdown(results, dd) {
  if (!dd) return;
  if (!results || results.length === 0) { hideDropdown(dd); return; }

  dd.innerHTML = results.map(r => {
    const badgeCls = CATEGORY_CLASS[r.category] || '';
    const sub = r.subtitle
      ? `<div class="sg-dd-sub">${escHtml(r.subtitle)}</div>`
      : '';
    return `<div class="sg-dd-item" data-category="${escAttr(r.category)}" data-value="${escAttr(r.value)}">
      <div class="sg-dd-row">
        <span class="sg-dd-label">${escHtml(r.label)}</span>
        <span class="search-badge ${badgeCls}">${escHtml(r.category)}</span>
      </div>
      ${sub}
    </div>`;
  }).join('');

  dd.querySelectorAll('.sg-dd-item').forEach(item => {
    item.addEventListener('mousedown', e => {
      e.preventDefault(); // prevent input blur before click
      const cat = item.dataset.category;
      const val = item.dataset.value;
      if      (cat === 'compound') window.location.href = `/compound/${encodeURIComponent(val)}`;
      else if (cat === 'target')   window.location.href = `/browse?target=${encodeURIComponent(val)}`;
      else if (cat === 'disease')  window.location.href = `/disease/${encodeURIComponent(val)}`;
    });
  });

  showDropdown(dd);
}

/**
 * Heuristic: input is treated as SMILES only when it contains
 * at least one character that is syntactically required by SMILES
 * but never appears in compound names, gene symbols, or disease terms.
 * The second (length-based) condition is removed because long
 * mixed-case strings (e.g. protein names) would otherwise be
 * misclassified.
 */
function looksLikeSmiles(q) {
    return /[=\#\(\)\[\]\/\\@\+]/.test(q);
}

/**
 * Show a "not found" guidance panel when a text query returns
 * zero results and does not look like a SMILES string.
 * Reuses the same #sg-tanimoto-panel element as the Tanimoto results.
 */
function showNotFoundGuidance(query) {
    let panel = document.getElementById('sg-tanimoto-panel');
    if (!panel) {
        panel = document.createElement('div');
        panel.id = 'sg-tanimoto-panel';
        panel.className = 'sg-tanimoto-panel mt-2';
        const activeInput =
            document.getElementById('hero-search') ||
            document.getElementById('global-search') ||
            document.querySelector('input[type="search"]');
        if (!activeInput) return;
        const wrapper = activeInput.closest('.position-relative, .input-group, form, div');
        if (wrapper) { wrapper.after(panel); } else { activeInput.after(panel); }
    }
    panel.style.display = 'block';
    panel.innerHTML = `
        <div class="alert alert-secondary mt-2 mb-0" style="font-size:0.88rem;">
          <strong>&ldquo;${escHtml(query)}&rdquo;</strong>
          was not found in SynerGPCR.<br>
          <span class="text-muted">
            If you know the SMILES for this compound, paste it into the
            search box to find structurally similar compounds in our database.
          </span>
          <div class="mt-2">
            <a href="/browse" class="btn btn-sm btn-outline-secondary me-1">
              Browse all GPCR targets
            </a>
            <a href="/download" class="btn btn-sm btn-outline-secondary">
              Download full compound list
            </a>
          </div>
        </div>`;
}

async function runTanimotoFallback(smiles) {
    // Always render Tanimoto results in a dedicated visible panel.
    // Create it below the active search input if it doesn't exist.
    let panel = document.getElementById('sg-tanimoto-panel');
    if (!panel) {
        panel = document.createElement('div');
        panel.id = 'sg-tanimoto-panel';
        panel.className = 'sg-tanimoto-panel mt-2';
        // Insert after whichever search input is active
        const activeInput =
            document.getElementById('hero-search') ||
            document.getElementById('global-search') ||
            document.querySelector('input[type="search"]');
        if (!activeInput) return;
        const wrapper = activeInput.closest('.position-relative, .input-group, form, div');
        if (wrapper) {
            wrapper.after(panel);
        } else {
            activeInput.after(panel);
        }
    }
    panel.style.display = 'block';

    panel.innerHTML = `
        <div class="sg-tanimoto-searching p-3 text-muted">
          <span class="spinner-border spinner-border-sm me-2"></span>
          Searching for structurally similar compounds
          (ECFP4 Tanimoto ≥ 0.5)…
        </div>`;

    try {
        const resp = await fetch('/api/tanimoto', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({ smiles, top_n: 10, threshold: 0.5 }),
        });

        if (!resp.ok) {
            const err = await resp.json().catch(() => ({}));
            panel.innerHTML = `
                <div class="alert alert-warning mt-2">
                  <strong>SMILES search failed:</strong>
                  ${err.detail || resp.statusText}
                </div>`;
            return;
        }

        const data = await resp.json();

        // ── Exact match → redirect straight to compound page ──
        if (data.exact_match && data.ikey) {
            window.location.href = `/compound/${data.ikey}`;
            return;
        }

        // ── No similar compounds found ────────────────────────
        if (!data.similar_compounds || data.similar_compounds.length === 0) {
            panel.innerHTML = `
                <div class="alert alert-secondary mt-2">
                  No compounds with Tanimoto similarity ≥ 0.5 were found
                  in SynerGPCR for the provided SMILES.
                  Try a different structure or browse GPCR targets directly.
                </div>`;
            return;
        }

        // ── Render similar compounds ──────────────────────────
        const rows = data.similar_compounds.map(c => {
            const badge = c.is_approved === 'True'
                ? '<span class="badge bg-success ms-1">Approved</span>'
                : `<span class="badge bg-secondary ms-1">${c.clinical_phase || 'Investigational'}</span>`;
            return `
                <a href="/compound/${c.ikey}"
                   class="list-group-item list-group-item-action d-flex
                          justify-content-between align-items-center">
                  <span>
                    <strong>${c.name || c.ikey}</strong>
                    ${badge}
                  </span>
                  <span class="text-muted" style="font-size:0.82rem;white-space:nowrap;">
                    Tanimoto ${c.similarity.toFixed(3)}
                  </span>
                </a>`;
        }).join('');

        panel.innerHTML = `
            <div class="mt-2">
              <div class="text-muted mb-1" style="font-size:0.82rem;">
                No exact match found. Showing structurally similar compounds
                in SynerGPCR (ECFP4 Tanimoto ≥ 0.5):
              </div>
              <div class="list-group list-group-flush">${rows}</div>
            </div>`;

    } catch (e) {
        panel.innerHTML = `
            <div class="alert alert-danger mt-2">
              Tanimoto search error: ${e.message}
            </div>`;
    }
}

function escHtml(str) {
  return String(str)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function escAttr(str) { return String(str).replace(/"/g,'&quot;'); }

// ── Wire a search input to a dropdown ────────────────────────────────────────

function wireSearch(inputId, dropdownId) {
  const input = document.getElementById(inputId);
  const dd    = document.getElementById(dropdownId);
  if (!input || !dd) return;

  let timer = null;

  input.addEventListener('input', function () {
    clearTimeout(timer);
    const q = this.value.trim();
    if (q.length < 2) { hideDropdown(dd); return; }

    timer = setTimeout(() => {
      fetch(`/api/search?q=${encodeURIComponent(q)}&type=all&limit=8`)
        .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .then(data => renderDropdown(data, dd))
        .catch(() => hideDropdown(dd));
    }, 300);
  });

  // Enter key → navigate to first result, or SMILES fallback
  input.addEventListener('keydown', function (e) {
    if (e.key === 'Enter') {
      const first = dd.querySelector('.sg-dd-item');
      if (first) {
        first.dispatchEvent(new MouseEvent('mousedown', { bubbles: true }));
      } else {
        // No autocomplete result — try SMILES fallback or show guidance
        const q = this.value.trim();
        if (q.length < 2) return;
        hideDropdown(dd);
        if (looksLikeSmiles(q)) {
          dd.style.display = 'none';
          runTanimotoFallback(q);
        } else {
          showNotFoundGuidance(q);
        }
      }
    }
    if (e.key === 'Escape') hideDropdown(dd);
  });

  // Close on blur (but mousedown fires first, so we use a small delay)
  input.addEventListener('blur', () => {
    setTimeout(() => hideDropdown(dd), 150);
  });

  // Close when clicking outside
  document.addEventListener('click', e => {
    if (!input.contains(e.target) && !dd.contains(e.target)) hideDropdown(dd);
  });
}

// ── Hero search button ────────────────────────────────────────────────────────

function wireHeroButton() {
  const btn = document.getElementById('hero-search-btn');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const input = document.getElementById('hero-search');
    if (!input) return;
    const q = input.value.trim();
    if (q.length < 2) return;
    const dd = document.getElementById('hero-search-dropdown');
    fetch(`/api/search?q=${encodeURIComponent(q)}&type=all&limit=8`)
      .then(r => r.json())
      .then(data => {
        // ── SMILES fallback on explicit button submit ──────────────────
        if (data.length === 0) {
          hideDropdown(dd);
          if (dd) dd.style.display = 'none';
          if (looksLikeSmiles(q)) {
            runTanimotoFallback(q);
          } else {
            showNotFoundGuidance(q);
          }
          return;
        }
        renderDropdown(data, dd);
      })
      .catch(() => {});
  });
}

// ── Target class badge helper (used by browse page table) ───────────────────

const CLASS_COLORS = {
  A: { bg: '#EFF6FF', color: '#1D4ED8' },  // Class A Rhodopsin — blue
  B: { bg: '#F0FDF4', color: '#15803D' },  // Class B Secretin  — green
  C: { bg: '#FFF7ED', color: '#C2410C' },  // Class C Glutamate  — orange
  F: { bg: '#FDF4FF', color: '#7E22CE' },  // Class F Frizzled   — purple
};

/**
 * Return a styled class badge HTML string for a GPCR class letter.
 * @param {string} letter - e.g. 'A', 'B', 'C', 'F'
 * @returns {string} HTML badge
 */
function gpcr_class_badge(letter) {  // eslint-disable-line camelcase
  const s = (letter || '?').toUpperCase();
  const c = CLASS_COLORS[s] || { bg: '#F8FAFC', color: '#64748B' };
  return `<span class="badge" style="background:${c.bg};color:${c.color};font-weight:700;">${s}</span>`;
}

// ── Init ──────────────────────────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', () => {
  wireSearch('global-search', 'search-dropdown');
  wireSearch('hero-search',   'hero-search-dropdown');
  wireHeroButton();
});
