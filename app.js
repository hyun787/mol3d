import { GROUPS } from './chem.js';

const $ = id => document.getElementById(id);
const empty = $('empty'), q = $('q');
const state = { style: 'ball', h: true, hl: -1 };
let cur = null, home = null, token = 0, seq = 0;

const showEmpty = html => { empty.innerHTML = html; empty.hidden = false; };
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

if (typeof $3Dmol === 'undefined') {
  showEmpty('<span class="err">The 3D viewer did not load. Check the connection and reload.</span>');
  throw new Error('3Dmol.js missing');
}
const viewer = $3Dmol.createViewer($('stage'), { backgroundColor: tok('--stage'), antialias: true });

// Chemistry runs in a worker. A newer request settles the older one as stale.
const worker = new Worker('worker.js', { type: 'module' });
let pending = null;
worker.onmessage = ({ data }) => { if (pending?.id === data.id) pending.resolve(data); };
worker.onerror = () => pending?.resolve({ error: 'The chemistry library did not load. Check the connection and reload.' });
const build = smiles => new Promise(resolve => {
  pending?.resolve({ stale: true });
  pending = { id: ++token, resolve };
  worker.postMessage({ id: token, smiles });
});

// Name → SMILES. OPSIN parses systematic (IUPAC) names; PubChem covers common and trade names.
async function resolveName(name) {
  try {
    const r = await fetch(`https://www.ebi.ac.uk/opsin/ws/${encodeURIComponent(name)}.json`);
    const j = r.ok ? await r.json() : null;
    if (j?.smiles) return { smiles: j.smiles, title: name, source: 'opsin' };
  } catch { /* OPSIN down or name rejected: try PubChem */ }
  const r = await fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/SMILES,IUPACName,Title/JSON',
    { method: 'POST', body: new URLSearchParams({ name }) });
  if (r.status === 404) return null;
  if (!r.ok) throw new Error(r.status === 503 ? 'PubChem is busy right now. Wait a few seconds and try again.' : `PubChem error (${r.status}).`);
  const p = (await r.json()).PropertyTable.Properties[0];
  return { smiles: p.SMILES, title: p.Title || name, iupac: p.IUPACName, cid: p.CID, source: 'pubchem' };
}

function styleFor(color) {
  const c = color ? { color } : { colorscheme: 'Jmol' };
  if (state.style === 'stick') return { stick: { radius: 0.2, ...c } };
  if (state.style === 'fill') return { sphere: { scale: 1.0, ...c } };
  return { stick: { radius: 0.14, ...c }, sphere: { scale: 0.25, ...c } };
}

function draw() {
  if (!cur) return;
  const g = cur.groups[state.hl];
  const set = g ? new Set(g.atoms) : null;
  const shown = a => state.h || a.elem !== 'H';
  viewer.setStyle({}, {});
  viewer.setStyle({ predicate: a => shown(a) && (!set || set.has(a.index)) }, styleFor(null));
  if (set) viewer.setStyle({ predicate: a => shown(a) && !set.has(a.index) }, styleFor(tok('--dim')));
  $('note').textContent = g ? g.info
    : cur.groups.length ? 'Tap a group to highlight it in the model.'
    : 'None of the 21 functional groups on our list are in this molecule.';
  document.querySelectorAll('#chips button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.hl)));
  viewer.render();
}

function buildChips() {
  const box = $('chips');
  box.textContent = '';
  [{ name: 'All' }, ...cur.groups].forEach((g, k) => {
    const b = document.createElement('button');
    b.textContent = g.name + (g.count > 1 ? ` ×${g.count}` : '');
    b.dataset.i = String(k - 1);
    b.addEventListener('click', () => { state.hl = k - 1; draw(); });
    box.appendChild(b);
  });
}

function formulaHTML(f, charge) {
  const sign = charge ? `<sup>${Math.abs(charge) > 1 ? Math.abs(charge) : ''}${charge > 0 ? '+' : '−'}</sup>` : '';
  return esc(f).replace(/(\d+)/g, '<sub>$1</sub>') + sign;
}

async function run(name) {
  name = name.trim();
  if (!name) return;
  q.value = name;
  const my = ++seq;
  $('go').disabled = true;
  $('panel').hidden = $('foot').hidden = true; // don't leave the previous molecule's info next to a new error
  $('name').textContent = name; $('formula').textContent = '';
  showEmpty('Looking up the name…');
  try {
    const hit = await resolveName(name);
    if (my !== seq) return; // a newer search started
    if (!hit) return showEmpty(`No structure found for “${esc(name)}”. Check the spelling, or try the IUPAC name.`);
    showEmpty('Building the 3D model…');
    const res = await build(hit.smiles);
    if (my !== seq) return;
    if (res.error) return showEmpty(`<span class="err">${esc(res.error)}</span>`);

    const info = Object.fromEntries(GROUPS.map(g => [g.id, g]));
    cur = { ...hit, ...res, groups: res.groups.map(g => ({ ...info[g.id], ...g })) };
    state.hl = -1;
    viewer.removeAllModels();
    viewer.addModel(res.molfile, 'sdf', { keepH: true });
    empty.hidden = true; $('panel').hidden = false; $('foot').hidden = false;
    document.title = `${cur.title} · Mol3D`;
    $('name').textContent = cur.title;
    $('formula').innerHTML = `<b>${formulaHTML(res.formula, res.charge)}</b>`
      + (cur.iupac && cur.iupac.toLowerCase() !== cur.title.toLowerCase() ? ' · ' + esc(cur.iupac) : '');
    $('f-source').innerHTML = cur.source === 'opsin'
      ? 'Structure from the name by <a href="https://www.ebi.ac.uk/opsin/" target="_blank" rel="noopener">OPSIN</a> (EMBL-EBI).'
      : `Structure from <a href="https://pubchem.ncbi.nlm.nih.gov/compound/${cur.cid}" target="_blank" rel="noopener">PubChem CID ${cur.cid}</a>.`;
    $('f-smiles').textContent = 'SMILES: ' + cur.smiles;
    $('f-note').textContent = '3D model built in your browser (OpenChemLib, MMFF94s+). It shows one low-energy shape; real molecules move between many.'
      + (res.stripped ? ' This is a salt or mixture, so only the largest part is shown.' : '')
      + (res.minimized ? '' : ' The force field does not cover every atom here, so bond lengths are approximate.');
    buildChips();
    viewer.resize(); viewer.zoomTo(); draw(); home = viewer.getView(); // resize: the stage may have changed size since the viewer was made
  } catch (e) {
    showEmpty(`<span class="err">${esc(e.message || 'The lookup failed. Check the connection and try again.')}</span>`);
  } finally {
    if (my === seq) $('go').disabled = false;
  }
}

// The URL hash holds the query, so a link like …/mol3d/#caffeine opens that molecule.
const fromHash = () => decodeURIComponent(location.hash.slice(1));
window.addEventListener('hashchange', () => run(fromHash()));
$('search').addEventListener('submit', e => {
  e.preventDefault();
  const name = q.value.trim();
  if (fromHash() === name) run(name); else location.hash = name;
});
document.querySelectorAll('[data-style]').forEach(b => b.addEventListener('click', () => {
  state.style = b.dataset.style;
  document.querySelectorAll('[data-style]').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
  draw();
}));
$('btn-h').addEventListener('click', e => { state.h = !state.h; e.currentTarget.setAttribute('aria-pressed', String(state.h)); draw(); });
$('btn-reset').addEventListener('click', () => { if (home) { viewer.setView(home); viewer.render(); } });
const retheme = () => { viewer.setBackgroundColor(tok('--stage')); draw(); };
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', retheme);
window.addEventListener('resize', () => { viewer.resize(); viewer.render(); });
if (location.hash) run(fromHash());
