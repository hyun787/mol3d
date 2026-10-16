import { GROUPS, KIT_PARTS, IDEAL } from './chem.js';
import { resolveName } from './lookup.js';

const $ = id => document.getElementById(id);
const empty = $('empty'), q = $('q');
// guess: the guess-first mode; tries[i] = { pick: 'R' | 'S', hint: true } for stereo row i in that mode
const state = { style: 'ball', h: true, tab: 'structure', hl: -1, atom: -1, st: -1, guess: false, tries: {}, mirror: false };
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

function styleFor(color) {
  const c = color ? { color } : { colorscheme: 'Jmol' };
  if (state.style === 'stick') return { stick: { radius: 0.2, ...c } };
  if (state.style === 'fill') return { sphere: { scale: 1.0, ...c } };
  return { stick: { radius: 0.14, ...c }, sphere: { scale: 0.25, ...c } };
}

// Atoms to keep in color for the open tab; everything else is dimmed.
function highlight() {
  if (state.tab === 'structure') return state.atom >= 0 ? [state.atom, ...cur.atoms[state.atom].nbrs] : cur.groups[state.hl]?.atoms;
  if (state.tab === 'stereo') return cur.stereo[state.st]?.atoms;
  return null;
}

function draw() {
  if (!cur) return;
  const g = cur.groups[state.hl], hl = highlight();
  const set = hl ? new Set(hl) : null;
  const shown = a => state.h || a.elem !== 'H';
  viewer.setStyle({}, {});
  viewer.setStyle({ predicate: a => shown(a) && (!set || set.has(a.index)) }, styleFor(null));
  if (set) viewer.setStyle({ predicate: a => shown(a) && !set.has(a.index) }, styleFor(tok('--dim')));
  if (state.tab === 'structure' && state.atom >= 0) // the tapped atom gets a bigger ball than its neighbors
    viewer.addStyle({ predicate: a => shown(a) && a.index === state.atom }, { sphere: { scale: state.style === 'fill' ? 1.15 : 0.42, colorscheme: 'Jmol' } });
  $('note').textContent = g ? g.info
    : cur.groups.length ? 'Tap a group to highlight it in the model.'
    : 'None of the 21 functional groups on our list are in this molecule.';
  document.querySelectorAll('#chips button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.hl)));
  document.querySelectorAll('#pane-stereo .rows button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.st)));
  showAtom();
  showCenter();
  addLabels();
  viewer.render();
}

const HYB = { sp3: 'sp³', sp2: 'sp²', sp: 'sp' };
const SUB = n => String(n).replace(/\d/g, d => '₀₁₂₃₄₅₆₇₈₉'[d]);
const chargeSign = q => (q ? (Math.abs(q) > 1 ? Math.abs(q) : '') + (q > 0 ? '⁺' : '⁻') : '');

function showAtom() {
  const box = $('atom');
  if (state.atom < 0) { box.innerHTML = '<p class="small">Tap an atom in the model to see its hybridization, lone pairs and bond angles.</p>'; return; }
  const a = cur.atoms[state.atom], el = i => cur.atoms[i].el;
  const lp = a.lp === null ? '—' : `${a.lp} lone pair${a.lp === 1 ? '' : 's'}`;
  let html = `<p class="facts"><b>${esc(a.el + chargeSign(a.charge))}</b> · ${a.hyb ? HYB[a.hyb] : a.el === 'H' ? 'no hybridization' : 'hybridization not shown'} · ${lp}</p>`;
  if (a.el === 'H') html += '<p class="small">H bonds with its 1s orbital.</p>';
  else if (!a.hyb) html += '<p class="small">Hybridization is shown for C, N and O only.</p>';
  if (a.resonance) html += `<p class="small">${a.lp === 1 ? 'Its lone pair' : 'One lone pair'} sits in a p orbital and is shared with the π bond next to it (resonance), so the atom counts as sp².</p>`;
  if (!a.angles.length) html += `<p class="small">${a.nbrs.length ? 'Only one bond, so' : 'No bonds, so'} there is no bond angle at this atom.</p>`;
  else {
    const rows = a.angles.map(([p, q, deg]) => [[el(p), el(q)].sort().join(`–${a.el}–`), deg])
      .sort((x, y) => x[0].localeCompare(y[0]) || x[1] - y[1]);
    if (a.hyb) html += `<p class="facts">Ideal angle (${HYB[a.hyb]}): <b>${IDEAL[a.hyb]}°</b></p>`;
    html += `<ul class="angles" aria-label="Bond angles in this model">${rows.map(([n, d]) => `<li><span>${n}</span><b>${d.toFixed(1)}°</b></li>`).join('')}</ul>`
      + '<p class="small">Angles measured on this 3D model. Lone pairs, rings and crowding push real angles away from the ideal.</p>';
  }
  box.innerHTML = html;
}

function buildKit() {
  const { n, check } = cur.kit;
  $('parts').innerHTML = KIT_PARTS.filter(p => n[p.id]).map(p => `<li><span>${esc(p.label)}`
    + (p.id === 'red' && n.redOne !== n.red ? ` (or ${n.redOne} with one piece per 2-bond O)` : '') + `</span><b>${n[p.id]}</b></li>`).join('')
    || '<li><span>No parts from the guide</span></li>';
  $('kit-check').innerHTML = (check.length ? '<p class="small">Not counted above:</p><ul class="check">'
    + check.map(c => `<li>${esc(c.what)}${c.count > 1 ? ` (×${c.count})` : ''}</li>`).join('') + '</ul>' : '')
    + (n.bend ? '<p class="small">The guide does not say how many bendable pieces each ring atom takes.</p>' : '');
}

// A ring substituent named by its first atom: CH₃, OH, Cl; "–" means the group goes on.
function groupName(g) {
  const a = cur.atoms[g], h = a.nbrs.filter(x => cur.atoms[x].el === 'H').length;
  return a.el + (h ? 'H' + (h > 1 ? SUB(h) : '') : '') + (a.nbrs.length - h > 1 ? '–' : '');
}

// A stereocenter named by its groups in atom order (not priority order, so the name gives nothing away).
function centerName(r) {
  const nb = cur.atoms[r.atom].nbrs, isH = x => cur.atoms[x].el === 'H';
  const names = [...nb.filter(x => !isH(x)), ...nb.filter(isH)].map(groupName);
  return `${cur.atoms[r.atom].el} with ${names.slice(0, -1).join(', ')} and ${names.at(-1)}`;
}

function buildStereo() {
  const centers = $('centers'), others = $('stereo');
  centers.textContent = others.textContent = '';
  cur.stereo.forEach((r, i) => {
    const b = document.createElement('button');
    const title = r.kind === 'center' ? `${centerName(r)}: ${revealed(i) ? r.rs : '?'}`
      : r.kind === 'bond' ? `${r.label} double bond: ${r.ez || 'E/Z not assigned'}${r.ct ? ` (${r.ct})` : ''}`
      : `${groupName(r.groups[0])} and ${groupName(r.groups[1])} on a ${r.size}-membered ring: ${r.ct}`;
    b.innerHTML = `<span>${esc(title)}</span>` + (r.given ? '' : '<span class="sub">Not set in the source structure, so the 3D build picked this one.</span>');
    b.dataset.i = String(i);
    b.addEventListener('click', () => { state.st = state.st === i ? -1 : i; draw(); });
    (r.kind === 'center' ? centers : others).appendChild(b);
  });
  if (!centers.children.length) centers.innerHTML = '<p class="small">No stereocenter (R/S) in this molecule.</p>';
  $('center-tools').hidden = !centers.querySelector('button');
  if (!others.children.length) others.innerHTML = '<p class="small">No double bond with E/Z isomers and no cis/trans ring pair in this molecule.</p>';
}

const selCenter = () => { const r = state.tab === 'stereo' && cur.stereo[state.st]; return r?.kind === 'center' ? r : null; };
const revealed = i => !state.guess || !!state.tries[i]?.pick;
const atomAt = (i, m = 0) => viewer.getModel(m).selectedAtoms({ index: i })[0];

// On the model: R/S beside the selected stereocenter, and 1–4 on its groups.
function addLabels() {
  viewer.removeAllLabels();
  const r = selCenter();
  if (!r) return;
  const show = revealed(state.st);
  if (!show && !state.tries[state.st]?.hint) return;
  for (const m of state.mirror ? [0, 1] : [0]) { // model 1 = the mirror image: same priorities, R and S swapped
    // R/S sits up and to the right of the center, so it stays readable when group 4 is right behind the center
    if (show) viewer.addLabel(m ? { R: 'S', S: 'R' }[r.rs] : r.rs, { position: atomAt(r.atom, m), alignment: 'bottomLeft',
      screenOffset: { x: 10, y: 10 }, fontSize: 18, bold: true, fontColor: tok('--bg'), backgroundColor: tok('--accent'), backgroundOpacity: 1, inFront: true });
    r.order.forEach((a, k) => viewer.addLabel(String(k + 1), { position: atomAt(a, m), alignment: 'center', fontSize: 15,
      bold: true, fontColor: tok('--bg'), backgroundColor: tok('--fg'), backgroundOpacity: 1, inFront: true }));
  }
}

// Below the stereocenter rows: the priorities of the selected one and the turn button.
// In guess mode an open center shows only R / S / Hint; picking either one opens the answer.
function showCenter() {
  const box = $('center-detail'), r = selCenter(), i = state.st, t = state.tries[i] || {};
  if (!r) { box.textContent = ''; return; }
  const prio = `<ol class="prio" aria-label="CIP priorities">${r.order.map((a, k) => `<li>${esc(cur.atoms[a].el + (r.sets[k].length ? ` (${r.sets[k].join(', ')})` : ''))}</li>`).join('')}</ol>`
    + '<p class="small">In parentheses: the atoms on that first atom. A double bond counts twice.</p>';
  if (!revealed(i)) {
    box.innerHTML = (t.hint ? prio + '<p class="small">Group 4 now points away from you. Follow 1 → 2 → 3.</p>' : '')
      + '<div class="row"><button data-pick="R">R</button><button data-pick="S">S</button><button id="btn-hint">Hint</button></div>';
    box.querySelectorAll('[data-pick]').forEach(b => b.addEventListener('click', () => {
      state.tries[i] = { ...t, pick: b.dataset.pick };
      buildStereo(); draw();
    }));
    $('btn-hint').addEventListener('click', () => { state.tries[i] = { ...t, hint: true }; draw(); lowestToBack(r); });
    return;
  }
  box.innerHTML = (state.guess ? `<p class="facts">${t.pick === r.rs ? 'Right' : 'Not quite'}: it is <b>${r.rs}</b>.</p>` : '')
    + prio + '<div class="row"><button id="btn-back">Turn so 4 points away</button></div>';
  $('btn-back').addEventListener('click', () => lowestToBack(r));
}

// ---- View rotation. 3Dmol: screen = rotation quaternion q applied to model coordinates; the camera sits on +z,
// so "away from the viewer" is −z. Quaternions are [x, y, z, w].
const qMul = ([ax, ay, az, aw], [bx, by, bz, bw]) => [aw * bx + ax * bw + ay * bz - az * by, aw * by - ax * bz + ay * bw + az * bx,
  aw * bz + ax * by - ay * bx + az * bw, aw * bw - ax * bx - ay * by - az * bz];
const qRot = (q, v) => qMul(qMul(q, [...v, 0]), [-q[0], -q[1], -q[2], q[3]]).slice(0, 3);
const unit = v => { const n = Math.hypot(...v); return v.map(x => x / n); };
// Shortest turn taking unit vector u onto unit vector t (when u is opposite to t, turn 180° about the screen's y axis).
function qArc(u, t) {
  const w = 1 + u[0] * t[0] + u[1] * t[1] + u[2] * t[2];
  return w < 1e-6 ? [0, 1, 0, 0] : unit([u[1] * t[2] - u[2] * t[1], u[2] * t[0] - u[0] * t[2], u[0] * t[1] - u[1] * t[0], w]);
}
let anim = 0;
function animateView(to, ms = 600, done = null) {
  const from = viewer.getView(), my = ++anim, t0 = performance.now();
  if (from[4] * to[4] + from[5] * to[5] + from[6] * to[6] + from[7] * to[7] < 0) for (let k = 4; k < 8; k++) to[k] = -to[k]; // same turn, shorter way
  const step = now => {
    if (my !== anim) return; // a newer turn took over
    const f = Math.min(1, (now - t0) / ms), e = f * f * (3 - 2 * f);
    const v = from.map((x, k) => x + (to[k] - x) * e), q = unit(v.slice(4, 8));
    viewer.setView([...v.slice(0, 4), ...q]);
    if (f < 1) requestAnimationFrame(step); else done?.();
  };
  requestAnimationFrame(step);
}
// Center the view on the stereocenter and turn it so group 4 points straight away from the viewer.
// With the mirror image on, the mirror is rebuilt for the new view afterwards, so both group 4s point away.
function lowestToBack(r) {
  const mirror = state.mirror;
  if (mirror) setMirror(false, false);
  const c = atomAt(r.atom), l = atomAt(r.order[3]), v = viewer.getView(), q0 = v.slice(4, 8);
  const d = qRot(q0, unit([l.x - c.x, l.y - c.y, l.z - c.z])); // where group 4 points on screen now
  animateView([-c.x, -c.y, -c.z, v[3], ...qMul(qArc(d, [0, 0, -1]), q0)], 600, mirror ? () => setMirror(true) : null);
}

// ---- Mirror image (Stereo tab): a reflected copy (model 1) to the right of the model. It is built from the current
// view, so the mirror plane stands upright on screen between the two; both then turn together.
function setMirror(on, zoom = true) {
  state.mirror = on;
  $('btn-mirror').setAttribute('aria-pressed', String(on));
  $('mirror-note').hidden = !on;
  const old = viewer.getModel(1); // 3Dmol returns the last model for a missing id, so check it is really model 1
  if (old?.getID() === 1) viewer.removeModel(old);
  viewer.removeAllShapes();
  if (on) {
    const q = viewer.getView().slice(4, 8), qi = [-q[0], -q[1], -q[2], q[3]];
    const [n, up, out] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(v => qRot(qi, v)); // screen axes in model coordinates
    const pts = viewer.getModel(0).selectedAtoms({}).map(a => [a.x, a.y, a.z]);
    const c = [0, 1, 2].map(k => pts.reduce((s, p) => s + p[k], 0) / pts.length);
    const along = (p, u) => (p[0] - c[0]) * u[0] + (p[1] - c[1]) * u[1] + (p[2] - c[2]) * u[2];
    const size = u => 2 * Math.max(...pts.map(p => Math.abs(along(p, u)))) + 3;
    const half = size(n) / 2; // the plane sits 1.5 Å past the model's right edge
    const reflect = p => { const d = 2 * (along(p, n) - half); return p.map((x, k) => x - d * n[k]); };
    const lines = cur.molfile.split('\n'), count = parseInt(lines[3].slice(0, 3));
    for (let i = 4; i < 4 + count; i++) // V2000 atom line: x, y, z in the first 30 characters
      lines[i] = reflect(lines[i].slice(0, 30).trim().split(/\s+/).map(Number)).map(x => x.toFixed(4).padStart(10)).join('') + lines[i].slice(30);
    viewer.addModel(lines.join('\n'), 'sdf', { keepH: true });
    viewer.addBox({ center: { x: c[0] + half * n[0], y: c[1] + half * n[1], z: c[2] + half * n[2] }, color: tok('--muted'), opacity: 0.18,
      dimensions: { w: { x: n[0] * 0.05, y: n[1] * 0.05, z: n[2] * 0.05 }, h: Object.fromEntries('xyz'.split('').map((k, i) => [k, up[i] * size(up)])),
        d: Object.fromEntries('xyz'.split('').map((k, i) => [k, out[i] * size(out)])) } });
    $('mirror-note').textContent = cur.superimposable
      ? 'Superimposable: turn the mirror image and it matches the original, so it is the same compound. A molecule with stereocenters that matches its mirror image is a meso compound. Build both with your kit and check.'
      : 'Not superimposable: no turn makes the mirror image match the original, so the two are enantiomers, and every R becomes S. Build both with your kit and try to stack them.';
  }
  if (zoom) viewer.zoomTo({}, 400);
  draw();
}

function showTab() {
  document.querySelectorAll('[data-tab]').forEach(t => {
    const on = t.dataset.tab === state.tab;
    t.setAttribute('aria-selected', String(on));
    t.tabIndex = on ? 0 : -1;
    $('pane-' + t.dataset.tab).hidden = !on;
  });
}

function buildChips() {
  const box = $('chips');
  box.textContent = '';
  [{ name: 'All' }, ...cur.groups].forEach((g, k) => {
    const b = document.createElement('button');
    b.textContent = g.name + (g.count > 1 ? ` ×${g.count}` : '');
    b.dataset.i = String(k - 1);
    b.addEventListener('click', () => { state.hl = k - 1; state.atom = -1; draw(); });
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
    state.hl = state.atom = state.st = -1; state.tries = {};
    if (state.mirror) setMirror(false, false);
    viewer.removeAllModels();
    viewer.addModel(res.molfile, 'sdf', { keepH: true });
    viewer.setClickable({}, true, atom => { // tapping an atom only means something on the Structure tab
      if (state.tab !== 'structure') return;
      state.atom = state.atom === atom.index ? -1 : atom.index;
      state.hl = -1;
      draw();
    });
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
    buildChips(); buildKit(); buildStereo();
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
const tabs = [...document.querySelectorAll('[data-tab]')];
tabs.forEach((t, i) => {
  t.addEventListener('click', () => {
    state.tab = t.dataset.tab; showTab();
    if (state.mirror && state.tab !== 'stereo') setMirror(false); else draw(); // the mirror image is a Stereo-tab view
  });
  t.addEventListener('keydown', e => { // arrow keys move between tabs, as in any tab list
    const d = { ArrowRight: 1, ArrowLeft: -1 }[e.key];
    if (!d) return;
    const next = tabs[(i + d + tabs.length) % tabs.length];
    next.focus(); next.click();
  });
});
showTab();
$('btn-h').addEventListener('click', e => {
  state.h = !state.h;
  if (!state.h && cur?.atoms[state.atom]?.el === 'H') state.atom = -1;
  e.currentTarget.setAttribute('aria-pressed', String(state.h));
  draw();
});
$('btn-guess').addEventListener('click', e => {
  state.guess = !state.guess; state.tries = {};
  e.currentTarget.setAttribute('aria-pressed', String(state.guess));
  buildStereo(); draw();
});
$('btn-reset').addEventListener('click', () => {
  if (!home) return;
  anim++; viewer.setView(home);
  if (state.mirror) setMirror(true); else viewer.render();
});
$('btn-mirror').addEventListener('click', () => setMirror(!state.mirror));
const retheme = () => { viewer.setBackgroundColor(tok('--stage')); draw(); };
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', retheme);
window.addEventListener('resize', () => { viewer.resize(); viewer.render(); });
if (location.hash) run(fromHash());
