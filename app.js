import { GROUPS, KIT_PARTS, IDEAL, conformer, stericPairs, diaxial, moreStable, huckel } from './chem.js';
import { resolveName } from './lookup.js';

const $ = id => document.getElementById(id);
const empty = $('empty'), q = $('q');
// guess: the guess-first mode; tries[i] = { pick: 'R' | 'S', hint: true } for stereo row i in that mode
// chair: selected chair row; flipped: the other chair is shown; cguess/ctries: guess-first for which chair is more stable
const state = { style: 'ball', h: true, tab: 'structure', hl: -1, atom: -1, st: -1, guess: false, tries: {}, mirror: false, bond: -1, phi: 0, sguess: false, stries: {}, ring: -1,
  chair: -1, flipped: false, cguess: false, ctries: {},
  // Structure-tab views, one at a time: pol (bond dipoles), ar (aromatic ring index), nmr (equivalent H labels);
  // pguess/ppick and nguess/npick: guess-first for polar/nonpolar and for the number of signals
  pol: false, ar: -1, nmr: false, pguess: false, ppick: null, nguess: false, npick: null };
let cur = null, home = null, home3d = [], token = 0, seq = 0; // home3d: atom positions as the 3D build made them

const showEmpty = html => { empty.innerHTML = html; empty.hidden = false; };
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

if (typeof $3Dmol === 'undefined') {
  showEmpty('<span class="err">The 3D viewer did not load. Check the connection and reload.</span>');
  throw new Error('3Dmol.js missing');
}
const viewer = $3Dmol.createViewer($('stage'), { backgroundColor: tok('--stage'), antialias: true });

// Chemistry runs in a worker: {smiles} builds a molecule, {flip: i} finds the other chair of chair i of the last one.
// A newer request settles the older one as stale.
const worker = new Worker('worker.js', { type: 'module' });
let pending = null;
worker.onmessage = ({ data }) => { if (pending?.id === data.id) pending.resolve(data); };
worker.onerror = () => pending?.resolve({ error: 'The chemistry library did not load. Check the connection and reload.' });
const ask = msg => new Promise(resolve => {
  pending?.resolve({ stale: true });
  pending = { id: ++token, resolve };
  worker.postMessage({ id: token, ...msg });
});

function styleFor(color) {
  const c = color ? { color } : { colorscheme: 'Jmol' };
  if (state.style === 'stick') return { stick: { radius: 0.2, ...c } };
  if (state.style === 'fill') return { sphere: { scale: 1.0, ...c } };
  return { stick: { radius: 0.14, ...c }, sphere: { scale: 0.25, ...c } };
}

// Atoms to keep in color for the open tab; everything else is dimmed.
function highlight() {
  if (state.tab === 'structure') return state.atom >= 0 ? [state.atom, ...cur.atoms[state.atom].nbrs] : state.ar >= 0 ? cur.arom.rings[state.ar].atoms : cur.groups[state.hl]?.atoms;
  if (state.tab === 'stereo') return cur.stereo[state.st]?.atoms;
  if (state.tab === 'shape') {
    const b = cur.bonds[state.bond], c = chairNow();
    return b ? [...b.atoms, ...b.front, ...b.back] : c ? [...c.atoms, ...c.subs.map(s => s.sub)] : cur.rings[state.ring]?.atoms;
  }
  return null;
}

// The selected chair as shown now: the flipped chair's axial/equatorial when it is flipped.
const chairNow = () => {
  const c = state.tab === 'shape' && cur.chairs[state.chair];
  return !c ? null : state.flipped ? { ...c, subs: cur.flips[state.chair].subs } : c;
};

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
  document.querySelectorAll('#bonds button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.bond)));
  document.querySelectorAll('#rings button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.ring)));
  document.querySelectorAll('#chairs button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.chair)));
  document.querySelectorAll('#arom button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.ar)));
  $('btn-pol').setAttribute('aria-pressed', String(state.pol));
  $('btn-nmr').setAttribute('aria-pressed', String(state.nmr));
  showAtom();
  showPolarity();
  showArom();
  showNMR();
  showCenter();
  showNewman();
  showChair();
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

// ---- Structure tab: polarity, aromatic rings, equivalent H. One of the three views is on at a time.
function setView(k, v) {
  state.pol = state.nmr = false; state.ar = -1;
  state[k] = v;
  if (k === 'ar' && v >= 0) state.hl = state.atom = -1; // the ring takes the highlight
  draw();
}

const NO_POLARITY = {
  ion: 'No verdict: this is an ion. An ion has a net charge, so polar or nonpolar does not apply in the usual way.',
  element: 'No verdict: this molecule has I, S, P or another element for which the electronegativity difference misleads (C–I and C–S bonds are polar in practice). Arrows show only bonds with a difference of 0.4 or more.',
};
function showPolarity() {
  const p = cur.pol, box = $('pol');
  $('btn-pguess').hidden = p.polar === null || !p.bonds.length;
  if (p.polar === null) { box.innerHTML = `<p class="small">${esc(NO_POLARITY[p.why])}</p>`; return; }
  if (!p.bonds.length) { box.innerHTML = '<p class="facts"><b>Nonpolar</b>: no bond here has an electronegativity difference of 0.4 or more.</p>'; return; }
  if (state.pguess && !state.ppick) {
    box.innerHTML = '<p class="facts">Polar or nonpolar? Look at the bond dipoles, then pick one.</p>'
      + '<div class="row"><button data-ppick="polar">Polar</button><button data-ppick="nonpolar">Nonpolar</button></div>';
    box.querySelectorAll('[data-ppick]').forEach(b => b.addEventListener('click', () => { state.ppick = b.dataset.ppick; setView('pol', true); }));
    return;
  }
  const v = p.polar ? 'polar' : 'nonpolar';
  box.innerHTML = (state.pguess ? `<p class="facts">${state.ppick === v ? 'Right' : 'Not quite'}: it is ${v}.</p>` : '')
    + (p.polar ? '<p class="facts"><b>Polar</b>: the bond dipoles do not cancel. The gold arrow points toward the δ− end of the molecule.</p>'
      : '<p class="facts"><b>Nonpolar</b>: the polar bonds point in directions that cancel.</p>');
}

const PI_WHY = {
  double: p => (p.planar ? 'double bond; its lone pair lies in the ring plane' : 'double bond'),
  lp: () => 'lone pair in the p orbital',
  empty: () => 'empty p orbital',
};
const piLine = pi => { const n = huckel(pi); return `${pi} π electrons · ${n === null ? 'not 4n + 2' : `4n + 2 with n = ${n}`}`; };

function buildArom() {
  const box = $('arom'), { rings, systems } = cur.arom;
  box.textContent = '';
  rings.forEach((r, i) => {
    const row = document.createElement('button'), sys = systems[r.system], het = [...new Set(r.atoms.map(a => cur.atoms[a].el).filter(e => e !== 'C'))];
    const fused = sys.rings.length > 1;
    row.innerHTML = `<span>${esc(`${r.atoms.length}-membered ring${het.length ? ` (${het.join(', ')})` : ''}: ${fused ? `${r.pi} π electrons` : piLine(r.pi)}`)}</span>`
      + (fused ? `<span class="sub">${esc(`Fused ring system of ${sys.rings.length} rings: ${sys.pi === null ? 'π electrons not counted' : piLine(sys.pi)}`)}</span>` : '');
    row.dataset.i = String(i);
    row.addEventListener('click', () => { setView('ar', state.ar === i ? -1 : i); if (state.ar === i) alignChair(r); }); // ring seen edge-on: lobes above and below
    box.appendChild(row);
  });
  if (!rings.length) box.innerHTML = '<p class="small">No aromatic ring in this molecule.</p>';
}

function showArom() {
  const r = cur.arom.rings[state.ar], box = $('arom-detail');
  if (!r) { box.textContent = ''; return; }
  const name = a => { const h = cur.atoms[a].nbrs.filter(x => cur.atoms[x].el === 'H').length; return cur.atoms[a].el + (h ? `–H` : ''); };
  box.innerHTML = `<ul class="angles" aria-label="π electrons from each ring atom">${r.parts.map((p, k) => `<li><span>${k + 1} · ${esc(name(p.atom))}: ${esc(PI_WHY[p.why](p))}</span><b>${p.e}</b></li>`).join('')}</ul>`
    + `<p class="facts">Total: <b>${r.pi}</b> π electrons in this ring.</p>`;
}

// Letter and color for equivalence class c: the letter is the main mark, since 8 colors repeat after h.
const NMR_COLORS = ['#1d4ed8', '#b91c1c', '#15803d', '#7e22ce', '#c2410c', '#0f766e', '#a21caf', '#4d7c0f'];
const letter = c => (c < 26 ? String.fromCharCode(97 + c) : String.fromCharCode(97 + Math.floor(c / 26) - 1) + String.fromCharCode(97 + c % 26));

function showNMR() {
  const { h, cls } = cur.nmr, box = $('nmr'), n = new Set(cls).size;
  $('btn-nguess').hidden = !h.length;
  if (!h.length) { box.innerHTML = '<p class="small">No H atoms in this molecule.</p>'; return; }
  if (state.nguess && state.npick === null) {
    box.innerHTML = '<div class="guess"><label class="facts" for="npick">How many ¹H NMR signals?</label><input id="npick" type="number" min="1" max="99" inputmode="numeric"><button id="btn-ncheck">Check</button></div>';
    const check = () => { const v = parseInt($('npick').value, 10); if (v > 0) { state.npick = v; setView('nmr', true); } };
    $('btn-ncheck').addEventListener('click', check);
    $('npick').addEventListener('keydown', e => { if (e.key === 'Enter') check(); });
    return;
  }
  const rows = [...new Set(cls)].map(c => {
    const hs = h.filter((_, i) => cls[i] === c), on = [...new Set(hs.map(a => plain(cur.atoms[a].nbrs[0])))];
    return `<li><span><i class="swatch" style="background:${NMR_COLORS[c % 8]}"></i>${letter(c)}</span><b>${hs.length} H on ${esc(on.join(', '))}</b></li>`;
  });
  box.innerHTML = (state.nguess ? `<p class="facts">${state.npick === n ? 'Right' : 'Not quite'}: ${n === 1 ? 'there is 1 signal' : `there are ${n} signals`}.</p>` : '')
    + `<p class="facts"><b>${n} signal${n === 1 ? '' : 's'}</b> (sets of equivalent H)</p>`
    + `<ul class="angles" aria-label="Sets of equivalent H">${rows.join('')}</ul>`
    + (state.h ? '' : '<p class="small">H atoms are hidden: turn on Show H to see the letters.</p>');
}

// On the model (Structure tab only): dipole arrows, p-orbital lobes, or H letters. Only these shapes are removed here.
let structShapes = [];
function drawStructure() {
  structShapes.forEach(x => viewer.removeShape(x));
  structShapes = [];
  if (state.tab !== 'structure' || !cur) return;
  const P = a => { const t = atomAt(a); return [t.x, t.y, t.z]; }, xyz = v => ({ x: v[0], y: v[1], z: v[2] });
  if (state.pol && cur.pol.bonds.length) {
    const all = home3d, c = [0, 1, 2].map(k => all.reduce((s, v) => s + v[k], 0) / all.length);
    for (const [p, q] of cur.pol.bonds.filter(([x, y]) => state.h || (cur.atoms[x].el !== 'H' && cur.atoms[y].el !== 'H'))) {
      const a = P(p), b = P(q), d = unit(b.map((x, k) => x - a[k])), m = a.map((x, k) => (x + b[k]) / 2);
      // shift the arrow 0.35 Å off the bond, away from the molecule's center, so the stick does not hide it
      let off = m.map((x, k) => x - c[k]);
      off = off.map((x, k) => x - d[k] * (off[0] * d[0] + off[1] * d[1] + off[2] * d[2]));
      if (Math.hypot(...off) < 1e-3) off = Math.abs(d[2]) < 0.9 ? [d[1], -d[0], 0] : [0, d[2], -d[1]];
      off = unit(off).map(x => x * 0.35);
      const len = Math.hypot(...b.map((x, k) => x - a[k]));
      structShapes.push(viewer.addArrow({ start: xyz(a.map((x, k) => x + off[k] + d[k] * len * 0.15)), end: xyz(b.map((x, k) => x + off[k] - d[k] * len * 0.15)),
        radius: 0.06, radiusRatio: 2.2, mid: 0.7, color: tok('--dipole'), opacity: 0.75 }));
    }
    const net = cur.pol.net, open = !state.pguess || state.ppick;
    if (cur.pol.polar && open) {
      const d = unit(net);
      structShapes.push(viewer.addArrow({ start: xyz(c.map((x, k) => x - d[k] * 1.2)), end: xyz(c.map((x, k) => x + d[k] * 1.8)),
        radius: 0.12, radiusRatio: 2, mid: 0.75, color: tok('--accent'), opacity: 0.85 }));
    }
  }
  const r = cur.arom.rings[state.ar];
  if (r) {
    const { n } = ringFrame(r.atoms.map(P));
    for (const a of r.atoms) for (const s of [1, -1]) // one see-through lobe above and one below the ring plane
      structShapes.push(viewer.addSphere({ center: xyz(P(a).map((x, k) => x + s * 0.75 * n[k])), radius: 0.5, color: tok('--lobe'), opacity: 0.6 }));
  }
  if (state.nmr && state.h && (!state.nguess || state.npick !== null)) {
    const { h, cls } = cur.nmr;
    h.forEach((a, i) => viewer.addLabel(letter(cls[i]), { position: atomAt(a), alignment: 'center', fontSize: 13, bold: true,
      fontColor: '#ffffff', backgroundColor: NMR_COLORS[cls[i] % 8], backgroundOpacity: 0.95, inFront: true }));
  }
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

// Shape tab: one row per C–C bond, named by its two carbons, with the groups on each.
function buildShape() {
  const box = $('bonds'), carbon = a => groupName(a).replace(/–$/, '');
  box.textContent = '';
  cur.bonds.forEach((b, i) => {
    const row = document.createElement('button'), [f, k] = b.atoms;
    row.innerHTML = `<span>${esc(`${carbon(f)}–${carbon(k)} bond`)}</span><span class="sub">${esc(`Front: ${b.front.map(groupName).join(', ')} · Back: ${b.back.map(groupName).join(', ')}`)}</span>`;
    row.dataset.i = String(i);
    row.addEventListener('click', () => setBond(state.bond === i ? -1 : i));
    box.appendChild(row);
  });
  if (!cur.bonds.length) box.innerHTML = '<p class="small">No C–C single bond outside a ring with two sp³ carbons in this molecule.</p>';
  buildChairs();
  buildRings();
}

// Chair rows: one per chair, with the ring's N or O named and the groups (not H) on it.
function buildChairs() {
  const box = $('chairs');
  box.textContent = '';
  cur.chairs.forEach((c, i) => {
    const row = document.createElement('button'), het = [...new Set(c.atoms.map(a => cur.atoms[a].el).filter(e => e !== 'C'))];
    const groups = c.subs.filter(s => !s.h).map(s => plain(s.sub));
    row.innerHTML = `<span>${esc(`Chair${het.length ? ` (${het.join(' and ')} in the ring)` : ''}`)}</span>`
      + `<span class="sub">${esc(groups.length ? `Groups on the ring: ${groups.join(', ')}` : 'Only H on the ring')}</span>`;
    row.dataset.i = String(i);
    row.addEventListener('click', () => setChair(state.chair === i ? -1 : i));
    box.appendChild(row);
  });
  if (!cur.chairs.length) box.innerHTML = '<p class="small">No six-membered ring of sp³ atoms in a chair shape in this molecule.</p>';
}

function setChair(i) {
  state.chair = i; state.flipped = false;
  if (i >= 0) state.bond = state.ring = -1;
  turn(null); // back to the 3D build: a Newman turn or the other chair goes
  draw();
  if (i >= 0) alignChair(cur.chairs[i]);
}

// Look at the chair (or any ring) from the side and a little above, as textbooks draw it: the ring's axis (normal) points up the
// screen, tilted 20° toward the viewer, so axial bonds run up and down and equatorial bonds fan out.
function alignChair(c) {
  const { ctr, n } = ringFrame(c.atoms.map(i => home3d[i]));
  const v = viewer.getView(), q0 = v.slice(4, 8), t = 20 * Math.PI / 180;
  animateView([-ctr[0], -ctr[1], -ctr[2], v[3], ...qMul(qArc(qRot(q0, unit(n)), [0, Math.cos(t), Math.sin(t)]), q0)]);
}

// Center and unit normal of a ring from its atom positions, in ring order (Newell: sum of cross products).
function ringFrame(p) {
  const ctr = [0, 1, 2].map(k => p.reduce((s, v) => s + v[k], 0) / p.length);
  let n = [0, 0, 0];
  p.forEach((v, i) => { const a = v.map((x, k) => x - ctr[k]), b = p[(i + 1) % p.length].map((x, k) => x - ctr[k]);
    n = n.map((x, k) => x + [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]][k]); });
  return { ctr, n: unit(n) };
}

// Why a chair does not flip here (chem.js flipBlock).
const NO_FLIP = {
  h: 'Only H atoms sit on this ring, so the two chairs match and a flip changes nothing.',
  bridged: 'This ring is part of a bridged ring system. The bridge holds it in this shape, so it cannot flip.',
  trans: 'This ring is trans-fused to another ring, as in trans-decalin. The other ring locks it, so it cannot flip.',
  cis: 'This ring is cis-fused to another ring, as in cis-decalin. It can flip, but only with the other ring flipping too, and this page does not show that.',
};
let flipping = -1; // chair index whose other chair the worker is looking for

// Below the chair rows: each group's axial/equatorial in the chair shown, and the ring flip.
function showChair() {
  const c = chairNow();
  $('chair').hidden = !c;
  if (!c) return;
  const g = c.subs.filter(s => !s.h), f = cur.flips[state.chair];
  $('axeq').innerHTML = g.length ? `<ul class="angles" aria-label="Axial and equatorial groups">${g.map(s => `<li><span>${esc(plain(s.sub))}</span><b>${s.ax ? 'axial' : 'equatorial'}</b></li>`).join('')}</ul>`
    : '<p class="small">Only H on this ring: each ring carbon holds one axial and one equatorial H.</p>';
  const isH = a => cur.atoms[a].el === 'H'; // pairs of two groups first: they crowd the most
  const pairs = diaxial(c).sort((p, q) => isH(p[1]) - isH(q[1])).map(([x, y]) => `${plain(x)} and ${plain(y)}`), seen = [...new Set(pairs)];
  $('diax').innerHTML = pairs.length ? 'Steric strain: 1,3-diaxial interactions between '
    + seen.map(t => { const n = pairs.filter(x => x === t).length; return esc(t) + (n > 1 ? ` (×${n})` : ''); }).join('; ') : '';
  const btn = $('btn-flip');
  btn.hidden = !!c.block || f === null;
  $('btn-cguess').hidden = !!c.block || f === null || !!c.judge; // nothing to guess without a verdict
  btn.disabled = flipping === state.chair;
  btn.setAttribute('aria-pressed', String(state.flipped));
  $('flip-out').textContent = c.block ? NO_FLIP[c.block]
    : flipping === state.chair ? 'Looking for the other chair…'
    : f === null ? 'Flip not available: the conformer search (OpenChemLib) did not build the other chair. This happens when a group larger than CH₃ or OH would turn axial, as with tert-butyl.'
    : f ? (state.flipped ? 'Showing the flipped chair: each axial group is now equatorial, and each equatorial group axial.' : 'Showing the chair as built.') : '';
  $('flip-note').textContent = f ? 'The animation slides each atom straight to its new place. A real ring flip passes through half-chair and twist-boat forms.' : '';
  $('stable').innerHTML = f ? stableHTML(c, f) : '';
  $('stable').querySelectorAll('[data-cpick]').forEach(b => b.addEventListener('click', () => { state.ctries[state.chair] = b.dataset.cpick; draw(); }));
}

// Which chair is more stable, once both are known. MMFF94s+ is trusted only for groups checked against A-values.
const NO_VERDICT = {
  halogen: 'No verdict on which chair is more stable. The force field used here (MMFF94s+) does not put F, Cl or Br lower when equatorial, though measured A-values favor equatorial, so the page gives no verdict for halogens.',
  nitrile: 'No verdict on which chair is more stable. The force field used here (MMFF94s+) puts C≡N lower when axial, though its measured A-value favors equatorial.',
  unchecked: 'No verdict on which chair is more stable. The force field used here (MMFF94s+) was checked against measured A-values only for alkyl groups, OH and NH₂.',
};
function stableHTML(c, f) {
  if (c.judge || f.dE === null) return `<p class="small">${esc(NO_VERDICT[c.judge?.why] || 'No verdict: the energy of the chair as built is missing.')}</p>`;
  const v = moreStable(f.dE), shown = state.flipped ? 'flipped' : 'built', pick = state.cguess && state.ctries[state.chair];
  if (state.cguess && !pick) return '<p class="facts">Which chair is more stable? Pick one to see the answer.</p>'
    + '<div class="row"><button data-cpick="built">Chair as built</button><button data-cpick="flipped">Flipped chair</button><button data-cpick="same">About the same</button></div>';
  return (pick ? `<p class="facts">${pick === v ? 'Right' : 'Not quite'}: ${v === 'same' ? 'the two chairs are about the same' : `the ${v === 'built' ? 'chair as built' : 'flipped chair'} is more stable`}.</p>` : '')
    + `<p class="facts"><b>${v === 'same' ? 'Both chairs are about the same' : v === shown ? 'This chair is more stable' : 'The other chair is more stable'}</b>.</p>`
    + '<p class="small">From the MMFF94s+ force field: the lowest energy found for each chair, compared. Less than 1 kJ/mol apart counts as the same. Groups other than H sit better equatorial, where they avoid 1,3-diaxial crowding.</p>';
}

// Ring flip: the first press asks the worker for the other chair (a conformer search, a second or more), then each press
// slides between the two.
$('btn-flip').addEventListener('click', async () => {
  const i = state.chair, mol = cur;
  if (!(i in mol.flips)) {
    flipping = i; draw();
    const res = await ask({ flip: i });
    if (res.stale) return; // a new molecule took over
    flipping = -1;
    mol.flips[i] = res.flip || null;
    if (state.chair !== i || state.tab !== 'shape' || !mol.flips[i]) return draw();
  }
  const f = cur.flips[i];
  state.flipped = !state.flipped;
  morph(state.flipped ? home3d : f.xyz, state.flipped ? f.xyz : home3d);
});

// Slide every atom in a straight line from one set of positions to another (labels and rods wait for the end).
let morphing = 0;
function morph(from, to, ms = 900) {
  const atoms = viewer.getModel(0).selectedAtoms({}), my = ++morphing, t0 = performance.now();
  const step = now => {
    if (my !== morphing) return; // turn() or a newer slide took over
    const k = Math.min(1, (now - t0) / ms), e = k * k * (3 - 2 * k);
    atoms.forEach((a, i) => { [a.x, a.y, a.z] = from[i].map((x, j) => x + (to[i][j] - x) * e); });
    if (k < 1) { state.sliding = true; draw(); requestAnimationFrame(step); } else { state.sliding = false; draw(); }
  };
  requestAnimationFrame(step);
}

// Ring-size notes, shown only for rings of sp³ atoms (sp² or aromatic rings get their own line).
const RING_NOTE = {
  3: 'Three ring atoms must lie flat in a triangle, so the ring angles are 60°, far below 109.5°. The bonds bend, all its C–H bonds are eclipsed, and three-membered rings hold the most strain.',
  4: 'A four-membered ring puckers a little, but its angles stay near 88°, still far below 109.5°. Angle strain is large.',
  5: 'A five-membered ring folds into an envelope. Its angles sit near 105°, close to 109.5°, so angle strain is small. Most of its strain comes from eclipsed bonds.',
  6: 'A six-membered ring folds into a chair. Its angles sit near 111° and every bond is staggered, so a ring of six sp³ atoms has almost no strain.',
  7: 'A seven-membered ring folds out of plane. Its angles sit a little above 109.5°, and crowding across the ring adds some strain.',
};

function buildRings() {
  const box = $('rings'), deg = x => `${x.toFixed(1)}°`;
  box.textContent = '';
  cur.rings.forEach((r, i) => {
    const row = document.createElement('button'), hybs = new Set(r.hyb);
    const d = r.ideal === null || hybs.size > 1 ? null : r.mean - r.ideal; // a mixed sp²/sp³ ring has no single ideal
    const note = r.aromatic ? `This ring is aromatic and flat, with sp² atoms (ideal 120°). A flat ring of ${r.size} atoms averages ${deg(180 * (r.size - 2) / r.size)} per angle.`
      : r.ideal === null ? 'Hybridization is shown for C, N and O only, so this ring has no ideal angle here.'
      : hybs.size === 1 && hybs.has('sp3') ? RING_NOTE[r.size]
      : 'This ring has sp² atoms (ideal 120°) as well as sp³ atoms (ideal 109.5°), so compare each angle with its own ideal.';
    row.innerHTML = `<span>${r.size}-membered ring${r.aromatic ? ' (aromatic)' : ''}: average angle ${deg(r.mean)}`
      + (d === null ? '' : ` · ideal ${deg(r.ideal)} · ${Math.abs(d) < 0.5 ? 'about the same' : `${deg(Math.abs(d))} ${d < 0 ? 'smaller' : 'larger'}`}`) + '</span>'
      + `<span class="sub mono">${esc(r.angles.map(deg).join(', '))}</span>` + (note ? `<span class="sub">${esc(note)}</span>` : '');
    row.dataset.i = String(i);
    row.addEventListener('click', () => { state.bond = state.chair = -1; state.flipped = false; turn(null); state.ring = state.ring === i ? -1 : i; draw(); });
    box.appendChild(row);
  });
  if (!cur.rings.length) box.innerHTML = '<p class="small">No ring in this molecule.</p>';
}

const plain = a => groupName(a).replace(/–$/, '');

function showNewman() {
  const b = cur.bonds[state.bond];
  $('newman').hidden = !b;
  if (!b) return;
  const c = conformer(state.phi, b.ref), [p, q] = b.x.map(plain);
  $('phi-out').textContent = `Turn the back carbon · ${p}–C–C–${q} dihedral: ${state.phi}°`;
  $('phi').value = String(state.phi);
  $('conf').innerHTML = !c.kind ? 'In between: neither staggered nor eclipsed.'
    : `<b>${c.kind === 'staggered' ? 'Staggered' : 'Eclipsed'}</b>` + (c.rel ? ` · ${esc(p)} and ${esc(q)} are <b>${c.rel}</b>` : '');
  const pairs = stericPairs(b, state.phi).map(([f, k, a]) => `${esc(plain(f))} and ${esc(plain(k))} are <b>${a <= 5 ? 'eclipsed' : Math.abs(a - 60) <= 5 ? 'gauche' : `${a}° apart`}</b>`);
  const steric = [...new Set(pairs)].map(t => { const n = pairs.filter(x => x === t).length; return n > 1 ? `${t} (×${n})` : t; });
  $('strain').innerHTML = steric.length ? 'Steric (van der Waals) strain: ' + steric.join('; ') : '';
  // Guess first: the curve stays hidden until the student turns to the form they think is least stable and taps.
  const guessing = state.sguess && b.levels, open = !guessing || state.bond in state.stries, top = b.levels && Math.max(...b.levels);
  $('btn-sguess').hidden = !b.levels;
  $('shape-guess').innerHTML = !guessing ? ''
    : !open ? '<p class="facts">Turn the back carbon to the least stable form, then tap the button.</p><div class="row"><button id="btn-pick">This one is least stable</button></div>'
    : state.stries[state.bond] ? `<p class="facts">Right: the least stable form is at <b>${anglesAt(b.levels, top)}</b>.</p>`
    : `<p class="facts">Not quite: the least stable form is at <b>${anglesAt(b.levels, top)}</b>, the highest point on the curve.</p>`;
  $('btn-pick')?.addEventListener('click', () => {
    const k = conformer(state.phi, b.ref);
    state.stries[state.bond] = !!k.kind && b.levels[k.at / 60] === top;
    draw();
  });
  $('curve').innerHTML = b.levels && open ? curveSVG(b, `${p}–C–C–${q}`) : '';
  $('curve-note').textContent = !open ? '' : b.levels
    ? 'Schematic: the curve shows which forms are higher in energy, not by how much. It follows the textbook rule: eclipsed forms sit above staggered ones, and each pair of carbon groups that is eclipsed or gauche raises the energy.'
      + (b.flat ? ' A C=C or benzene-ring group turns aside, so its gauche form costs almost nothing (under 1 kJ/mol by MMFF94s+), and only its eclipsed forms count.' : '')
      + ' Tap the curve to turn to that angle.'
    : 'No curve for this bond. Groups other than H, sp³ carbon, and C=C or benzene-ring carbon (Cl, O, N, F, or a C=O carbon) bring in polarity and other effects that can change which form is lowest: 1,2-difluoroethane, for example, prefers gauche. Staggered forms are still usually lower than eclipsed ones.';
}

// The angles among 0°, 60° … 300° whose form has order r, as an English list.
const anglesAt = (lv, r) => new Intl.ListFormat('en').format([0, 60, 120, 180, 240, 300].filter((_, i) => lv[i] === r).map(d => d + '°'));

// Schematic energy vs. dihedral (no numbers): the six forms at 0°, 60° … 300° placed by their order, staggered in the low
// band and eclipsed in the high band, joined by half-cosine curves. The dot is the slider angle.
const CW = 280, CH = 140, CL = 26, CR = 14, CT = 10, CB = 24;
function curveSVG(b, name) {
  const w = CW - CL - CR, h = CH - CT - CB, lv = b.levels;
  const band = idx => [...new Set(idx.map(i => lv[i]))].sort((x, y) => x - y);
  const st = band([1, 3, 5]), ec = band([0, 2, 4]);
  const y6 = lv.map(r => st.includes(r) ? 0.3 * st.indexOf(r) / Math.max(1, st.length - 1)
    : ec.length === 1 ? 1 : 0.65 + 0.35 * ec.indexOf(r) / (ec.length - 1));
  const E = phi => { const i = Math.floor(phi / 60) % 6, t = (phi % 60) / 60, a = y6[i], c = y6[(i + 1) % 6]; return a + (c - a) * (1 - Math.cos(Math.PI * t)) / 2; };
  const X = phi => CL + w * phi / 360, Y = e => CT + h * (1 - e);
  const line = Array.from({ length: 73 }, (_, k) => `${X(5 * k).toFixed(1)},${Y(E(5 * k)).toFixed(1)}`).join(' ');
  const ticks = [0, 60, 120, 180, 240, 300, 360].map(d => `<line x1="${X(d)}" x2="${X(d)}" y1="${CT + h}" y2="${CT + h + 4}" stroke="var(--muted)"/>`
    + `<text x="${X(d)}" y="${CH - 6}" text-anchor="middle">${d}°</text>`).join('');
  return `<svg viewBox="0 0 ${CW} ${CH}" role="img" aria-label="${esc(`Schematic energy curve for the ${name} dihedral. Least stable at ${anglesAt(lv, Math.max(...lv))}; most stable at ${anglesAt(lv, 0)}.`)}">
    <g font-family="var(--mono)" font-size="10" fill="var(--muted)">${ticks}
      <text transform="translate(10 ${CT + h / 2}) rotate(-90)" text-anchor="middle">Energy →</text></g>
    <line x1="${CL}" x2="${CL + w}" y1="${CT + h}" y2="${CT + h}" stroke="var(--line)"/>
    <line x1="${CL}" x2="${CL}" y1="${CT}" y2="${CT + h}" stroke="var(--line)"/>
    <polyline points="${line}" fill="none" stroke="var(--fg)" stroke-width="2" stroke-linejoin="round"/>
    <circle cx="${X(state.phi)}" cy="${Y(E(state.phi))}" r="5" fill="var(--accent)" stroke="var(--bg)" stroke-width="2"/></svg>`;
}

// Turn the back side of bond b so its x-pair dihedral is state.phi. b = null puts every atom back where the 3D build put it,
// so the other tabs (and their bond angles) always see the original model.
function turn(b) {
  morphing++; state.sliding = false; // stop a ring-flip slide
  const atoms = viewer.getModel(0).selectedAtoms({});
  atoms.forEach((a, i) => { [a.x, a.y, a.z] = home3d[i]; });
  if (!b) return;
  const [F, B] = b.atoms.map(i => home3d[i]), k = unit(B.map((v, j) => v - F[j])), h = (state.phi - b.phi0) * Math.PI / 360;
  const q = [...k.map(v => v * Math.sin(h)), Math.cos(h)];
  for (const i of b.side) {
    const v = qRot(q, home3d[i].map((x, j) => x - F[j]));
    [atoms[i].x, atoms[i].y, atoms[i].z] = v.map((x, j) => x + F[j]);
  }
}

// Look straight down bond b: front carbon toward the viewer, its x group pointing up, as in a textbook Newman projection.
function alignBond(b) {
  const [F, B] = b.atoms.map(i => home3d[i]), x = home3d[b.x[0]], v = viewer.getView(), q0 = v.slice(4, 8);
  const q1 = qMul(qArc(qRot(q0, unit(B.map((c, j) => c - F[j]))), [0, 0, -1]), q0);
  const p = qRot(q1, x.map((c, j) => c - F[j])), h = (Math.PI / 2 - Math.atan2(p[1], p[0])) / 2;
  const c = F.map((s, j) => (s + B[j]) / 2);
  animateView([-c[0], -c[1], -c[2], v[3], ...qMul([0, 0, Math.sin(h), Math.cos(h)], q1)]);
}

function setBond(i, zoom = false) {
  state.bond = i; state.flipped = false;
  if (i >= 0) state.ring = state.chair = -1;
  const b = cur.bonds[i];
  if (b) state.phi = Math.round(b.phi0 / 5) * 5 % 360;
  turn(b || null);
  draw();
  if (b) alignBond(b); else if (zoom) viewer.zoomTo({}, 400);
}

const selCenter = () => { const r = state.tab === 'stereo' && cur.stereo[state.st]; return r?.kind === 'center' ? r : null; };
const revealed = i => !state.guess || !!state.tries[i]?.pick;
const atomAt = (i, m = 0) => viewer.getModel(m).selectedAtoms({ index: i })[0];

// On the model: R/S beside the selected stereocenter, and 1–4 on its groups.
function addLabels() {
  viewer.removeAllLabels();
  drawStrain();
  drawStructure();
  const c = !state.sliding && chairNow(); // Shape tab: ax/eq beside each group (not H) on the selected chair
  if (c) for (const s of c.subs.filter(x => !x.h)) viewer.addLabel(s.ax ? 'ax' : 'eq', { position: atomAt(s.sub), alignment: 'bottomLeft',
    screenOffset: { x: 6, y: 6 }, fontSize: 12, bold: true, fontColor: tok('--bg'), backgroundColor: tok(s.ax ? '--accent' : '--fg'), backgroundOpacity: 1, inFront: true });
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

// Shape tab: a see-through pink rod with a label between each pair of groups in steric strain: "vdW strain" on a Newman
// bond, "1,3-diaxial" on a chair. Only these rods are removed here, so the mirror plane (also a shape) is never touched.
let strainShapes = [];
function drawStrain() {
  strainShapes.forEach(s => viewer.removeShape(s));
  strainShapes = [];
  if (state.tab !== 'shape' || state.sliding) return;
  const b = cur.bonds[state.bond], c = chairNow();
  const pairs = b ? stericPairs(b, state.phi).map(([f, k]) => [f, k, 'vdW strain']) : c ? diaxial(c).map(([x, y]) => [x, y, '1,3-diaxial']) : [];
  for (const [f, k, text] of pairs) {
    const [p, q] = [atomAt(f), atomAt(k)].map(a => ({ x: a.x, y: a.y, z: a.z }));
    strainShapes.push(viewer.addCylinder({ start: p, end: q, radius: 0.15, color: tok('--cut'), opacity: 0.6 }));
    // the label sits just above the rod's middle so it does not hide the rod
    viewer.addLabel(text, { position: { x: (p.x + q.x) / 2, y: (p.y + q.y) / 2, z: (p.z + q.z) / 2 }, alignment: 'bottomCenter', screenOffset: { x: 0, y: 10 },
      fontSize: 12, bold: true, fontColor: tok('--bg'), backgroundColor: tok('--cut'), backgroundOpacity: 0.85, inFront: true });
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
    b.addEventListener('click', () => { state.hl = k - 1; state.atom = -1; state.ar = -1; draw(); });
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
    const res = await ask({ smiles: hit.smiles });
    if (my !== seq) return;
    if (res.error) return showEmpty(`<span class="err">${esc(res.error)}</span>`);

    const info = Object.fromEntries(GROUPS.map(g => [g.id, g]));
    cur = { ...hit, ...res, groups: res.groups.map(g => ({ ...info[g.id], ...g })), flips: {} }; // flips[i]: the other chair of chair i, null if none
    state.hl = state.atom = state.st = state.bond = state.ring = state.chair = -1; state.tries = {}; state.stries = {}; state.ctries = {}; state.flipped = state.sliding = false;
    state.ar = -1; state.ppick = state.npick = null;
    morphing++; flipping = -1; // a ring-flip slide or search for the old molecule is over
    if (state.mirror) setMirror(false, false);
    viewer.removeAllModels();
    viewer.addModel(res.molfile, 'sdf', { keepH: true });
    home3d = viewer.getModel(0).selectedAtoms({}).map(a => [a.x, a.y, a.z]);
    viewer.setClickable({}, true, atom => { // tapping an atom only means something on the Structure tab
      if (state.tab !== 'structure') return;
      state.atom = state.atom === atom.index ? -1 : atom.index;
      state.hl = -1; state.ar = -1;
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
    buildChips(); buildKit(); buildStereo(); buildShape(); buildArom();
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
    const left = state.tab; state.tab = t.dataset.tab; showTab();
    if (left === 'shape' && state.tab !== 'shape' && (state.bond >= 0 || state.flipped)) setBond(-1, true); // a turned model or the other chair is a Shape-tab view
    else if (state.mirror && state.tab !== 'stereo') setMirror(false); else draw(); // the mirror image is a Stereo-tab view
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
$('btn-sguess').addEventListener('click', e => {
  state.sguess = !state.sguess; state.stries = {};
  e.currentTarget.setAttribute('aria-pressed', String(state.sguess));
  draw();
});
$('btn-cguess').addEventListener('click', e => {
  state.cguess = !state.cguess; state.ctries = {};
  e.currentTarget.setAttribute('aria-pressed', String(state.cguess));
  draw();
});
$('btn-reset').addEventListener('click', () => {
  if (!home) return;
  anim++; viewer.setView(home);
  if (state.mirror) setMirror(true); else viewer.render();
});
$('btn-mirror').addEventListener('click', () => setMirror(!state.mirror));
$('btn-pol').addEventListener('click', () => setView('pol', !state.pol));
$('btn-nmr').addEventListener('click', () => setView('nmr', !state.nmr));
$('btn-pguess').addEventListener('click', e => {
  state.pguess = !state.pguess; state.ppick = null;
  e.currentTarget.setAttribute('aria-pressed', String(state.pguess));
  draw();
});
$('btn-nguess').addEventListener('click', e => {
  state.nguess = !state.nguess; state.npick = null;
  e.currentTarget.setAttribute('aria-pressed', String(state.nguess));
  if (state.nguess) state.nmr = false; // the letters would give the count away
  draw();
});
const setPhi = phi => { state.phi = phi; turn(cur.bonds[state.bond]); draw(); };
$('phi').addEventListener('input', e => setPhi(+e.target.value));
$('curve').addEventListener('click', e => { // tap the curve: turn to that angle (5° steps, 360° = 0°)
  const r = e.currentTarget.getBoundingClientRect(), d = ((e.clientX - r.left) / r.width * CW - CL) / (CW - CL - CR) * 360;
  setPhi(Math.round(Math.min(360, Math.max(0, d)) / 5) * 5 % 360);
});
const retheme = () => { viewer.setBackgroundColor(tok('--stage')); draw(); };
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', retheme);
window.addEventListener('resize', () => { viewer.resize(); viewer.render(); });
if (location.hash) run(fromHash());
