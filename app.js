import { GROUPS, KIT_PARTS, IDEAL, conformer, stericPairs, diaxial, moreStable, huckel, CONE, sn2Rate, sn1Rate } from './chem.js';
import { resolveName } from './lookup.js';

const $ = id => document.getElementById(id);
const empty = $('empty'), q = $('q');
// guess: the guess-first mode; tries[i] = { pick: 'R' | 'S', hint: true } for stereo row i in that mode
// chair: selected chair row; flipped: the other chair is shown; cguess/ctries: guess-first for which chair is more stable
const state = { style: 'ball', h: true, tab: 'structure', hl: -1, atom: -1, st: -1, guess: false, tries: {}, mirror: false, bond: -1, phi: 0, sguess: false, stries: {}, ring: -1,
  chair: -1, flipped: false, cguess: false, ctries: {},
  // Structure-tab views, one at a time: pol (bond dipoles), ar (aromatic ring index), nmr (equivalent H labels);
  // pguess/ppick and nguess/npick: guess-first for polar/nonpolar and for the number of signals
  pol: false, ar: -1, nmr: false, pguess: false, ppick: null, nguess: false, npick: null,
  // Stereo tab: cmp = the compared structure and its relation; side = it is shown beside the model (model 1)
  cmp: null, side: false,
  // plane: the meso mirror plane is shown (cur.mplane: undefined = not looked for yet, null = none on reachable shapes)
  plane: false,
  // guess first: rguess/rpick for the relation, fguess with fdl (D/L pick) or fpicks[row] (R/S picks) for the Fischer projection
  rguess: false, rpick: null, fguess: false, fdl: null, fpicks: {},
  // Reactions tab: rp = principle (-1 none, 0 SN2, 1 SN1, 2 E2, 3 addition), rx = C–X row; xguess/xpicks: guess first, picks by result key
  rp: -1, rx: 0, xguess: false, xpicks: {},
  sn1v: 'cation', // SN1: show the carbocation or the two products
  rh: -1, // E2: the β-H row picked
  rdb: 0, rgt: 'br2' }; // addition: the C=C row and the reagent
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
  if (state.tab === 'reactions') return rxAtoms();
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
  // the compared structure and the Reactions tab's intermediates and products are never dimmed
  const lit = a => !set || set.has(a.index) || (a.model > 0 && (state.side || state.tab === 'reactions'));
  viewer.setStyle({}, {});
  viewer.setStyle({ predicate: a => shown(a) && lit(a) }, styleFor(null));
  if (set) viewer.setStyle({ predicate: a => shown(a) && !lit(a) }, styleFor(tok('--dim')));
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
  showPlane();
  showFischer();
  showRel();
  showNewman();
  showChair();
  showReact();
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
  drawPlane();
  drawRx();
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
  const mirror = state.mirror, side = state.side;
  if (mirror) setMirror(false, false);
  if (side) setSide(false, false);
  const c = atomAt(r.atom), l = atomAt(r.order[3]), v = viewer.getView(), q0 = v.slice(4, 8);
  const d = qRot(q0, unit([l.x - c.x, l.y - c.y, l.z - c.z])); // where group 4 points on screen now
  animateView([-c.x, -c.y, -c.z, v[3], ...qMul(qArc(d, [0, 0, -1]), q0)], 600, mirror ? () => setMirror(true) : side ? () => setSide(true) : null);
}

// ---- Mirror image (Stereo tab): a reflected copy (model 1) to the right of the model. It is built from the current
// view, so the mirror plane stands upright on screen between the two; both then turn together.
// The screen's right, up and toward-the-viewer directions in model coordinates.
function screenAxes() {
  const q = viewer.getView().slice(4, 8), qi = [-q[0], -q[1], -q[2], q[3]];
  return [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(v => qRot(qi, v));
}

function setMirror(on, zoom = true) {
  if (on && state.plane) planeOff(); // the mirror image is built from the model as built
  state.mirror = on;
  if (on) { state.side = false; $('btn-side').setAttribute('aria-pressed', 'false'); } // the compared structure is model 1 too
  $('btn-mirror').setAttribute('aria-pressed', String(on));
  $('mirror-note').hidden = !on;
  const old = viewer.getModel(1); // 3Dmol returns the last model for a missing id, so check it is really model 1
  if (old?.getID() === 1) viewer.removeModel(old);
  viewer.removeAllShapes();
  if (on) {
    const [n, up, out] = screenAxes();
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

// ---- Compare two structures (Stereo tab): the relation, and the second structure beside the model as model 1 (never
// together with the mirror image, which is model 1 as well).
const REL = {
  same: ['Same compound', 'the two structures match, stereo included. Two models that differ only by turning about single bonds or by a ring flip (conformations) are the same compound.'],
  enantiomers: ['Enantiomers', 'mirror images that do not superimpose. Each stereocenter has the opposite configuration (R ↔ S), while double bonds keep E/Z.'],
  diastereomers: ['Diastereomers', 'stereoisomers (atoms joined the same way) that are not mirror images of each other.'],
  constitutional: ['Constitutional isomers', f => `the same formula (${formulaHTML(f[0], 0)}), with the atoms joined in a different order.`],
  'not isomers': ['Not isomers', f => `the formulas differ (${formulaHTML(f[0], 0)} and ${formulaHTML(f[1], 0)}).`],
};
let cmpSeq = 0;

function showRel() {
  const c = state.cmp, box = $('rel');
  $('btn-side').hidden = !c;
  if (!c) return; // lookup messages stay until a result arrives
  const r = c.rel, [title, text] = REL[r.kind] || [];
  const which = r.open ? (r.open[0] && r.open[1] ? 'either structure' : r.open[0] ? `the first structure (${esc(cur.title)})` : `the second structure (${esc(c.title)})`) : '';
  const tail = `<p class="small">${state.side ? `Left: ${esc(cur.title)} · Right: ${esc(c.title)}` : `Compared with ${esc(c.title)}`}${c.stripped ? ' (a salt or mixture: only its largest part)' : ''}.</p>`;
  if (state.rguess && !state.rpick && r.kind !== 'open') { // nothing to guess when the stereo is not set
    box.innerHTML = '<p class="facts">How are the two structures related? Pick one.</p><div class="row">'
      + Object.entries(REL).map(([k, [t]]) => `<button data-rpick="${k}">${t}</button>`).join('') + '</div>' + tail;
    box.querySelectorAll('[data-rpick]').forEach(b => b.addEventListener('click', () => { state.rpick = b.dataset.rpick; draw(); }));
    return;
  }
  box.innerHTML = (state.rguess && state.rpick && r.kind !== 'open' ? `<p class="facts">${state.rpick === r.kind ? 'Right' : 'Not quite'}.</p>` : '')
    + (r.kind === 'open'
    ? `<p class="facts"><b>Cannot tell</b>: the stereo (R/S, E/Z or cis/trans) of ${which} is not set in the name, so the relationship cannot be told. Add it to the name, as in (R)-2-butanol or cis-2-butene.</p>`
    : `<p class="facts"><b>${title}</b>: ${typeof text === 'function' ? text(r.formulas) : text}</p>`) + tail;
}

// Show the compared structure to the right of the model, as the 3D build made it (not turned to match).
function setSide(on, zoom = true) {
  if (on && state.mirror) setMirror(false, false);
  if (on && state.plane) planeOff();
  state.side = on && !!state.cmp;
  $('btn-side').setAttribute('aria-pressed', String(state.side));
  const old = viewer.getModel(1); // 3Dmol returns the last model for a missing id, so check it is really model 1
  if (old?.getID() === 1) viewer.removeModel(old);
  if (state.side) {
    const [n] = screenAxes(), dot3 = (p, u) => p[0] * u[0] + p[1] * u[1] + p[2] * u[2];
    const lines = state.cmp.molfile.split('\n'), count = parseInt(lines[3].slice(0, 3));
    const mean = P => [0, 1, 2].map(k => P.reduce((t, p) => t + p[k], 0) / P.length);
    const a = viewer.getModel(0).selectedAtoms({}).map(x => [x.x, x.y, x.z]);
    const b = lines.slice(4, 4 + count).map(l => l.slice(0, 30).trim().split(/\s+/).map(Number)); // V2000 atom lines: x, y, z first
    const ca = mean(a), cb = mean(b);
    const gap = Math.max(...a.map(p => dot3(p, n) - dot3(ca, n))) - Math.min(...b.map(p => dot3(p, n) - dot3(cb, n))) + 2.5;
    const shift = [0, 1, 2].map(k => ca[k] + gap * n[k] - cb[k]);
    for (let i = 0; i < count; i++) lines[4 + i] = b[i].map((x, k) => (x + shift[k]).toFixed(4).padStart(10)).join('') + lines[4 + i].slice(30);
    viewer.addModel(lines.join('\n'), 'sdf', { keepH: true });
  }
  if (zoom) viewer.zoomTo({}, 400);
  draw();
}

$('cmp-form').addEventListener('submit', async e => {
  e.preventDefault();
  const name = $('q2').value.trim(), mol = cur, my = ++cmpSeq, box = $('rel');
  if (!name || !mol) return;
  $('cmp-go').disabled = true;
  state.cmp = null; if (state.side) setSide(false, false);
  box.innerHTML = '<p class="small">Looking up the name…</p>';
  try {
    const hit = await resolveName(name);
    if (my !== cmpSeq || mol !== cur) return; // a newer comparison or molecule took over
    if (!hit) { box.innerHTML = `<p class="small">No structure found for “${esc(name)}”. Check the spelling, or try the IUPAC name.</p>`; return; }
    box.innerHTML = '<p class="small">Building the 3D model…</p>';
    const res = await ask({ compare: hit.smiles });
    if (res.stale || my !== cmpSeq || mol !== cur) return;
    if (res.error) { box.innerHTML = `<p class="small err">${esc(res.error)}</p>`; return; }
    state.cmp = { title: hit.title, ...res.compare }; state.rpick = null;
    setSide(true);
  } catch (err) {
    if (my === cmpSeq) box.innerHTML = `<p class="small err">${esc(err.message || 'The lookup failed. Check the connection and try again.')}</p>`;
  } finally {
    if (my === cmpSeq) $('cmp-go').disabled = false;
  }
});
$('btn-side').addEventListener('click', () => setSide(!state.side));
$('btn-rguess').addEventListener('click', e => {
  state.rguess = !state.rguess; state.rpick = null;
  e.currentTarget.setAttribute('aria-pressed', String(state.rguess));
  draw();
});

// ---- Fischer projection (Stereo tab): drawn from chem.js fischer() as SVG, with R/S beside each stereocenter.
const FISCHER_WHY = {
  none: 'No stereocenter, so there is no Fischer projection to draw.',
  ring: 'No Fischer projection: the stereocenters are in a ring.',
  chain: 'No Fischer projection: the stereocenters are not all on one carbon chain.',
  end: 'No Fischer projection here: a stereocenter needs a carbon above and below it on the chain.',
};
// A group on the left of a cross reads toward the center: OH → HO, NH₂ → H₂N, CHO → OHC, COOH → HOOC, CH₂OH → HOH₂C.
function leftLabel(t) {
  const units = t.match(/[A-Z][a-z]?(?:H\d*)?[⁺⁻]?/g) || [];
  if (t === 'H' || units.join('') !== t) return t; // parentheses, =, or a ring formula: as is
  const flip = u => u.replace(/^([A-Z][a-z]?)(H\d*)([⁺⁻]?)$/, '$2$1$3');
  const r = units.reverse();
  r[0] = flip(r[0]); if (r.length > 1) r[r.length - 1] = flip(r[r.length - 1]);
  return r.join('');
}
const RH = 46, FP = 24, FX = 122, FW = 300; // row height, padding, the vertical line's x, width
function fischerSVG(f, tag) {
  const n = f.rows.length, H = 2 * FP + (n - 1) * RH, Y = i => FP + i * RH;
  const text = (x, y, s, anchor = 'middle') => `<text x="${x}" y="${y + 5}" text-anchor="${anchor}">${esc(SUB(s))}</text>`;
  let lines = '', words = '', tags = '';
  if (f.mirror !== null) lines += `<line class="mp" x1="8" x2="${FW - 44}" y1="${FP + f.mirror * RH}" y2="${FP + f.mirror * RH}"/>`;
  for (let i = 0; i + 1 < n; i++) lines += `<line x1="${FX}" x2="${FX}" y1="${Y(i) + (f.rows[i].cross ? 0 : 12)}" y2="${Y(i + 1) - (f.rows[i + 1].cross ? 0 : 12)}"/>`;
  f.rows.forEach((r, i) => {
    if (!r.cross) { words += text(FX, Y(i), r.label); return; }
    lines += `<line x1="${FX - 44}" x2="${FX + 44}" y1="${Y(i)}" y2="${Y(i)}"/>`;
    words += text(FX - 50, Y(i), leftLabel(r.left.label), 'end') + text(FX + 50, Y(i), r.right.label, 'start');
    const t = tag(i), q = t?.q ? ' class="q"' : '';
    if (t) tags += `<rect${q} x="${FW - 34}" y="${Y(i) - 12}" width="24" height="24" rx="5"/><text${q} x="${FW - 22}" y="${Y(i) + 5}" text-anchor="middle">${esc(t.text)}</text>`;
  });
  const said = f.rows.map((r, i) => (r.cross ? `${r.left.label} left, ${r.right.label} right${tag(i) && !tag(i).q ? ` (${r.rs})` : ''}` : r.label)).join('; ');
  return `<svg viewBox="0 0 ${FW} ${H}" role="img" aria-label="${esc(`Fischer projection, top to bottom: ${said}`)}">
    <g stroke="var(--fg)" stroke-width="1.6" stroke-linecap="round">${lines}</g>
    <g font-family="var(--sans)" font-size="15" fill="var(--fg)" stroke="var(--bg)" stroke-width="5" paint-order="stroke" stroke-linejoin="round">${words}</g>
    <g class="tags" font-family="var(--sans)" font-size="14" font-weight="600">${tags}</g></svg>`;
}

function showFischer() {
  const f = cur.fischer, box = $('fischer'), out = $('fischer-out');
  const crosses = f.rows ? f.rows.flatMap((r, i) => (r.cross ? [i] : [])) : [], named = crosses.filter(i => f.rows[i].rs);
  $('fischer-tools').hidden = !!f.why;
  $('btn-fview').hidden = crosses.length !== 1;
  $('btn-fguess').hidden = !f.dl && !named.length; // nothing to guess
  if (f.why) { box.innerHTML = `<p class="facts">${esc(FISCHER_WHY[f.why])}</p>`; out.textContent = ''; return; }
  // Guess first: D/L when the drawing has one (its pick opens R/S too), otherwise R/S at each named center, numbered from the top.
  const g = state.fguess, dlOpen = !g || !f.dl || !!state.fdl;
  const tag = i => {
    const r = f.rows[i];
    if (!r.rs || (state.guess && !revealed(r.row))) return null; // the stereocenter rows' own guess still hides it
    if (!g || (f.dl ? state.fdl : state.fpicks[i])) return { text: r.rs };
    return f.dl ? null : { text: String(named.indexOf(i) + 1), q: true };
  };
  box.innerHTML = fischerSVG(f, tag);
  const notes = [];
  if (g && f.dl && !state.fdl) notes.push(`<p class="facts">D or L? ${f.dl.kind === 'sugar' ? 'Look at the OH on the lowest stereocenter.' : 'Look at the NH₂ on the α-carbon.'}</p><div class="row"><button data-fdl="D">D</button><button data-fdl="L">L</button></div>`);
  if (g && f.dl && state.fdl) notes.push(`<p class="facts">${state.fdl === f.dl.dl ? 'Right' : 'Not quite'}: it is ${f.dl.dl}.</p>`);
  if (g && !f.dl) named.forEach((i, k) => {
    const pick = state.fpicks[i], rs = f.rows[i].rs;
    notes.push(pick ? `<p class="facts">Stereocenter ${k + 1}: ${pick === rs ? 'Right' : 'Not quite'}, it is <b>${rs}</b>.</p>`
      : `<p class="facts">Stereocenter ${k + 1} (counting from the top): R or S?</p><div class="row"><button data-fpick="${i}" data-v="R">R</button><button data-fpick="${i}" data-v="S">S</button></div>`);
  });
  if (f.dl && dlOpen) notes.push(`<p class="facts"><b>${f.dl.dl}</b>-${f.dl.kind === 'sugar' ? 'sugar: the OH on the lowest stereocenter points' : 'amino acid: the NH₂ on the α-carbon points'} ${f.dl.dl === 'D' ? 'right' : 'left'}.</p>`);
  if (f.mirror !== null) notes.push('<p class="facts">The top half mirrors the bottom half across the dashed line: a meso compound.</p>');
  if (f.ringLeft) notes.push('<p class="small">Stereocenters in a ring are left out.</p>');
  if (f.open) notes.push('<p class="small">Some stereocenters are not set in the source structure, so this shows the ones the 3D build picked.</p>');
  out.innerHTML = notes.join('');
  out.querySelectorAll('[data-fdl]').forEach(b => b.addEventListener('click', () => { state.fdl = b.dataset.fdl; draw(); }));
  out.querySelectorAll('[data-fpick]').forEach(b => b.addEventListener('click', () => { state.fpicks[b.dataset.fpick] = b.dataset.v; draw(); }));
}
$('btn-fguess').addEventListener('click', e => {
  state.fguess = !state.fguess; state.fdl = null; state.fpicks = {};
  e.currentTarget.setAttribute('aria-pressed', String(state.fguess));
  draw();
});

// One stereocenter: turn the model so its chain neighbors point up-and-back and down-and-back and its other two groups
// come toward the viewer, as the cross draws it. Only the view turns; the coordinates stay.
function fischerView() {
  const r = cur.fischer.rows.find(x => x.cross), P = i => { const a = atomAt(i); return [a.x, a.y, a.z]; };
  const c = P(r.atom), u = P(r.up).map((x, k) => x - c[k]), d = P(r.down).map((x, k) => x - c[k]);
  const y = unit(u.map((x, k) => x - d[k])), b = u.map((x, k) => -(x + d[k])), t = b[0] * y[0] + b[1] * y[1] + b[2] * y[2];
  const zv = unit(b.map((x, k) => x - t * y[k])); // toward the viewer: away from where the two vertical bonds point
  const mirror = state.mirror, side = state.side;
  if (mirror) setMirror(false, false);
  if (side) setSide(false, false);
  const v = viewer.getView(), q1 = qMul(qArc(qRot(v.slice(4, 8), y), [0, 1, 0]), v.slice(4, 8)), q2 = qMul(qArc(qRot(q1, zv), [0, 0, 1]), q1);
  animateView([-c[0], -c[1], -c[2], v[3], ...q2], 600, mirror ? () => setMirror(true) : side ? () => setSide(true) : null);
}
$('btn-fview').addEventListener('click', fischerView);

// ---- Meso mirror plane (Stereo tab): found in the worker on real coordinates. When the model as built has none, the
// worker turns single bonds (or tries other conformers) until one shows, and the model slides there; off slides it back.
let planeBusy = false;
const PLANE_HOW = {
  model: '',
  turned: ' To show it, the model was turned about single bonds, as on the Shape tab. Turning changes the conformation, not the compound.',
  conformer: ' This is another conformation from the conformer search, not the one built first.',
};
function showPlane() {
  const btn = $('btn-plane'), box = $('plane-note'), p = cur.mplane;
  btn.hidden = !cur.meso;
  btn.disabled = planeBusy;
  btn.setAttribute('aria-pressed', String(state.plane));
  box.hidden = !planeBusy && !state.plane;
  box.textContent = planeBusy ? 'Looking for the mirror plane…'
    : !state.plane ? ''
    : p ? 'Mirror plane: the see-through plate cuts the molecule into two halves that are mirror images of each other, so the compound is achiral (meso) even though it has stereocenters.'
      + PLANE_HOW[p.how] + (cur.mring ? ' Textbooks check a ring drawn flat, its average shape; the ring in this model has the plane too.' : '')
    : 'This compound has a mirror plane in another conformation.'
      + (cur.mring ? ' For a ring, that is the flat ring the textbook draws: the average of the shapes the ring flips between. No shape this page can reach has it, since ring bonds do not turn.' : '');
}

$('btn-plane').addEventListener('click', async () => {
  if (state.plane) return setPlane(false);
  const mol = cur;
  if (mol.mplane === undefined) {
    planeBusy = true; draw();
    const res = await ask({ plane: true });
    planeBusy = false;
    if (res.stale || mol !== cur) return draw(); // another request took over: try again later
    if (res.error) { mol.mplane = null; mol.mring = false; } else { mol.mplane = res.plane; mol.mring = res.ring; }
  }
  setPlane(true);
});

function setPlane(on) {
  if (on && state.mirror) setMirror(false, false);
  if (on && state.side) setSide(false, false);
  const p = cur.mplane, moved = p && p.how !== 'model';
  state.plane = on;
  let view = null;
  if (on && p) { // look at the plane almost edge-on (turned 30°), so one half sits left of it and its mirror half right
    const pts = moved ? p.xyz : home3d, g = [0, 1, 2].map(k => pts.reduce((t, v) => t + v[k], 0) / pts.length);
    const v = viewer.getView(), q0 = v.slice(4, 8), t = 30 * Math.PI / 180, q = qMul(qArc(qRot(q0, p.n), [Math.cos(t), 0, Math.sin(t)]), q0);
    planeUp = qRot([-q[0], -q[1], -q[2], q[3]], [0, 1, 0]); // the screen's up once the turn ends, for the plate's edges
    view = [-g[0], -g[1], -g[2], v[3], ...q];
  }
  if (moved) morph(on ? home3d : p.xyz, on ? p.xyz : home3d); else draw();
  if (view) animateView(view);
}
// Off at once, without a slide (another view takes over, or the tab changes).
function planeOff() { state.plane = false; turn(null); }

// The plate: through the plane's point, as wide as the model in the plane plus 1 Å each side. Drawn when the slide ends.
let planeShapes = [], planeUp = [0, 1, 0];
function drawPlane() {
  planeShapes.forEach(x => viewer.removeShape(x));
  planeShapes = [];
  const p = cur?.mplane;
  if (state.tab !== 'stereo' || !state.plane || !p || state.sliding) return;
  const n = p.n, d3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const pts = viewer.getModel(0).selectedAtoms({}).map(a => [a.x, a.y, a.z]);
  const g = [0, 1, 2].map(k => pts.reduce((t, v) => t + v[k], 0) / pts.length), off = d3(g.map((x, k) => x - p.c[k]), n);
  const ctr = g.map((x, k) => x - off * n[k]);
  const up = planeUp, k = d3(up, n), u = Math.abs(k) < 0.99 ? unit(up.map((x, i) => x - k * n[i])) : unit(cross3(n, [1, 0, 0]));
  const w = cross3(n, u); // two directions in the plane, one as close to the screen's up as it gets
  const size = e => 2 * Math.max(...pts.map(v => Math.abs(d3(v.map((x, k) => x - ctr[k]), e)))) + 2;
  const vec = (e, len) => ({ x: e[0] * len, y: e[1] * len, z: e[2] * len });
  planeShapes.push(viewer.addBox({ center: { x: ctr[0], y: ctr[1], z: ctr[2] }, dimensions: { w: vec(n, 0.04), h: vec(u, size(u)), d: vec(w, size(w)) },
    color: tok('--accent'), opacity: 0.4 }));
}
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

// Leaving the Stereo tab takes its model views away: the mirror image, the compared structure, the mirror plane.
function leaveStereo() {
  if (state.plane) planeOff();
  if (state.mirror) setMirror(false); else if (state.side) setSide(false); else draw();
}

// ---- Reactions tab: one principle at a time (state.rp): 0 backside attack (SN2), 1 flat intermediate (SN1), 2 anti-periplanar
// H and X (E2), 3 faces of a double bond. The worker builds intermediates and products on request (chem.js); cur.rxr keeps each
// result by request (null = it could not be built). Models 1 and 2 (intermediate, products) sit to the right of the model.
const P3 = (i, m = 0) => { const a = atomAt(i, m); return [a.x, a.y, a.z]; };
const xyz3 = v => ({ x: v[0], y: v[1], z: v[2] });
const DEG = ['methyl', '1°', '2°', '3°'];
const degName = h => (h.deg || h.plain ? `${DEG[h.deg]} carbon` : 'carbon with no other carbon on it');
const centerRS = a => cur.stereo.find(r => r.kind === 'center' && r.atom === a) || null;

// The worker request for the current pick, its cache key, and its result (undefined = not asked yet).
function rxReq() {
  const h = cur.hal.rows[state.rx];
  if (state.rp === 0 && h) return { kind: 'sn2', i: state.rx };
  if (state.rp === 1 && h && sn1Rate(h) !== 'none') return { kind: 'sn1', i: state.rx };
  if (state.rp === 2 && h) return { kind: 'e2', i: state.rx };
  if (state.rp === 3 && cur.dbl.rows[state.rdb]) return { kind: 'add', i: state.rdb, reagent: state.rgt };
  return null;
}
const rxKey = () => { const q = rxReq(); return q && JSON.stringify(q); };
const rxRes = () => { const k = rxKey(); return k ? cur.rxr[k] : undefined; };

// Atoms kept in color on the model.
function rxAtoms() {
  const h = cur.hal.rows[state.rx];
  if (state.rp === 0 && h) return [h.c, h.x, ...h.cone];
  if (state.rp === 1 && h) return [h.x, ...cur.atoms[h.c].nbrs];
  if (state.rp === 2 && h) {
    const r = rxRes(), p = e2Pick();
    if (p) return [h.c, h.x, r.betas[p.row.beta], ...p.row.hs.map(i => r.hs[i].h)];
    return [h.c, h.x, ...(r?.betas || []).flatMap(b => [b, ...cur.atoms[b].nbrs.filter(k => cur.atoms[k].el === 'H')])];
  }
  const d = state.rp === 3 && cur.dbl.rows[state.rdb];
  if (d) return d.atoms.flatMap(a => [a, ...cur.atoms[a].nbrs]);
  return null;
}

let rxBusy = null; // cache key of the request the worker is on
async function rxView() {
  const k = state.tab === 'reactions' && rxKey(), mol = cur;
  applyRx(true); // the view turns now; models come when the worker is done
  if (!k || k in mol.rxr) return;
  rxBusy = k; draw();
  const res = await ask({ react: rxReq() });
  if (res.stale || mol !== cur) return; // a newer request took over
  rxBusy = null;
  mol.rxr[k] = res.error ? null : res.react;
  if (rxKey() === k) applyRx(false); else draw();
}

// Turn the view so model direction u points along screen direction t; returns the new rotation.
const viewTurn = (u, t) => { const q0 = viewer.getView().slice(4, 8); return qMul(qArc(qRot(q0, u), t), q0); };
// The rotation that puts model direction `right` along the screen's x and `up` (made perpendicular to it) along y, then
// tips the top `tilt` degrees toward you.
function frameTurn(right, up, tilt = 0) {
  const x = unit(right), d = up[0] * x[0] + up[1] * x[1] + up[2] * x[2], y = unit(up.map((v, k) => v - d * x[k]));
  const z = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
  const [m00, m01, m02] = x, [m10, m11, m12] = y, [m20, m21, m22] = z, tr = m00 + m11 + m22; // rows: model → screen
  let q;
  if (tr > 0) { const k = 2 * Math.sqrt(tr + 1); q = [(m21 - m12) / k, (m02 - m20) / k, (m10 - m01) / k, k / 4]; }
  else if (m00 > m11 && m00 > m22) { const k = 2 * Math.sqrt(1 + m00 - m11 - m22); q = [k / 4, (m01 + m10) / k, (m02 + m20) / k, (m21 - m12) / k]; }
  else if (m11 > m22) { const k = 2 * Math.sqrt(1 + m11 - m00 - m22); q = [(m01 + m10) / k, k / 4, (m12 + m21) / k, (m02 - m20) / k]; }
  else { const k = 2 * Math.sqrt(1 + m22 - m00 - m11); q = [(m02 + m20) / k, (m12 + m21) / k, k / 4, (m10 - m01) / k]; }
  const t = -tilt * Math.PI / 360; // half angle; negative turns the top toward the viewer (+z)
  return qMul([Math.sin(t), 0, 0, Math.cos(t)], unit(q));
}

// Put the model back as built, add the intermediate or products beside it, and turn the view for the pick.
function applyRx(animate = true, slide = animate) {
  if (state.tab !== 'reactions') return;
  const r = rxRes(), h = cur.hal.rows[state.rx];
  turn(null);
  let models = [], q = null, shape = null;
  if (state.rp === 0 && h) {
    const C = home3d[h.c], X = home3d[h.x];
    q = viewTurn(unit(C.map((x, k) => x - X[k])), [-1, 0, 0]); // the back side to the left: the nucleophile comes in from there
    if (r && sn2Rate(h) !== 'none') models = [r.product.molfile];
  }
  if (state.rp === 1 && h && r) {
    const t = 20 * Math.PI / 180;
    if (r.cation) q = viewTurn(r.cation.n, [0, Math.cos(t), Math.sin(t)]); // the p orbital's axis up, tilted 20° toward you
    models = sn1View() === 'cation' ? [r.cation.molfile] : sn1View() === 'products' ? r.products.map(p => p.molfile) : [];
  }
  const p = state.rp === 2 && h && e2Pick();
  if (p) {
    const b = r.betas[p.row.beta], ci = cur.chairs.findIndex(c => c.atoms.includes(h.c) && c.atoms.includes(b));
    if (!p.q.ring) { // down Cα→Cβ with X up, so the anti H points down; turned 25° about the screen's up so the bond shows
      const t = 12.5 * Math.PI / 180;
      shape = p.q.xyz; q = qMul([0, Math.sin(t), 0, Math.cos(t)], lookDown(home3d[h.c], home3d[b], home3d[h.x]));
    }
    else if (ci >= 0) {
      const t = 20 * Math.PI / 180;
      q = viewTurn(ringFrame(cur.chairs[ci].atoms.map(i => home3d[i])).n, [0, Math.cos(t), Math.sin(t)]); // chair from the side
      if (p.q.ok && !p.q.now) { if (!(ci in cur.flips)) e2Flip(ci); else if (cur.flips[ci]) shape = cur.flips[ci].xyz; } // X axial in the other chair
    }
    if (p.q.ok && p.q.molfile && e2Open()) models = [p.q.molfile];
  }
  const d = state.rp === 3 && cur.dbl.rows[state.rdb];
  if (d) {
    q = frameTurn(home3d[d.atoms[1]].map((x, k) => x - home3d[d.atoms[0]][k]), addNormal(d), 20); // C=C across, p orbitals up and down
    if (r?.products && xOpen(rxKey())) models = r.products.map(x => x.molfile);
  }
  if (animate || !rxQ) rxQ = q; // the turn the models are placed for, kept while the view moves there
  if (shape) { if (slide) morph(home3d, shape); else viewer.getModel(0).selectedAtoms({}).forEach((a, i) => { [a.x, a.y, a.z] = shape[i]; }); }
  rxModels(models, rxQ);
  draw();
  rxFrame(rxQ);
}

// The normal of a C=C's plane at its first carbon, from the model as built (the same as chem.js addition, in the viewer's frame).
function addNormal(d) {
  const [p, q] = d.atoms, r = cur.atoms[p].nbrs.find(a => a !== q), P = home3d[p];
  const u = home3d[q].map((x, k) => x - P[k]), v = home3d[r].map((x, k) => x - P[k]);
  return unit([u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]);
}

// E2: the H rows. H on one β carbon that give the same alkene with the same verdict share a row (a CH₃'s three H).
function e2Rows(r) {
  const rows = [];
  r.hs.forEach((q, i) => {
    const g = rows.find(x => x.beta === q.beta && x.alkene === q.alkene && x.ok === q.ok && x.now === q.now);
    if (g) g.hs.push(i); else rows.push({ beta: q.beta, alkene: q.alkene, ok: q.ok, now: q.now, hs: [i] });
  });
  return rows;
}
// The picked row, its result and its first H (the one shown), or null.
function e2Pick() {
  const r = state.rp === 2 && rxRes(), row = r && e2Rows(r)[state.rh];
  return row ? { r, row, q: r.hs[row.hs[0]] } : null;
}
const e2Key = () => `${rxKey()}:${state.rh}`;
const e2Open = () => xOpen(e2Key());
// The other chair for an H that is trans to X but equatorial with it now (the Shape tab's ring flip, shared cache).
async function e2Flip(ci) {
  if (flipping === ci) return;
  const mol = cur;
  flipping = ci;
  const res = await ask({ flip: ci });
  if (flipping === ci) flipping = -1;
  if (res.stale || mol !== cur) return;
  mol.flips[ci] = res.flip || null;
  applyRx(false, true);
}
// Looking from F's side down the F→B bond, with F→U pointing up the screen (as alignBond does for the Newman view).
function lookDown(F, B, U) {
  const q0 = viewer.getView().slice(4, 8), q1 = qMul(qArc(qRot(q0, unit(B.map((c, j) => c - F[j]))), [0, 0, -1]), q0);
  const v = qRot(q1, U.map((c, j) => c - F[j])), h = (Math.PI / 2 - Math.atan2(v[1], v[0])) / 2;
  return qMul([0, 0, Math.sin(h), Math.cos(h)], q1);
}
let rxQ = null;

// Move the view to rotation q (null: keep the rotation) with everything in view, as one smooth move. 3Dmol's animated
// zoomTo adds steps to the view, so two at once (or one during our turn) overshoot; its instant zoomTo only finds the target.
function rxFrame(q) {
  const v0 = viewer.getView();
  if (q) viewer.setView([...v0.slice(0, 4), ...q]);
  viewer.zoomTo({});
  viewer.zoom(0.9); // a margin: zoomTo fits the models to the edge
  const to = viewer.getView();
  viewer.setView(v0);
  animateView(to, 600);
}

// Models 1, 2 …: each molfile in a row to the right (screen right of rotation q, or of the view now), 2.5 Å apart.
// Only shifted, so a product built onto the reactant keeps its turn.
function rxModels(list, q) {
  for (let m = viewer.getModel(); m && m.getID() > 0; m = viewer.getModel()) viewer.removeModel(m); // highest id first: 3Dmol then reuses the ids
  if (!list.length) return;
  const w = q || viewer.getView().slice(4, 8), r = qRot([-w[0], -w[1], -w[2], w[3]], [1, 0, 0]);
  const along = p => p[0] * r[0] + p[1] * r[1] + p[2] * r[2];
  let edge = Math.max(...viewer.getModel(0).selectedAtoms({}).map(a => along([a.x, a.y, a.z])));
  for (const mf of list) {
    const lines = mf.split('\n'), count = parseInt(lines[3].slice(0, 3));
    const b = lines.slice(4, 4 + count).map(l => l.slice(0, 30).trim().split(/\s+/).map(Number)); // V2000 atom lines: x, y, z first
    const s = edge + 2.5 - Math.min(...b.map(along));
    for (let i = 0; i < count; i++) lines[4 + i] = b[i].map((x, k) => (x + s * r[k]).toFixed(4).padStart(10)).join('') + lines[4 + i].slice(30);
    viewer.addModel(lines.join('\n'), 'sdf', { keepH: true });
    edge = Math.max(...b.map(p => along(p) + s));
  }
}
const rxShown = () => { const m = viewer.getModel(); return m ? m.getID() : 0; }; // how many models sit beside the model

function leaveReactions() { rxModels([]); turn(null); draw(); }

// Shapes and tags on the models. Only these shapes are removed here.
let rxShapes = [];
const rxTag = (text, at, color = '--accent', offset = { x: 10, y: 10 }, alignment = 'bottomLeft') => viewer.addLabel(text, { position: xyz3(at), alignment,
  screenOffset: offset, fontSize: 15, bold: true, fontColor: tok('--bg'), backgroundColor: tok(color), backgroundOpacity: 1, inFront: true });
function drawRx() {
  rxShapes.forEach(s => viewer.removeShape(s));
  rxShapes = [];
  if (state.tab !== 'reactions' || !cur || state.sliding) return;
  const h = cur.hal.rows[state.rx], r = rxRes();
  if (state.rp === 0 && h) {
    const C = P3(h.c), X = P3(h.x), d = unit(C.map((x, k) => x - X[k])), at = s => C.map((x, k) => x + s * d[k]);
    rxShapes.push(viewer.addArrow({ start: xyz3(at(4.4)), end: xyz3(at(1.9)), radius: 0.08, radiusRatio: 2.4, mid: 0.72, color: tok('--accent'), opacity: 0.9 }));
    rxTag('OH⁻', at(4.6), '--accent', { x: -14, y: 6 });
    for (const a of h.cone) rxShapes.push(viewer.addSphere({ center: xyz3(P3(a)), radius: 0.5, color: tok('--cut'), opacity: 0.45 }));
    const was = centerRS(h.c);
    if (was) rxTag(was.rs, C);
    if (rxShown() && r.product.rs && sn2Open()) rxTag(r.product.rs, P3(h.c, 1));
  }
  if (state.rp === 1 && h && rxShown() && sn1View() === 'cation') { // the empty p orbital, and water coming to each face
    const C = P3(r.cation.center, 1), n = r.cation.n, at = s => C.map((x, k) => x + s * n[k]);
    for (const s of [1, -1]) {
      rxShapes.push(viewer.addSphere({ center: xyz3(at(0.95 * s)), radius: 0.62, color: tok('--lobe'), opacity: 0.6 }));
      rxShapes.push(viewer.addArrow({ start: xyz3(at(3.4 * s)), end: xyz3(at(1.7 * s)), radius: 0.07, radiusRatio: 2.4, mid: 0.72, color: tok('--accent'), opacity: 0.9 }));
      rxTag('H₂O', at(3.7 * s), '--accent', { x: -16, y: s > 0 ? 6 : -10 });
    }
  }
  if (state.rp === 1 && h && rxShown() && sn1View() === 'products')
    r.products.forEach((p, i) => { if (p.rs) rxTag(p.rs, P3(r.cation.center, i + 1)); });
  if (state.rp === 2 && h && r) { // α and β1, β2 … on the model; a rod joins the picked H and X once they sit anti
    // α up and to the left, β down and to the right: seen down the bond, β sits right behind α
    rxTag('α', P3(h.c), '--fg', { x: -8, y: 10 }, 'bottomRight');
    r.betas.forEach((b, i) => rxTag(`β${i + 1}`, P3(b), '--fg', { x: 8, y: -10 }, 'topLeft'));
    const p = e2Pick();
    const anti = p && (!p.q.ring ? true : p.q.now || (p.q.ok && cur.flips[cur.chairs.findIndex(c => c.atoms.includes(h.c) && c.atoms.includes(r.betas[p.row.beta]))]));
    if (anti) rxShapes.push(viewer.addCylinder({ start: xyz3(P3(p.q.h)), end: xyz3(P3(h.x)), radius: 0.06, color: tok('--accent'), opacity: 0.6, dashed: true }));
  }
  const d = state.rp === 3 && cur.dbl.rows[state.rdb];
  if (d) { // the π bond's p orbitals on both faces, and where the reagent comes from
    const n = addNormal(d), [A, B] = d.atoms.map(a => P3(a)), M = A.map((x, k) => (x + B[k]) / 2), at = (o, s) => o.map((x, k) => x + s * n[k]);
    for (const C of [A, B]) for (const s of [1, -1]) rxShapes.push(viewer.addSphere({ center: xyz3(at(C, 0.85 * s)), radius: 0.5, color: tok('--lobe'), opacity: 0.55 }));
    const label = { h2: ['H₂ / Pd', 'H₂ / Pd'], oso4: ['OsO₄', 'OsO₄'], br2: ['Br₂', 'Br⁻'] }[state.rgt];
    [1, -1].forEach((s, i) => {
      rxShapes.push(viewer.addArrow({ start: xyz3(at(M, 3.4 * s)), end: xyz3(at(M, 1.6 * s)), radius: 0.07, radiusRatio: 2.4, mid: 0.72, color: tok('--accent'), opacity: 0.9 }));
      rxTag(label[i], at(M, 3.7 * s), '--accent', { x: -20, y: s > 0 ? 6 : -10 });
    });
    if (rxShown() && r?.products) r.products.forEach((x, i) => d.atoms.forEach(a => { if (x.rs[a]) rxTag(x.rs[a], P3(a, i + 1)); }));
  }
}

// SN1 view: 'cation' (flat ones only), 'products' (once guessed, in guess mode), or null (nothing to show yet).
function sn1View() {
  const r = rxRes();
  if (!r?.cation) return null;
  const products = r.products && xOpen(rxKey()) ? 'products' : null;
  return r.cation.flat && (state.sn1v === 'cation' || !products) ? 'cation' : products;
}

// Guess first: picks by request key; a result is open once picked (or when the mode is off).
const xOpen = k => !state.xguess || !!state.xpicks[k];
const sn2Open = () => xOpen(rxKey());
const verdictLine = (pick, right, answer) => `<p class="facts">${pick === right ? 'Right' : 'Not quite'}: ${answer}.</p>`;

const SN2_RATE = {
  fast: h => `<b>${h.deg ? '1°' : 'Methyl'}: fast.</b> Little blocks the back side of the carbon.`,
  slow: () => '<b>2°: slow.</b> Two carbon groups crowd the back side.',
  none: h => (h.neo ? '<b>Practically no S<sub>N</sub>2.</b> The carbon is 1°, but its neighbor carries three carbons (neopentyl), and one of them sits right behind it.'
    : '<b>3°: no S<sub>N</sub>2.</b> Three carbon groups block the back side.'),
};
const NO_CX = {
  F: 'C–F: fluoride leaves poorly, so alkyl fluorides rarely react this way.',
  OH: 'C–OH: OH⁻ leaves poorly. An acid first turns –OH into –OH₂⁺, which leaves as water.',
  sp2: 'A halogen on a C=C or ring carbon (vinyl or aryl halide) does not do S<sub>N</sub>1 or S<sub>N</sub>2.',
};

function buildReact() {
  const box = $('rx-rows');
  box.textContent = '';
  cur.hal.rows.forEach((h, i) => {
    const row = document.createElement('button'), others = cur.atoms[h.c].nbrs.filter(k => k !== h.x);
    const names = [...others.filter(k => cur.atoms[k].el !== 'H'), ...others.filter(k => cur.atoms[k].el === 'H')].map(plain);
    row.innerHTML = `<span>${esc(`C–${h.el} on a ${degName(h)}`)}</span><span class="sub">${esc(`On the carbon: ${names.join(', ')}`)}</span>`;
    row.dataset.i = String(i);
    row.addEventListener('click', () => { state.rx = i; state.rh = -1; rxView(); });
    box.appendChild(row);
  });
  const dbl = $('rx-dbl');
  dbl.textContent = '';
  cur.dbl.rows.forEach((d, i) => {
    const row = document.createElement('button'), side = a => { const o = cur.atoms[a].nbrs.filter(k => !d.atoms.includes(k)); return [...o.filter(k => cur.atoms[k].el !== 'H'), ...o.filter(k => cur.atoms[k].el === 'H')].map(plain).join(', '); };
    row.innerHTML = `<span>C=C${cur.rings.some(r => d.atoms.every(a => r.atoms.includes(a))) ? ' in a ring' : ''}</span><span class="sub">${esc(`On its carbons: ${side(d.atoms[0])} | ${side(d.atoms[1])}`)}</span>`;
    row.dataset.i = String(i);
    row.addEventListener('click', () => { state.rdb = i; rxView(); });
    dbl.appendChild(row);
  });
}

function showReact() {
  if (state.tab !== 'reactions') return;
  document.querySelectorAll('#rprin button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.rp === state.rp)));
  const cx = state.rp >= 0 && state.rp < 3, out = $('rx-out');
  $('rx-pick').hidden = cx ? !cur.hal.rows.length : state.rp !== 3 || !cur.dbl.rows.length;
  $('rx-pick-label').textContent = cx ? 'C–X bond' : 'C=C bond';
  $('rx-rows').hidden = !cx; $('rx-dbl').hidden = cx;
  document.querySelectorAll('#rx-rows button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.rx)));
  document.querySelectorAll('#rx-dbl button').forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.i === state.rdb)));
  if (state.rp < 0) { out.innerHTML = '<p class="small">Pick a principle to see it on this molecule.</p>'; return; }
  if (cx && !cur.hal.rows.length) {
    out.innerHTML = '<p class="facts">No C–Cl, C–Br or C–I on an sp³ carbon in this molecule.</p>'
      + cur.hal.why.map(w => `<p class="small">${NO_CX[w]}</p>`).join('');
    return;
  }
  const h = cur.hal.rows[state.rx], r = rxRes(), busy = rxBusy === rxKey();
  if (state.rp === 0) out.innerHTML = sn2HTML(h, r, busy);
  if (state.rp === 1) out.innerHTML = sn1HTML(h, r, busy);
  if (state.rp === 2) out.innerHTML = e2HTML(h, r, busy);
  if (state.rp === 3) out.innerHTML = addHTML(r, busy);
  if (cur.stereo.some(x => !x.given)) // the results hold for the isomer the build picked
    out.insertAdjacentHTML('afterbegin', '<p class="small">The name does not set every stereocenter or stereo double bond, so the 3D build picked one (see the Stereo tab); the results here are for that isomer.</p>');
  out.querySelectorAll('[data-rgt]').forEach(b => b.addEventListener('click', () => { state.rgt = b.dataset.rgt; rxView(); }));
  out.querySelectorAll('[data-rh]').forEach(b => b.addEventListener('click', () => { const i = +b.dataset.rh; state.rh = state.rh === i ? -1 : i; applyRx(true); }));
  out.querySelectorAll('[data-sn1v]').forEach(b => b.addEventListener('click', () => { state.sn1v = b.dataset.sn1v; applyRx(false); }));
  out.querySelectorAll('[data-xpick]').forEach(b => b.addEventListener('click', () => { // a pick can open models (products, the alkene)
    state.xpicks[state.rp === 2 ? e2Key() : rxKey()] = b.dataset.xpick;
    if (state.rp) applyRx(false, false); else draw();
  }));
}

function sn2HTML(h, r, busy) {
  const rate = sn2Rate(h), n = h.cone.length, was = centerRS(h.c);
  let s = `<p class="facts">${rate ? SN2_RATE[rate](h) : '<b>No verdict on the rate.</b> The textbook rule here covers carbons whose other groups are only C and H.'}</p>`
    + `<p class="facts">${n ? `${n} atom${n === 1 ? '' : 's'} (pink) sit${n === 1 ? 's' : ''} behind the carbon in this model, in the path of the nucleophile.` : 'Nothing sits behind the carbon in this model: its groups point to the sides.'}</p>`
    + `<p class="small">The nucleophile, here OH⁻, comes in from the back, 180° from the leaving ${h.el}; any nucleophile comes in the same way. As the C–O bond forms and the C–${h.el} bond breaks, the other three groups flip like an umbrella in the wind, so the configuration is inverted. Pink: atoms within ${CONE.deg}° of that line and ${CONE.r.toFixed(1)} Å of the carbon.</p>`;
  if (rate === 'none') return s + '<p class="small">No product is shown, since S<sub>N</sub>2 does not happen here.</p>';
  if (busy || r === undefined) return s + '<p class="small">Building the product…</p>';
  if (!r) return s + '<p class="small">The product could not be built in 3D here.</p>';
  const rs = r.product.rs, k = rxKey();
  s += '<p class="facts">Product, on the right: the alcohol, its OH where the nucleophile came in.</p>';
  if (!rs) return s + '<p class="small">This carbon is not a stereocenter in the product, so the inversion leaves no R/S to compare.</p>';
  if (!sn2Open()) return s + '<p class="facts">Is the product R or S at this carbon?</p><div class="row"><button data-xpick="R">R</button><button data-xpick="S">S</button></div>';
  return s + (state.xguess ? verdictLine(state.xpicks[k], rs, `it is <b>${rs}</b>`) : '')
    + `<p class="facts">${was ? `The carbon was <b>${was.rs}</b> and is <b>${rs}</b> in the product.` : `The carbon is <b>${rs}</b> in the product.`}</p>`
    + (was ? '<p class="small">Inversion turns the arrangement inside out. The R/S letter comes from the priorities, so it can stay the same when OH ranks in a different place than the leaving group did.</p>' : '');
}

const SN1_RATE = {
  good: () => '<b>3°: works well.</b> Three carbon groups hold up the + of the carbocation.',
  slow: () => '<b>2°: slow.</b> Two carbon groups hold up the + less well.',
  none: h => `<b>${h.deg ? '1°' : 'Methyl'}: no S<sub>N</sub>1.</b> The carbocation would be too unstable to form, so none is shown.`,
  resonance: h => `<b>${h.benzylic ? 'Benzylic' : 'Allylic'}: works.</b> The ${h.benzylic ? 'benzene ring' : 'C=C next to it'} shares the + by resonance.`,
};
const REL_FACES = {
  same: '<b>One product</b>: water gives the same compound from either face.',
  enantiomers: '<b>A pair of enantiomers</b>: racemic, 50:50. The two faces are mirror images, so water reaches each one equally.',
  diastereomers: '<b>Diastereomers</b>: not equal amounts. The rest of the molecule makes the two faces different.',
};
// Guess first for a product relationship: the choices that apply, then the verdict.
const REL_SAID = { same: 'one product', one: 'one product', meso: 'one meso compound', enantiomers: 'a pair of enantiomers', diastereomers: 'a pair of diastereomers' };
function relGuess(choices) {
  return '<p class="facts">How are the products from the two faces related? Pick one.</p><div class="row">'
    + choices.map(([k, t]) => `<button data-xpick="${k}">${t}</button>`).join('') + '</div>';
}

function sn1HTML(h, r, busy) {
  const rate = sn1Rate(h);
  let s = `<p class="facts">${rate ? SN1_RATE[rate](h) : '<b>No verdict.</b> The textbook rule here covers carbons whose other groups are only C and H.'}</p>`
    + `<p class="small">First the C–${h.el} bond breaks and ${h.el}⁻ leaves on its own. The carbon keeps three groups and gets a + charge: a carbocation, sp² and flat, with an empty p orbital above and below the plane. Water then bonds to either lobe and loses H⁺. E1 starts the same way, but the carbocation loses an H⁺ from a neighboring carbon and gives an alkene.</p>`;
  if (rate === 'none') return s;
  if (rate !== 'resonance' && h.shift) s += '<p class="facts">A hydride or methyl shift can move the + to a more stable carbon; products from that are not shown.</p>';
  if (h.allylic && !h.benzylic) s += '<p class="small">In an allylic carbocation the + is shared with the far end of the C=C, so water can also add there; those products are not shown.</p>';
  if (busy || r === undefined) return s + '<p class="small">Building the carbocation and its products…</p>';
  if (!r?.cation) return s + '<p class="small">The carbocation could not be built in 3D here.</p>';
  const c = r.cation, k = rxKey(), view = sn1View();
  if (!c.flat) s += `<p class="small">The force field (MMFF94s+) builds this carbocation slightly pyramidal (its three bond angles add to ${c.sum.toFixed(1)}°; flat is 360°), so it is not drawn. The textbook carbocation is flat.</p>`;
  if (!r.products) return s + '<p class="small">The products could not be built in 3D here.</p>';
  if (c.flat && xOpen(k)) s += `<div class="row"><button data-sn1v="cation" aria-pressed="${view === 'cation'}">Carbocation</button><button data-sn1v="products" aria-pressed="${view === 'products'}">Products</button></div>`;
  if (view === 'cation') s += `<p class="facts">Carbocation, on the right: its three bond angles add to ${c.sum.toFixed(1)}°, so it is flat (360°). Green: the empty p orbital. Arrows: water coming to the top face and to the bottom face.</p>`;
  if (!xOpen(k)) return s + relGuess([['same', 'One product'], ['enantiomers', 'A pair of enantiomers'], ['diastereomers', 'Diastereomers']]);
  if (view === 'products') s += `<p class="facts">Middle: water added to the top face${r.products[0].rs ? ` (${r.products[0].rs})` : ''}. Right: to the bottom face${r.products[1].rs ? ` (${r.products[1].rs})` : ''}.</p>`;
  return s + (state.xguess ? verdictLine(state.xpicks[k], r.rel, `it is ${REL_SAID[r.rel]}`) : '') + `<p class="facts">${REL_FACES[r.rel] || esc(r.rel)}</p>`
    + '<p class="small">A carbon radical is flat (or nearly) too, so a reaction through a radical at this carbon reaches both faces the same way.</p>';
}

const ADD_WHY = {
  h2: '<b>H₂ on Pd: syn.</b> Both H atoms come from the metal surface, so they add to the face of the C=C that lies on the metal.',
  br2: '<b>Br₂: anti.</b> Br⁺ bonds to one face first as a bromonium ion, which blocks that face, so Br⁻ comes in from the opposite face.',
  oso4: '<b>OsO₄: syn.</b> Both O atoms add at once from one face, as a cyclic osmate ester; workup turns them into two OH on the same side.',
};
const ADD_KIND = {
  none: '<b>One product</b> (no new stereocenter): both faces give the same compound.',
  meso: '<b>One meso compound</b>: both faces give the same compound, achiral though it has stereocenters.',
  same: '<b>One product</b>: both faces give the same compound.',
  enantiomers: '<b>A pair of enantiomers</b>: racemic, 50:50. The two faces are mirror images, so the reagent reaches each one equally.',
  diastereomers: '<b>Diastereomers</b>: not equal amounts. The rest of the molecule makes the two faces different.',
};
const addAnswer = k => (k === 'none' || k === 'same' ? 'one' : k);

function addHTML(r, busy) {
  if (!cur.dbl.rows.length) return '<p class="facts">No C=C outside an aromatic ring in this molecule.</p>'
    + (cur.dbl.aromatic ? '<p class="small">A benzene ring’s C=C bonds do not add H₂, Br₂ or OsO₄ like this: adding would break the aromatic ring.</p>' : '');
  const d = cur.dbl.rows[state.rdb], rs = x => d.atoms.map(a => x.rs[a]).filter(Boolean).join(', '), k = rxKey();
  let s = `<div class="row" role="group" aria-label="Reagent">${[['h2', 'H₂, Pd'], ['br2', 'Br₂'], ['oso4', 'OsO₄']].map(([v, t]) => `<button data-rgt="${v}" aria-pressed="${state.rgt === v}">${t}</button>`).join('')}</div>`
    + `<p class="facts">${ADD_WHY[state.rgt]}</p>`
    + '<p class="small">Green: the p orbitals of the π bond, above and below the plane of the C=C, where a reagent can reach it from either face. Only additions of the same atom to both carbons are shown; HBr, hydroboration and halohydrins need a Markovnikov choice and are left out here.</p>';
  if (busy || r === undefined) return s + '<p class="small">Building the products…</p>';
  if (!r?.products) return s + '<p class="small">The products could not be built in 3D here.</p>';
  if (!xOpen(k)) return s + relGuess([['one', 'One product'], ['meso', 'One meso compound'], ['enantiomers', 'A pair of enantiomers'], ['diastereomers', 'Diastereomers']]);
  const [a, b] = r.products, syn = state.rgt !== 'br2';
  s += `<p class="facts">${syn ? `Middle: both atoms added to the top face${rs(a) ? ` (${rs(a)})` : ''}. Right: both to the bottom face${rs(b) ? ` (${rs(b)})` : ''}.`
    : `Middle: Br⁺ on the top face, Br⁻ from the bottom${rs(a) ? ` (${rs(a)})` : ''}. Right: the other way round${rs(b) ? ` (${rs(b)})` : ''}.`}</p>`;
  return s + (state.xguess ? verdictLine(state.xpicks[k], addAnswer(r.kind), `it is ${REL_SAID[addAnswer(r.kind)]}`) : '') + `<p class="facts">${ADD_KIND[r.kind] || esc(r.kind)}</p>`;
}

const SUBST = ['unsubstituted', 'monosubstituted', 'disubstituted', 'trisubstituted', 'tetrasubstituted'];
const alkeneName = a => `C=C to β${a.beta + 1}${a.ez ? ` (${a.ez})` : ''}`;
const E2_WHY = {
  open: (x) => `Turn the C–C bond until this H is anti to ${x}`,
  now: (x) => `Axial, anti to the axial ${x}: ready now`,
  flip: (x) => `Trans to ${x}, both equatorial now: a ring flip puts both axial`,
  cis: (x) => `Cis to ${x}: never anti, so no E2 with this H`,
  chair: () => 'In a ring that is not a chair in this model: no verdict',
};
const e2Why = q => (!q.ring ? 'open' : q.ok === null ? 'chair' : !q.ok ? 'cis' : q.now ? 'now' : 'flip');

function e2HTML(h, r, busy) {
  let s = `<p class="small">E2 is one step: a strong base takes an H from a β carbon (next to the C–${h.el} carbon α) as ${h.el}⁻ leaves. The C–H and C–${h.el} bonds must be anti-periplanar: in one plane and pointing opposite ways (H–C–C–${h.el} 180°), so the electrons of the C–H bond push ${h.el} out from the back. E2 needs a strong base, and it competes with S<sub>N</sub>2.</p>`;
  if (busy || r === undefined) return s + '<p class="small">Building the alkenes…</p>';
  if (!r) return s + '<p class="small">The alkenes could not be built in 3D here.</p>';
  if (!r.hs.length) return s + '<p class="facts">No H on a β carbon (an sp³ carbon next to the C–X carbon), so no E2.</p>';
  const rows = e2Rows(r), name = b => plain(r.betas[b]);
  s += '<div class="rows" role="group" aria-label="Pick a β H">' + rows.map((row, i) => {
    const same = rows.filter(x => x.beta === row.beta), n = row.hs.length, j = same.indexOf(row) + 1;
    const which = n > 1 ? ` · any of its ${n} H` : same.length > 1 ? ` · H ${j} of ${same.length}` : '';
    return `<button data-rh="${i}" aria-pressed="${state.rh === i}"><span>${esc(`H on β${row.beta + 1} (${name(row.beta)})${which}`)}</span><span class="sub">${esc(E2_WHY[e2Why(r.hs[row.hs[0]])](h.el))}</span></button>`;
  }).join('') + '</div>';
  const p = e2Pick(), k = e2Key();
  if (p) {
    const q = p.q, why = e2Why(q), al = q.alkene !== null ? r.alkenes[q.alkene] : null;
    const ci = cur.chairs.findIndex(c => c.atoms.includes(h.c) && c.atoms.includes(r.betas[p.row.beta])), f = cur.flips[ci];
    if (why === 'open') s += `<p class="facts">Turned about the Cα–Cβ bond so this H sits anti to ${h.el} (H–C–C–${h.el} 180°), as in a Newman view down the bond. The groups that sit gauche end up cis on the C=C.</p>`;
    if (why === 'now') s += `<p class="facts">This H and ${h.el} are both axial, on opposite faces: anti-periplanar in this chair.</p>`;
    if (why === 'flip') s += `<p class="facts">This H and ${h.el} are trans, but both equatorial in this chair: a ring flip puts ${h.el} axial, and this H with it.</p>`
      + (flipping === ci ? '<p class="small">Looking for the other chair…</p>' : f ? '<p class="small">Showing the flipped chair. The slide is a straight line; a real ring flip passes through half-chair and twist-boat forms.</p>'
        : f === null ? '<p class="small">Flip not available: the conformer search (OpenChemLib) did not build the other chair, as happens when a group larger than CH₃ or OH would turn axial. The chair is not shown.</p>' : '');
    if (why === 'cis') s += `<p class="facts">This H is cis to ${h.el}, on the same face of the ring. In a chair, cis neighbors are one axial and one equatorial, so they are never anti: no E2 with this H.</p>`;
    if (why === 'chair') s += '<p class="facts">This bond is in a ring that is not a chair in this model, so the page gives no verdict for it.</p>';
    if (state.xguess && !xOpen(k)) {
      const opts = r.alkenes.map((a, i) => a && [String(i), alkeneName(a)]).filter(Boolean);
      if (r.hs.some(x => x.ok === false)) opts.push(['none', 'No alkene: cannot be anti']);
      return s + '<p class="facts">Which alkene does this H give? Pick one.</p><div class="row">' + opts.map(([v, t]) => `<button data-xpick="${v}">${esc(t)}</button>`).join('') + '</div>';
    }
    const answer = q.ok && al ? String(q.alkene) : 'none';
    if (state.xguess) s += verdictLine(state.xpicks[k], answer, answer === 'none' ? 'this H cannot be anti, so it gives no alkene' : `it gives the ${esc(alkeneName(al))}`);
    if (q.ok && al) s += `<p class="facts">Alkene, on the right: <b>${esc(alkeneName(al))}</b>, ${SUBST[al.subst]}${al.ez ? '' : ', no E/Z'}.</p>`;
    if (q.ok && !al) s += '<p class="small">This alkene could not be built in 3D here.</p>';
  }
  const can = r.alkenes.filter(a => a?.possible), z = r.alkenes.filter(a => a?.zaitsev);
  s += '<p class="facts">Alkenes E2 can give here:</p><ul class="angles">' + (can.length ? can.map(a => `<li><span>${esc(alkeneName(a))}, ${SUBST[a.subst]}</span><b>${a.zaitsev ? 'Zaitsev' : ''}</b></li>`).join('') : '<li><span>None: no β H can sit anti</span></li>') + '</ul>';
  if (z.some(a => a.possible)) s += '<p class="small">Zaitsev: the more substituted alkene is usually major with a small base (OH⁻, CH₃O⁻). With a bulky base such as tert-butoxide, the less substituted one is usually major (Hofmann).</p>';
  else if (z.length) s += `<p class="facts">The more substituted alkene (${esc(z.map(a => alkeneName(a)).join(', '))}) cannot form: no H on that carbon can sit anti to ${h.el}. So E2 gives the less substituted alkene here, against Zaitsev's rule.</p>`;
  return s;
}

document.querySelectorAll('#rprin button').forEach(b => b.addEventListener('click', () => {
  const i = +b.dataset.rp;
  state.rp = state.rp === i ? -1 : i; state.rh = -1;
  rxView();
}));
$('btn-xguess').addEventListener('click', e => {
  state.xguess = !state.xguess; state.xpicks = {};
  e.currentTarget.setAttribute('aria-pressed', String(state.xguess));
  draw();
});

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
    cur = { ...hit, ...res, groups: res.groups.map(g => ({ ...info[g.id], ...g })), flips: {}, rxr: {} }; // flips[i]: the other chair of chair i, null if none
    state.rx = state.rdb = 0; state.rh = -1; state.xpicks = {}; rxBusy = null; rxQ = null;
    state.hl = state.atom = state.st = state.bond = state.ring = state.chair = -1; state.tries = {}; state.stries = {}; state.ctries = {}; state.flipped = state.sliding = false;
    state.ar = -1; state.ppick = state.npick = null;
    state.plane = false; planeBusy = false; state.rpick = state.fdl = null; state.fpicks = {};
    state.cmp = null; state.side = false; cmpSeq++; $('rel').textContent = ''; $('cmp-go').disabled = false; $('btn-side').setAttribute('aria-pressed', 'false');
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
    buildChips(); buildKit(); buildStereo(); buildShape(); buildArom(); buildReact();
    viewer.resize(); viewer.zoomTo(); draw(); home = viewer.getView(); // resize: the stage may have changed size since the viewer was made
    if (state.tab === 'reactions') rxView();
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
    if (!cur) return;
    if (left === 'shape' && state.tab !== 'shape' && (state.bond >= 0 || state.flipped)) setBond(-1, true); // a turned model or the other chair is a Shape-tab view
    else if (left === 'stereo' && state.tab !== 'stereo') leaveStereo();
    else if (left === 'reactions' && state.tab !== 'reactions') leaveReactions(); else draw();
    if (state.tab === 'reactions' && left !== 'reactions') rxView();
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
  if (state.mirror) setMirror(true); else if (state.side) setSide(true);
  else if (state.tab === 'reactions' && rxShown()) { rxQ = null; applyRx(false); } else viewer.render();
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
