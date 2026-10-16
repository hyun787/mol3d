// Name → SMILES. OPSIN parses systematic (IUPAC) names; PubChem covers common and trade names.
const PUBCHEM = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/SMILES,IUPACName,Title/JSON';
const sleep = ms => new Promise(r => setTimeout(r, ms));

export async function resolveName(name, wait = sleep) {
  try {
    const r = await fetch(`https://www.ebi.ac.uk/opsin/ws/${encodeURIComponent(name)}.json`);
    const j = r.ok ? await r.json() : null;
    if (j?.smiles) return { smiles: j.smiles, title: name, source: 'opsin' };
  } catch { /* OPSIN down or name rejected: try PubChem */ }

  // PubChem throttles by request source (about 5 per second). A class behind one school IP can hit that
  // together, so on "busy" wait a random 1–3 s, then 2–6 s, then 3–9 s and retry: the randomness spreads
  // the students out. A network error counts as busy too, in case the throttle reply lacks CORS headers.
  for (let attempt = 1; ; attempt++) {
    let r = null;
    try { r = await fetch(PUBCHEM, { method: 'POST', body: new URLSearchParams({ name }) }); } catch { /* retry below */ }
    const busy = !r || r.status === 503 || r.status === 429;
    if (busy && attempt <= 3) { await wait(attempt * (1000 + Math.random() * 2000)); continue; }
    if (busy) throw new Error('PubChem is busy right now. Wait a few seconds and try again.');
    if (r.status === 404) return null;
    if (!r.ok) throw new Error(`PubChem error (${r.status}).`);
    const p = (await r.json()).PropertyTable.Properties[0];
    return { smiles: p.SMILES, title: p.Title || name, iupac: p.IUPACName, cid: p.CID, source: 'pubchem' };
  }
}
