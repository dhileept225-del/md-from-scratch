// Runs real LAMMPS (WebAssembly build from the lammps.js package) in the reader's browser.
// The ~3 MB engine is downloaded (into a Web Worker) only when the reader first presses Run.
let worker = null, nextId = 1;
const waiting = new Map();

function getWorker() {
  if (!worker) {
    worker = new Worker(new URL('./lmp-worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = e => { const w = waiting.get(e.data.id); if (w) { waiting.delete(e.data.id); w(e.data); } };
    worker.onerror = e => { for (const [id, w] of waiting) w({ id, lines: [], error: 'LAMMPS could not start in this browser (' + (e.message || 'worker error') + ')', seconds: 0 }); waiting.clear(); worker = null; };
  }
  return worker;
}

// Returns { lines, version, thermo: { cols, rows }, error, seconds }
export function runLammps(script) {
  return new Promise(resolve => {
    const id = nextId++;
    waiting.set(id, res => {
      const lines = res.lines || [];
      const errLine = lines.find(l => /^ERROR/.test(l.trim()));
      resolve({ lines, error: errLine ? errLine.trim() : res.error, seconds: res.seconds, version: res.banner || lines.find(l => /^LAMMPS \(/.test(l)) || '', thermo: parseThermo(lines) });
    });
    getWorker().postMessage({ id, script });
  });
}

export function parseThermo(lines) {
  const i = lines.findIndex(l => /^\s*Step\s/.test(l));
  if (i < 0) return { cols: [], rows: [] };
  const cols = lines[i].trim().split(/\s+/), rows = [];
  for (let k = i + 1; k < lines.length; k++) {
    const w = lines[k].trim().split(/\s+/);
    if (w.length !== cols.length || !/^\d+$/.test(w[0])) break;
    rows.push(w);
  }
  return { cols, rows };
}

// LAMMPS prints thermo numbers with 8 significant digits, like C's %.8g
export const g8 = x => {
  if (!isFinite(x)) return String(x);
  if (x === 0) return '0';
  const ex = Math.floor(Math.log10(Math.abs(Number(x.toPrecision(8)))));
  if (ex < -4 || ex >= 8) {
    let [m, e] = x.toExponential(7).split('e');
    m = m.replace(/\.?0+$/, '');
    const n = Math.abs(+e);
    return `${m}e${+e < 0 ? '-' : '+'}${n < 10 ? '0' + n : n}`;
  }
  return Number(x.toPrecision(8)).toString();
};

// Mark how many leading characters of two printed numbers agree
export function matchMarkup(hand, lmp) {
  let n = 0;
  while (n < hand.length && n < lmp.length && hand[n] === lmp[n]) n++;
  const full = hand === lmp;
  return { full, html: `<span class="m-ok">${hand.slice(0, n)}</span><span class="m-bad">${hand.slice(n)}</span>` };
}
