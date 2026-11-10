// Web Worker: real LAMMPS (WebAssembly, from the lammps.js package) runs here so the page never freezes.
import createModule from 'lammps.js/wasm';

const out = [];
let lmp = null, banner = '';
const ready = createModule({ print: m => out.push(String(m)), printErr: m => out.push(String(m)) }).then(M => {
  lmp = new M.LAMMPSWeb();
  lmp.start();
  banner = out.find(l => /^LAMMPS \(/.test(l)) || '';
});

onmessage = async e => {
  const { id, script } = e.data;
  try { await ready; } catch (err) { postMessage({ id, lines: [], error: 'LAMMPS could not start: ' + String(err && err.message || err), seconds: 0 }); return; }
  out.length = 0;
  let error = null;
  const t0 = performance.now();
  try {
    lmp.runScript('clear');
    lmp.runScript(script.endsWith('\n') ? script : script + '\n');
  } catch (err) {
    const msg = lmp.getLastErrorMessage ? lmp.getLastErrorMessage() : '';
    error = msg || String(err && err.message || err);
  }
  postMessage({ id, lines: out.slice(), error, banner, seconds: (performance.now() - t0) / 1000 });
};
