// Web Worker: real LAMMPS (WebAssembly, from the lammps.js package) runs here so the page never freezes.
import createModule from 'lammps.js/wasm';

const out = [];
let lmp = null, banner = '', mod = null;
const ready = createModule({ print: m => out.push(String(m)), printErr: m => out.push(String(m)) }).then(M => {
  mod = M;
  lmp = new M.LAMMPSWeb();
  lmp.start();
  banner = out.find(l => /^LAMMPS \(/.test(l)) || '';
});

onmessage = async e => {
  const { id, script, files = {}, read = [] } = e.data;
  try { await ready; } catch (err) { postMessage({ id, lines: [], error: 'LAMMPS could not start: ' + String(err && err.message || err), seconds: 0 }); return; }
  for (const [name, text] of Object.entries(files)) mod.FS.writeFile(name, text);
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
  const got = {};
  for (const name of read) { try { got[name] = mod.FS.readFile(name, { encoding: 'utf8' }); mod.FS.unlink(name); } catch { got[name] = null; } }
  postMessage({ id, lines: out.slice(), error, banner, files: got, seconds: (performance.now() - t0) / 1000 });
};
