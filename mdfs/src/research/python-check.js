// Research · E.12: check from Python. LAMMPS built as a shared library and driven from Python: numpy views of x and f,
// an independent numpy recomputation of lj/mine's forces and energy, and the derivative test (force = −dE/dx by central
// differences) that catches a force/energy mismatch no reference file would. Every finite difference rechecked here.
import { NATIVE } from './python-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const BUILD = `cd lammps/build
cmake -D BUILD_SHARED_LIBS=on .     # liblammps.so instead of liblammps.a: every file recompiles once (${NATIVE.shared[0]} s here)
make -j2
export LD_LIBRARY_PATH=$PWD                 # where liblammps.so is
export PYTHONPATH=$PWD/../python            # the lammps Python package (pure Python, ctypes)
python3 check.py`;
// parse the finite-difference lines and the analytic force
function parseFD(text) {
  const fk = Number(text.match(/analytic f\[7\]\[0\] = (\S+)/)[1]);
  const rows = [...text.matchAll(/h (\S+)\s+E\+ (\S+)\s+E- (\S+)\s+fd (\S+)\s+err (\S+)/g)].map(m => ({ h: Number(m[1]), ep: Number(m[2]), em: Number(m[3]), fd: Number(m[4]), err: Number(m[5]) }));
  return { fk, rows };
}

export default {
  id: 'python-check',
  title: 'Check from Python',
  minutes: 18,
  summary: 'Build LAMMPS as a shared library and every array inside it becomes a numpy array in Python, without copying. That turns two checks into a few lines: recompute the forces of pair lj/mine (page E.3) independently in numpy, and test that the force really is minus the derivative of the energy, by moving one atom by ±h and differencing the energies LAMMPS returns. The error falls as h² down to 10⁻⁸, then rounding takes over, exactly as numerical analysis predicts. Then a planted bug, a wrong coefficient in the force only, which no YAML reference made with the same code would catch: the derivative test flags it at every h.',
  learn: [
    'The library interface: liblammps.so, the C API, the lammps Python module',
    'numpy views of x and f, and writing into LAMMPS\'s own arrays',
    'The derivative test: O(h²) truncation, O(ε/h) rounding, the best h',
    'A check that needs no second implementation',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: LAMMPS as a library</h2>
    <div class="prose"><p>The lmp program is a thin main() around a library: create a LAMMPS object, feed it commands, destroy it. ${gh('src/library.h')} exposes that as C functions (lammps_open, lammps_command, lammps_extract_atom, lammps_extract_compute, …), and the <code>lammps</code> Python package calls them with ctypes. <code>lmp.numpy.extract_atom("x")</code> returns a numpy array that points at LAMMPS's own x: reading it costs nothing, writing to it moves atoms. It needs the library as a shared object, which is one more full rebuild:</p></div>
    <pre class="lmp-mini">${esc(BUILD)}</pre>
    <div class="prose"><p>After that, changing a style costs a recompile of that file and a relink of liblammps.so (2 s here); the Python script picks up the new library at its next start.</p></div>
  </section>

  <section class="block">
    <h2>2 · The script</h2>
    <pre class="lmp-mini">${esc(NATIVE.script)}</pre>
    <div class="prose"><p>Two details from earlier pages: <code>atom_modify sort 0 0.0</code>, because each <code>run 0</code> would otherwise re-sort the atoms in memory (pages C.6 and D.18) and row 7 of x would be a different atom after the first run; and <code>post no</code> to keep the output short. 108 atoms of a disordered fcc crystal, Mie 14–6.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: numpy against LAMMPS</h2>
    <pre class="lmp-mini">${esc(NATIVE.good.split('\n').slice(0, 4).join('\n'))}</pre>
    <div class="prose"><p>The numpy version is a different algorithm (all 108 × 107 ordered pairs, nearest image by rounding, forces summed per row) using the same formula: forces agree to 2 × 10⁻¹³ on values up to 565 (4 × 10⁻¹⁶ relative), the energy to 1.4 × 10⁻¹⁵. Not bit for bit, because the sums run in a different order; close enough that any real error would stand out by ten orders of magnitude.</p></div>
  </section>

  <section class="block">
    <h2>4 · Observe: is the force the derivative of the energy?</h2>
    <div class="prose"><p>Atom 7 is moved by +h and −h along x; after each <code>run 0</code> LAMMPS returns the total energy, and the central difference −(E₊ − E₋)/2h should approach the force LAMMPS computed, −1.6212.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Recheck the differences</button></div>
    <div class="tbl py-a"></div>
    <div class="prose py-asum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: a bug only this test catches</h2>
    <div class="prose"><p>One coefficient wrong in the force, the energy untouched:</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.bugdiff.split('\n').filter(l => !/^(\+\+\+|---)/.test(l)).join('\n'))}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Show the same test with the bug</button></div>
    <div class="tbl py-b"></div>
    <div class="prose py-bsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why does the error grow again below h = 10⁻⁵?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>E₊ and E₋ are each correct to about one unit in the last place of 280, i.e. ±3 × 10⁻¹⁴. Their difference is divided by 2h, so the rounding error of the derivative is about 3 × 10⁻¹⁴ / h: 3 × 10⁻⁷ at h = 10⁻⁷. The truncation error is (h²/6)·|E‴|; with |E‴|/6 ≈ 258 from the table, the sum of the two is smallest near h = (3 × 10⁻¹⁴ / (2 × 258))^(1/3) ≈ 4 × 10⁻⁶, between the table's two best rows.</div></div>
    <div class="ex"><p><b>2.</b> Would E.11's unit test have caught the bug of section 5?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only with a reference from an independent implementation, as E.11 used. A reference generated by the buggy code itself (-g) stores the wrong forces and passes forever. The derivative test needs no second implementation: it checks the code against itself, energy against force.</div></div>
    <div class="ex"><p><b>3.</b> What else would you do from Python?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Anything that needs a loop around LAMMPS: scan a parameter (variable or pair_coeff) and fit, drive umbrella windows or replica exchange by hand, compute an analysis on x, v, f every N steps without dump files, or couple LAMMPS with another code. The same calls exist in C, C++ and Fortran.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const table = (text, sel, sumSel, buggy) => {
      const { fk, rows } = parseFD(text); let ok = 0;
      const tr = rows.map(r => { const fd = -(r.ep - r.em) / (2 * r.h), err = Math.abs(fd - fk), same = fd === r.fd && Math.abs(err - r.err) <= 5e-4 * err; ok += same; return row([r.h.toExponential(0), `<span class="mono">${s17(r.ep)}</span>`, `<span class="mono">${s17(r.em)}</span>`, `<span class="mono">${s17(fd)}</span>`, err.toExponential(3), (err / r.h ** 2).toExponential(2), same ? OK : BAD]); });
      $(sel).innerHTML = `<table>${head(['h', 'E₊', 'E₋', '−(E₊ − E₋)/2h', '|error|', 'error / h²', 'recomputed'])}${tr.join('')}</table>`;
      return { fk, rows, ok };
    };
    $('[data-k=a]').onclick = () => {
      const { fk, rows, ok } = table(NATIVE.good, '.py-a');
      $('.py-asum').innerHTML = `<p>${ok} of ${rows.length} differences recomputed here from the printed energies, identical to Python's. From h = 10⁻² to 10⁻⁴ the error falls 100-fold per decade and error/h² is constant (${(rows[2].err / rows[2].h ** 2).toFixed(3)}): the central difference is second order, and that constant is |E‴|/6 along this coordinate. At h = 10⁻⁵ the error is ${rows[4].err.toExponential(1)}, 10⁻⁸ of the force; below it, rounding in E₊ − E₋ (two numbers near −280 that agree in their first eight or nine digits) takes over and the error grows like 1/h. At h = 0.1 the step is too large for the expansion. The force LAMMPS computes (${s17(fk)}) is minus the derivative of the energy it computes, to as many digits as this test can resolve.</p>`;
    };
    $('[data-k=b]').onclick = () => {
      const { fk, rows } = table(NATIVE.bug, '.py-b');
      const head4 = NATIVE.bug.split('\n')[1];
      $('.py-bsum').innerHTML = `<p>The energies are the same as before (the bug does not touch them), the force is now ${s17(fk)}, and the error stalls at ${rows[5].err.toFixed(3)} for every h from 10⁻³ down: it is not converging to 0, so force and energy belong to different potentials. The numpy comparison sees it too (<code>${esc(head4)}</code>), but only because numpy was written independently; the derivative test needed nothing but LAMMPS. Unit tests check that a style still does what it did; this checks that it is self-consistent. Restored and relinked afterwards; the script then printed exactly the output of section 3 again.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
