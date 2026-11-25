// Research · F.1: Python and the library interface. Beyond E.12: a force written in Python and called by LAMMPS every
// step (fix external, pf/callback), byte-identical to the built-in fix spring/self and to your browser's run; the same
// library from 20 lines of C; what each call returns.
import { runLammps } from './lammps-runner.js';
import { NATIVE } from './library-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const BROWSER = `units           lj
atom_style      atomic
atom_modify     map array
lattice         fcc 0.8442
region          box block 0 8 0 8 0 8
region          blob block 2.5 5.5 2.5 5.5 2.5 5.5
create_box      1 box
create_atoms    1 region blob
mass            1 1.0
velocity        all create 1.0 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
fix             2 all spring/self 5.0
fix_modify      2 energy yes
dump            d all custom 50 b.dump id x y z vx vy vz
dump_modify     d format float %.17g sort id
run             200`;
const COMPILE = `mpicc -I lammps/src capi.c -L lammps/build -llammps -o capi     # 0.1 s
LD_LIBRARY_PATH=lammps/build ./capi`;
const f64 = s => new Float64Array(Uint8Array.from(atob(s), c => c.charCodeAt(0)).buffer);

export default {
  id: 'python-interface',
  title: 'Python and the library interface',
  minutes: 15,
  summary: 'Page E.12 read LAMMPS\'s arrays from Python. This page goes the other way: Python computes a force and LAMMPS calls it every step. fix external hands the callback the positions and atom IDs; the callback fills the force array with a harmonic tether written in four lines of numpy. 200 steps of 171 atoms give a dump byte-identical to LAMMPS\'s built-in fix spring/self, and to the same run with spring/self in your browser. Then the same library from 20 lines of C, with every number checked against the lmp program.',
  learn: [
    'fix external pf/callback: LAMMPS calls your function in post_force',
    'Atom IDs vs local indices, and why reference data must be keyed by ID',
    'Reporting energy from a callback so thermo includes it',
    'The C API in one file; what lammps_get_thermo and lammps_extract_compute return',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: calls in both directions</h2>
    <div class="prose"><p>The library interface (${gh('src/library.cpp')}) lets a program call LAMMPS: commands, extract, gather/scatter. ${gh('src/fix_external.cpp', 150)} lets LAMMPS call the program back. With <code>fix ID group external pf/callback N 1</code>, every N steps, in post_force (stage 5 on page C.2), LAMMPS calls a registered function with the step, the number of local atoms, their IDs, their positions and an array fexternal; whatever the function writes there is added to the forces. The Python module wraps both pointers as numpy arrays, so the function is plain numpy.</p>
<p>The local index of an atom changes whenever LAMMPS sorts or migrates atoms; its ID does not. So the reference positions are stored by ID, <code>x0[tag]</code>, never by row.</p></div>
  </section>

  <section class="block">
    <h2>2 · A force written in Python</h2>
    <pre class="lmp-mini">${esc(NATIVE.py)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.pyout)}</pre>
    <div class="prose"><p>Why identical and not just close: fix spring/self computes dx = x − x_ref and does <code>f -= k*dx</code>; the callback writes <code>-(K*dx)</code> and fix external does <code>f += fext</code>. Negation is exact, so f − k·dx and f + (−(k·dx)) are the same double. Unwrapping adds 0·L for atoms that never crossed a boundary (the 171 atoms sit in the middle of a larger box), which is exact too. The energy reported through fix_external_set_energy_global and the spring energy also agree to the last digit, so thermo's PotEng is identical.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: the same run in your browser</h2>
    <div class="prose"><p>Your browser has no Python inside LAMMPS, but it has fix spring/self. The page compares its dump with the dump of the Python run above (positions and velocities, steps 0 to 200):</p></div>
    <pre class="lmp-mini">${esc(BROWSER)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run and compare</button><span class="lmp-status muted pi-sta"></span></div>
    <div class="prose pi-asum"></div>
  </section>

  <section class="block">
    <h2>4 · The same library from C</h2>
    <pre class="lmp-mini">${esc(NATIVE.c)}</pre>
    <pre class="lmp-mini">${esc(COMPILE)}\n\n${esc(NATIVE.cout)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Check the C output</button><span class="lmp-status muted pi-stb"></span></div>
    <div class="prose pi-bsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> In the callback, why x0[tag] and not x0[:nlocal]?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The rows of x are in LAMMPS's memory order, which atom sorting (every 1000 steps and at each setup) and migration between processors change. tag holds each row's atom ID, so x0[tag] lines up the right reference with each row on any step and any number of processors.</div></div>
    <div class="ex"><p><b>2.</b> What does fix_modify 2 energy yes do, and what if it is left out?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>It adds the fix's energy (here the value the callback reports) to thermo's PotEng. Without it the forces are still applied but the energy is missing from PotEng and TotEng, so the total energy appears not to be conserved.</div></div>
    <div class="ex"><p><b>3.</b> The callback costs a Python call per step. When does that matter?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each call has a fixed overhead of microseconds plus whatever numpy does, done on one core per MPI rank. For a few hundred atoms that is nothing (0.016 s against 0.011 s for 200 steps here); for millions of atoms or a cheap step it can dominate. Calling every N steps (pf/callback N) or moving hot code into a C++ fix are the usual answers.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.pi-sta'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(BROWSER, { read: ['b.dump'] }); b.disabled = false;
      if (res.error || !res.files['b.dump']) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump')}</span>`; return; } st.textContent = '';
      const D = f64(NATIVE.dump), L = res.files['b.dump'].split('\n'); let ok = 0, n = 0, f = 0;
      for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const na = +L[i + 3]; for (let a = 0; a < na; a++) { const w = L[i + 9 + a].trim().split(/\s+/).map(Number); for (let c = 0; c < 6; c++) { n++; if (w[1 + c] === D[(f * NATIVE.n + a) * 6 + c]) ok++; } } f++; i += 8 + na; }
      $('.pi-asum').innerHTML = `<p>${ok} of ${n} numbers identical ${ok === n ? OK : BAD}, over ${f} frames: a force computed by numpy inside a native LAMMPS, called back 201 times, and a force computed by LAMMPS's own C++ in WebAssembly in your browser give the same trajectory to the last bit.</p>`;
    };
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.pi-stb'); b.disabled = true; st.textContent = 'Running…';
      const sc = `units lj\natom_style atomic\nlattice fcc 0.8442\nregion box block 0 4 0 4 0 4\ncreate_box 1 box\ncreate_atoms 1 box\nmass 1 1.0\nvelocity all create 1.44 87287 loop geom\npair_style lj/cut 2.5\npair_coeff 1 1 1.0 1.0 2.5\nfix 1 all nve\nthermo_style custom step pe\nthermo_modify format float %.17g\nrun 100\nprint "TOT $(pe*atoms:%.17g)"`;
      const res = await runLammps(sc); b.disabled = false; st.textContent = '';
      const m = NATIVE.cout.match(/thermo\) (\S+), total \(compute\) (\S+)/), per = Number(m[1]), tot = Number(m[2]);
      const lp = Number(res.thermo.rows[res.thermo.rows.length - 1][1]);
      $('.pi-bsum').innerHTML = `<p>PotEng at step 100 in your browser: ${s17(lp)} ${lp === per ? OK : BAD} the C program's lammps_get_thermo("pe"). lammps_extract_compute("thermo_pe") returns the compute's own scalar, the total energy, ${s17(tot)}; thermo divides it by the 256 atoms because LJ units default to thermo_modify norm yes: ${s17(tot / 256)} ${tot / 256 === per ? OK : BAD}. Same library, three front ends (lmp, Python, C) and the browser build: one number.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
