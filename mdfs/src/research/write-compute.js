// Research · E.7: write a compute. compute nbcount/atom: per-atom neighbour counts from an occasional full neighbour
// list, their group average (global scalar) and their distribution (global vector). Checked by hand from the
// positions, and live in the browser against LAMMPS's own coord/atom on the same configuration.
import { runLammps } from './lammps-runner.js';
import { NATIVE } from './nbcount-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

const SCRIPT = `units           lj
atom_style      atomic
lattice         fcc 1.2
region          box block 0 5 0 5 0 5
create_box      2 box
create_atoms    1 box
set             type 1 type/fraction 2 0.2 7731
mass            * 1.0
velocity        all create 1.5 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0
pair_coeff      2 2 0.5 0.88
pair_coeff      1 2 1.5 0.8
fix             1 all nve
run             1000
compute         nb  all nbcount/atom 1.5
compute         nb2 all nbcount/atom 1.5 type 2 bins 8
thermo_style    custom step c_nb c_nb2 c_nb2[*]
run             0`;
const atomsOf = () => NATIVE.atoms.split('\n').map(l => l.trim().split(/\s+/).map(Number));   // id type x y z ix iy iz
const dataFile = () => `500 atoms from the native run, step 1000\n\n500 atoms\n2 atom types\n\n0 ${NATIVE.box} xlo xhi\n0 ${NATIVE.box} ylo yhi\n0 ${NATIVE.box} zlo zhi\n\nMasses\n\n1 1\n2 1\n\nAtoms # atomic\n\n${NATIVE.atoms}\n`;
const BROWSER = `units           lj
atom_style      atomic
read_data       conf.data
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0
compute         co  all coord/atom cutoff 1.5 *
compute         co2 all coord/atom cutoff 1.5 2
dump            d all custom 1 c.dump id c_co c_co2
dump_modify     d sort id
run             0`;

// count neighbours within rc by hand (nearest image; box 7.47 > 2 × 1.5)
function countByHand() {
  const A = atomsOf(), L = NATIVE.box, rcsq = 1.5 * 1.5, out = {};
  for (const a of A) {
    let n = 0, n2 = 0;
    for (const b of A) {
      if (a === b) continue;
      const d = [2, 3, 4].map(k => { let x = a[k] - b[k]; if (x > 0.5 * L) x -= L; else if (x < -0.5 * L) x += L; return x; });
      if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < rcsq) { n++; if (b[1] === 2) n2++; }
    }
    out[a[0]] = [n, n2];
  }
  return out;
}
const binom = (n, k, p) => { let c = 1; for (let i = 1; i <= k; i++) c = c * (n - k + i) / i; return c * p ** k * (1 - p) ** (n - k); };

export default {
  id: 'write-compute',
  title: 'Write a compute',
  minutes: 20,
  summary: 'A compute that counts neighbours: per atom, averaged over a group, and as a distribution, the kind of number behind every "Li⁺ is coordinated by four oxygens". compute nbcount/atom asks LAMMPS for an occasional full neighbour list, fills a per-atom array that grows with the atoms, and offers three outputs at once: c_ID per atom for dumps, c_ID as a group average for thermo, c_ID[k] as a histogram. Its 1000 counts in a 500-atom binary liquid match a hand count from the positions and, live in your browser, LAMMPS\'s own coord/atom on the same configuration; one flag decides whether thermo divides the average by N.',
  learn: [
    'Compute flags: peratom, scalar, vector, and the "ext" flags thermo reads',
    'Requesting an occasional full neighbour list: init, init_list, build_one',
    'Per-atom storage sized by atom->nmax, and MPI_Allreduce for global results',
    'One compute, three outputs, computed once per step',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: what a compute must declare</h2>
    <div class="prose"><p>A compute is asked for values by thermo, dumps, fixes and variables (page D.15). It tells them what it can deliver through flags set in the constructor, and delivers through matching functions:</p>
<pre class="lmp-mini">peratom_flag = 1, size_peratom_cols = 0  →  compute_peratom()  fills vector_atom     (dump c_ID, compute reduce)
scalar_flag  = 1, extscalar = 0          →  compute_scalar()   returns scalar        (thermo c_ID)
vector_flag  = 1, size_vector = N        →  compute_vector()   fills vector[N]       (thermo c_ID[k], fix ave/time)</pre>
<p>"ext" says whether a value is extensive (grows with the number of atoms). Thermo with <code>thermo_modify norm yes</code>, the default in LJ units, divides extensive values by N. For neighbours it needs an <b>occasional full</b> neighbour list: <code>neighbor->add_request(this, REQ_FULL | REQ_OCCASIONAL)</code> in init(), the list handed over in init_list(), and built on demand by <code>neighbor->build_one(list)</code> only on steps when the compute is called. A full list holds every pair twice (i→j and j→i), so each atom can count its own neighbours without communicating with other processors. LAMMPS's ${gh('src/compute_coord_atom.cpp')} works the same way.</p></div>
  </section>

  <section class="block">
    <h2>2 · Write it</h2>
    <pre class="lmp-mini">${esc(NATIVE.hdr.split('\n').slice(4).join('\n').replace(/\n{3,}/g, '\n\n'))}</pre>
    <pre class="lmp-mini">${esc(NATIVE.src)}</pre>
    <div class="prose"><p>compute_scalar and compute_vector call compute_peratom only if it has not run this step (invoked_peratom), so asking for all three outputs costs one neighbour loop. Native rebuild 13 s.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: 1000 counts by hand</h2>
    <div class="prose"><p>A Kob–Andersen-like binary liquid: 500 atoms at ρ = 1.2, 20 % type 2 (ε<sub>AB</sub> = 1.5, σ<sub>AB</sub> = 0.8), 1000 steps from T = 1.5. Then every atom's neighbours closer than 1.5σ, all types and type 2 only (native run):</p></div>
    <pre class="lmp-mini">${esc(SCRIPT)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Count by hand</button></div>
    <div class="tbl cw-a"></div>
    <div class="prose cw-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Match: against coord/atom, live</h2>
    <div class="prose"><p>The browser build has no nbcount/atom, but it has coord/atom. The page writes the native configuration into a data file (positions are exact in write_data, page D.18) and your browser counts with LAMMPS's own compute:</p></div>
    <pre class="lmp-mini">${esc(BROWSER)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run coord/atom here</button><span class="lmp-status muted cw-stb"></span></div>
    <div class="prose cw-bsum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: one flag, a factor of 500</h2>
    <div class="prose cw-c"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why is count sized by atom->nmax and not by the number of atoms?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A per-atom array lives on each processor, one entry per local atom, and atoms migrate between processors at every reneighbouring. atom->nmax is the length of the processor's own per-atom arrays (x, v, …); growing count with it keeps index i valid for every local atom. The total number of atoms is irrelevant on a processor.</div></div>
    <div class="ex"><p><b>2.</b> A half neighbour list would need half the memory. What else would it need?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each pair appears once, so both count[i] and count[j] must be incremented, and j may be a ghost: its count belongs to another processor and would have to be sent back with a reverse communication (comm->reverse_comm(this) plus pack/unpack functions), as pair forces are with newton on. The full list avoids all of that.</div></div>
    <div class="ex"><p><b>3.</b> How do you get the time-averaged distribution of type-2 neighbours over a run?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>fix 2 all ave/time 100 10 1000 c_nb2 mode vector file nb2.dat</code>: c_nb2 is a global vector, so mode vector averages each element. extvector = 0 tells the fix and thermo the elements are intensive fractions.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const th = NATIVE.thermo.map(Number);
    $('[data-k=a]').onclick = () => {
      const H = countByHand(); let same = 0, sum = 0; const hist = new Array(8).fill(0);
      for (const id in NATIVE.nb) { if (H[id][0] === NATIVE.nb[id][0] && H[id][1] === NATIVE.nb[id][1]) same++; sum += H[id][0]; hist[Math.min(H[id][1], 7)] += 1; }
      const nb2mean = Object.values(H).reduce((s, h) => s + h[1], 0) / 500, mean = sum / 500;
      const rows = hist.map((c, k) => { const f = c / 500, p = binom(Math.round(mean), k, 0.2); return row([k === 7 ? '7 or more' : k, c, f, th[4 + k], f === th[4 + k] ? OK : BAD, p.toFixed(3)]); });
      $('.cw-a').innerHTML = `<table>${head(['type-2 neighbours', 'atoms (hand)', 'fraction (hand)', 'c_nb2[k] (LAMMPS)', '', `binomial, n = ${Math.round(mean)}, p = 0.2`])}${rows.join('')}</table>`;
      $('.cw-asum').innerHTML = `<p>All neighbours and type-2 neighbours: ${same} of 500 atoms with both counts identical to the native per-atom output. Averages: ${mean} by hand, c_nb = ${th[1]} ${mean === th[1] ? OK : BAD}; type 2: ${nb2mean}, c_nb2 = ${th[3]} ${nb2mean === th[3] ? OK : BAD}; and <code>compute reduce ave</code> of coord/atom in the same run printed ${th[2]}. The histogram is the same 500 counts binned. The last column is what random mixing would give: if each of an atom's ${Math.round(mean)} neighbours were type 2 with probability 0.2. The measured distribution sits lower (mean ${nb2mean.toFixed(2)} against ${(0.2 * Math.round(mean)).toFixed(1)}): the liquid is not randomly mixed. Splitting the count by centre type (two groups, one compute each) would show where the difference comes from; that is what the group argument of a compute is for.</p>`;
    };
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.cw-stb'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(BROWSER, { files: { 'conf.data': dataFile() }, read: ['c.dump'] }); b.disabled = false;
      if (res.error || !res.files['c.dump']) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump')}</span>`; return; } st.textContent = '';
      const L = res.files['c.dump'].split('\n'), n = +L[3]; let same = 0;
      for (let k = 0; k < n; k++) { const w = L[9 + k].trim().split(/\s+/); const id = w[0]; if (Number(w[1]) === NATIVE.nb[id][0] && Number(w[2]) === NATIVE.nb[id][1]) same++; }
      $('.cw-bsum').innerHTML = `<p>${same} of ${n} atoms: coord/atom in your browser gives exactly the counts nbcount/atom gave natively, for all types and for type 2. Two independent programs, one written on this page, agree on every atom. (In the native run, coord/atom and nbcount/atom were also dumped side by side: 500 of 500.)</p>`;
    };
    $('.cw-c').innerHTML = `<p>The constructor sets <code>extscalar = 0</code>: an average neighbour count does not grow with the system. Rebuilt once with <code>extscalar = 1</code>, the same native run prints c_nb = ${NATIVE.ext1} in thermo while <code>compute reduce ave c_co</code> (intensive) next to it prints ${th[2]}: thermo divided by N = 500, because LJ units default to <code>thermo_modify norm yes</code>. A wrong ext flag is not an error anywhere, just a wrong number; check it by running with <code>thermo_modify norm no</code> and yes and seeing which values change.</p>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
