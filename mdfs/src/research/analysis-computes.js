// Research · D.15: computes for analysis. compute rdf (histogram + normalisation, special-bond exclusion) and compute msd
// (unwrapping, com, average) reproduced bit for bit from dumped positions, and the triclinic "average yes" bug (compute_msd.cpp:233-234).
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

// ---------- dump and fix ave/time readers ----------
function frames(text) {
  const L = text.split('\n'), out = [];
  for (let i = 0; i < L.length; i++) {
    if (!L[i].startsWith('ITEM: TIMESTEP')) continue;
    const step = +L[i + 1], n = +L[i + 3], A = [];
    for (let k = 0; k < n; k++) { const w = L[i + 9 + k].trim().split(/\s+/); A.push({ id: +w[0], x: [+w[1], +w[2], +w[3]], im: w.length > 4 ? [+w[4], +w[5], +w[6]] : [0, 0, 0] }); }
    out.push({ step, box: [0, 1, 2].map(k => L[i + 5 + k].trim().split(/\s+/).map(Number)), A }); i += 8 + n;
  }
  return out;
}
const aveRows = text => text.split('\n').filter(l => l.trim() && !l.startsWith('#')).slice(1).map(l => l.trim().split(/\s+/).map(Number));

// ---------- compute rdf by hand (ComputeRDF::compute_array, one type) ----------
function rdfHand(A, prd, nbin, cut, skip) {
  const N = A.length, delr = cut / nbin, delrinv = 1.0 / delr, hist = new Array(nbin).fill(0);
  for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
    if (i === j || (skip && skip(A[i].id, A[j].id))) continue;
    const d = [0, 1, 2].map(k => { let xj = A[j].x[k]; const dk = A[i].x[k] - xj; if (dk > 0.5 * prd) xj = xj + prd; else if (dk < -0.5 * prd) xj = xj - prd; return A[i].x[k] - xj; });   // j as a ghost: x + prd
    const r = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), b = Math.trunc(r * delrinv);
    if (b < nbin) hist[b] += 1.0;
  }
  const constant = 4.0 * Math.PI / (3.0 * prd * prd * prd), normfac = N - N / N;   // jcount − duplicates/icount
  let nc = 0.0; const g = [], ncs = [], r = [];
  for (let b = 0; b < nbin; b++) {
    const rl = b * delr, ru = (b + 1) * delr, vfrac = constant * (ru * ru * ru - rl * rl * rl);
    let gr = 0.0; if (vfrac * normfac !== 0.0) { gr = hist[b] / (vfrac * normfac * N); nc += gr * vfrac * normfac; }
    g.push(gr); ncs.push(nc); r.push((b + 0.5) * delr);
  }
  return { hist, g, nc: ncs, r };
}

// ---------- compute msd by hand (ComputeMSD constructor + compute_vector), atoms in memory order ----------
function msdHand(F0, fr, H, { com, average, asWritten = true }) {
  const tri = H[3] !== 0 || H[4] !== 0 || H[5] !== 0;
  const un = (x, im) => tri ? [x[0] + H[0] * im[0] + H[5] * im[1] + H[4] * im[2], x[1] + H[1] * im[1] + H[3] * im[2], x[2] + H[2] * im[2]]
    : [x[0] + im[0] * H[0], x[1] + im[1] * H[1], x[2] + im[2] * H[2]];
  const xcm = A => { const c = [0.0, 0.0, 0.0]; for (const a of A) { const u = un(a.x, a.im); for (let k = 0; k < 3; k++) c[k] += u[k] * 1.0; } return c.map(v => v / A.length); };
  const orig = {}, A0 = [...F0.A].sort((a, b) => a.id - b.id);           // at definition: creation order = id order
  for (const a of A0) orig[a.id] = un(a.x, a.im);
  if (com) { const c = xcm(A0); for (const a of A0) for (let k = 0; k < 3; k++) orig[a.id][k] -= c[k]; }
  let nav = 0; const out = [];
  for (const f of fr) {
    const cm = com ? xcm(f.A) : [0.0, 0.0, 0.0]; let navfac = 0;
    if (average) { nav++; navfac = 1.0 / (nav + 1); }
    const m = [0.0, 0.0, 0.0, 0.0];
    for (const a of f.A) {
      const u = un(a.x, a.im), t = [u[0] - cm[0], u[1] - cm[1], u[2] - cm[2]], o = orig[a.id];
      if (average) {
        if (tri && asWritten) { o[0] = (o[0] * nav + t[0]) * navfac; o[1] = (o[0] * nav + t[0]) * navfac; o[2] = (o[0] * nav + t[0]) * navfac; }  // lines 232-234 as written
        else { o[0] = (o[0] * nav + t[0]) * navfac; o[1] = (o[1] * nav + t[1]) * navfac; o[2] = (o[2] * nav + t[2]) * navfac; }
      }
      const dx = t[0] - o[0], dy = t[1] - o[1], dz = t[2] - o[2];
      m[0] += dx * dx; m[1] += dy * dy; m[2] += dz * dz; m[3] += dx * dx + dy * dy + dz * dz;
    }
    out.push({ step: f.step, v: m.map(s => s / f.A.length) });
  }
  return out;
}

// ---------- scripts ----------
const LIQ = `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 5 0 5 0 5
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
run             500
compute         r all rdf 50
fix             2 all ave/time 1 1 1 c_r[*] file rdf.out mode vector format " %.17g"
write_dump      all custom x.dump id x y z modify format float %.17g
run             0`;

function dimerData() {
  const n = 8, L = Math.cbrt(512 / 0.8), a = L / n, at = [], bo = [];
  let id = 0;
  for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { id++; at.push(`${id} ${(id + 1) >> 1} 1 ${(i + 0.5) * a} ${(j + 0.5) * a} ${(k + 0.5) * a}`); if (id % 2 === 0) bo.push(`${bo.length + 1} 1 ${id - 1} ${id}`); }
  return `dimers: 256 two-atom molecules\n\n512 atoms\n256 bonds\n1 atom types\n1 bond types\n\n0 ${L} xlo xhi\n0 ${L} ylo yhi\n0 ${L} zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms # bond\n\n${at.join('\n')}\n\nBonds\n\n${bo.join('\n')}\n`;
}
const DIM = `units           lj
atom_style      bond
read_data       dim.data
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
bond_style      harmonic
bond_coeff      1 100.0 1.0
velocity        all create 1.0 4928459 loop geom
fix             1 all nve
run             1000
write_dump      all custom x.dump id mol x y z modify format float %.17g
compute         r all rdf 50
fix             2 all ave/time 1 1 1 c_r[*] file rdf0.out mode vector format " %.17g"
run             0                  # special_bonds lj 0.0 0.0 0.0 (default)
unfix           2
special_bonds   lj 1.0e-50 0.0 0.0
fix             2 all ave/time 1 1 1 c_r[*] file rdf1.out mode vector format " %.17g"
run             0`;

const MSD = ({ tri, com, average }) => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box ${tri ? 'prism 0 5 0 5 0 5 1 0 0' : 'block 0 5 0 5 0 5'}
create_box      1 box
print           "HBOX $(lx:%.17g) $(ly:%.17g) $(lz:%.17g) $(yz:%.17g) $(xz:%.17g) $(xy:%.17g)"
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
compute         m all msd com ${com ? 'yes' : 'no'} average ${average ? 'yes' : 'no'}
write_dump      all custom x0.dump id x y z ix iy iz modify format float %.17g
thermo_style    custom step c_m[1] c_m[2] c_m[3] c_m[4]
thermo_modify   format float %.17g
thermo          50
dump            d all custom 50 x.dump id x y z ix iy iz
dump_modify     d format float %.17g
run             500`;

async function runMsd(o) {
  const res = await runLammps(MSD(o), { read: ['x0.dump', 'x.dump'] });
  if (res.error || !res.files['x.dump']) return { error: res.error || 'no dump' };
  const H = res.lines.find(l => l.startsWith('HBOX')).split(/\s+/).slice(1).map(Number);
  const th = res.thermo.rows.map(r => ({ step: +r[0], v: r.slice(1, 5).map(Number) }));
  return { H, th, F0: frames(res.files['x0.dump'])[0], fr: frames(res.files['x.dump']) };
}
const sameRows = (h, th) => h.reduce((n, r, k) => n + (th[k] && r.step === th[k].step && r.v.every((v, c) => v === th[k].v[c]) ? 1 : 0), 0);

export default {
  id: 'analysis-computes',
  title: 'Computes for analysis',
  minutes: 20,
  summary: 'A compute does nothing until something asks for its value. You reproduce two of the most used ones from dumped positions, bit for bit: compute rdf (a histogram from an occasional neighbour list, divided by an ideal-gas shell with N − 1 partners) and compute msd (unwrapped with image flags, measured from positions stored when the compute was defined, optionally about the centre of mass or a running average). Along the way: why the default special_bonds makes the bond peak vanish from a molecular g(r), and a bug in compute msd with average yes in a triclinic box that turns 0.11 σ² into 24 σ².',
  learn: [
    'When a compute runs, and what scalar, vector, array and per-atom outputs are',
    'compute rdf: binning, the N − 1 normalisation, the running coordination number',
    'Excluded special pairs are missing from g(r), and how to bring them back',
    'compute msd: image flags, the reference positions, com and average',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a compute answers when asked</h2>
    <div class="prose"><p>A <code>compute</code> line only creates an object. It calculates when thermo, a fix, a dump or a variable references it (<code>c_ID</code>, <code>c_ID[2]</code>, <code>c_ID[*]</code>), at most once per step: the caller checks the <code>invoked_*</code> step stamp first. The result is a scalar, a vector, an array or one value per atom.</p>
<p><b>compute rdf</b> (${gh('src/compute_rdf.cpp', 277)}) builds an <i>occasional</i> neighbour list, loops over pairs and adds 1 to bin int(r/Δr), Δr = r_c/nbin. It then divides by what an ideal gas of the same density would put in each shell:</p>
<pre class="lmp-mini">if (factor_lj == 0.0 &amp;&amp; factor_coul == 0.0) continue;        // line 347: excluded special pairs
vfrac   = 4π/(3V) · (r_upper³ − r_lower³);                      // fraction of the box in the shell
normfac = jcount − duplicates/icount;                          // one type: N − N/N = N − 1
g(r)    = hist / (vfrac · normfac · icount);
ncoord += g · vfrac · normfac;                                 // running coordination number</pre>
<p><b>compute msd</b> (${gh('src/compute_msd.cpp', 142)}) stores every atom's unwrapped position when the compute is <i>defined</i> (constructor, line 84, kept by an internal fix STORE that travels with the atoms), then at each call:</p>
<pre class="lmp-mini">x_u = x + ix·Lx (+ xy·iy + xz·iz in a triclinic box) − x_cm   // x_cm only with com yes
if (average) x_ref = (x_ref·n + x_u)/(n + 1);                  // running mean, n = number of calls
msd = Σ (x_u − x_ref)² / N                                     // per component, then the sum</pre></div>
  </section>

  <section class="block">
    <h2>2 · Match: g(r) of an LJ liquid</h2>
    <div class="prose"><p>500 atoms at ρ = 0.8442, 500 steps after starting the fcc lattice at T = 1.44, 50 bins to the pair cutoff 2.5σ. The positions at the moment of the compute are written with 17 digits; your browser counts the pairs (nearest image) and normalises.</p></div>
    <pre class="lmp-mini">${esc(LIQ)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run and count by hand</button><span class="lmp-status muted ac-sta"></span></div>
    <div class="tbl ac-a"></div>
    <div class="prose ac-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: the bond that g(r) cannot see</h2>
    <div class="prose"><p>256 dimers (bond 1.0σ, harmonic), same density. With the default <code>special_bonds lj 0.0 0.0 0.0</code> a bonded pair has both weights 0, so it is left out of the neighbour list and skipped at line 347. Setting the 1-2 weight to 10⁻⁵⁰ puts it back (its LJ energy is then 10⁻⁵⁰ of a normal pair: nothing).</p></div>
    <pre class="lmp-mini">${esc(DIM)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run both settings</button><span class="lmp-status muted ac-stb"></span></div>
    <div class="tbl ac-b"></div>
    <div class="prose ac-bsum"></div>
  </section>

  <section class="block">
    <h2>4 · Match: mean-squared displacement, four ways</h2>
    <div class="prose"><p>The same liquid from the lattice, <code>compute msd</code> defined before the run, thermo and a dump (x, y, z and image flags, unsorted so the atoms come in memory order) every 50 steps. The hand code needs the positions at definition too (<code>write_dump</code> right after the compute): the reference is taken then, not at step 0.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run com no/yes × average no/yes</button><span class="lmp-status muted ac-stc"></span></div>
    <div class="tbl ac-c"></div>
    <div class="prose ac-csum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: average yes in a triclinic box</h2>
    <div class="prose"><p>The triclinic branch of the running average, ${gh('src/compute_msd.cpp', 232)}:</p>
<pre class="lmp-mini">xoriginal[i][0] = (xoriginal[i][0] * naverage + xtmp) * navfac;
xoriginal[i][1] = (xoriginal[i][0] * naverage + xtmp) * navfac;   // [0] and xtmp: should be [1] and ytmp
xoriginal[i][2] = (xoriginal[i][0] * naverage + xtmp) * navfac;   // should be [2] and ztmp</pre>
<p>The orthogonal branch a few lines up is correct. Predict: the y and z references become the running mean of <i>x</i>, so their MSD measures the distance between an atom's y (or z) and its own x, several σ. The box: the same fcc liquid with tilt xy = 1 lattice unit (an exact lattice vector, so the crystal fits).</p></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Run triclinic, average no and yes</button><span class="lmp-status muted ac-std"></span></div>
    <div class="tbl ac-d"></div>
    <div class="prose ac-dsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why N − 1 and not N in the normalisation?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each of the N atoms has N − 1 partners, and the ideal-gas reference spreads those N − 1 over the box. With N, g(r) would sit at (N − 1)/N below 1 at long range: 0.998 for 500 atoms, 0.9 for 10.</div></div>
    <div class="ex"><p><b>2.</b> g(r) beyond the pair cutoff?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>compute rdf 100 cutoff 4.0</code>, plus <code>comm_modify cutoff 4.3</code> (cutoff + skin) so that ghost atoms reach that far; init() stops with an error otherwise. Keep it below half the box length, or one pair is counted at several images.</div></div>
    <div class="ex"><p><b>3.</b> In section 3, how much did the coordination number at 1.3σ grow when the bonds were put back, and why exactly that?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>By 1 (to round-off): every atom has exactly one bonded partner, and all bonds are shorter than 1.3σ. The ncoord column is Σ hist / N, so 512 extra counts (each bond twice) over 512 atoms.</div></div>
    <div class="ex"><p><b>4.</b> What is <code>average yes</code> for, if not diffusion?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>In a solid the running mean converges to the lattice site, so the MSD becomes the mean-square vibration amplitude ⟨u²⟩ (Lindemann ratio, Debye–Waller factor) even if the crystal drifts or the sites are not where the atoms started. In a liquid it has no clear meaning.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.ac-sta'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(LIQ, { read: ['rdf.out', 'x.dump'] }); b.disabled = false;
      if (res.error || !res.files['rdf.out']) return err(st, res.error); st.textContent = '';
      const f = frames(res.files['x.dump'])[0], prd = f.box[0][1] - f.box[0][0], R = aveRows(res.files['rdf.out']), H = rdfHand(f.A, prd, 50, 2.5);
      let ok = 0; R.forEach((r, k) => { if (r[1] === H.r[k] && r[2] === H.g[k] && r[3] === H.nc[k]) ok++; });
      const pk = H.g.indexOf(Math.max(...H.g)); let mn = pk; while (mn + 1 < 50 && H.g[mn + 1] < H.g[mn]) mn++;
      const show = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 35, 40, 45, 49];
      $('.ac-a').innerHTML = `<table>${head(['r (bin centre)', 'pairs ×2', 'g(r) LAMMPS', 'ncoord LAMMPS', 'hand'])}${show.map(k => row([H.r[k].toFixed(3), H.hist[k], `<span class="mono">${s17(R[k][2])}</span>`, `<span class="mono">${s17(R[k][3])}</span>`, R[k][2] === H.g[k] && R[k][3] === H.nc[k] ? OK : BAD])).join('')}</table>`;
      const pairs = H.hist.reduce((s, v) => s + v, 0) / 2;
      $('.ac-asum').innerHTML = `<p>${ok} of 50 rows identical to the last bit (bin centre, g and ncoord). ${pairs} pairs closer than 2.5σ, each counted from both ends. Box ${s17(prd)}σ, so 4π/(3V) = ${s17(4 * Math.PI / (3 * prd ** 3))} and normfac = 499. The first peak is at r = ${H.r[pk].toFixed(3)}σ (g = ${H.g[pk].toFixed(3)}); the first minimum is at ${H.r[mn].toFixed(3)}σ, where ncoord = ${H.nc[mn].toFixed(2)}: each atom has about ${Math.round(H.nc[mn])} nearest neighbours, as in a dense simple liquid (12 in the fcc crystal it melted from).</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.ac-stb'); b.disabled = true; st.textContent = 'Running 1000 steps…';
      const res = await runLammps(DIM, { files: { 'dim.data': dimerData() }, read: ['rdf0.out', 'rdf1.out', 'x.dump'] }); b.disabled = false;
      if (res.error || !res.files['rdf1.out']) return err(st, res.error); st.textContent = '';
      const L = res.files['x.dump'].split('\n'), n = +L[3], prd = (() => { const w = L[5].trim().split(/\s+/).map(Number); return w[1] - w[0]; })();
      const A = []; for (let k = 0; k < n; k++) { const w = L[9 + k].trim().split(/\s+/).map(Number); A.push({ id: w[0], mol: w[1], x: [w[2], w[3], w[4]] }); }
      const mol = {}; for (const a of A) mol[a.id] = a.mol;
      const R0 = aveRows(res.files['rdf0.out']), R1 = aveRows(res.files['rdf1.out']);
      const H0 = rdfHand(A, prd, 50, 2.5, (i, j) => mol[i] === mol[j]), H1 = rdfHand(A, prd, 50, 2.5);
      let ok = 0; for (let k = 0; k < 50; k++) { if (R0[k][2] === H0.g[k] && R0[k][3] === H0.nc[k]) ok++; if (R1[k][2] === H1.g[k] && R1[k][3] === H1.nc[k]) ok++; }
      const show = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 30, 40];
      $('.ac-b').innerHTML = `<table>${head(['r', 'g, default', 'g, 1-2 = 10⁻⁵⁰', 'ncoord, default', 'ncoord, 1-2 = 10⁻⁵⁰', 'difference'])}${show.map(k => row([H0.r[k].toFixed(3), R0[k][2].toFixed(4), R1[k][2].toFixed(4), R0[k][3].toFixed(4), R1[k][3].toFixed(4), (R1[k][3] - R0[k][3]).toFixed(6)])).join('')}</table>`;
      const intra = H1.hist.map((v, k) => v - H0.hist[k]), lo = intra.findIndex(v => v > 0), hi = 49 - [...intra].reverse().findIndex(v => v > 0);
      $('.ac-bsum').innerHTML = `<p>Hand code, once skipping pairs in the same molecule and once counting all: ${ok} of 100 rows identical to LAMMPS. The bonded pairs fall in the bins from ${H0.r[lo].toFixed(3)} to ${H0.r[hi].toFixed(3)}σ, right under the first intermolecular peak, so with the default they are simply missing and nothing in the curve tells you. Past them the coordination number is higher by exactly ${(R1[30][3] - R0[30][3]).toFixed(6)}: one bonded partner per atom.</p><p>For a molecular liquid (the anion–cation g(r) of an electrolyte, say) the default is usually what you want, intermolecular structure only. For the intramolecular part, either the 10⁻⁵⁰ trick (it changes the neighbour lists, so do it in a separate <code>rerun</code> or after the production run, as here), or a separate analysis of the bond lengths.</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.ac-stc'); b.disabled = true; const rows = []; let tot = 0, all = 0;
      for (const com of [false, true]) for (const average of [false, true]) {
        st.textContent = `Running com ${com ? 'yes' : 'no'}, average ${average ? 'yes' : 'no'}…`;
        const r = await runMsd({ tri: false, com, average }); if (r.error) { b.disabled = false; return err(st, r.error); }
        const h = msdHand(r.F0, r.fr, r.H, { com, average }), ok = sameRows(h, r.th); tot += ok; all += r.th.length;
        const last = r.th[r.th.length - 1];
        rows.push(row([`com ${com ? 'yes' : 'no'}, average ${average ? 'yes' : 'no'}`, `<span class="mono">${s17(last.v[3])}</span>`, `<span class="mono">${s17(r.th[1].v[3])}</span>`, `${ok === r.th.length ? OK : BAD} ${ok} of ${r.th.length}`]));
      }
      b.disabled = false; st.textContent = '';
      $('.ac-c').innerHTML = `<table>${head(['compute msd', 'MSD at step 500 (σ²)', 'at step 50', 'thermo rows × 4 values'])}${rows.join('')}</table>`;
      $('.ac-csum').innerHTML = `<p>${tot} of ${all} thermo rows identical in all four columns. Three things the hand code had to get right: the atoms are summed in memory order, which the sort at setup changed, while each atom keeps its own reference (keyed by ID here, carried by fix STORE in LAMMPS); the reference comes from the positions at the <code>compute</code> line; and the centre of mass is summed over unwrapped positions (Group::xcm). com yes changes only the last digits here, because <code>loop geom</code> velocities have zero total momentum and the centre of mass does not move; with a drifting system it removes the drift. average yes gives a smaller number because the reference follows the atom.</p>`;
    };

    $('[data-k=d]').onclick = async e => {
      const b = e.currentTarget, st = $('.ac-std'); b.disabled = true; const rows = []; let sum = '';
      for (const average of [false, true]) {
        st.textContent = `Running average ${average ? 'yes' : 'no'}…`;
        const r = await runMsd({ tri: true, com: false, average }); if (r.error) { b.disabled = false; return err(st, r.error); }
        const last = r.th[r.th.length - 1], hw = msdHand(r.F0, r.fr, r.H, { com: false, average }), okw = sameRows(hw, r.th);
        rows.push(row([`average ${average ? 'yes' : 'no'}, LAMMPS`, ...last.v.map(v => v.toFixed(4)), '']));
        rows.push(row([average ? 'hand, lines 232–234 as written' : 'hand', ...hw[hw.length - 1].v.map(v => v.toFixed(4)), `${okw === r.th.length ? OK : BAD} ${okw} of ${r.th.length} rows bitwise`]));
        if (average) {
          const hc = msdHand(r.F0, r.fr, r.H, { com: false, average, asWritten: false });
          rows.push(row(['hand, corrected', ...hc[hc.length - 1].v.map(v => v.toFixed(4)), '']));
          sum = `<p>The hand code with the three lines as written reproduces LAMMPS to the last bit; with [1]/ytmp and [2]/ztmp it gives ${hc[hc.length - 1].v[3].toFixed(3)} σ², below the average-no value, as in the cubic box (section 4). LAMMPS reports ${last.v[3].toFixed(2)} σ²: the x component is right, y and z are not displacements at all. Box ${r.H[0].toFixed(3)}σ with tilt ${r.H[5].toFixed(3)}σ; an atom's y and its x are unrelated numbers between 0 and about ${(r.H[0] + r.H[5]).toFixed(1)}σ, so (y − mean x)² is several σ², and that is what ends up in c_m[2] and c_m[3].</p><p>The branch is in the current develop code as well, and I found no issue about it: a candidate for a one-line bug report. Until it is fixed, compute msd with average yes in a triclinic box (any box with a tilt factor, even 0 after <code>change_box … triclinic</code>) gives wrong y, z and total; c_m[1] is still correct.</p>`;
        }
      }
      b.disabled = false; st.textContent = '';
      $('.ac-d').innerHTML = `<table>${head(['', 'x²', 'y²', 'z²', 'total', ''])}${rows.join('')}</table>`;
      $('.ac-dsum').innerHTML = sum;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
