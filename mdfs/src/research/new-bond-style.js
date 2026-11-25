// Research · E.9: a new bond style. bond_style quartic/mine (harmonic plus cubic and quartic terms, a drop-in extension
// of harmonic) written from bond_harmonic; a C-H diatomic reproduced by hand for 600 steps; identical to CLASS2's
// bond class2 and, with K3 = K4 = 0, to bond harmonic run live in the browser; anharmonic period and bond stretching.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const FTM2V = 1.0 / 48.88821291 / 48.88821291, M = [12.011, 1.008], K = 345.0, R0 = 1.101, K3 = -691.892, K4 = 844.5998, DT = 0.1;
const DATA = x2 => `C-H diatomic\n\n2 atoms\n1 bonds\n2 atom types\n1 bond types\n\n-10 10 xlo xhi\n-10 10 ylo yhi\n-10 10 zlo zhi\n\nMasses\n\n1 12.011\n2 1.008\n\nAtoms # bond\n\n1 1 1 0.0 0.0 0.0\n2 1 2 ${x2} 0.0 0.0\n\nBonds\n\n1 1 1 2\n`;
const SCRIPT = (style, coeff) => `units           real
atom_style      bond
read_data       ch.data
pair_style      zero 5.0
pair_coeff      * *
bond_style      ${style}
bond_coeff      1 ${coeff}
fix             1 all nve
timestep        0.1
dump            d all custom 1 x.dump id x vx
dump_modify     d format float %.17g sort id
thermo_style    custom step pe ke etotal
thermo_modify   format float %.17g
thermo          1
run             600`;
const f64 = s => new Float64Array(Uint8Array.from(atob(s), c => c.charCodeAt(0)).buffer);

// BondQuarticMine::compute + FixNVE for the diatomic along x (bond list 1-2)
function handRun(x2, k3, k4, n = 600) {
  const x = [0.0, x2], v = [0.0, 0.0], dtf = 0.5 * DT * FTM2V, out = [];
  const force = () => {
    const delx = x[0] - x[1], rsq = delx * delx + 0.0 * 0.0 + 0.0 * 0.0, r = Math.sqrt(rsq), dr = r - R0, dr2 = dr * dr, dr3 = dr2 * dr, dr4 = dr3 * dr;
    const de = 2.0 * K * dr + 3.0 * k3 * dr2 + 4.0 * k4 * dr3, fb = r > 0.0 ? -de / r : 0.0;
    return { f: [0.0 + delx * fb, 0.0 - delx * fb], e: K * dr2 + k3 * dr3 + k4 * dr4 };
  };
  let { f, e } = force(); out.push({ x: [...x], v: [...v], e });
  for (let s = 1; s <= n; s++) {
    for (let i = 0; i < 2; i++) { const m = dtf / M[i]; v[i] += m * f[i]; x[i] += DT * v[i]; }
    ({ f, e } = force());
    for (let i = 0; i < 2; i++) { const m = dtf / M[i]; v[i] += m * f[i]; }
    out.push({ x: [...x], v: [...v], e });
  }
  return out;
}
// period from downward zero crossings of r − r0, and the mean of r − r0
function analyse(r) {
  const d = Array.from(r, q => q - R0), cr = [];
  for (let k = 1; k < d.length; k++) if (d[k - 1] > 0 && d[k] <= 0) cr.push(k - 1 + d[k - 1] / (d[k - 1] - d[k]));
  const P = cr.slice(1).map((c, j) => (c - cr[j]) * DT);
  return { T: P.reduce((s, p) => s + p, 0) / P.length, mean: d.reduce((s, q) => s + q, 0) / d.length, lo: Math.min(...d), hi: Math.max(...d) };
}

export default {
  id: 'new-bond-style',
  title: 'A new bond style',
  minutes: 20,
  summary: 'Real bonds are not springs: they stiffen when compressed and soften when stretched. bond_style quartic/mine adds cubic and quartic terms to harmonic, E = K Δr² + K₃ Δr³ + K₄ Δr⁴, as an optional extension, so old inputs still work. You write it from bond_harmonic, follow a C–H bond with COMPASS-like constants for 600 steps of 0.1 fs by hand to the last bit, find it identical to CLASS2\'s class2 bond and, with K₃ = K₄ = 0, to bond harmonic run live in your browser; then measure what the anharmonic terms do: the vibration slows with amplitude, as perturbation theory predicts for small amplitudes, and the average bond length grows.',
  learn: [
    'The anatomy of a bond style: bondlist, compute, coeff, restart, data, single',
    'Extending a style without breaking its old input',
    'Arithmetic that makes two styles agree bit for bit (and where it cannot)',
    'Anharmonic vibration: frequency shift and bond stretching',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a bond style</h2>
    <div class="prose"><p>Bond styles are smaller than pair styles. The neighbour code already lists every bond owned by this processor (<code>neighbor->bondlist</code>: atom i₁, atom i₂, type), so compute() is one loop: distance, force, energy, add ±F to the two atoms (${gh('src/MOLECULE/bond_harmonic.cpp', 48)}). The rest is bookkeeping: coeff() reads <code>bond_coeff</code>, equilibrium_distance() serves SHAKE and the special-bonds machinery, write/read_restart and write_data store the coefficients, single() and born_matrix() give energy and derivatives of one bond on request, extract() lets fix adapt reach the arrays.</p>
<p>The new style: <code>bond_coeff N K r0 [K3 K4]</code>. With two numbers it is bond harmonic; with four, E = KΔr² + K₃Δr³ + K₄Δr⁴, the CLASS2 (COMPASS) bond, which LAMMPS has in the CLASS2 package with the order r0 K2 K3 K4.</p></div>
  </section>

  <section class="block">
    <h2>2 · Write it</h2>
    <div class="prose"><p>Copy bond_harmonic.{h,cpp} with the E.2 renames (class BondQuarticMine, style quartic/mine), then:</p></div>
    <pre class="lmp-mini gd bs-diff">loading…</pre>
    <div class="prose"><p>The force is written as 2.0·K·Δr + 3.0·K₃·Δr² + 4.0·K₄·Δr³, in class2's order. With K₃ = K₄ = 0 the two extra terms are exactly 0, and 2.0·K·Δr equals harmonic's 2.0·(K·Δr) because multiplying by 2 is exact; so forces and trajectories agree with bond harmonic to the last bit. The energy K·(Δr·Δr) and harmonic's (K·Δr)·Δr round differently, so energies may differ in the last digit. Native rebuild 11 s; CLASS2's bond was cloned as class2/ref for comparison (page E.3's trick).</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: a C–H bond, 600 steps by hand</h2>
    <div class="prose"><p>C (12.011) and H (1.008), K = 345 kcal/mol/Å², r₀ = 1.101 Å, K₃ = −691.892, K₄ = 844.5998 (COMPASS-like), started at rest 0.1 Å stretched, Δt = 0.1 fs, 60 fs (native run):</p></div>
    <pre class="lmp-mini">${esc(SCRIPT('quartic/mine', `${K} ${R0} ${K3} ${K4}`))}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Integrate by hand</button></div>
    <div class="prose bs-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Match: two other implementations</h2>
    <div class="prose bs-b"></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run bond harmonic here</button><span class="lmp-status muted bs-stc"></span></div>
    <div class="prose bs-csum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: what the extra terms do</h2>
    <div class="prose"><p>Native runs started at Δr = 0.05, 0.1 and 0.2 Å, with and without K₃, K₄. Period from the downward zero crossings of Δr, averaged over the run.</p></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Analyse</button></div>
    <div class="tbl bs-d"></div>
    <div class="prose bs-dsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> The harmonic C–H period by hand?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>E = KΔr² means a spring constant 2K. μ = 12.011 × 1.008 / 13.019 = 0.92998 g/mol; ω² = 2K·ftm2v/μ with ftm2v = 1/48.888² converting kcal/mol/Å/(g/mol) to Å/fs²: ω = 0.5572 fs⁻¹, T = 11.277 fs, about 2960 cm⁻¹, a C–H stretch.</div></div>
    <div class="ex"><p><b>2.</b> Why is the harmonic run's measured period 11.2755 fs, not 11.2769?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Velocity Verlet on a harmonic oscillator oscillates exactly, but at a shifted frequency: sin(ω′Δt/2) = ωΔt/2, i.e. T′ = πΔt / arcsin(ωΔt/2), slightly shorter. With ωΔt/2 = 0.027859 that is 11.2755 fs, 0.0015 fs shorter: exactly what the run measured.</div></div>
    <div class="ex"><p><b>3.</b> Which other places would need the new terms for SHAKE or fix bond/break to work?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>SHAKE uses only equilibrium_distance() (r₀), unchanged. Anything that evaluates a single bond uses single(), updated here; born_matrix() too. The restart and data files now carry four numbers per type, so a data file written by quartic/mine cannot be read by bond_style harmonic.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const h0 = x => x.toFixed(4);
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    import('./bond-native.js').then(({ NATIVE }) => {
      $('.bs-diff').textContent = NATIVE.diff.split('\n').filter(l => !/^(\+\+\+|---)/.test(l)).join('\n');
      const Q = f64(NATIVE['q_0.1']), H = f64(NATIVE['h_0.1']);   // per frame: atom1 x vx, atom2 x vx
      $('[data-k=a]').onclick = () => {
        const hq = handRun(Number('1.201'), K3, K4); let ok = 0, n = 0, pe = 0;
        hq.forEach((h, s) => { for (let i = 0; i < 2; i++) { n += 2; ok += (Q[s * 4 + i * 2] === h.x[i]) + (Q[s * 4 + i * 2 + 1] === h.v[i]); } if (Number(NATIVE.pe_q[s]) === 0.0 + 0.0 + 0.0 + h.e) pe++; });
        const last = hq[600];
        $('.bs-asum').innerHTML = `<p>${ok} of ${n} positions and velocities identical to the last bit, and ${pe} of 601 bond energies (PotEng). After 60 fs the H atom is at x = ${s17(last.x[1])} Å. The hand code is the compute() loop above plus fix nve, with ftm2v = 1/48.88821291² converting kcal/mol/Å and g/mol to Å/fs²; pair_style zero contributes nothing.</p>`;
      };
      const ts = NATIVE.thermoSame, hh = NATIVE.hhPe;
      $('.bs-b').innerHTML = `<p><b>CLASS2's bond class2</b>, cloned as class2/ref, with <code>bond_coeff 1 1.101 345.0 -691.892 844.5998</code>: the same three runs (Δr₀ = 0.05, 0.1, 0.2 Å), every thermo line identical: ${Object.entries(ts).map(([a, v]) => `${v[0]}/${v[1]}`).join(', ')}, and the dumps byte for byte.</p><p><b>bond harmonic</b> with <code>bond_coeff 1 345.0 1.101</code> against quartic/mine with the same two numbers: dumps (x, v every step) identical; PotEng differs in the last digit on ${Object.values(hh).map(v => v[0]).join(', ')} of 601 lines, as predicted. Your browser has bond harmonic, so check the trajectory yourself:</p>`;
      $('[data-k=c]').onclick = async e => {
        const b = e.currentTarget, st = $('.bs-stc'); b.disabled = true; st.textContent = 'Running…';
        const res = await runLammps(SCRIPT('harmonic', `${K} ${R0}`), { files: { 'ch.data': DATA('1.201') }, read: ['x.dump'] }); b.disabled = false;
        if (res.error || !res.files['x.dump']) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump')}</span>`; return; } st.textContent = '';
        const L = res.files['x.dump'].split('\n'); let ok = 0, n = 0, f = 0;
        for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { for (let a = 0; a < 2; a++) { const w = L[i + 9 + a].trim().split(/\s+/).map(Number); n += 2; ok += (w[1] === H[f * 4 + a * 2]) + (w[2] === H[f * 4 + a * 2 + 1]); } f++; i += 10; }
        $('.bs-csum').innerHTML = `<p>${ok} of ${n} positions and velocities from bond harmonic in your browser equal those of quartic/mine with K₃ = K₄ = 0 natively, to the last bit, over all ${f} frames. The extension really is an extension.</p>`;
      };
      $('[data-k=d]').onclick = () => {
        const mu = M[0] * M[1] / (M[0] + M[1]), w0 = Math.sqrt(2 * K * FTM2V / mu), T0 = 2 * Math.PI / w0, Tv = Math.PI * DT / Math.asin(w0 * DT / 2);
        const e2 = 3 * K3 * FTM2V / mu, e3 = 4 * K4 * FTM2V / mu, rows = [];
        for (const a of ['0.05', '0.1', '0.2']) {
          const q = analyse(f64(NATIVE[`r_q_${a}`])), h = analyse(f64(NATIVE[`r_h_${a}`])), A = (q.hi - q.lo) / 2, wp = w0 + (3 * e3 / (8 * w0) - 5 * e2 * e2 / (12 * w0 ** 3)) * A * A;
          rows.push(row([a, h.T.toFixed(4), q.T.toFixed(4), (2 * Math.PI / wp).toFixed(4), `${q.lo.toFixed(3)} … ${q.hi.toFixed(3)}`, q.mean.toFixed(4)]));
        }
        $('.bs-d').innerHTML = `<table>${head(['start Δr (Å)', 'T harmonic (fs)', 'T quartic (fs)', 'T, 2nd-order theory', 'Δr range (quartic)', '⟨Δr⟩ (Å)'])}${rows.join('')}</table>`;
        $('.bs-dsum').innerHTML = `<p>Harmonic: the same period at every amplitude, ${h0(T0)} fs analytically, ${h0(Tv)} fs for Verlet with Δt = 0.1 fs (exercise 2), which is what the runs give. Quartic: the period grows with amplitude, because the negative K₃ softens the stretched side more than K₄ stiffens it. Lindstedt–Poincaré perturbation theory, ω ≈ ω₀ + [3ε₃/(8ω₀) − 5ε₂²/(12ω₀³)]A² with ε₂ = 3K₃/μ, ε₃ = 4K₄/μ and A half the peak-to-peak range, matches the small-amplitude runs and fails at 0.2 Å, where higher orders matter. The bond also spends more time stretched than compressed: ⟨Δr⟩ is positive and grows roughly like A², the microscopic origin of thermal expansion.</p>`;
      };
    });
    return () => {};
  },
};
