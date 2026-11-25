// Research · D.1: units and constants. Update::set_units read line by line, every constant re-derived from today's exact SI values, and five of them recovered bit for bit from a live run.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

// update.cpp, 22 Jul 2025 (first line of each block)
const LMP = {
  real: { line: 170, boltz: 0.0019872067, hplanck: 95.306976368, mvv2e: 48.88821291 * 48.88821291, ftm2v: 1.0 / 48.88821291 / 48.88821291, nktv2p: 68568.415, qqr2e: 332.06371 },
  metal: { line: 191, boltz: 8.617343e-5, hplanck: 4.135667403e-3, mvv2e: 1.0364269e-4, ftm2v: 1.0 / 1.0364269e-4, nktv2p: 1.6021765e6, qqr2e: 14.399645 },
  si: { line: 212, boltz: 1.3806504e-23, hplanck: 6.62606896e-34, mvv2e: 1, ftm2v: 1, nktv2p: 1, qqr2e: 8.9876e9 },
  cgs: { line: 233, boltz: 1.3806504e-16, hplanck: 6.62606896e-27, mvv2e: 1, ftm2v: 1, nktv2p: 1, qqr2e: 1 },
  electron: { line: 254, boltz: 3.16681534e-6, hplanck: 0.1519829846, mvv2e: 1.06657236, ftm2v: 0.937582899, nktv2p: 2.94210108e13, qqr2e: 1 },
  micro: { line: 275, boltz: 1.3806504e-8, hplanck: 6.62606896e-13, mvv2e: 1, ftm2v: 1, nktv2p: 1, qqr2e: 8.987556e6 },
  nano: { line: 296, boltz: 0.013806504, hplanck: 6.62606896e-4, mvv2e: 1, ftm2v: 1, nktv2p: 1, qqr2e: 230.7078669 },
};
// SI 2019 exact + CODATA 2018
const NA = 6.02214076e23, E = 1.602176634e-19, KB = 1.380649e-23, H = 6.62607015e-34, EPS0 = 8.8541878128e-12, CAL = 4184, ATM = 101325;
const K4 = 1 / (4 * Math.PI * EPS0), AMU = 1.66053906660e-27, BOHR = 0.529177210903e-10, EH = 4.3597447222071e-18;
const DER = {
  real: { units: 'mass g/mol · length Å · time fs · energy kcal/mol · pressure atm · charge e', boltz: ['k_B N_A / 4184', KB * NA / CAL], hplanck: ['h N_A / 4184 / 10⁻¹⁵', H * NA / CAL / 1e-15], mvv2e: ['(10⁻³/N_A)(10⁻¹⁰/10⁻¹⁵)² / (4184/N_A) = 10⁷/4184', 1e7 / CAL], ftm2v: ['4184/10⁷', CAL / 1e7], nktv2p: ['(4184/N_A) / 10⁻³⁰ / 101325', CAL / NA / 1e-30 / ATM], qqr2e: ['e² / (4πε₀ · 10⁻¹⁰) · N_A / 4184', E * E * K4 / 1e-10 * NA / CAL] },
  metal: { units: 'mass g/mol · length Å · time ps · energy eV · pressure bar · charge e', boltz: ['k_B / e', KB / E], hplanck: ['h / e / 10⁻¹²', H / E / 1e-12], mvv2e: ['(10⁻³/N_A)(10⁻¹⁰/10⁻¹²)² / e = 10 / (N_A e)', 10 / (NA * E)], ftm2v: ['N_A e / 10', NA * E / 10], nktv2p: ['e / 10⁻³⁰ / 10⁵', E / 1e-30 / 1e5], qqr2e: ['e / (4πε₀ · 10⁻¹⁰)', E * K4 / 1e-10] },
  si: { units: 'kg · m · s · J · Pa · C', boltz: ['k_B', KB], hplanck: ['h', H], mvv2e: ['1', 1], ftm2v: ['1', 1], nktv2p: ['1', 1], qqr2e: ['1 / (4πε₀)', K4] },
  cgs: { units: 'g · cm · s · erg · dyne/cm² · statcoulomb', boltz: ['k_B · 10⁷', KB * 1e7], hplanck: ['h · 10⁷', H * 1e7], mvv2e: ['1', 1], ftm2v: ['1', 1], nktv2p: ['1', 1], qqr2e: ['1 (Gaussian)', 1] },
  electron: { units: 'mass amu · length bohr · time fs · energy hartree · pressure Pa · charge e', boltz: ['k_B / E_h', KB / EH], hplanck: ['h / E_h / 10⁻¹⁵', H / EH / 1e-15], mvv2e: ['u (a₀/10⁻¹⁵)² / E_h', AMU * (BOHR / 1e-15) ** 2 / EH], ftm2v: ['E_h / (u a₀²) · 10⁻³⁰', EH / (AMU * BOHR * BOHR) * 1e-30], nktv2p: ['E_h / a₀³', EH / BOHR ** 3], qqr2e: ['1 (atomic units)', 1] },
  micro: { units: 'pg · µm · µs · pg µm²/µs² · picopascal · pC', boltz: ['k_B / 10⁻¹⁵', KB / 1e-15], hplanck: ['h / 10⁻¹⁵ / 10⁻⁶', H / 1e-15 / 1e-6], mvv2e: ['1', 1], ftm2v: ['1', 1], nktv2p: ['1', 1], qqr2e: ['(10⁻¹²)² / (4πε₀ 10⁻⁶) / 10⁻¹⁵', 1e-24 * K4 / 1e-6 / 1e-15] },
  nano: { units: 'ag · nm · ns · ag nm²/ns² · attopascal · e', boltz: ['k_B / 10⁻²¹', KB / 1e-21], hplanck: ['h / 10⁻²¹ / 10⁻⁹', H / 1e-21 / 1e-9], mvv2e: ['1', 1], ftm2v: ['1', 1], nktv2p: ['1', 1], qqr2e: ['e² / (4πε₀ 10⁻⁹) / 10⁻²¹', E * E * K4 / 1e-9 / 1e-21] },
};
const KEYS = ['boltz', 'hplanck', 'mvv2e', 'ftm2v', 'nktv2p', 'qqr2e'];

// live: NaCl pair at run 0, and one atom pushed by a constant force for one step
const VEL = { real: [0.005, -0.003], metal: [5.0, -3.0] }, DT = { real: 1.0, metal: 0.001 };
const scriptA = u => `units           ${u}
atom_style      charge
region          box block 0 20 0 20 0 20 units box
create_box      2 box
create_atoms    1 single 5.0 5.0 5.0 units box
create_atoms    2 single 8.0 5.0 5.0 units box
mass            1 22.99
mass            2 35.45
set             type 1 charge 1.0
set             type 2 charge -1.0
set             atom 1 vx ${VEL[u][0]}
set             atom 2 vx ${VEL[u][1]}
pair_style      coul/cut 10.0
pair_coeff      * *
compute         k all ke/atom
dump            d all custom 1 a.dump id c_k
dump_modify     d format float %.17g sort id
thermo_style    custom step temp ke ecoul press
thermo_modify   norm no format float %.17g
run             0`;
const scriptB = u => `units           ${u}
atom_style      atomic
region          box block 0 20 0 20 0 20 units box
create_box      1 box
create_atoms    1 single 5.0 5.0 5.0 units box
mass            1 22.99
fix             1 all nve
fix             2 all addforce 1.0 0.0 0.0
dump            d all custom 1 b.dump id x vx
dump_modify     d format float %.17g
run             1`;
function hand(u) {
  const c = LMP[u], m = [22.99, 35.45], v = VEL[u], x = [5.0, 8.0], q = [1.0, -1.0];
  const ke = [0, 1].map(i => 0.5 * c.mvv2e * m[i] * (v[i] * v[i] + 0 + 0));
  let t = 0; for (let i = 0; i < 2; i++) t += (v[i] * v[i] + 0 + 0) * m[i];
  const dof = 3.0, T = t * (c.mvv2e / (dof * c.boltz)); let KE = T; KE *= 0.5 * dof * c.boltz;
  const qq = c.qqr2e / 1.0, delx = x[0] - x[1], rsq = delx * delx, r2inv = 1.0 / rsq, rinv = Math.sqrt(r2inv);
  const fpair = 1.0 * (qq * 1.0 * q[0] * q[1] * rinv) * r2inv, ec = 1.0 * qq * 1.0 * q[0] * q[1] * rinv;
  const f = [0.0 + delx * fpair, 0.0 - delx * fpair]; let W = 0; for (let i = 0; i < 2; i++) W += f[i] * x[i];
  const P = (dof * c.boltz * T + W + 0 + 0) / 3.0 * (1.0 / 8000.0) * c.nktv2p;
  const dtfm = 0.5 * DT[u] * c.ftm2v / 22.99; let vb = 0.0 + dtfm * 1.0; const xb = 5.0 + DT[u] * vb; vb += dtfm * 1.0;
  return { ke, T, KE, ec: 0.0 + (0.0 + ec), P, xb, vb, W };
}

export default {
  id: 'units-constants',
  title: 'Units and constants',
  minutes: 20,
  summary: 'Every unit style is 16 numbers typed into Update::set_units. You re-derive them from today\'s exact SI constants and find most agree to 6–8 digits, not more: they are older CODATA values, si units use 8.9876e9 for 1/(4πε₀), and real\'s Planck constant uses a different calorie and is off by 0.07 %. Then you recover boltz, mvv2e, ftm2v, nktv2p and qqr2e from a live run in real and metal units, bit for bit, and see the same atoms come out at temperatures differing in the 8th digit.',
  learn: [
    'Where each constant is used: kinetic energy, the kick, temperature, pressure, Coulomb',
    'Derive every constant of 7 unit styles from SI, and the size of the disagreement',
    'Five constants recovered from one run, matched to the last bit',
    'Why real and metal units do not give exactly the same physics',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: 16 typed numbers per style</h2>
    <div class="prose"><p>${gh('src/update.cpp', 141)}. The <code>units</code> command just fills the <code>force</code> object; every other file multiplies by these. Real units:</p>
<pre class="lmp-mini">force-&gt;boltz  = 0.0019872067;                     // 171  k_B in kcal/mol/K
force-&gt;mvv2e  = 48.88821291 * 48.88821291;        // 173  (g/mol)(Å/fs)² → kcal/mol
force-&gt;ftm2v  = 1.0 / 48.88821291 / 48.88821291;  // 174  (kcal/mol/Å)·fs/(g/mol) → Å/fs
force-&gt;nktv2p = 68568.415;                        // 176  (kcal/mol)/Å³ → atm
force-&gt;qqr2e  = 332.06371;                        // 177  e²/Å → kcal/mol</pre>
<p>Where they are used, on pages you have met:</p></div>
    <div class="tbl"><table>${head(['constant', 'formula in the code', 'where'])}
      ${row(['mvv2e', 'ke = ½ · mvv2e · m · v²', gh('src/compute_ke_atom.cpp', 91)])}
      ${row(['ftm2v', 'dtf = ½ · dt · ftm2v; v += dtf · F / m', gh('src/fix_nve.cpp', 59)])}
      ${row(['boltz', 'T = mvv2e · Σ m v² / (dof · boltz)', gh('src/compute_temp.cpp', 65)])}
      ${row(['nktv2p', 'P = (dof · boltz · T + W) / 3V · nktv2p', gh('src/compute_pressure.cpp', 256)])}
      ${row(['qqr2e', 'E = qqrd2e · qᵢ qⱼ / r, qqrd2e = qqr2e / dielectric', gh('src/pair_coul_cut.cpp', 116)])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Derive: every constant from SI</h2>
    <div class="prose"><p>Since 2019, N_A, e, k_B and h are exact; ε₀, the bohr, the hartree and the atomic mass unit come from CODATA 2018; 1 cal = 4.184 J, 1 atm = 101 325 Pa.</p></div>
    <div class="row"><label class="lc-in">units<select class="lc-u uc-u">${Object.keys(LMP).map(k => `<option${k === 'real' ? ' selected' : ''}>${k}</option>`).join('')}</select></label></div>
    <p class="muted uc-desc"></p>
    <div class="tbl uc-der"></div>
  </section>

  <section class="block">
    <h2>3 · Match: recover the constants from a run</h2>
    <div class="prose"><p>Two runs. (A) Na⁺ at x = 5, Cl⁻ at x = 8 Å, moving apart along x, <code>run 0</code>: kinetic energy per atom, T, KE, E_coul, Press. (B) one Na atom at rest pushed by <code>fix addforce 1.0 0 0</code> for one step: x and v. Your hand values use only the constants of section 1, in the order the C++ uses them.</p></div>
    <div class="row"><label class="lc-in">units<select class="lc-u uc-m"><option>real</option><option>metal</option></select></label><button class="btn sm primary" data-k="run">▶ Run A and B</button><span class="lmp-status muted uc-st"></span></div>
    <div class="tbl uc-out"></div>
    <div class="prose uc-sum"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe</h2>
    <div class="tbl uc-cmp"></div>
    <ul>
      <li><b>Same atoms, different temperature.</b> 0.005 Å/fs is exactly 5 Å/ps, yet T differs between real and metal in the 8th digit: each style's mvv2e/boltz was typed separately, and the two ratios disagree at 3.5 × 10⁻⁸.</li>
      <li><b>The constants are old.</b> metal's boltz 8.617343e-5 and si's 1.3806504e-23 are the CODATA 2006 values; today's exact k_B gives 8.617333e-5, a difference of 1.1 × 10⁻⁶. si's qqr2e = 8.9876e9 is 1/(4πε₀) cut to 5 digits: 5 × 10⁻⁶ too large.</li>
      <li><b>One constant is off by 0.07 %.</b> real's hplanck, 95.306976368, is h·N_A divided by 4186.8 J/kcal (the International Table calorie), while everything else in the style, and the comment at the top of set_units, uses 4184 J. From SI it should be 95.37076. It enters wherever ħ does in real units: the de Broglie wavelength of fix gcmc (${gh('src/MC/fix_gcmc.cpp', 704)}), the bead springs of fix pimd/nvt (${gh('src/REPLICA/fix_pimd_nvt.cpp', 220)}), fix qtb. It is already reported upstream as <a href="https://github.com/lammps/lammps/issues/5295" target="_blank" rel="noopener">issue #5295</a> (open as of November 2026).</li>
      <li><b>CHARMM changes qqr2e behind your back.</b> ${gh('src/force.cpp', 64)} keeps two values, 332.06371 and 332.0716; the newer CHARMM pair styles switch to the second (with a warning), so Coulomb energies in real units move by 2.4 × 10⁻⁵.</li>
      <li><b>None of this matters for most physics</b>: 10⁻⁶ in a constant is far below statistical error. It matters when you compare codes to many digits, or convert results between unit styles and expect exact agreement.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Derive metal's nktv2p: 1 eV/Å³ in bar.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>1.602176634 × 10⁻¹⁹ J / 10⁻³⁰ m³ = 1.602176634 × 10¹¹ Pa = <b>1.602176634 × 10⁶ bar</b>. LAMMPS has 1.6021765e6.</div></div>
    <div class="ex"><p><b>2.</b> In real units, k_B T at 300 K in kcal/mol, and the same in metal units in eV?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>0.0019872067 × 300 = <b>0.59616 kcal/mol</b>; 8.617343e-5 × 300 = <b>0.025852 eV</b>.</div></div>
    <div class="ex"><p><b>3.</b> You switch a real-units script to metal. Which numbers in the input must you change by hand?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Every energy (pair_coeff ε, bond K, …) from kcal/mol to eV (÷ 23.0605), the timestep (fs → ps, ÷ 1000), damping times, and pressures (atm → bar, × 1.01325). Lengths, masses, charges and temperatures stay. LAMMPS also resets the default timestep and skin when units change (${gh('src/update.cpp', 325)}).</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const der = () => {
      const u = $('.uc-u').value, L = LMP[u], D = DER[u];
      $('.uc-desc').innerHTML = `${esc(D.units)} · ${gh('src/update.cpp', L.line)}`;
      $('.uc-der').innerHTML = `<table>${head(['constant', 'in LAMMPS', 'from SI', '= ', 'relative difference'])}${KEYS.map(k => { const d = D[k][1], l = L[k], rd = (l - d) / d; return row([k, `<span class="mono">${l.toPrecision(11).replace(/\.?0+(e|$)/, '$1')}</span>`, esc(D[k][0]), `<span class="mono">${d.toPrecision(11).replace(/\.?0+(e|$)/, '$1')}</span>`, Math.abs(rd) < 1e-12 ? 'exact' : `<b>${rd.toExponential(1)}</b>`]); }).join('')}</table>`;
    };
    $('.uc-u').onchange = der; der();

    const R = hand('real'), M = hand('metal');
    $('.uc-cmp').innerHTML = `<table>${head(['', 'real', 'metal', 'ratio', 'expected'])}
      ${row(['T (K)', s17(R.T), s17(M.T), (R.T / M.T).toFixed(12), '1'])}
      ${row(['E_coul', `${s17(R.ec)} kcal/mol`, `${s17(M.ec)} eV`, (R.ec / M.ec).toFixed(8), `${(NA * E / CAL).toFixed(8)} kcal/mol per eV (SI)`])}
      ${row(['Press', `${s17(R.P)} atm`, `${s17(M.P)} bar`, (M.P / R.P).toFixed(8), '1.01325 bar per atm'])}</table><p class="muted">Computed by your hand code, which section 3 shows to be LAMMPS, bit for bit.</p>`;

    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.uc-st'), u = $('.uc-m').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const ra = await runLammps(scriptA(u), { read: ['a.dump'] }), rb = await runLammps(scriptB(u), { read: ['b.dump'] }); b.disabled = false;
      const da = ra.files?.['a.dump'], db = rb.files?.['b.dump'];
      if (ra.error || rb.error || !da || !db || !ra.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(ra.error || rb.error || 'no output')}</span>`; return; }
      st.textContent = '';
      const h = hand(u), th = ra.thermo.rows[0], c = ra.thermo.cols, g = k => parseFloat(th[c.indexOf(k)]);
      const ka = da.split('ITEM: ATOMS')[1].trim().split('\n').slice(1).map(l => parseFloat(l.trim().split(/\s+/)[1]));
      const bb = db.split('ITEM: TIMESTEP').pop().split('ITEM: ATOMS')[1].trim().split('\n').slice(1)[0].trim().split(/\s+/).map(Number);
      const U = u === 'real' ? ['kcal/mol', 'Å', 'Å/fs', 'atm'] : ['eV', 'Å', 'Å/ps', 'bar'];
      const rows = [
        ['ke/atom Na', 'mvv2e', '½ · mvv2e · 22.99 · v₁²', h.ke[0], ka[0], U[0]], ['ke/atom Cl', 'mvv2e', '½ · mvv2e · 35.45 · v₂²', h.ke[1], ka[1], U[0]],
        ['Temp', 'mvv2e/boltz', 'mvv2e · Σmv² / (3 · boltz)', h.T, g('Temp'), 'K'], ['KinEng', 'boltz', '½ · 3 · boltz · T', h.KE, g('KinEng'), U[0]],
        ['E_coul', 'qqr2e', 'qqr2e · (+1)(−1) / 3', h.ec, g('E_coul'), U[0]], ['Press', 'nktv2p', '(3 boltz T + W) / 3V · nktv2p', h.P, g('Press'), U[3]],
        ['x after 1 step', 'ftm2v', '5 + dt · (½ dt ftm2v / m)', h.xb, bb[1], U[1]], ['vx after 1 step', 'ftm2v', '2 · (½ dt ftm2v / m) · 1', h.vb, bb[2], U[2]]];
      let ok = 0;
      $('.uc-out').innerHTML = `<table>${head(['', 'constant', 'by hand', '', 'LAMMPS', ''])}${rows.map(([n, k, f, hv, lv, un]) => { const s = hv === lv; if (s) ok++; return row([n, k, esc(f), `<span class="mono">${s17(hv)}</span> ${un}`, `<span class="mono">${s ? '✓ identical' : s17(lv)}</span>`, s ? OK : BAD]); }).join('')}</table>`;
      $('.uc-sum').innerHTML = `<p>${ok} of ${rows.length} identical, in ${u} units. W = Σ x·f = ${s17(h.W)}.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
