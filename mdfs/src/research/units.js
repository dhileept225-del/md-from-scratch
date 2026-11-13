// Research · A.4: units. The same melt in LJ and in real units (argon): every number converted by hand and matched step by step.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/src/update.cpp';
const SIG = 3.405, EPS = 0.2381, KB = 0.0019872067, NKTV2P = 68568.415, N = 4000;

const LJ = IN_MELT.replace('loop geom', 'loop all');
const REAL = `# melt in real units: argon, sigma = 3.405 A, epsilon = 0.2381 kcal/mol, m = 39.948 g/mol
units           real

variable        sig equal 3.405
variable        eps equal 0.2381
atom_style      atomic

lattice         fcc $((4/0.8442)^(1/3)*v_sig)
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 39.948

velocity        all create $(3*v_eps/0.0019872067) 87287 loop all

pair_style      lj/cut $(2.5*v_sig)
pair_coeff      1 1 \${eps} \${sig} $(2.5*v_sig)

neighbor        $(0.3*v_sig) bin
neigh_modify    every 20 delay 0 check no
timestep        $(0.005*v_sig*sqrt(39.948/v_eps)*48.88821291)

fix             1 all nve
thermo          50
run             250
`;

export default {
  id: 'units',
  title: 'Units',
  minutes: 25,
  summary: 'The units command decides what every number in a script means. You see the conversion constants inside LAMMPS, derive one by hand, then run melt twice: in Lennard-Jones units and as real argon in Å, femtoseconds and kcal/mol. Converted by hand, the two runs agree at every printed step.',
  learn: [
    'What each unit style means for mass, length, time, energy, temperature and pressure',
    'The constants LAMMPS uses inside (Boltzmann, force-to-acceleration, pressure) and where one comes from',
    'Convert melt to argon in real units, predict its output by hand and match it step by step',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: the same number means different things</h2>
    <p>Set <code>units</code> before creating the simulation box or defining force-field settings, so LAMMPS interprets later values correctly. This example places it near the top. <code>timestep 1.0</code> is one femtosecond in real units and one picosecond in metal units.</p>
    <div class="tbl wide"><table>${head(['', 'lj', 'real', 'metal', 'si'])}
${row(['mass', 'm (set to 1)', 'g/mol', 'g/mol', 'kg'])}
${row(['distance', 'σ', 'Å', 'Å', 'm'])}
${row(['time', 'τ = σ√(m/ε)', 'fs', 'ps', 's'])}
${row(['energy', 'ε', 'kcal/mol', 'eV', 'J'])}
${row(['temperature', 'ε/k<sub>B</sub>', 'K', 'K', 'K'])}
${row(['pressure', 'ε/σ³', 'atm', 'bar', 'Pa'])}
${row(['default time step', '0.005 τ', '1 fs', '0.001 ps', '10⁻⁸ s'])}
${row(['default neighbour skin', '0.3 σ', '2 Å', '2 Å', '0.001 m'])}
${row(['thermo energies are', 'per atom', 'totals', 'totals', 'totals'])}</table></div>
    <p>The last row matters when comparing numbers: in LJ units thermo divides energies by the number of atoms; in all other styles it prints the total.</p>
  </section>

  <section class="block prose">
    <h2>2 · The constants inside LAMMPS</h2>
    <p>Equations of motion need a few conversion factors. LAMMPS sets them in <a href="${GH}#L141" target="_blank" rel="noopener">Update::set_units()</a>. For real units (lines 171–176 and nearby):</p>
    <div class="calc">force->boltz  = 0.0019872067;                     <b>k_B in kcal/mol/K</b>
force->mvv2e  = 48.88821291 * 48.88821291;        <b>mass × velocity² → energy</b>
force->ftm2v  = 1.0 / 48.88821291 / 48.88821291;  <b>force / mass → acceleration</b>
force->nktv2p = 68568.415;                        <b>kcal/mol/Å³ → atm</b>
force->qqr2e  = 332.06371;                        <b>charge²/distance → kcal/mol</b></div>
    <p>In LJ units all five are 1. Where does 48.88821291 come from? A force of 1 kcal/mol/Å on a mass of 1 g/mol gives an acceleration of:</p>
    <div class="calc">a = F/m = (4184 J per mol / 10⁻¹⁰ m) / (10⁻³ kg per mol)
  = 4.184 × 10¹⁶ m/s²
  = 4.184 × 10¹⁶ × 10¹⁰ Å / (10¹⁵ fs)² = <b>4.184 × 10⁻⁴ Å/fs²</b>

ftm2v = 1/48.88821291² = 1/2390.0574 = <b>4.1840 × 10⁻⁴</b>   ✓</div>
  </section>

  <section class="block prose">
    <h2>3 · Modify: melt as real argon</h2>
    <p>Argon is often modelled with σ = 3.405 Å, ε = 0.2381 kcal/mol (0.996 kJ/mol) and m = 39.948 g/mol. Every LJ number in melt becomes a real number by multiplying with the right combination:</p>
    <div class="calc">${esc(REAL)}</div>
    <p>Both scripts use <code>loop all</code> instead of <code>loop geom</code> for the velocities. With <code>loop geom</code> the random velocity of each atom is seeded from its coordinates, and the coordinates are different numbers in Å and in σ, so the two runs would start with different (equally valid) velocities.</p>
    <h3>Predict by hand</h3>
    <div class="calc">lattice constant  a  = 1.6795962 σ × 3.405 Å/σ             = <b>5.719025 Å</b>
box               L  = 10 a                                  = <b>57.19025 Å</b>
temperature       T  = 3 ε/k_B = 3 × 0.2381 / 0.0019872067    = <b>359.44927 K</b>
time unit         τ  = σ √(m/ε) × 48.88821291 fs
                     = 3.405 × 12.952924 × 48.88821291      = 2156.2003 fs
time step         Δt = 0.005 τ                               = <b>10.781001 fs</b>
energies are totals here: multiply the per-atom LJ value by ε × N = 0.2381 × 4000 = 952.4
E_pair = −6.77336805 × 952.4                                 = <b>−6450.9557 kcal/mol</b>
TotEng = −2.27449305 × 952.4                                 = <b>−2166.2272 kcal/mol</b>
pressure: ε/σ³ × nktv2p = 0.2381/39.477655 × 68568.415      = 413.55393 atm
Press  = −3.7033504 × 413.55393                              = <b>−1531.5351 atm</b></div>
  </section>

  <section class="block">
    <h2>4 · Observe and match: run both</h2>
    <p class="muted">Each button runs the real LAMMPS (250 steps, about 2 seconds). After both runs, every real-unit number is converted back to LJ units by hand and compared with the LJ run.</p>
    <div class="row"><button class="btn sm primary" data-k="lj">▶ Run in LJ units</button><button class="btn sm primary" data-k="real">▶ Run in real units</button><span class="lmp-status muted ust"></span></div>
    <div class="u-step0 lmp-compare"></div>
    <div class="u-cmp lmp-compare"></div>
  </section>

  <section class="block prose">
    <h2>What the two runs show</h2>
    <ul>
      <li>Converted, the two runs agree at every printed step to 6–8 digits: units change the numbers, never the physics. The last digits differ because the computer rounds different numbers (Å instead of σ).</li>
      <li>Every quantity needs its own conversion factor: k<sub>B</sub> for temperature, ε × N for energies (and the per-atom normalisation), ε/σ³ × nktv2p for pressure.</li>
      <li>LJ results are therefore general: the same run describes argon, krypton or any LJ fluid, once you choose σ, ε and m.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> A real-units run prints Temp = 200.79919 K at step 50. What is that in LJ units?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>T* = T k<sub>B</sub>/ε = 200.79919 × 0.0019872067 / 0.2381 = <b>1.6758903</b>, the LJ run's value at step 50.</div></div>
    <div class="ex"><p><b>2.</b> You copy <code>timestep 2.0</code> from a real-units water script into a metal-units script. What time step do you get?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>2 ps instead of 2 fs: a thousand times too large. The simulation will blow up within a few steps. Always change the time step when you change units.</div></div>
    <div class="ex"><p><b>3.</b> Why do the real-units energies look 952 times larger than the LJ ones?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Two effects multiply: the energy unit (1 ε = 0.2381 kcal/mol, a factor 0.2381) and the normalisation (LJ prints per atom, real prints the total for 4000 atoms). 0.2381 × 4000 = 952.4.</div></div>
  </section>`,

  mount(root) {
    const st = root.querySelector('.ust'); const got = {};
    const conv = { Temp: v => v * KB / EPS, E_pair: v => v / (EPS * N), TotEng: v => v / (EPS * N), Press: v => v / (EPS / SIG ** 3 * NKTV2P) };
    const cell = (h, l) => { const m = matchMarkup(g8(h), l); return `<td class="mono">${m.html}</td>`; };
    function draw() {
      if (got.real) {
        const r = got.real.thermo.rows.find(w => w[0] === '0'), c = got.real.thermo.cols;
        const H = { Temp: 359.44927, E_pair: -6450.9557, TotEng: -2166.2272, Press: -1531.5351 };
        root.querySelector('.u-step0').innerHTML = `<div class="tbl"><table><tr><th>real units, step 0</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${Object.entries(H).map(([k, h]) => { const lv = r[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
      }
      if (got.real && got.lj) {
        const R = got.real.thermo, L = got.lj.thermo, keys = ['Temp', 'E_pair', 'TotEng', 'Press'];
        const rows = R.rows.map(w => { const l = L.rows.find(x => x[0] === w[0]); if (!l) return ''; return `<tr><td>${w[0]}</td>${keys.map(k => cell(conv[k](+w[R.cols.indexOf(k)]), l[L.cols.indexOf(k)]) + `<td class="mono">${l[L.cols.indexOf(k)]}</td>`).join('')}</tr>`; }).join('');
        root.querySelector('.u-cmp').innerHTML = `<p><b>Real-unit output converted by hand</b> (green: digits that agree with the LJ run, printed next to it)</p><div class="tbl"><table><tr><th>step</th>${keys.map(k => `<th>${k} real→LJ</th><th>LJ run</th>`).join('')}</tr>${rows}</table></div>`;
      }
    }
    const run = async (which, btn) => {
      btn.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(which === 'lj' ? LJ : REAL); btn.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      got[which] = res; st.textContent = `${which === 'lj' ? 'LJ' : 'Real'} run done in ${res.seconds.toFixed(1)} s.${got.lj && got.real ? '' : ' Now run the other one.'}`; draw();
    };
    root.querySelector('[data-k=lj]').onclick = e => run('lj', e.currentTarget);
    root.querySelector('[data-k=real]').onclick = e => run('real', e.currentTarget);
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
