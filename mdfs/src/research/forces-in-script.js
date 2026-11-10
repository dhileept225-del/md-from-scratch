// Research · A.7: forces and neighbours in the script: pair_style, pair_coeff, pair_modify (shift, tail, mix), neighbor, neigh_modify.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const run0 = s => s.replace('run             250', 'run             0');

const TWO_TYPES = mix => `units           lj
atom_style      atomic
atom_modify     map yes
region          box block -5 5 -5 5 -5 5
create_box      2 box
create_atoms    1 single 0.0 0.0 0.0
create_atoms    2 single 1.2 0.0 0.0
mass            * 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0
pair_coeff      2 2 0.5 1.2
pair_modify     mix ${mix}
variable        f equal fx[1]
thermo_style    custom step pe v_f
thermo_modify   norm no format float %.8g
run             0
`;

const EXP = [
  { key: 'shift', title: 'Shift the energy to zero at the cutoff', line: 'pair_modify     shift yes',
    script: () => run0(IN_MELT.replace('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 1 1.0 1.0 2.5\npair_modify     shift yes')),
    calc: `without the shift, u jumps from u(2.5) to 0 at the cutoff:
u(2.5) = 4(2.5⁻¹² − 2.5⁻⁶) = 4(0.0000168 − 0.004096) = −0.01631689
shift yes subtracts u(2.5) from every pair inside the cutoff:
each atom has 54 neighbours within 2.5, each pair counted half:
ΔE_pair = −½ × 54 × (−0.01631689) = +0.44055606
E_pair = −6.77336805 + 0.44055606  = <b>−6.332812</b>
TotEng = −6.33281199 + 4.49887500 = <b>−1.833937</b>
Press  = <b>−3.7033504</b>  unchanged: a constant shift does not change the force`,
    hand: { E_pair: -6.332812, TotEng: -1.833937, Press: -3.7033504 } },
  { key: 'tail', title: 'Add the missing tail beyond the cutoff', line: 'pair_modify     tail yes',
    script: () => run0(IN_MELT.replace('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 1 1.0 1.0 2.5\npair_modify     tail yes')),
    calc: `beyond 2.5 the liquid is assumed uniform (g = 1):
E_tail = (8/3)πρ [ (1/3) r_c⁻⁹ − r_c⁻³ ]
       = 7.0723534 × (0.0000874 − 0.064) = −0.45201262
P_tail = (16/3)πρ² [ (2/3) r_c⁻⁹ − r_c⁻³ ]
       = 11.940961 × (0.0001748 − 0.064) = −0.76213470
E_pair = −6.77336805 − 0.45201262  = <b>−7.2253807</b>
TotEng = −2.27449305 − 0.45201262  = <b>−2.7265057</b>
Press  = −3.70335042 − 0.76213470  = <b>−4.4654851</b>`,
    hand: { E_pair: -7.2253807, TotEng: -2.7265057, Press: -4.4654851 } },
];

const NEIGH = [
  ['every 20 delay 0 check no', 'melt as shipped'],
  ['every 20 delay 0 check yes', 'same schedule, but count the dangerous builds'],
  ['every 10 delay 0 check yes', 'look twice as often'],
  ['every 1 delay 0 check yes', 'look every step, rebuild only when needed'],
];

export default {
  id: 'forces-in-script',
  title: 'Forces and neighbours in the script',
  minutes: 25,
  summary: 'Five commands control how forces are computed: pair_style, pair_coeff, pair_modify, neighbor and neigh_modify. You predict by hand what shifting the energy, adding a tail correction and mixing two atom types do, run each live, and then discover that the neighbour settings of the official melt example are not safe.',
  learn: [
    'What pair_style, pair_coeff and pair_modify (shift, tail, mix) change, predicted and checked',
    'How many neighbours each atom has in the list, counted by hand',
    'How neigh_modify trades speed against missed interactions ("dangerous builds")',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: five commands</h2>
    <div class="tbl wide"><table>${head(['Command', 'Sets', 'In melt'])}
${row(['<code>pair_style lj/cut 2.5</code>', 'the kind of force, and the default cutoff', 'Lennard-Jones, cut at 2.5'])}
${row(['<code>pair_coeff 1 1 1.0 1.0 2.5</code>', 'ε, σ (and optionally the cutoff) for one pair of atom types', 'types 1–1: ε = 1, σ = 1'])}
${row(['<code>pair_modify …</code>', '<code>shift yes</code>: energy 0 at the cutoff; <code>tail yes</code>: add the part beyond the cutoff; <code>mix</code>: rule for unlike types', 'not used'])}
${row(['<code>neighbor 0.3 bin</code>', 'the skin: neighbour lists reach cutoff + 0.3 = 2.8; bin: how they are built (Lesson 2.6)', '2.8'])}
${row(['<code>neigh_modify every 20 delay 0 check no</code>', 'when to rebuild the lists', 'every 20 steps, no checking'])}</table></div>
  </section>

  <section class="block prose">
    <h2>2 · Modify, predict, observe, match</h2>
    ${EXP.map(e => `
    <div class="exp" data-e="${e.key}">
      <h3>${e.title}</h3>
      <div class="calc">${e.line}</div>
      <p><b>Predict by hand:</b></p>
      <div class="calc">${e.calc}</div>
      <div class="row"><button class="btn sm primary" data-run="${e.key}">▶ Run it</button><span class="lmp-status muted"></span></div>
      <div class="lmp-compare"></div>
    </div>`).join('')}
    <div class="exp" data-e="mix">
      <h3>Two atom types: mixing rules</h3>
      <p>Give type 2 its own ε = 0.5 and σ = 1.2 and say nothing about the pair 1–2. LAMMPS invents it with a mixing rule; <code>lj/cut</code> uses the geometric rule unless told otherwise. Two atoms, type 1 at x = 0 and type 2 at x = 1.2:</p>
      <div class="calc">geometric : ε₁₂ = √(1.0 × 0.5) = 0.70710678   σ₁₂ = √(1.0 × 1.2) = 1.0954451
arithmetic: ε₁₂ = √(1.0 × 0.5) = 0.70710678   σ₁₂ = (1.0 + 1.2)/2 = 1.1
E  = 4ε₁₂[(σ₁₂/r)¹² − (σ₁₂/r)⁶]          fx on atom 1 = −24ε₁₂[2(σ₁₂/r)¹² − (σ₁₂/r)⁶]/r
geometric : (σ/r)⁶ = 0.57870370  E = <b>−0.68958673</b>   fx = <b>−1.2882389</b>
arithmetic: (σ/r)⁶ = 0.59329219  E = <b>−0.68248975</b>   fx = <b>−1.5655211</b></div>
      <div class="row"><button class="btn sm primary" data-mix="geometric">▶ Run geometric</button><button class="btn sm primary" data-mix="arithmetic">▶ Run arithmetic</button><span class="lmp-status muted"></span></div>
      <div class="lmp-compare"></div>
      <p class="muted small">Force fields are built for one rule: OPLS uses geometric, AMBER and CHARMM use arithmetic (Lorentz–Berthelot). Using the wrong one changes every unlike interaction.</p>
    </div>
  </section>

  <section class="block prose">
    <h2>3 · The neighbour list, counted by hand</h2>
    <p>The list holds every pair closer than cutoff + skin = 2.8. In the perfect crystal at step 0 the shells inside 2.8 are 12 + 6 + 24 + 12 + 24 = 78 atoms (the fifth, at 2.6557, is inside the skin but outside the cutoff). Each pair is stored once (Newton's third law), so:</p>
    <div class="calc">neighbours per atom = 78 / 2 = <b>39</b>
total in the list    = 4000 × 39 = <b>156000</b></div>
    <div class="row"><button class="btn sm primary" data-k="nb">▶ Run step 0 only</button><span class="lmp-status muted nb-st"></span></div>
    <div class="lmp-compare nb-cmp"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe: are melt's neighbour settings safe?</h2>
    <p>The list is valid only until some atom has moved more than half the skin (0.15). Then a pair that was outside 2.8 may already be inside 2.5 and is missed. A rebuild that happens later than that is a <b>dangerous build</b>. Estimate: after melting T ≈ 1.7, typical speed √(3T) ≈ 2.3, the fastest atoms about three times that. In 20 steps (0.1 τ) a fast atom moves about 0.7, far more than 0.15.</p>
    <div class="tbl wide"><table class="nb-table">${head(['neigh_modify', 'idea', 'list builds', 'dangerous', 'TotEng at 250', ''])}
${NEIGH.map(([s, idea], i) => `<tr data-n="${i}"><td><code>${s}</code></td><td>${idea}</td><td class="mono b"></td><td class="mono d"></td><td class="mono te"></td><td><button class="btn sm ghost" data-nrun="${i}">Run</button></td></tr>`).join('')}</table></div>
    <p>Our runs gave: melt as shipped rebuilt 12 times and every one of them was dangerous (it simply did not check). Looking every step, LAMMPS rebuilt 42 times (about every 6 steps) with no dangerous builds. The total energy at step 250 also moves: −2.2812174 as shipped, −2.2797308 when safe, against −2.2744931 at step 0. The safe run conserves energy better. The melt example is a speed demonstration, not a model for real work.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> With <code>pair_modify shift yes</code>, why does the pressure not change at all?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The shift adds the same constant to the energy of every pair inside the cutoff. The force is minus the slope of the energy, and a constant has no slope. Pressure comes from forces and velocities only, so it is identical to every digit.</div></div>
    <div class="ex"><p><b>2.</b> With a skin of 0.5 instead of 0.3, how many neighbours per atom does the step-0 list hold?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Cutoff + skin = 3.0. The sixth shell (8 atoms at a√3 = 2.9091) is now inside: 12 + 6 + 24 + 12 + 24 + 8 = 86 atoms, so <b>43</b> per atom in the half list, 172000 in total. Bigger lists cost more per step but need rebuilding less often.</div></div>
    <div class="ex"><p><b>3.</b> The tail correction lowers the pressure by 0.76. Is that small?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No: the whole pressure is −3.70, so the tail is 20% of it. Energies change by 7%. When you compare a cut potential with experiment or with a paper, the cutoff and the tail setting matter.</div></div>
  </section>`,

  mount(root) {
    const table = (H, r, c) => `<div class="tbl"><table><tr><th>step 0</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${Object.entries(H).map(([k, h]) => { const lv = r[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
    root.querySelectorAll('[data-run]').forEach(b => b.onclick = async () => {
      const e = EXP.find(x => x.key === b.dataset.run), box = b.closest('.exp'), st = box.querySelector('.lmp-status');
      b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(e.script()); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      box.querySelector('.lmp-compare').innerHTML = table(e.hand, res.thermo.rows.find(w => w[0] === '0'), res.thermo.cols);
    });
    const MIX = { geometric: { PotEng: -0.68958673, v_f: -1.2882389 }, arithmetic: { PotEng: -0.68248975, v_f: -1.5655211 } };
    root.querySelectorAll('[data-mix]').forEach(b => b.onclick = async () => {
      const box = b.closest('.exp'), st = box.querySelector('.lmp-status'), mix = b.dataset.mix;
      b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(TWO_TYPES(mix)); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = (res.lines.find(l => /mixed pair_coeff/.test(l)) || '').trim();
      const r = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols;
      box.querySelector('.lmp-compare').innerHTML = `<p><b>${mix}</b></p>` + table({ PotEng: MIX[mix].PotEng, v_f: MIX[mix].v_f }, r, c).replace('<td>v_f</td>', '<td>fx on atom 1</td>').replace('<td>PotEng</td>', '<td>energy</td>');
    });
    root.querySelector('[data-k=nb]').onclick = async e => {
      const b = e.currentTarget, st = root.querySelector('.nb-st'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(run0(IN_MELT)); b.disabled = false; st.textContent = '';
      const tot = (res.lines.find(l => l.startsWith('Total # of neighbors')) || '').split('=')[1]?.trim(), ave = (res.lines.find(l => l.startsWith('Ave neighs/atom')) || '').split('=')[1]?.trim();
      root.querySelector('.nb-cmp').innerHTML = `<div class="tbl"><table><tr><th></th><th>by hand</th><th>LAMMPS</th><th></th></tr><tr><td>neighbours per atom</td><td>39</td><td class="mono">${ave}</td><td>${ave === '39' ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr><tr><td>total</td><td>156000</td><td class="mono">${tot}</td><td>${tot === '156000' ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr></table></div>`;
    };
    root.querySelectorAll('[data-nrun]').forEach(b => b.onclick = async () => {
      const i = +b.dataset.nrun, tr = root.querySelector(`tr[data-n="${i}"]`); b.disabled = true; b.textContent = '…';
      const res = await runLammps(IN_MELT.replace('neigh_modify    every 20 delay 0 check no', `neigh_modify    ${NEIGH[i][0]}`)); b.disabled = false; b.textContent = 'Run';
      const get = p => (res.lines.find(l => l.startsWith(p)) || '').split('=')[1]?.trim() ?? (res.lines.find(l => l.startsWith(p)) ? 'not checked' : '');
      tr.querySelector('.b').textContent = get('Neighbor list builds');
      tr.querySelector('.d').textContent = res.lines.some(l => l.startsWith('Dangerous builds not checked')) ? 'not checked' : get('Dangerous builds');
      const r = res.thermo.rows.find(w => w[0] === '250'); tr.querySelector('.te').textContent = r ? r[res.thermo.cols.indexOf('TotEng')] : '';
    });
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
