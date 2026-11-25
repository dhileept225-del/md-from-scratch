// Research · C.6: setup, multiple runs and rerun. What Verlet::setup redoes, why run 50 + run 50 is not run 100, ramps restarting, and rerun's fresh neighbour lists.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const rowsOf = res => res.thermo.rows.map(r => Object.fromEntries(res.thermo.cols.map((c, i) => [c, r[i]])));
const L = 10 * 1.6795961913825073;

const NB = { melt: 'neigh_modify    every 20 delay 0 check no', safe: 'neigh_modify    every 1 delay 0 check yes' };
const base = nb => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 3.0 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
${NB[nb]}
fix             1 all nve
`;
const SPLIT = { one: 'run             100', two: 'run             50\nrun             50', pre: 'run             50\nrun             50 pre no' };
const scriptA = (nb, sp) => `${base(nb)}thermo_style    custom step pe etotal
thermo_modify   norm no format float %.17g
thermo          50
${SPLIT[sp]}
write_dump      all custom f.dump id x y z vx vy vz modify format float %.17g sort id`;
const RAMP = { one: 'run             100', two: 'run             50\nrun             50', ss: 'run             50 start 0 stop 100\nrun             50 start 0 stop 100' };
const scriptB = r => `${base('melt')}fix             2 all temp/rescale 10 1.0 2.0 0.0 1.0
thermo_style    custom step temp
thermo_modify   format float %.17g
thermo          10
${RAMP[r]}`;
const scriptC = nb => `${base(nb)}dump            d all custom 10 r.dump id x y z
dump_modify     d format float %.17g
thermo_style    custom step temp pe press
thermo_modify   norm no format float %.17g
thermo          10
run             50
undump          d
rerun           r.dump dump x y z`;

const parseDump = t => { const m = {}; const ls = t.split('ITEM: ATOMS')[1].trim().split('\n').slice(1); for (const l of ls) { const w = l.trim().split(/\s+/).map(Number); m[w[0]] = w.slice(1); } return m; };

export default {
  id: 'multiple-runs',
  title: 'Setup, multiple runs and rerun',
  minutes: 20,
  summary: 'Every run command starts with Verlet::setup: wrap, sort, rebuild the neighbour list, recompute forces. So run 50 + run 50 is not run 100. You measure how different it is with the melt settings (a jump in the energy, from pairs the old list had missed) and with safe settings (round-off only), predict by hand how a temperature ramp restarts at the second run, and use rerun to find the missed pairs in a finished trajectory.',
  learn: [
    'What setup redoes and what carries over between runs',
    'run 50 + run 50 vs run 100 vs pre no: bitwise, round-off or real difference',
    'beginstep and endstep: a ramp restarts unless you give start/stop',
    'rerun recomputes every frame with a fresh list; its pressure uses leftover velocities',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: what a new run redoes</h2>
    <div class="prose"><p>${gh('src/run.cpp', 163)}: each <code>run</code> sets the ramp window, then, unless <code>pre no</code>, calls <code>lmp->init()</code> and <code>Verlet::setup()</code> (${gh('src/verlet.cpp', 93)}):</p>
<pre class="lmp-mini">163  if (startflag) update->beginstep = start; else update->beginstep = update->firststep;
165  if (stopflag) update->endstep = stop;     else update->endstep = update->laststep;
168  if (preflag || update->first_update == 0) {
169    lmp->init();
170    update->integrate->setup(1);
         ↓ verlet.cpp
118      domain->pbc();                 // wrap every atom into the box
123      if (atom->sortfreq > 0) atom->sort();   // reorder atoms in memory
129      neighbor->build(1);            // a fresh list, whatever the last one was
140      force->pair->compute(eflag,vflag);
159      modify->setup(vflag);          // fixes' setup()
160      output->setup(flag);           // the thermo line at the first step, again</pre></div>
    <div class="tbl"><table>${head(['carries over', 'is redone or reset'])}${row(['x, v, image flags, the step counter', 'neighbour list, atom order in memory, forces'])}${row(['fix state: Nosé–Hoover ξ, ecouple, random-number streams, fix ave sums', 'ramps of fixes (temp/rescale, nvt Tstart→Tstop): beginstep, endstep'])}${row(['compute msd reference positions', 'thermo header line and the first-step output'])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Modify: one run of 100, or two of 50?</h2>
    <div class="row">
      <label class="lc-in">neighbours<select class="lc-u mr-nb"><option value="melt">every 20 check no (melt)</option><option value="safe">every 1 check yes (default)</option></select></label>
      <label class="lc-in">runs<select class="lc-u mr-sp"><option value="two">run 50 + run 50</option><option value="pre">run 50 + run 50 pre no</option></select></label>
      <button class="btn sm primary" data-k="a">▶ Run both ways</button><span class="lmp-status muted mr-sta"></span></div>
    <div class="tbl mr-outa"></div>
    <div class="prose mr-suma"></div>
  </section>

  <section class="block">
    <h2>3 · Match: a ramp and its restart</h2>
    <div class="prose"><p><code>fix temp/rescale 10 1.0 2.0 0.0 1.0</code> sets T exactly to the target every 10 steps (window 0, fraction 1). The target, ${gh('src/fix_temp_rescale.cpp', 144)}:</p>
<pre class="lmp-mini">delta = (step − beginstep) / (endstep − beginstep)
t_target = t_start + delta · (t_stop − t_start)</pre>
<p>Predict the printed T at steps 10 … 100, then run.</p></div>
    <div class="row"><label class="lc-in">runs<select class="lc-u mr-r"><option value="one">run 100</option><option value="two">run 50 + run 50</option><option value="ss">run 50 + run 50, both start 0 stop 100</option></select></label><button class="btn sm primary" data-k="b">▶ Run</button><span class="lmp-status muted mr-stb"></span></div>
    <div class="tbl mr-outb"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: rerun a finished trajectory</h2>
    <div class="prose"><p><code>rerun</code> reads each frame of a dump and computes energy and pressure for it, with no time step. It calls <code>setup_minimal(1)</code> per frame (${gh('src/rerun.cpp', 170)}), so every frame gets a fresh neighbour list. Run 50 steps with a dump every 10, then rerun it:</p></div>
    <div class="row"><label class="lc-in">neighbours<select class="lc-u mr-nbc"><option value="melt">every 20 check no (melt)</option><option value="safe">every 1 check yes (default)</option></select></label><button class="btn sm primary" data-k="c">▶ Run, then rerun</button><span class="lmp-status muted mr-stc"></span></div>
    <div class="tbl mr-outc"></div>
    <div class="prose mr-sumc"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You equilibrate with <code>fix nvt temp 300 300 100</code> for one run, then want a 300 → 400 K ramp over the next three runs of 10 000 steps each. What do you write?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>fix 1 all nvt temp 300 400 100</code> once, then <code>run 10000 start S stop E</code> three times, with S = the step at the start of the first and E = S + 30000. Without start/stop, each run would ramp 300 → 400 on its own.</div></div>
    <div class="ex"><p><b>2.</b> A loop in your script does <code>run 100</code> 1000 times to print something between runs. Why is it slow, and what fixes it?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each run does setup: init, sort, a full neighbour build and a force call, 1000 times. <code>run 100 pre no post no</code> skips it after the first, and the trajectory stays bitwise the same as one long run (section 2).</div></div>
    <div class="ex"><p><b>3.</b> In section 4, why is the rerun's pressure at step 50 exactly the original, but not at step 0?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The dump has only x, y, z. The kinetic part of the pressure uses the velocities still in memory, which are those of step 50. At step 50 they are right; at step 0 the rerun uses T(50) instead of T(0) = 3. Dump vx vy vz and rerun with <code>dump x y z vx vy vz</code> to fix it.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.mr-sta'), nb = $('.mr-nb').value, sp = $('.mr-sp').value; b.disabled = true;
      st.textContent = 'Running run 100…'; const r1 = await runLammps(scriptA(nb, 'one'), { read: ['f.dump'] });
      if (r1.error || !r1.files?.['f.dump']) { b.disabled = false; return err(st, r1); }
      st.textContent = 'Running the split version…'; const r2 = await runLammps(scriptA(nb, sp), { read: ['f.dump'] }); b.disabled = false;
      if (r2.error || !r2.files?.['f.dump']) return err(st, r2);
      st.textContent = '';
      const t1 = rowsOf(r1), t2 = [], th = r2.lines.join('\n').split(/^\s*Step.*$/m).slice(1).map(s => s.trim().split('\n').filter(l => /^\s*\d+\s/.test(l)).map(l => l.trim().split(/\s+/)));
      th.forEach((blk, k) => blk.forEach(w => t2.push({ run: k + 1, Step: w[0], PotEng: w[1], TotEng: w[2] })));
      const a = parseDump(r1.files['f.dump']), c = parseDump(r2.files['f.dump']); let n = 0, d = 0, m = 0;
      for (const id in a) for (let k = 0; k < 6; k++) { n++; const u = a[id][k], v = c[id][k]; if (u !== v) d++; let dd = u - v; if (k < 3) dd = ((dd + L / 2) % L + L) % L - L / 2; m = Math.max(m, Math.abs(dd)); }
      const pe1 = Object.fromEntries(t1.map(r => [r.Step, r.PotEng]));
      $('.mr-outa').innerHTML = `<table>${head(['', 'step', 'PotEng (split)', 'PotEng (run 100)', 'same?'])}${t2.map(r => row([`run ${r.run}`, r.Step, `<span class="mono">${r.PotEng}</span>`, `<span class="mono">${pe1[r.Step]}</span>`, r.PotEng === pe1[r.Step] ? '<span class="m-ok">✓ bit for bit</span>' : `Δ = ${(parseFloat(r.PotEng) - parseFloat(pe1[r.Step])).toExponential(3)}`])).join('')}</table>`;
      const msg = d === 0 ? 'Every one of the 24 000 final positions and velocities is identical: with <code>pre no</code> there is no setup, so the second run just continues the loop.'
        : nb === 'melt' ? `At the start of run 2 (step 50), setup builds a fresh list and the energy jumps by ${(parseFloat(t2.find(r => r.run === 2).PotEng) - parseFloat(pe1['50'])).toFixed(4)}: those are pairs that came inside 2.5 after the list of step 40 was made, and that run 100 never saw. The trajectories then part: final positions differ by up to ${m.toExponential(2)}.`
        : `With a list that is always safe, the fresh list only changes the order of the sums (sorted atoms, new list): PotEng at step 50 differs in the last digits, and after 50 more steps positions differ by up to ${m.toExponential(2)}. Chaos will grow that, but it is round-off, not physics.`;
      $('.mr-suma').innerHTML = `<p>${d} of ${n} final numbers differ (positions compared modulo the box). ${msg}</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.mr-stb'), r = $('.mr-r').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptB(r)); b.disabled = false;
      if (res.error) return err(st, res); st.textContent = '';
      const blocks = res.lines.join('\n').split(/^\s*Step.*$/m).slice(1).map(s => s.trim().split('\n').filter(l => /^\s*\d+\s/.test(l)).map(l => l.trim().split(/\s+/)));
      const rows = [];
      blocks.forEach((blk, k) => {
        const begin = r === 'ss' ? 0 : k * (r === 'one' ? 100 : 50), end = r === 'ss' ? 100 : begin + (r === 'one' ? 100 : 50);
        blk.forEach(([s, T], j) => { s = +s; if (s === 0) return; if (j === 0 && k > 0) { rows.push(row([`run ${k + 1}`, s, 'setup line, no rescale', `<span class="mono">${T}</span>`, ''])); return; }
          let delta = s - begin; delta /= end - begin; const tt = 1.0 + delta * (2.0 - 1.0);
          rows.push(row([`run ${k + 1}`, s, `1 + (${s} − ${begin}) / ${end - begin} = ${+tt.toFixed(12)}`, `<span class="mono">${T}</span>`, Math.abs(tt - T) < 1e-12 ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'])); });
      });
      $('.mr-outb').innerHTML = `<table>${head(['', 'step', 'target by hand', 'LAMMPS T', ''])}${rows.join('')}</table><p class="muted">${r === 'two' ? 'Run 2 starts the ramp again: beginstep = 50, endstep = 100, so the target at step 60 is 1.2, not 1.6.' : r === 'ss' ? 'start 0 stop 100 gives both runs the window of one long run.' : 'One window from 0 to 100.'} The match is to 12 digits: T after rescaling is the target only up to round-off (page C.3).</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.mr-stc'), nb = $('.mr-nbc').value; b.disabled = true; st.textContent = 'Running, then rerunning…';
      const res = await runLammps(scriptC(nb)); b.disabled = false;
      if (res.error) return err(st, res); st.textContent = '';
      const blocks = res.lines.join('\n').split(/^\s*Step.*$/m).slice(1).map(s => s.trim().split('\n').filter(l => /^\s*\d+\s/.test(l)).map(l => l.trim().split(/\s+/)));
      const [run, rer] = blocks; if (!rer) return err(st, { error: 'no rerun output' });
      const V = L * L * L, T50 = parseFloat(run[run.length - 1][1]);
      $('.mr-outc').innerHTML = `<table>${head(['step', 'PotEng (run)', 'ΔPotEng (rerun − run)', 'ΔPress', 'ΔPress by hand: 11997 (T₅₀ − T)/3V'])}${run.map((w, i) => { const v = rer[i], dpe = parseFloat(v[2]) - parseFloat(w[2]), dp = parseFloat(v[3]) - parseFloat(w[3]), hand = 11997 * (T50 - parseFloat(w[1])) / 3 / V; const big = Math.abs(dpe) > 1e-6;
        return row([w[0], `<span class="mono">${w[2]}</span>`, big ? `<b>${dpe.toFixed(6)}</b>` : dpe === 0 ? '0' : dpe.toExponential(1), dp === 0 ? '0' : dp.toFixed(10), `${hand === 0 ? '0' : hand.toFixed(10)} ${Math.abs(hand - dp) < 1e-9 ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗ pairs missing</span>'}`]); }).join('')}</table>`;
      $('.mr-sumc').innerHTML = nb === 'melt'
        ? '<p>At step 50 the rerun finds the energy the run lost: the list from step 40 had already let some pairs slip inside the cutoff unseen. Steps 10–40 agree to round-off. A rerun with the same input is a cheap test of your neighbour settings.</p>'
        : '<p>With safe lists, the energies agree to round-off at every frame. The pressure differences are all kinetic: the rerun uses the velocities of step 50 for every frame, and your hand column explains them.</p>';
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
