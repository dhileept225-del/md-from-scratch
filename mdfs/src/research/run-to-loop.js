// Research · C.1: from `run 250` to the loop. What the run command does before the first step, predicted (thermo lines, neighbour builds) and matched for two runs in a row.
import { runLammps } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;

const TRACE = [
  { what: 'input.cpp finds the word run and calls the Run command', where: ['src/run.cpp', 37], code: `37  void Run::command(int narg, char **arg)
64      if (strcmp(arg[iarg],"upto") == 0) {      // optional keywords
77      } else if (strcmp(arg[iarg],"pre") == 0) {
108     nsteps = static_cast<int> (nsteps_input);  // run 250 → 250
113     nsteps = static_cast<int> (delta);         // run 400 upto → 400 − current step` },
  { what: 'Write down where this run starts and ends', where: ['src/run.cpp', 157], code: `157    update->nsteps = nsteps;
158    update->firststep = update->ntimestep;
159    update->laststep = update->ntimestep + nsteps;` },
  { what: 'Initialise every member (page B.2): cutoffs, neighbour settings, unit checks, fixes', where: ['src/run.cpp', 168], code: `168    if (preflag || update->first_update == 0) {
169      lmp->init();
170      update->integrate->setup(1);
171    } else output->setup(0);          // "pre no": skip both` },
  { what: 'Setup, before step 1: wrap atoms, share them, build the list, compute forces once', where: ['src/verlet.cpp', 109], code: `118  domain->pbc();               // put every atom back inside the box
122  comm->exchange();            // hand atoms to their owner (page B.7)
124  comm->borders();             // make ghosts (page B.6)
129  neighbor->build(1);          // first neighbour list
131  neighbor->ncalls = 0;        // this build is not counted
137  force_clear();
140  force->pair->compute(eflag,vflag);   // forces at the starting step
160  output->setup(flag);         // prints the header and the first thermo line` },
  { what: 'The loop itself: N steps (page C.2)', where: ['src/run.cpp', 175], code: `173    timer->init();
175    update->integrate->run(nsteps);` },
  { what: 'Clean up and print the summary (page A.9)', where: ['src/run.cpp', 178], code: `178    update->integrate->cleanup();
180    Finish finish(lmp);          // Loop time, timing table, neighbour statistics` },
];

// the rules the code follows, written out
function predict(N1, T, E, N2, mode) {
  const out = []; let step = 0, ago = 0, printed = 0;
  for (const [Nin, m] of [[N1, 'normal'], [N2, mode]]) {
    const N = m === 'upto' ? Nin - step : Nin, first = step, last = step + N, pre = m !== 'preno';
    if (N < 0) return { error: `run ${Nin} upto: ${Nin} is before the current step ${step}` };
    let builds = 0; if (pre) { ago = 0; printed = 0; }
    for (let s = first + 1; s <= last; s++) { ago++; if (ago % E === 0) { builds++; ago = 0; } }
    printed += builds;
    const th = [first]; for (let n = Math.floor(first / T) * T + T; n < last; n += T) th.push(n); if (last !== first) th.push(last);
    out.push({ N, first, last, pre, th, builds, printed, forces: N + (pre ? 1 : 0) });
    step = last;
  }
  return { runs: out };
}

export default {
  id: 'run-to-loop',
  title: 'From run 250 to the loop',
  minutes: 15,
  summary: 'Typing run 250 does much more than 250 steps: LAMMPS initialises every member, wraps and shares the atoms, builds the first neighbour list and computes the forces once before step 1, and prints a summary after the last. You follow the six stops in the code, then predict for two runs in a row which steps get a thermo line and how many neighbour builds the log reports, including the surprises of pre no and upto, and check them against LAMMPS.',
  learn: [
    'The path from the run command to Verlet::run(): parse, init, setup, loop, finish',
    'Why the first thermo line comes from setup, and forces are computed N + 1 times',
    'Predict thermo lines and neighbour builds for consecutive runs',
    'What run … upto and run … pre no change',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: six stops from run to the loop</h2>
    <div class="stepper st3"><div class="row"><button class="btn sm ghost" data-t="-1">◀</button><span class="st3-n mono"></span><button class="btn sm ghost" data-t="1">▶</button></div>
      <p class="st3-when"></p><pre class="lmp-mini st3-code"></pre></div>
    <div class="prose"><p>Three consequences you can see in every log:</p>
    <ul>
      <li>The thermo line at the starting step is printed by <b>setup</b>, not by the loop. A run of N steps computes the forces N + 1 times.</li>
      <li>The neighbour list built in setup is not counted: "Neighbor list builds" counts only the loop's builds.</li>
      <li>After the loop, thermo always prints the last step, whatever the <code>thermo</code> interval (${gh('src/output.cpp', 215)}).</li>
    </ul></div>
  </section>

  <section class="block">
    <h2>2 · Predict: two runs in a row</h2>
    <div class="prose"><p>in.melt with your settings, then a second <code>run</code>. The rules from the code: thermo prints the first step, every multiple of the interval, and the last step; the list is rebuilt when the steps since the last build reach a multiple of <code>every</code> (${gh('src/neighbor.cpp', 2373)}: <code>ago++; if (ago % every == 0)</code>).</p></div>
    <div class="lc-in">
      <label>first run N₁ <input type="number" class="lmp-in sm" data-c="N1" value="250" min="1" max="500"></label>
      <label>thermo every <input type="number" class="lmp-in sm" data-c="T" value="50" min="1" max="500"></label>
      <label>neigh_modify every <input type="number" class="lmp-in sm" data-c="E" value="20" min="1" max="100"></label>
      <label>second run N₂ <input type="number" class="lmp-in sm" data-c="N2" value="100" min="1" max="800"></label>
      <label>second run is <select class="lc-u" data-c="mode"><option value="normal">run N₂</option><option value="upto">run N₂ upto</option><option value="preno">run N₂ pre no</option></select></label>
    </div>
    <pre class="lmp-mini c1-script"></pre>
    <div class="tbl c1-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run both with LAMMPS</button><span class="lmp-status muted c1-st"></span></div>
    <div class="tbl c1-cmp"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe: what upto and pre no change</h2>
    <ul>
      <li><b>run N upto</b>: N is the step to stop at, not the number of steps (${gh('src/run.cpp', 113)}). After run 250, <code>run 330 upto</code> runs 80 steps. Useful in scripts that are restarted: the end step stays the same.</li>
      <li><b>run N pre no</b>: skips <code>lmp->init()</code> and setup (${gh('src/run.cpp', 168)}). No new neighbour list, no force recomputation, and the build counter is not reset, so the log prints the builds of <b>both</b> runs added together. It saves time when a script runs many short runs in a loop, but only if nothing was changed between them: a changed pair_coeff or fix would be ignored without init.</li>
      <li>Between two normal runs the step and the time carry on (page B.2), and the thermo line of the shared step is printed twice: once as the last line of run 1, once by setup of run 2.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> <code>thermo 100</code>, <code>run 250</code>. Which steps are printed?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><b>0, 100, 200, 250</b>: 0 by setup, 100 and 200 as multiples, 250 because the last step is always printed.</div></div>
    <div class="ex"><p><b>2.</b> <code>neigh_modify every 20 check no</code>, <code>run 250</code>, then <code>run 30 pre no</code>. What does the second run's log say for "Neighbor list builds"?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Run 1 builds at steps 20, 40 … 240: 12, and the last build was 10 steps before its end. Run 2 continues the count: steps 260 and 280 are 20 and 40 steps after 240, so 2 more. Without a reset the log says <b>14</b>. Check it above.</div></div>
    <div class="ex"><p><b>3.</b> How many times are forces computed by <code>run 0</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><b>Once</b>, in setup. That is why <code>run 0</code> is the standard way to get the energy and pressure of a structure without moving it (pages A.5, A.8).</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let t = 0;
    const showT = () => { const T = TRACE[t]; $('.st3-n').textContent = `${t + 1} / ${TRACE.length}`; $('.st3-when').innerHTML = `${esc(T.what)} · ${gh(...T.where)}`; $('.st3-code').textContent = T.code; };
    root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => { t = (t + +b.dataset.t + TRACE.length) % TRACE.length; showT(); }); showT();

    let C = null;
    const v = k => root.querySelector(`[data-c=${k}]`).value;
    const int = (k, d, lo, hi) => { const x = parseInt(v(k)); return isFinite(x) ? Math.min(hi, Math.max(lo, x)) : d; };
    const script = c => IN_MELT.replace('thermo          50', `thermo          ${c.T}`).replace('every 20', `every ${c.E}`).replace('run             250', `run             ${c.N1}\nrun             ${c.N2}${c.mode === 'upto' ? ' upto' : c.mode === 'preno' ? ' pre no' : ''}`);
    const draw = () => {
      C = { N1: int('N1', 250, 1, 500), T: int('T', 50, 1, 500), E: int('E', 20, 1, 100), N2: int('N2', 100, 1, 800), mode: v('mode') };
      $('.c1-script').textContent = script(C).split('\n').filter(l => /^(thermo|neigh_modify|run)/.test(l)).join('\n');
      const P = predict(C.N1, C.T, C.E, C.N2, C.mode); C.P = P;
      if (P.error) { $('.c1-hand').innerHTML = `<p><span class="m-bad">${esc(P.error)}</span>: LAMMPS will stop with "Invalid run command upto value".</p>`; $('.c1-cmp').innerHTML = ''; return; }
      $('.c1-hand').innerHTML = `<table>${head(['by hand', 'steps', 'thermo lines at steps', 'forces computed', 'log: Neighbor list builds'])}${P.runs.map((r, i) => row([`run ${i + 1}${r.pre ? '' : ' (pre no)'}`, `${r.first} → ${r.last}`, r.th.join(', '), r.forces, r.pre ? r.builds : `${r.printed} (${r.printed - r.builds} + ${r.builds}: not reset)`])).join('')}</table>`;
      $('.c1-cmp').innerHTML = ''; $('.c1-st').textContent = '';
    };
    root.querySelectorAll('[data-c]').forEach(i => { i.oninput = draw; i.onchange = draw; }); draw();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.c1-st'), c = C; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(script(c)); b.disabled = false; st.textContent = '';
      const got = []; let th = null;
      for (const l of res.lines) {
        if (/^\s+Step\s/.test(l)) { th = []; continue; }
        if (th && /^\s+\d+\s/.test(l)) { th.push(+l.trim().split(/\s+/)[0]); continue; }
        const m = l.match(/^Neighbor list builds = (\d+)/); if (m) { got.push({ th: th || [], builds: +m[1] }); th = null; }
      }
      const errLine = res.lines.find(l => /^ERROR/.test(l));
      if (c.P.error) { $('.c1-cmp').innerHTML = `<p class="mono">${esc(errLine || res.error || 'no error?')}</p>`; return; }
      const ok = x => x ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>';
      $('.c1-cmp').innerHTML = `<table>${head(['LAMMPS', 'thermo lines at steps', '', 'Neighbor list builds', ''])}${c.P.runs.map((r, i) => { const g = got[i] || { th: [], builds: '—' }; return row([`run ${i + 1}`, g.th.join(', '), ok(g.th.join() === r.th.join()), g.builds, ok(g.builds === r.printed)]); }).join('')}</table>${errLine ? `<p class="mono">${esc(errLine)}</p>` : ''}`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
