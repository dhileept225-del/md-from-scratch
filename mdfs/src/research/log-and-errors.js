// Research · A.9: reading the log and the error messages. Every derived number in the end-of-run summary recomputed by hand; ten common mistakes reproduced live.
import { runLammps } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

const swap = (a, b) => s => s.replace(a, b);
const ERRS = [
  { key: 'e1', what: 'A typo in a command name', old: 'pair_style      lj/cut 2.5', now: 'pair_stlye      lj/cut 2.5', f: swap('pair_style      lj/cut 2.5', 'pair_stlye      lj/cut 2.5'),
    why: 'LAMMPS looked the first word up in its list of commands and did not find it. The message repeats the whole line, so the typo is easy to see.', fix: 'Correct the spelling.' },
  { key: 'e2', what: 'A style from a package that is not built in', old: 'pair_style      lj/cut 2.5', now: 'pair_style      tersoff', f: swap('pair_style      lj/cut 2.5', 'pair_style      tersoff'),
    why: 'The style exists, but in the MANYBODY package, which this binary was compiled without (page 0.4).', fix: 'Use a LAMMPS built with that package (cmake -D PKG_MANYBODY=on), or a different style.' },
  { key: 'e3', what: 'A timestep 10 times too large', old: '(default timestep 0.005)', now: 'timestep        0.05', f: swap('fix             1 all nve', 'fix             1 all nve\ntimestep        0.05'),
    why: 'Atoms move so far in one step that two of them overlap, the LJ repulsion throws them across the box, and they leave the region this processor tracks. LAMMPS counts the atoms and finds most missing.', fix: 'Reduce the timestep. "Lost atoms" almost always means the dynamics blew up: too large dt, overlapping atoms, or a wrong force field.' },
  { key: 'e4', what: 'pair_coeff missing', old: 'pair_coeff      1 1 1.0 1.0 2.5', now: '(deleted)', f: swap('pair_coeff      1 1 1.0 1.0 2.5\n', ''),
    why: 'Nothing is wrong when the line is read. Only at run, when the pair style checks its table of ε and σ, does it find an empty entry. The line after it adds a status note; the message itself, "All pair coeffs are not set", is the clue.', fix: 'Give a pair_coeff for every pair of types (or for i i, and let mixing fill i j).' },
  { key: 'e5', what: 'A word where a number is expected', old: 'pair_coeff      1 1 1.0 1.0 2.5', now: 'pair_coeff      1 1 1.0 one 2.5', f: swap('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 1 1.0 one 2.5'),
    why: 'The 4th word must be σ, a number. The message names the file and line of the C++ code that checked it.', fix: 'Write the number.' },
  { key: 'e6', what: 'A command with too few arguments', old: 'fix             1 all nve', now: 'fix             1 all nvt', f: swap('fix             1 all nve', 'fix             1 all nvt'),
    why: 'fix nvt needs a target temperature and a damping time (page D.10); nve needs nothing. "Illegal … command" means the arguments do not match the documentation page of that command.', fix: 'fix 1 all nvt temp 3.0 3.0 0.5' },
  { key: 'e7', what: 'A file that is not there', old: 'lattice … create_atoms (4 lines)', now: 'read_data       melt.data', f: s => s.replace(/lattice[\s\S]*?create_atoms    1 box/, 'read_data       melt.data'),
    why: 'LAMMPS looks for the file in the directory where it was started. On a real computer this is the commonest error of all: running lmp from the wrong folder.', fix: 'cd into the folder that holds the file (page 0.2), or give the full path.' },
  { key: 'e8', what: 'A warning, not an error', old: 'neigh_modify    every 20 delay 0 check no', now: 'neigh_modify    every 20 delay 0 check yes', f: swap('check no', 'check yes'),
    why: 'The run finishes. But at the end, every rebuild of the neighbour list was "dangerous": an atom had already moved more than half the skin before the list was rebuilt, so some pairs may have been missed (page A.7). LAMMPS does not stop you; you must read it.', fix: 'every 1 delay 0 check yes (the safe default) or a larger skin.', warn: /Dangerous builds = \d+/ },
  { key: 'e9', what: 'Commands in the wrong order', old: 'create_box before create_atoms', now: 'create_atoms moved above create_box', f: swap('create_box      1 box\ncreate_atoms    1 box', 'create_atoms    1 box\ncreate_box      1 box'),
    why: 'Atoms can only be placed in a box that exists (block 2 of page A.3 has an order).', fix: 'Put create_box first.' },
  { key: 'e10', what: 'A type that does not exist', old: 'pair_coeff      1 1 1.0 1.0 2.5', now: 'pair_coeff      1 2 1.0 1.0 2.5', f: swap('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 2 1.0 1.0 2.5'),
    why: 'create_box 1 made one atom type, so the only valid type number is 1.', fix: 'Use type numbers between 1 and the number given to create_box (or in the data file header).' },
  { key: 'e11', what: 'A mass missing', old: 'mass            1 1.0', now: '(deleted)', f: swap('mass            1 1.0\n', ''),
    why: 'Without a mass, F = ma cannot give an acceleration, and velocity create cannot compute a kinetic energy. So the error comes at the velocity line, the first command that needs masses.', fix: 'Give a mass for every type (or a Masses section in the data file).' },
];

function summary(lines) {
  const get = re => { const l = lines.find(x => re.test(x)); return l ? l.match(re) : null; };
  const loop = get(/^Loop time of ([\d.eE+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms/);
  const perf = get(/^Performance: ([\d.]+) tau\/day, ([\d.]+) timesteps\/s, ([\d.]+) Matom-step\/s/);
  const sec = {}; lines.forEach(l => { const m = l.match(/^(Pair|Neigh|Comm|Output|Modify|Other)\s*\|(.*)$/); if (m) { const c = m[2].split('|').map(x => x.trim()); sec[m[1]] = { avg: +c[1], pct: c[4] }; } });
  const neighs = get(/^Neighs:\s+(\d+) ave/), ave = get(/^Ave neighs\/atom = ([\d.]+)/), nloc = get(/^Nlocal:\s+(\d+) ave/), nghost = get(/^Nghost:\s+(\d+) ave/);
  if (!loop || !perf) return null;
  return { t: +loop[1], procs: +loop[2], steps: +loop[3], N: +loop[4], tauday: perf[1], tps: perf[2], mas: perf[3], sec, neighs: neighs && +neighs[1], ave: ave && ave[1], nloc: nloc && +nloc[1], nghost: nghost && +nghost[1] };
}
const fix = (x, d) => x.toFixed(d);
const ok = (a, b) => a === b ? '<span class="m-ok">✓</span>' : `<span class="m-bad">✗</span>`;

export default {
  id: 'log-and-errors',
  title: 'Reading the log and the error messages',
  minutes: 20,
  summary: 'Everything LAMMPS prints after a run, explained and recomputed by hand from the loop time and the neighbour count; then eleven common mistakes in the melt script, each run live with its real error message, what it means and how to fix it.',
  learn: [
    'The parts of a log: echo, setup, thermo table, timing summary, neighbour statistics',
    'Recompute the performance line, the timing percentages and the average neighbour count by hand',
    'Recognise the most common errors and warnings and know what each one means',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: what the log contains</h2>
    <p>Everything LAMMPS prints goes to the screen and to <code>log.lammps</code> in the directory where it ran. In order:</p>
    <div class="tbl wide"><table>${head(['Part', 'What it shows', 'Page'])}
${row(['banner', 'the version: <code>LAMMPS (22 Jul 2025 - Update 6)</code>', '0.3'])}
${row(['echo of commands', 'what each command created: lattice spacing, box, "Created 4000 atoms"', 'A.5'])}
${row(['setup', 'neighbour list settings, units, timestep, memory', 'A.7'])}
${row(['thermo table', 'one line every <code>thermo N</code> steps', 'A.3'])}
${row(['<code>Loop time</code>', 'wall-clock seconds for the run, and speed', 'below'])}
${row(['timing breakdown', 'where the time went: Pair, Neigh, Comm, Output, Modify, Other', 'C.8'])}
${row(['neighbour statistics', 'owned atoms, ghost atoms, pairs in the list, rebuilds', 'A.7, D.3'])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Observe: run melt and read the end of the log</h2>
    <p class="muted">The times are measured on your device, so your numbers will differ from anyone else's. The hand checks use your numbers.</p>
    <div class="row"><button class="btn sm primary" data-k="melt">▶ Run in.melt</button><span class="lmp-status muted lg-st"></span></div>
    <pre class="lmp-mini lg-tail" hidden></pre>
  </section>

  <section class="block prose">
    <h2>3 · Match: every derived number, by hand</h2>
    <p>Only three things are measured: the loop time, the time in each section, and the number of pairs. Everything else is arithmetic.</p>
    <div class="tbl wide"><table>${head(['Printed', 'Formula'])}
${row(['timesteps/s', 'steps ÷ loop time'])}
${row(['tau/day', 'timesteps/s × timestep × 86400 (seconds in a day). In real units: ns/day, with dt in fs'])}
${row(['Matom-step/s', 'timesteps/s × atoms ÷ 10⁶'])}
${row(['%total', '100 × section time ÷ loop time'])}
${row(['Ave neighs/atom', 'Neighs ÷ Nlocal'])}</table></div>
    <div class="lmp-compare lg-cmp"><p class="muted">Run melt above to fill in this table.</p></div>
    <p>Ave neighs/atom is 37.9, not the 39 counted around a perfect lattice on page A.7: the list was built at the last rebuild, when the crystal had already melted. Nghost (5506) is the number of copies of atoms near the box faces that this processor keeps so that every atom sees all its neighbours across the periodic boundary (page D.2).</p>
  </section>

  <section class="block">
    <h2>4 · Modify: eleven mistakes, run live</h2>
    <p class="muted">Each changes one line of in.melt. Predict what LAMMPS says, then run it.</p>
    ${ERRS.map(e => `<div class="exp" data-e="${e.key}">
      <p><b>${e.what}</b></p>
      <div class="diff"><span class="del">− ${esc(e.old)}</span>
<span class="add">+ ${esc(e.now)}</span></div>
      <div class="row"><button class="btn sm primary" data-err="${e.key}">▶ Run it</button><span class="lmp-status muted"></span></div>
      <div class="er-out" hidden><pre class="lmp-mini"></pre><p><b>What it means.</b> ${e.why}</p><p><b>Fix.</b> ${esc(e.fix)}</p></div>
    </div>`).join('')}
  </section>

  <section class="block prose">
    <h2>How to read any error</h2>
    <ul>
      <li><b>Read the first ERROR line, not the last line.</b> It names the command, and often the C++ file and line that raised it, e.g. <code>(src/fix_nh.cpp:61)</code>. That file is where the check lives (page B.8).</li>
      <li><b>Look at "Last input line".</b> It is the command LAMMPS was executing, which may be later than the line you changed: a missing mass is only noticed at the velocity line.</li>
      <li><b>Some errors appear only at <code>run</code>.</b> Missing pair coefficients are checked when the run starts, not when the script is read.</li>
      <li><b>Warnings do not stop the run.</b> Search the log for <code>WARNING</code> and for "Dangerous builds" after every new script.</li>
      <li><b>Follow the link.</b> Common errors end with "For more information see https://docs.lammps.org/err…", a manual page that explains that exact message. For the rest, search the message in the manual and the LAMMPS forum.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> A run of 1000 steps with 32000 atoms and timestep 0.005 takes 40 s. What does the Performance line say?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>1000/40 = <b>25 timesteps/s</b>; 25 × 0.005 × 86400 = <b>10800 tau/day</b>; 25 × 32000/10⁶ = <b>0.8 Matom-step/s</b>.</div></div>
    <div class="ex"><p><b>2.</b> Pair takes 85 % of the time in melt. If you change the cutoff from 2.5 to 5.0, roughly how much longer does the run take?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The number of neighbours grows with the volume of the cutoff sphere, (5.3/2.8)³ ≈ 6.8 times with the 0.3 skin. Pair time grows about that much, so the run is about 0.85 × 6.8 + 0.15 ≈ <b>6 times</b> longer. Try it on the Start page.</div></div>
    <div class="ex"><p><b>3.</b> A run stops with "Lost atoms" at step 3000, after running fine. What do you check first?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The thermo lines just before it: a temperature or pressure that rises sharply shows the blow-up. Then the timestep (too large for the hottest moment), close contacts in the starting structure, and the force field parameters.</div></div>
  </section>`,

  mount(root) {
    root.querySelector('[data-k=melt]').onclick = async e => {
      const b = e.currentTarget, st = root.querySelector('.lg-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(IN_MELT); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      const k = res.lines.findIndex(l => /^Loop time/.test(l)), tail = root.querySelector('.lg-tail');
      tail.textContent = res.lines.slice(k, res.lines.findIndex(l => /^Dangerous builds|^Neighbor list builds/.test(l)) + 2).join('\n'); tail.hidden = false;
      st.textContent = 'Done. The end of the log:';
      const S = summary(res.lines); if (!S) return;
      const tps = S.steps / S.t, dt = 0.005;
      const r = [
        ['timesteps/s', `${S.steps} ÷ ${S.t}`, fix(tps, 3), S.tps],
        ['tau/day', `${fix(tps, 3)} × ${dt} × 86400`, fix(tps * dt * 86400, 3), S.tauday],
        ['Matom-step/s', `${fix(tps, 3)} × ${S.N} ÷ 10⁶`, fix(tps * S.N / 1e6, 3), S.mas],
        ...['Pair', 'Neigh'].filter(n => S.sec[n]).map(n => [`${n} %total`, `100 × ${S.sec[n].avg} ÷ ${S.t}`, fix(100 * S.sec[n].avg / S.t, 2), S.sec[n].pct]),
        ...(S.neighs ? [['Ave neighs/atom', `${S.neighs} ÷ ${S.nloc}`, String(+(S.neighs / S.nloc).toFixed(3)), S.ave]] : []),
      ];
      root.querySelector('.lg-cmp').innerHTML = `<div class="tbl"><table><tr><th></th><th>by hand</th><th></th><th>LAMMPS</th><th></th></tr>${r.map(([n, f, h, l]) => `<tr><td>${n}</td><td class="mono">${f}</td><td class="mono">${h}</td><td class="mono">${l}</td><td>${ok(h, l)}</td></tr>`).join('')}</table></div>
        <p class="muted">The loop time is printed with 6 digits, so a hand value can differ from LAMMPS in its last digit: LAMMPS divides by the unrounded time.</p>`;
    };
    root.querySelectorAll('[data-err]').forEach(b => b.onclick = async () => {
      const E = ERRS.find(x => x.key === b.dataset.err), box = b.closest('.exp'), st = box.querySelector('.lmp-status'), out = box.querySelector('.er-out');
      b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(E.f(IN_MELT)); b.disabled = false; st.textContent = '';
      let msg;
      if (E.warn) msg = (res.error || 'No error: the run finished.') + '\n' + res.lines.filter(x => /^Neighbor list builds|^Dangerous builds/.test(x)).join('\n');
      else { const i = res.lines.findIndex(x => /^ERROR/.test(x)); if (i < 0) msg = res.error || 'No error.'; else { const j = res.lines.findIndex((x, k) => k > i && !x.trim()); msg = res.lines.slice(i, j < 0 ? i + 4 : Math.min(j, i + 4)).join('\n'); } }
      out.querySelector('pre').textContent = msg.trim(); out.hidden = false;
    });
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
