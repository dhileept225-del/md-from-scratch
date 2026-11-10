// Research · Part 0.2: the terminal in 10 commands, with a practice terminal that runs the real LAMMPS.
import { mountTerminal, IN_MELT } from './terminal.js';
import { handStep0 } from './lammps-start.js';
import { g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

const TASKS = [
  { cmd: 'pwd', what: 'Where am I? <b>pwd</b> prints the folder you are in.', predict: '/home/you (your home folder)' },
  { cmd: 'ls', what: 'What is here? <b>ls</b> lists files and folders. Folders end with /.', predict: 'one folder: lammps-course/' },
  { cmd: 'cd lammps-course/melt', what: 'Go into a folder with <b>cd</b>. It prints nothing; the prompt changes. Then type <b>ls</b> again.', predict: 'the prompt ends in ~/lammps-course/melt; ls shows in.melt' },
  { cmd: 'cat in.melt', what: 'Print a whole file with <b>cat</b>.', predict: 'the 33 lines of the melt script' },
  { cmd: 'grep -n pair in.melt', what: 'Find lines containing a word with <b>grep</b>; <b>-n</b> adds line numbers.', predict: 'two lines, numbers 14 and 15' },
  { cmd: 'less in.melt', what: '<b>less</b> shows a long file one screen at a time (q quits).', predict: 'the same text as cat' },
  { cmd: 'mkdir cold', what: 'Make a new folder with <b>mkdir</b>, copy the script into it with <b>cp</b>, and move in.', more: ['cp in.melt cold/', 'cd cold'], predict: 'no output; ls inside cold shows in.melt' },
  { cmd: 'nano in.melt', what: 'Edit with <b>nano</b>: change <code>create 3.0</code> to <code>create 1.0</code> on line 12, then Save and exit.', predict: 'the copy in cold/ now starts at T = 1.0; the original is untouched' },
  { cmd: 'lmp -in in.melt', what: 'Run LAMMPS on the edited copy. On your computer it may be <b>./lmp</b> (the program in this folder), <b>lmp</b>, <b>lmp_serial</b> or <b>lmp.exe</b>.', predict: 'step 0: E_pair −6.7733681 (positions unchanged), TotEng −5.2737431, Press −5.3913283' },
  { cmd: 'mpirun -np 4 lmp -in in.melt', what: 'Run in parallel with <b>mpirun</b>: <b>-np 4</b> starts 4 copies that share the atoms.', predict: 'the same numbers: splitting the work does not change the physics' },
];

export default {
  id: 'terminal-basics',
  title: 'The terminal in 10 commands',
  minutes: 20,
  summary: 'LAMMPS is used by typing commands, on your own computer and on every supercomputer. Ten commands are enough for everything in this course. You practise them here in a terminal that runs the real LAMMPS, edit the melt script, and check the result by hand.',
  learn: [
    'What a terminal is, and how a command is built from a name and its arguments',
    'The ten commands you need: ls, cd, cat, less, grep, mkdir, cp, nano, lmp and mpirun',
    'Edit an input script, run LAMMPS on it, and match the new numbers by hand',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: talking to the computer in sentences</h2>
    <p>A terminal is a window where you type instructions instead of clicking. LAMMPS has no buttons of its own (LAMMPS-GUI adds some), and supercomputers offer nothing but a terminal, so this is how almost everyone runs it.</p>
    <p>Every instruction is a short sentence: the <b>name of a program</b>, then <b>options</b> (starting with -), then what it should work on. Press Enter and the program runs.</p>
    <div class="calc">grep   -n        pair       in.melt
 │      │          │          └ the file to search
 │      │          └ what to look for
 │      └ option: show line numbers
 └ the program</div>
  </section>

  <section class="block prose">
    <h2>2 · The ten commands</h2>
    <div class="tbl wide"><table>${head(['Command', 'What it does', 'Example'])}
${row(['<code>ls</code>', 'list what is in the current folder', '<code>ls</code>'])}
${row(['<code>cd</code>', 'change folder; <code>cd ..</code> goes up one, <code>cd ~</code> goes home', '<code>cd lammps-course/melt</code>'])}
${row(['<code>cat</code>', 'print a whole file', '<code>cat in.melt</code>'])}
${row(['<code>less</code>', 'read a long file page by page (q quits)', '<code>less log.lammps</code>'])}
${row(['<code>grep</code>', 'show only the lines that contain a word', '<code>grep -n pair in.melt</code>'])}
${row(['<code>mkdir</code>', 'make a new folder', '<code>mkdir cold</code>'])}
${row(['<code>cp</code>', 'copy a file', '<code>cp in.melt cold/</code>'])}
${row(['<code>nano</code>', 'edit a text file', '<code>nano in.melt</code>'])}
${row(['<code>./lmp</code>', 'run a program; <code>./</code> means "the one in this folder"', '<code>./lmp -in in.melt</code>'])}
${row(['<code>mpirun</code>', 'run a program on several processors', '<code>mpirun -np 4 lmp -in in.melt</code>'])}</table></div>
    <p>One more is handy: <code>pwd</code> prints where you are. Folder paths work like addresses: <code>~</code> is your home folder, <code>..</code> the folder above, and <code>cold/in.melt</code> means "the file in.melt inside the folder cold".</p>
  </section>

  <section class="block">
    <h2>3 · Modify and observe: practise here</h2>
    <p class="muted">This terminal keeps its files in your browser only. The lmp command runs the real LAMMPS. Tap a command to type it, or type it yourself.</p>
    <div class="term-wrap"><div id="term"></div>
      <ol class="tasks">${TASKS.map((t, i) => `<li><p>${t.what}</p><div class="row">${[t.cmd, ...(t.more || [])].map(c => `<button class="btn sm ghost tcmd" data-cmd="${c}">${c}</button>`).join('')}</div><p class="muted small"><b>Predict:</b> ${t.predict}</p></li>`).join('')}</ol>
    </div>
  </section>

  <section class="block prose">
    <h2>4 · Match: the edited run by hand</h2>
    <p>Task 8 changed only the starting temperature. At step 0 the atoms still sit on the crystal, so the pair energy cannot change; only the kinetic parts do (the Start page, change D):</p>
    <div class="calc">E_pair = <b>−6.7733681</b>                      (positions unchanged)
KE/N   = 1.5 × T × (N−1)/N = 1.5 × 1.0 × 0.99975 = 1.4996250
TotEng = −6.7733681 + 1.4996250             = <b>−5.2737431</b>
Press  = 0.8442 × 1.0 × 0.99975 − 6.2353173 = <b>−5.3913283</b></div>
    <div class="term-match"><p class="muted">Run LAMMPS in the terminal (task 9) and its step-0 numbers appear here next to the hand values.</p></div>
  </section>

  <section class="block prose">
    <h2>5 · On your own computer</h2>
    <div class="tbl wide"><table>${head(['System', 'Where the terminal is', 'Note'])}
${row(['Linux', 'the Terminal application (Ctrl+Alt+T on Ubuntu)', 'all ten commands work as here'])}
${row(['macOS', 'Applications → Utilities → Terminal', 'all ten commands work as here'])}
${row(['Windows', 'PowerShell or Windows Terminal; or Ubuntu inside Windows (WSL)', 'in PowerShell use <code>type</code> instead of cat and <code>Select-String</code> instead of grep; WSL behaves like Linux'])}</table></div>
    <p>Two habits save hours: press <b>Tab</b> to complete a file name, and press <b>↑</b> to bring back the last command. ↑ also works in the practice terminal above.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You are in <code>~/lammps-course/melt/cold</code>. Where does <code>cd ../..</code> take you?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Up two folders: <code>~/lammps-course</code>. Check with <code>pwd</code>.</div></div>
    <div class="ex"><p><b>2.</b> How many lines of in.melt contain "create", and what does <code>grep -c create in.melt</code> print?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Three: create_box (line 8), create_atoms (line 9) and velocity all create (line 12). <code>-c</code> prints only the count: <b>3</b>.</div></div>
    <div class="ex"><p><b>3.</b> You type <code>lmp -in in.melt</code> in your home folder and get "Cannot open input script in.melt". Why, and what are two fixes?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>LAMMPS looks for the file in the current folder, and in.melt is not in your home folder. Either <code>cd lammps-course/melt</code> first, or give the path: <code>lmp -in lammps-course/melt/in.melt</code>.</div></div>
  </section>`,

  mount(root) {
    const match = root.querySelector('.term-match');
    const term = mountTerminal(root.querySelector('#term'), {
      files: { 'lammps-course/melt/in.melt': IN_MELT },
      onRun: ({ script, res }) => {
        const r = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols, h = handStep0(script);
        if (!r || !h.ok) { match.innerHTML = `<p class="muted">${h.ok ? 'No step-0 line in the output.' : h.why}</p>`; return; }
        const L = k => r[c.indexOf(k)];
        const line = (k, hv) => { const m = matchMarkup(g8(hv), L(k)); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${L(k)}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; };
        match.innerHTML = `<div class="tbl"><table><tr><th>step 0, T = ${h.T}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${line('E_pair', h.E)}${line('TotEng', h.TotEng)}${line('Press', h.Press)}</table></div>`;
      },
    });
    root.querySelectorAll('.tcmd').forEach(b => b.onclick = () => term.type(b.dataset.cmd));
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
