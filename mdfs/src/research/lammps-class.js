// Research · B.2: the conductor. One LAMMPS object owns 14 members; every script line writes into one of them; their bookkeeping checked by hand.
import { runLammps } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const g15 = x => String(Number(x.toPrecision(15)));

const MEMBERS = [
  ['memory', 'Memory', 'allocates and frees all big arrays'],
  ['error', 'Error', 'prints ERROR and WARNING messages and stops'],
  ['universe', 'Universe', 'all processors, and partitions for multi-replica runs'],
  ['input', 'Input', 'reads the script line by line and calls the right code'],
  ['atom', 'Atom', 'every per-atom array: x, v, f, type, mass, charge …'],
  ['update', 'Update', 'timestep, step number, units, the integrator (run)'],
  ['neighbor', 'Neighbor', 'neighbour lists and the skin'],
  ['comm', 'Comm', 'sends atoms between processors; ghost atoms'],
  ['domain', 'Domain', 'the box, periodic boundaries, lattice, regions'],
  ['force', 'Force', 'pair, bond, angle, dihedral, improper, kspace styles; unit constants'],
  ['modify', 'Modify', 'the lists of fixes and computes'],
  ['group', 'Group', 'named groups of atoms ("all" exists from the start)'],
  ['output', 'Output', 'thermo, dumps and restart files'],
  ['timer', 'Timer', 'stopwatches for the timing table (page A.9)'],
];
// which member(s) each line of in.melt writes into
const WHO = [
  ['units lj', [['update', 'unit_style = "lj", dt = 0.005'], ['force', 'boltz = 1, mvv2e = 1 … (page A.4)'], ['neighbor', 'skin = 0.3']]],
  ['atom_style atomic', [['atom', 'avec = AtomVecAtomic: which per-atom arrays exist']]],
  ['lattice fcc 0.8442', [['domain', 'lattice: a = (4/0.8442)^⅓ = 1.6796']]],
  ['region box block 0 10 0 10 0 10', [['domain', 'a region called "box"']]],
  ['create_box 1 box', [['domain', 'box 0 … 16.795962 in x, y, z; box_exist = 1'], ['atom', 'ntypes = 1']]],
  ['create_atoms 1 box', [['atom', 'x[i] for 4000 atoms, natoms = 4000']]],
  ['mass 1 1.0', [['atom', 'mass[1] = 1.0']]],
  ['velocity all create 3.0 87287 loop geom', [['atom', 'v[i] for every atom'], ['group', 'uses group "all"']]],
  ['pair_style lj/cut 2.5', [['force', 'pair = new PairLJCut, cut_global = 2.5']]],
  ['pair_coeff 1 1 1.0 1.0 2.5', [['force', 'pair->epsilon[1][1] = 1, sigma[1][1] = 1']]],
  ['neighbor 0.3 bin', [['neighbor', 'skin = 0.3, binning method']]],
  ['neigh_modify every 20 delay 0 check no', [['neighbor', 'every = 20, dist_check = 0']]],
  ['fix 1 all nve', [['modify', 'fix list: [FixNVE "1"]']]],
  ['thermo 50', [['output', 'thermo_every = 50']]],
  ['run 250', [['update', 'nsteps = 250; integrate->run() conducts every other member each step (page C.2)']]],
];

const UNITS = [['lj', 0.005, 0.3, 167], ['real', 1.0, 2.0, 188], ['metal', 0.001, 2.0, 209], ['si', 1e-8, 0.001, 230], ['cgs', 1e-8, 0.1, 251], ['electron', 0.001, 2.0, 272], ['micro', 2.0, 0.1, 293], ['nano', 0.00045, 0.1, 314]];

const BOOK = p => IN_MELT.replace('neighbor        0.3 bin', `neighbor        ${p.skin} bin`).replace('pair_style      lj/cut 2.5', `pair_style      lj/cut ${p.rc}`).replace('pair_coeff      1 1 1.0 1.0 2.5', `pair_coeff      1 1 1.0 1.0 ${p.rc}`)
  .replace('thermo          50', 'thermo          1000').replace('run             250', `run             ${p.n1}
timestep        ${p.dt2}
run             ${p.n2}
print           "time = $(time:%.15g)"
info            system communication`);

export default {
  id: 'lammps-class',
  title: 'The conductor: the LAMMPS class',
  minutes: 20,
  summary: 'The whole program is one object of class LAMMPS. It owns fourteen members (atom, update, domain, force, neighbor …), and every line of an input script writes into one of them. You see which member each line of in.melt changes, then predict and check the bookkeeping that the members keep: step number, timestep, simulation time, ghost cutoff, and what clear throws away.',
  learn: [
    'main() is three lines: create LAMMPS, read the script, delete LAMMPS',
    'The fourteen members, what each one keeps, and the order they are created in',
    'Which member each line of a script writes into',
    'Predict step, dt, simulation time and the communication cutoff by hand, and what clear resets',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: the whole program in three lines</h2>
    <p>When you type <code>lmp -in in.melt</code>, the function <code>main()</code> in ${gh('src/main.cpp', 69)} does:</p>
<pre class="lmp-mini">69  auto *lammps = new LAMMPS(argc, argv, lammps_comm);   // build the conductor
70  lammps->input->file();                                // read the script, line by line
71  delete lammps;                                        // clean up and finish</pre>
    <p>Everything else happens inside that one object. Its constructor builds the members in a fixed order (${gh('src/lammps.cpp', 844)}), and the comments in the code say why the order matters:</p>
<pre class="lmp-mini">852  comm = new CommBrick(this);       // Comm before Atom, so threads are known
855  neighbor = new Neighbor(this);
861  domain = new Domain(this);
865  atom = new Atom(this);
873  group = new Group(this);
875  force = new Force(this);          // after group, to create temperature
878  modify = new Modify(this);
880  output = new Output(this);        // after group, so "all" exists
882  update = new Update(this);        // after output, force, neighbor
883  timer = new Timer(this);</pre>
    <p><code>memory</code>, <code>error</code>, <code>universe</code> and <code>input</code> are made just before. Every class in LAMMPS inherits from <code>Pointers</code> (${gh('src/pointers.h', 110)}), which gives it a short name for each member. That is why code anywhere can simply write <code>update->dt</code> or <code>atom->x</code>.</p>
    <div class="tbl wide"><table>${head(['Member', 'Class', 'Keeps'])}${MEMBERS.map(([m, c, k]) => row([`<code>${m}</code>`, c, k])).join('')}</table></div>
  </section>

  <section class="block">
    <h2>2 · Which member does each line change?</h2>
    <p class="muted">Tap a line of in.melt. The members it writes into light up.</p>
    <div class="lc-grid">
      <div class="lc-lines">${WHO.map(([l], i) => `<button class="lc-line" data-w="${i}"><code>${esc(l)}</code></button>`).join('')}</div>
      <div><div class="lc-mem">${MEMBERS.map(([m]) => `<div class="lc-m" data-m="${m}">${m}</div>`).join('')}</div><div class="lc-what"></div></div>
    </div>
  </section>

  <section class="block">
    <h2>3 · Modify: the bookkeeping, predicted by hand</h2>
    <div class="prose">
      <p>The members keep running totals. Two of them are easy to get wrong. The <b>simulation time</b> is not step × dt: <code>update</code> adds each stretch of steps with the timestep that was in force, every time the timestep changes and at the end of every run (${gh('src/update.cpp', 552)}):</p>
<pre class="lmp-mini">552  void Update::update_time()
554    atime += (ntimestep - atimestep) * dt;   // time of the steps since the last update
555    atimestep = ntimestep;</pre>
      <p>The <b>communication cutoff</b> is how far beyond its own region a processor keeps copies (ghosts) of atoms: the pair cutoff plus the skin, because the neighbour list reaches that far (${gh('src/neighbor.cpp', 379)}: <code>cut = cutoff + delta</code>, delta = skin).</p>
    </div>
    <div class="lc-in">
      <label>first run <input type="number" class="lmp-in sm" data-p="n1" value="100" min="0" max="300" step="10"></label>
      <label>then timestep <input type="number" class="lmp-in sm" data-p="dt2" value="0.01" min="0.001" max="0.01" step="0.001"></label>
      <label>second run <input type="number" class="lmp-in sm" data-p="n2" value="50" min="0" max="300" step="10"></label>
      <label>skin <input type="number" class="lmp-in sm" data-p="skin" value="0.3" min="0" max="1" step="0.1"></label>
      <label>pair cutoff <input type="number" class="lmp-in sm" data-p="rc" value="2.5" min="1.5" max="3.5" step="0.1"></label>
    </div>
    <pre class="lmp-mini lc-script"></pre>
    <div class="calc lc-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run and ask the members</button><span class="lmp-status muted lc-st"></span></div>
    <div class="lmp-compare lc-cmp"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe: clear rebuilds the conductor</h2>
    <p><code>clear</code> is three lines of ${gh('src/input.cpp', 886)}: <code>lmp->destroy(); lmp->create(); lmp->post_create();</code>. Every member is thrown away and built fresh, as if LAMMPS had just started. (This page runs <code>clear</code> before every script for exactly that reason.) The new <code>update</code> sets the defaults of the chosen units (${gh('src/update.cpp', 149)}):</p>
    <div class="tbl"><table>${head(['units', 'default dt', 'default skin', 'code'])}${UNITS.map(([u, dt, sk, l]) => row([`<code>${u}</code>`, dt, sk, gh('src/update.cpp', l)])).join('')}</table></div>
    <div class="row"><label class="lab" style="margin:0">after clear, units</label><select class="lc-u">${UNITS.map(u => `<option>${u[0]}</option>`).join('')}</select><button class="btn sm primary" data-k="clear">▶ Run</button></div>
    <pre class="lmp-mini lc-clear-script"></pre>
    <div class="calc lc-clear-hand"></div>
    <pre class="lmp-mini lc-clear-out" hidden></pre>
  </section>

  <section class="block prose">
    <h2>What the conductor shows</h2>
    <ul>
      <li>An input script is a list of instructions to fourteen members. Knowing which member a command talks to tells you where its code is (page B.1) and what it can change.</li>
      <li>Members remember: step count and simulation time carry over from one run to the next, until <code>reset_timestep</code> or <code>clear</code>.</li>
      <li>Some commands check other members first. <code>units</code> asks <code>domain</code> whether a box exists and refuses if it does (${gh('src/input.cpp', 2011)}): the box was measured in the old units.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> run 200 (dt 0.005), then <code>timestep 0.002</code>, run 500. What are step and time at the end?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Step = 700. Time = 200 × 0.005 + 500 × 0.002 = 1 + 1 = <b>2</b>, not 700 × 0.002 = 1.4.</div></div>
    <div class="ex"><p><b>2.</b> After that, <code>reset_timestep 0</code>. What is the time now?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><b>0</b>. When the new step is smaller than the old one, <code>Update::reset_timestep</code> sets <code>atime = 0</code> (${gh('src/update.cpp', 508)}). Moving the step forward (e.g. reset_timestep 1000) would instead add the missing steps at the current dt.</div></div>
    <div class="ex"><p><b>3.</b> Which member answers <code>thermo_style custom vol</code>, and which one answers <code>temp</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>vol</code> is <code>domain->xprd × yprd × zprd</code>. <code>temp</code> is a compute (compute temp, made by output for thermo), kept in <code>modify</code>, which reads <code>atom->v</code> and <code>atom->mass</code>.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    root.querySelectorAll('[data-w]').forEach(b => b.onclick = () => {
      const w = WHO[+b.dataset.w][1], on = new Set(w.map(x => x[0]));
      root.querySelectorAll('[data-w]').forEach(x => x.classList.toggle('on', x === b));
      root.querySelectorAll('.lc-m').forEach(m => m.classList.toggle('on', on.has(m.dataset.m)));
      $('.lc-what').innerHTML = w.map(([m, t]) => `<p><code>${m}</code>: ${esc(t)}</p>`).join('');
    });
    root.querySelector('[data-w="0"]').click();

    const P = () => { const p = {}; root.querySelectorAll('[data-p]').forEach(i => { p[i.dataset.p] = i.value.trim(); }); return p; };
    let H = null;
    const draw = () => {
      const p = P(), n1 = parseInt(p.n1) || 0, n2 = parseInt(p.n2) || 0, dt2 = parseFloat(p.dt2), skin = parseFloat(p.skin), rc = parseFloat(p.rc);
      $('.lc-script').textContent = BOOK(p).split('\n').filter(l => /^(neighbor |pair_style|run|timestep|print|info)/.test(l)).join('\n');
      let atime = 0; atime += n1 * 0.005; atime += n2 * dt2;
      const cut = Math.sqrt(rc * rc) + skin;
      H = { step: String(n1 + n2), dt: String(dt2), time: g15(atime), cut: String(cut), p };
      $('.lc-hand').innerHTML = esc(`step   = ${n1} + ${n2} = `) + `<b>${H.step}</b>` + esc(`
dt     = the last timestep = `) + `<b>${H.dt}</b>` + esc(`
time   = ${n1} × 0.005 + ${n2} × ${dt2} = `) + `<b>${H.time}</b>` + esc(`
         (step × dt would give ${g15((n1 + n2) * dt2)})
ghosts = cutoff + skin = ${rc} + ${skin} = `) + `<b>${H.cut}</b>` + esc(`
atoms  = 4 × 10³ = `) + `<b>4000</b>`;
      $('.lc-cmp').innerHTML = '';
    };
    root.querySelectorAll('[data-p]').forEach(i => { i.oninput = draw; }); draw();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.lc-st'), p = H.p; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(BOOK(p)); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const get = re => { const l = res.lines.find(x => re.test(x)); return l ? l.match(re)[1] : '—'; };
      const L = { step: get(/^Current timestep number = (\S+)/), dt: get(/^Current timestep size = (\S+)/), time: get(/^time = (\S+)/), cut: get(/^Communication cutoff = (\S+)/), atoms: get(/^Atoms\s+=\s+(\d+)/) };
      const rows = [['update->ntimestep', H.step, L.step], ['update->dt', H.dt, L.dt], ['update->atime', H.time, L.time], ['comm cutoff', H.cut, L.cut], ['atom->natoms', '4000', L.atoms]];
      $('.lc-cmp').innerHTML = `<div class="tbl"><table><tr><th>member</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${rows.map(([k, h, l]) => `<tr><td><code>${k}</code></td><td class="mono">${h}</td><td class="mono">${l}</td><td>${h === l ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`).join('')}</table></div>`;
    };
    const cs = u => `run             100\nclear\nunits           ${u}\ninfo            system`;
    const drawClear = () => {
      const u = $('.lc-u').value, U = UNITS.find(x => x[0] === u);
      $('.lc-clear-script').textContent = '(in.melt up to run, then)\n' + cs(u);
      $('.lc-clear-hand').innerHTML = esc(`Prediction: step = `) + '<b>0</b>' + esc(`, dt = `) + `<b>${U[1]}</b>` + esc(`, atoms = `) + '<b>0</b>' + esc(`, and no box.`);
      $('.lc-clear-out').hidden = true;
    };
    $('.lc-u').onchange = drawClear; drawClear();
    $('[data-k=clear]').onclick = async e => {
      const b = e.currentTarget; b.disabled = true;
      const res = await runLammps(IN_MELT.replace('run             250', cs($('.lc-u').value))); b.disabled = false;
      const o = $('.lc-clear-out'), k = res.lines.findIndex(x => /^System information/.test(x));
      o.textContent = res.error && k < 0 ? res.error : res.lines.slice(k, k + 14).filter(x => x.trim() && !/^Info-Info/.test(x)).join('\n'); o.hidden = false;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
