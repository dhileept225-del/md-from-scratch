// Research · C.2: the loop, line by line. Verlet::run() in six stages; the two half-kicks and the drift of fix nve recomputed for every atom and matched bit for bit.
import { runLammps } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const L = 10 * 1.6795961913825073;

const STAGES = [
  { name: '1 · Next step', where: ['src/verlet.cpp', 252], does: 'Add 1 to the step counter. ev_set decides whether this step must also add up energies and the virial (only on thermo steps).', code: `246  for (int i = 0; i < n; i++) {
252    ntimestep = ++update->ntimestep;
253    ev_set(ntimestep);` },
  { name: '2 · First half: kick and drift', where: ['src/verlet.cpp', 258], does: 'Every fix with initial_integrate runs. For fix nve (fix_nve.cpp:100): half a kick with the old force, then a whole drift.', code: `258    modify->initial_integrate(vflag);
         ↓ fix_nve.cpp
100        dtfm = dtf / mass[type[i]];       // dtf = ½ · dt · ftm2v
101        v[i][0] += dtfm * f[i][0];        // v(t + ½dt) = v(t) + ½ dt F(t)/m
104        x[i][0] += dtv * v[i][0];         // x(t + dt) = x(t) + dt v(t + ½dt)` },
  { name: '3 · Neighbours: reuse or rebuild', where: ['src/verlet.cpp', 264], does: 'Most steps only send the new positions of ghosts (forward_comm). On a rebuild step: wrap atoms into the box, hand them to their owners, make ghosts, build the list.', code: `264    nflag = neighbor->decide();             // page C.1
266    if (nflag == 0) {
268      comm->forward_comm();                 // ghosts get the new x
270    } else {
277      domain->pbc();                        // wrap into the box (only here!)
284      comm->exchange();
285      if (sortflag && ntimestep >= atom->nextsort) atom->sort();   // page B.6
286      comm->borders();
293      neighbor->build(1);
299    }` },
  { name: '4 · Forces', where: ['src/verlet.cpp', 306], does: 'Zero all forces, then each force style adds its share: pair first, then bonded, then long-range. With newton on, forces on ghosts go back to their owners.', code: `306    force_clear();
316      force->pair->compute(eflag,vflag);
321      if (force->bond) force->bond->compute(eflag,vflag);
322      if (force->angle) force->angle->compute(eflag,vflag);
323      if (force->dihedral) force->dihedral->compute(eflag,vflag);
324      if (force->improper) force->improper->compute(eflag,vflag);
329      force->kspace->compute(eflag,vflag);
341      comm->reverse_comm();` },
  { name: '5 · Second half: kick', where: ['src/verlet.cpp', 347], does: 'Fixes that change forces go first (post_force), then the second half-kick with the new force, then fixes that act at the end of a step (thermostats like fix temp/rescale, averages).', code: `347    if (n_post_force_any) modify->post_force(vflag);
348    modify->final_integrate();
         ↓ fix_nve.cpp
141        v[i][0] += dtfm * f[i][0];        // v(t + dt) = v(t + ½dt) + ½ dt F(t + dt)/m
349    if (n_end_of_step) modify->end_of_step();` },
  { name: '6 · Output', where: ['src/verlet.cpp', 354], does: 'Only on steps where thermo, a dump or a restart is due.', code: `354    if (ntimestep == output->next) {
356      output->write(ntimestep);` },
];

const SCRIPT = `${IN_MELT.replace('run             250', 'dump            1 all custom 1 s.dump id x y z vx vy vz fx fy fz\ndump_modify     1 format float %.17g sort id\nrun             2')}`;

export default {
  id: 'verlet-loop',
  title: 'The loop, line by line',
  minutes: 20,
  summary: 'Every timestep of every LAMMPS run is the same 115 lines of Verlet::run(). You walk through them in six stages, then take the positions, velocities and forces LAMMPS dumps at steps 0, 1 and 2, apply the two half-kicks and the drift of fix nve yourself, and get every one of the 24 000 new numbers exactly, bit for bit. On the way you see why atoms can sit outside the box for up to 20 steps.',
  learn: [
    'The six stages of one step: count, first half, neighbours, forces, second half, output',
    'Velocity Verlet as LAMMPS splits it: half kick and drift in initial_integrate, half kick in final_integrate',
    'Recompute x and v for every atom from the dumped numbers, exactly',
    'Why coordinates leave the box between neighbour-list rebuilds',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: one step in six stages</h2>
    <div class="prose"><p>${gh('src/verlet.cpp', 229)}. Tap a stage.</p></div>
    <div class="chips vl-chips">${STAGES.map((s, i) => `<button class="chip" data-s="${i}">${esc(s.name)}</button>`).join('')}</div>
    <p class="vl-does"></p>
    <pre class="lmp-mini vl-code"></pre>
    <div class="prose"><p>This is the velocity Verlet of lesson 1.5, split in two so that fixes can act in between: stage 2 does the first half-kick and the drift, stage 4 computes the new force, stage 5 does the second half-kick. With <code>dt = 0.005</code> and m = 1 in LJ units, each half-kick multiplies the force by <code>dtfm = ½ · 0.005 · 1 / 1 = 0.0025</code>.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: two steps for every atom, by hand</h2>
    <div class="prose"><p>Run melt for 2 steps with a dump of x, v and f at steps 0, 1 and 2, sorted by ID. The new force F(t + dt) comes from LAMMPS (page C.4 computes one by hand); everything else is the three lines of fix nve:</p>
<pre class="lmp-mini">v½   = v(t) + 0.0025 · F(t)
x(t+dt) = x(t) + 0.005 · v½
v(t+dt) = v½ + 0.0025 · F(t+dt)</pre></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run 2 steps and read the dump</button><span class="lmp-status muted vl-st"></span></div>
    <div class="row vl-pick" hidden><label class="lab" style="margin:0">atom ID</label><input type="number" class="lmp-in sm vl-id" value="1" min="1" max="4000"><label class="lab" style="margin:0">step</label><select class="lc-u vl-step"><option value="0">0 → 1</option><option value="1">1 → 2</option></select></div>
    <div class="calc vl-hand" hidden></div>
    <div class="tbl vl-cmp"></div>
    <div class="prose vl-all"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe: atoms outside the box</h2>
    <p>Atom 1 starts at (0, 0, 0). After one step its x is negative: it has left the box, and LAMMPS has not put it back. Wrapping (<code>domain->pbc()</code>, line 277) happens only on a rebuild step, every 20 steps in melt. In between, an atom may drift a little outside the box. That is allowed: ghosts and the neighbour list were built with the skin as a margin for exactly this movement (page A.7). Dump <code>x</code> on a non-rebuild step and you may see coordinates slightly below 0 or above L = ${L.toFixed(4)}; dump with <code>xu</code> or wait for a rebuild step if you need them inside.</p>
    <div class="vl-out"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> In real units, m = 39.948 g/mol, dt = 1 fs. What is dtfm?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>dtf = ½ · 1 · ftm2v, with ftm2v = 1/48.88821291² = 4.184 × 10⁻⁴ (page A.4), so dtf = 2.092 × 10⁻⁴ and dtfm = dtf / 39.948 = <b>5.237 × 10⁻⁶</b> (Å/fs per kcal/mol/Å).</div></div>
    <div class="ex"><p><b>2.</b> A thermostat fix wants to change velocities after the new forces are known. Which stage, and which function of the fix?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Stage 5: <code>end_of_step()</code> (fix temp/rescale) or <code>post_force()</code> (fix langevin adds random and friction forces before the second kick). Fix nvt instead replaces nve's own kicks with scaled ones (page C.3).</div></div>
    <div class="ex"><p><b>3.</b> If the dump wrote x before stage 2 instead of after stage 5, which numbers in section 2 would change?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>None of the arithmetic, but the dump at "step 1" would show x(0). LAMMPS writes output in stage 6, after the whole step, so x, v and f in a dump all belong to the same time t.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const show = i => { root.querySelectorAll('[data-s]').forEach(b => b.classList.toggle('on', +b.dataset.s === i)); const S = STAGES[i]; $('.vl-does').innerHTML = `${esc(S.does)} · ${gh(...S.where)}`; $('.vl-code').textContent = S.code; };
    root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => show(+b.dataset.s)); show(1);
    let D = null;
    const dtfm = (0.5 * 0.005 * 1.0) / 1.0, dtv = 0.005;
    const s17 = x => x.toPrecision(17).replace(/\.?0+(e|$)/, '$1');
    const pick = () => {
      if (!D) return;
      let id = parseInt($('.vl-id').value); if (!(id >= 1 && id <= 4000)) id = 1;
      const st = +$('.vl-step').value, a = D[st][id], b = D[st + 1][id];
      const rows = [], lines = [];
      for (let d = 0; d < 3; d++) {
        const ax = 'xyz'[d]; let v = a[3 + d]; v += dtfm * a[6 + d]; const x = a[d] + dtv * v; const vn = v + dtfm * b[6 + d];
        lines.push(`${ax}:  v½ = ${s17(a[3 + d])} + 0.0025 × ${s17(a[6 + d])} = ${s17(v)}`);
        lines.push(`    x  = ${s17(a[d])} + 0.005 × v½ = ${s17(x)}`);
        lines.push(`    v  = v½ + 0.0025 × ${s17(b[6 + d])} = ${s17(vn)}`);
        rows.push([ax, x, b[d]]); rows.push([`v${ax}`, vn, b[3 + d]]);
      }
      $('.vl-hand').hidden = false; $('.vl-hand').textContent = `atom ${id}, step ${st} → ${st + 1}\n` + lines.join('\n');
      $('.vl-cmp').innerHTML = `<table>${head(['', 'by hand', `LAMMPS step ${st + 1}`, ''])}${rows.map(([k, h, l]) => row([k, `<span class="mono">${s17(h)}</span>`, `<span class="mono">${s17(l)}</span>`, h === l ? '<span class="m-ok">✓ identical</span>' : '<span class="m-bad">✗</span>'])).join('')}</table>`;
    };
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.vl-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(SCRIPT, { read: ['s.dump'] }); b.disabled = false;
      const txt = res.files && res.files['s.dump'];
      if (res.error || !txt) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump')}</span>`; return; }
      st.textContent = '';
      D = {};
      for (const f of txt.split('ITEM: TIMESTEP').slice(1)) { const ls = f.trim().split('\n'), stp = +ls[0], m = {}; for (const l of ls.slice(8)) { const w = l.trim().split(/\s+/).map(Number); m[w[0]] = w.slice(1); } D[stp] = m; }
      let checked = 0, same = 0; const out = [0, 0, 0];
      for (const stp of [0, 1]) for (let id = 1; id <= 4000; id++) {
        const a = D[stp][id], c = D[stp + 1][id];
        for (let d = 0; d < 3; d++) { let v = a[3 + d]; v += dtfm * a[6 + d]; const x = a[d] + dtv * v, vn = v + dtfm * c[6 + d]; checked += 2; if (x === c[d]) same++; if (vn === c[3 + d]) same++; }
      }
      for (const stp of [0, 1, 2]) for (let id = 1; id <= 4000; id++) for (let d = 0; d < 3; d++) { const x = D[stp][id][d]; if (x < 0 || x >= L) out[stp]++; }
      $('.vl-pick').hidden = false; pick();
      $('.vl-all').innerHTML = `<p>All atoms, both steps: <b>${same} of ${checked}</b> new positions and velocities are identical to LAMMPS's, to the last bit. ${same === checked ? '<span class="m-ok">✓</span>' : ''}</p>`;
      $('.vl-out').innerHTML = `<div class="tbl"><table>${head(['step', 'coordinates outside 0 … L'])}${out.map((n, i) => row([i, n])).join('')}</table></div>`;
    };
    $('.vl-id').oninput = pick; $('.vl-step').onchange = pick;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
