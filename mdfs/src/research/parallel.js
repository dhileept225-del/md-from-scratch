// Research · B.7: running in parallel. Domain decomposition: processor grid, atoms and ghosts per processor predicted by hand and matched with real mpirun runs; speed-up from the timing table.
import { g8 } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;

// melt: 10 fcc cells per side; the lattice constant exactly as LAMMPS computes it (pow(4/0.8442, 1/3))
const A = 1.6795961913825073, L = 10 * A - 0 * A, H = A / 2, CUT = 2.5 + 0.3, EPS = 1.0e-6;
const X = []; for (let i = 0; i < 10; i++) for (const b of [0, 0.5]) X.push((i + b) * A);   // the 20 planes along one axis

// measured: in.melt with run 0, LAMMPS 22 Jul 2025, MPICH 5.0.2, mpirun -np P on a 2-core Linux computer
const RUNS = {
  1: ['1 by 1 by 1', [4000, 4000, 4000], [5841, 5841, 5841]],
  2: ['1 by 1 by 2', [2000, 2000, 2000], [4196, 4196, 4196]],
  3: ['1 by 1 by 3', [1333.33, 1400, 1200], [3648.33, 3903, 3339]],
  4: ['1 by 2 by 2', [1000, 1000, 1000], [2901, 2901, 2901]],
  5: ['1 by 1 by 5', [800, 1000, 600], [3209, 3409, 3009]],
  6: ['1 by 2 by 3', [666.667, 700, 600], [2470, 2613, 2284]],
  7: ['1 by 1 by 7', [571.429, 600, 400], [2917.29, 3245, 2680]],
  8: ['2 by 2 by 2', [500, 500, 500], [1956, 1956, 1956]],
  9: ['1 by 3 by 3', [444.444, 490, 360], [2077, 2286, 1791]],
  12: ['2 by 2 by 3', [333.333, 350, 300], [1641.67, 1723, 1529]],
  16: ['2 by 2 by 4', [250, 300, 200], [1484, 1534, 1434]],
};
const TIMING = [  // 32 000 atoms (20 cells per side), 500 steps, same computer
  { P: 1, loop: 6.10393, pair: 80.58, comm: 1.26 },
  { P: 2, loop: 3.16732, pair: 76.05, comm: 7.31 },
  { P: 4, loop: 4.01285, pair: 41.52, comm: 37.66 },
];

// all factorizations px·py·pz = P in the order procmap.cpp makes them; the first with the smallest surface wins
function grids(P) {
  const f = [];
  for (let i = 1; i <= P; i++) { if (P % i) continue; const n = P / i; for (let j = 1; j <= n; j++) if (n % j === 0) f.push([i, j, n / j]); }
  return f.map(g => ({ g, surf: 1 / (g[0] * g[1]) + 1 / (g[0] * g[2]) + 1 / (g[1] * g[2]) }));
}
function bestGrid(P) { let best = null; for (const x of grids(P)) if (!best || x.surf < best.surf - 1e-12) best = x; return best.g; }

// along one axis split n ways: owned planes and ghost planes (as images) of brick k
function axis(n) {
  return Array.from({ length: n }, (_, k) => {
    const lo = 0 + L * (k * 1.0 / n), hi = k < n - 1 ? 0 + L * ((k + 1) * 1.0 / n) : L;
    const clo = lo - (k === 0 ? EPS : 0), chi = hi - (k === n - 1 ? 2 * EPS : 0);   // create_atoms trims the box edges
    const own = []; X.forEach((x, p) => { if (x >= clo && x < chi) own.push(p); });
    const all = own.map(p => p); let below = 0, above = 0;
    X.forEach((x, p) => { for (const m of [-2, -1, 0, 1, 2]) { const y = x + m * L; if (m === 0 && own.includes(p)) continue; if (y >= lo - CUT && y < lo) { all.push(p); below++; } else if (y >= hi && y <= hi + CUT) { all.push(p); above++; } } });
    return { lo, hi, own, all, below, above };
  });
}
const sites = (a, b, c) => { let n = 0; for (const i of a) for (const j of b) for (const k of c) if ((i + j + k) % 2 === 0) n++; return n; };
function predict(P) {
  const g = bestGrid(P), ax = g.map(axis), nl = [], ng = [];
  for (const a of ax[0]) for (const b of ax[1]) for (const c of ax[2]) { const o = sites(a.own, b.own, c.own); nl.push(o); ng.push(sites(a.all, b.all, c.all) - o); }
  const st = v => [+(v.reduce((s, x) => s + x, 0) / v.length).toFixed(3), Math.max(...v), Math.min(...v)];
  return { g, ax, nl: st(nl), ng: st(ng) };
}
const fmtAvg = x => String(+x.toFixed(2));

export default {
  id: 'parallel',
  title: 'Running in parallel',
  minutes: 20,
  summary: 'LAMMPS runs on many processors by cutting the box into bricks, one per processor. You predict by hand which grid of bricks LAMMPS chooses, how many atoms and ghost atoms each processor holds, and why five processors get uneven shares, and check every number against real mpirun runs. Then you read the speed-up and the cost of communication from the timing table.',
  learn: [
    'Domain decomposition: each processor owns a brick and its atoms',
    'How LAMMPS chooses the processor grid (smallest surface)',
    'Atoms and ghosts per processor, predicted plane by plane',
    'Speed-up, efficiency, and why more processors can be slower',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: one box, many bricks</h2>
    <p>Run <code>mpirun -np 4 lmp -in in.melt</code> and four copies of LAMMPS start. Each one gets a <b>brick</b> of the box (${gh('src/domain.cpp', 371)}: <code>sublo = boxlo + xprd·xsplit</code>), stores only the atoms inside it plus ghosts around it (page B.6), and computes their forces. Every step they talk through MPI (${gh('src/verlet.cpp', 264)}):</p>
    <div class="tbl wide"><table>${head(['When', 'Call', 'What moves between processors'])}
${row(['every step', `<code>comm->forward_comm()</code> (${gh('src/verlet.cpp', 268)})`, 'new positions of the atoms other processors hold as ghosts'])}
${row(['every step', `<code>comm->reverse_comm()</code> (${gh('src/verlet.cpp', 341)})`, 'forces computed on ghosts, sent back to their owners (newton on)'])}
${row(['when the neighbour list is rebuilt', `<code>comm->exchange()</code> (${gh('src/verlet.cpp', 284)})`, 'atoms that left a brick, handed to the new owner'])}
${row(['same', `<code>comm->borders()</code> (${gh('src/verlet.cpp', 286)})`, 'the new list of ghosts, as on page B.6'])}</table></div>
    <p>The force loop itself is the same code as in serial: each processor just loops over fewer atoms. The browser runs one processor only, so the parallel numbers on this page were measured with real <code>mpirun</code> on a 2-core Linux computer (LAMMPS 22 Jul 2025, MPICH 5.0.2) and are matched by hand below.</p>
  </section>

  <section class="block">
    <h2>2 · Predict the processor grid</h2>
    <div class="prose"><p>For P processors LAMMPS tries every px × py × pz = P and keeps the one where the bricks have the <b>smallest surface</b>, because the surface is where ghosts live and communication happens (${gh('src/procmap.cpp', 875)}). For a cube of side L the surface per brick, divided by L², is 1/(px·py) + 1/(px·pz) + 1/(py·pz). A tie keeps the first one found (<code>surf &lt; bestsurf</code>, ${gh('src/procmap.cpp', 878)}).</p></div>
    <div class="row"><label class="lab" style="margin:0">processors P</label><select class="lc-u pp">${Object.keys(RUNS).map(p => `<option${p === '4' ? ' selected' : ''}>${p}</option>`).join('')}</select></div>
    <div class="tbl pg-grid"></div>
  </section>

  <section class="block">
    <h2>3 · Match: atoms and ghosts on each processor</h2>
    <div class="prose"><p>At step 0 the atoms sit on 20 planes per direction, ${g8(H)} apart, 200 atoms per plane. Along each split direction, a brick owns the planes between its two walls; its ghosts are the planes within the communication cutoff (2.5 + 0.3 = 2.8) beyond each wall, taken from the neighbouring bricks (or periodic images). Count the FCC sites in the grown brick, subtract the owned atoms: that is Nghost (page B.6).</p></div>
    <div class="tbl pg-axes"></div>
    <div class="tbl pg-cmp"></div>
    <div class="prose pg-note"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe: is it faster?</h2>
    <p>32 000 atoms (20 cells per side), 500 steps, on the same 2-core computer:</p>
    <div class="tbl"><table>${head(['P', 'loop time (s)', 'speed-up = t₁/t_P', 'efficiency = speed-up/P', 'Pair %', 'Comm %'])}${TIMING.map(t => row([t.P, t.loop, (TIMING[0].loop / t.loop).toFixed(2), `${(100 * TIMING[0].loop / t.loop / t.P).toFixed(0)} %`, t.pair, t.comm])).join('')}</table></div>
    <ul>
      <li><b>2 processors</b>: 6.104 / 3.167 = 1.93 times faster, 96 % efficient. Comm grows from 1.3 % to 7.3 % of the time: the bricks now have faces that need messages.</li>
      <li><b>4 processors on 2 cores</b>: slower than 2. Two processes share each core, and every exchange of positions waits for the process that is not running. Never ask for more MPI processes than physical cores.</li>
      <li>Ghosts grow relative to owned atoms as bricks shrink: 5841 / 4000 = 1.5 ghosts per atom with 1 processor, 1956 / 500 = 3.9 with 8. Each brick does less work but carries more surface. Below a few hundred atoms per processor, communication starts to dominate.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>On your own computer</h2>
<pre class="lmp-mini">mpirun -np 4 lmp -in in.melt                 # LAMMPS picks the grid
# force a grid in the script, before create_box:
processors      4 1 1</pre>
    <p>The log prints the grid (<code>1 by 2 by 2 MPI processor grid</code>), and Nlocal and Nghost with ave, max and min over the processors. A large gap between max and min means uneven work: the slowest processor sets the speed (the <code>balance</code> command moves the walls to even it out; page F.5 covers performance).</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Which grid does LAMMPS choose for 6 processors, and what is its surface?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>1 × 2 × 3: 1/2 + 1/3 + 1/6 = <b>1</b>. 1 × 1 × 6 gives 1 + 1/6 + 1/6 = 1.33. The orders 1 × 3 × 2, 2 × 1 × 3 … give the same 1, but 1 × 2 × 3 is found first.</div></div>
    <div class="ex"><p><b>2.</b> With 3 processors, why do two get 1400 atoms and one 1200?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The 20 planes along z are cut at 20/3 = 6.67 and 13.33 plane spacings: planes 0–6 (7 planes), 7–13 (7), 14–19 (6). 7 × 200 = 1400, 6 × 200 = 1200.</div></div>
    <div class="ex"><p><b>3.</b> You have a 32-core node and 4000 atoms. Is <code>mpirun -np 32</code> a good idea?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Probably not: 125 atoms per processor, each with over 1000 ghosts (try P = 16 above: 250 atoms, about 1500 ghosts). Most time would go into communication. Test 4, 8, 16 and keep the fastest; for small systems, fewer processors (or OpenMP threads) win.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const draw = () => {
      const P = +$('.pp').value, R = RUNS[P], pr = predict(P), best = pr.g;
      const all = grids(P), minS = Math.min(...all.map(x => x.surf));
      $('.pg-grid').innerHTML = `<table>${head(['px × py × pz', 'surface / L²', ''])}${all.map(x => row([x.g.join(' × '), x.surf.toFixed(4), x.g === best ? '<b>chosen</b> (first with the smallest)' : Math.abs(x.surf - minS) < 1e-12 ? 'same surface, found later' : ''])).join('')}</table>
        <p>Hand: <b>${best.join(' by ')}</b> · LAMMPS: <b>${R[0]}</b> ${best.join(' by ') === R[0] ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</p>`;
      const names = ['x', 'y', 'z'];
      $('.pg-axes').innerHTML = `<table>${head(['axis', 'brick', 'walls (σ)', 'owned planes', 'ghost planes below + above'])}${pr.ax.map((bricks, d) => bricks.map((b, k) => row([names[d], `${k + 1} of ${bricks.length}`, `${g8(b.lo)} – ${g8(b.hi)}`, b.own.length ? `${b.own[0]}–${b.own.at(-1)} (${b.own.length})` : '0', `${b.below} + ${b.above}`])).join('')).join('')}</table>`;
      const c = (h, l) => `<td class="mono">${h}</td><td class="mono">${l}</td><td>${String(h) === String(l) ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td>`;
      $('.pg-cmp').innerHTML = `<table><tr><th>P = ${P}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>
        <tr><td>Nlocal ave</td>${c(fmtAvg(pr.nl[0]), fmtAvg(R[1][0]))}</tr><tr><td>Nlocal max</td>${c(pr.nl[1], R[1][1])}</tr><tr><td>Nlocal min</td>${c(pr.nl[2], R[1][2])}</tr>
        <tr><td>Nghost ave</td>${c(fmtAvg(pr.ng[0]), fmtAvg(R[2][0]))}</tr><tr><td>Nghost max</td>${c(pr.ng[1], R[2][1])}</tr><tr><td>Nghost min</td>${c(pr.ng[2], R[2][2])}</tr></table>`;
      const onWall = pr.ax.some(bricks => bricks.length > 1 && bricks.some((b, k) => k > 0 && Math.abs(b.lo / H - Math.round(b.lo / H)) < 1e-9));
      $('.pg-note').innerHTML = P === 5 || P === 16 ? `<p><b>Why uneven?</b> With ${P} processors some walls fall exactly on a plane of atoms (20/${P === 5 ? 5 : 4} = ${P === 5 ? 4 : 5} spacings per brick). Whether that plane belongs to the brick below or above is decided by the last bit of two floating-point numbers: the wall L·(k/n) and the atom position (i + ½)·a. One brick gets the plane, its neighbour loses it: 1000 and 600 atoms instead of 800 each for P = 5. The hand calculation reproduces it only because it does the same floating-point arithmetic as LAMMPS.</p>`
        : onWall ? '<p>Here the walls fall exactly on planes too, but the rounding happens to split them evenly.</p>' : '';
    };
    $('.pp').onchange = draw; draw();
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
