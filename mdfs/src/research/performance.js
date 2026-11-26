// Research · F.5: performance. Reading the timing table, the neighbor skin trade-off measured and modelled,
// O(N) scaling, and what two MPI ranks buy; every printed performance number and the ghost and rebuild counts
// recomputed by hand, plus an exact neighbor count and a skin scan in the browser.
import { NATIVE } from './performance-native.js';
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const RHO = 0.8442, A = Math.cbrt(4 / RHO), RC = 2.5, DT = 0.005;
const min = (a, k) => Math.min(...a.map(r => r[k]));

const INPUT = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 \${L} 0 \${L} 0 \${L}     # L = 10: 4000 atoms
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.44 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
neighbor \${s} bin                          # the skin
neigh_modify every 1 delay 0 check yes
fix 1 all nve
thermo 500
run 500                                    # melt; T settles near 0.70
run \${n}                                   # timed

lmp -in in.sk -var L 10 -var s 0.3 -var n 2000`;
const SCAN = (s, L) => `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 ${L} 0 ${L} 0 ${L}
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.44 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
neighbor ${s} bin
neigh_modify every 1 delay 0 check yes
fix 1 all nve
thermo 500
run 300
run 600`;
const EXACT = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 6 0 6 0 6
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.44 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
neighbor 0.3 bin
fix 1 all nve
run 300
dump d all custom 1 e.dump id x y z
dump_modify d sort id format float %.17g
run 0                 # setup rebuilds the list from these positions`;
const TIPS = [
  ['neighbor skin', 'Default 0.3 σ (lj) or 2 Å (real, metal) is usually near optimal; the minimum is flat. Check "Dangerous builds = 0" whenever you change it.'],
  ['neigh_modify every N', 'every 1 checks displacements each step (cheap). every 10 check no would skip the check but risks missed pairs.'],
  ['newton on', 'Default. Each pair is computed once and the ghost forces are sent back; off doubles the pair work to save one communication.'],
  ['output', 'thermo and dumps every step cost real time (page D.18: %g formatting). Here Output is 0.00 % at thermo 500.'],
  ['processors / balance', 'Choose a processor grid that minimizes surface; balance or fix balance when the density is uneven (interfaces, vacuum).'],
  ['accelerators', '-sf omp (OPENMP), -sf opt (OPT), -sf intel, -sf gpu, -k on -sf kk (KOKKOS): same input, different implementations of the same styles. Check one against another with a short run, as on page F.3.'],
];

export default {
  id: 'performance',
  title: 'Performance',
  minutes: 15,
  summary: 'Where the time goes in a LAMMPS run, and how to make it go less. The timing table at the end of every run, the neighbor skin trade-off measured over ten skins and explained with a two-parameter cost model, linear scaling with system size, and what a second MPI rank buys and costs. Every number in the "Performance:" line, the ghost-atom counts, and the number of neighbor rebuilds are recomputed by hand, and the neighbor count in your browser is checked exactly.',
  learn: [
    'Reading Performance:, the section table and Nghost/Neighs',
    'The skin trade-off: Pair grows with the skin, Neigh shrinks',
    'Predicting rebuilds from the fastest atom of a Maxwell distribution',
    'O(N) scaling, and the price of a second MPI rank',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the timing table</h2>
    <div class="prose"><p>Every run ends with this block (${gh('src/finish.cpp')}), here for 4000 LJ atoms, 2000 steps, skin 0.3 σ, Release build on a 2-core cloud machine:</p></div>
    <pre class="lmp-mini">${esc(INPUT)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.blk)}</pre>
    <div class="prose"><ul>
    <li><b>Pair</b>: the force loop of page C.5, two-thirds of the time. <b>Neigh</b>: building the lists of page C.4 and checking whether a rebuild is due. <b>Comm</b>: ghost atoms (in serial, periodic copies of the atoms themselves). <b>Modify</b>: fixes, here only nve. <b>Other</b>: what is left of the loop time.</li>
    <li>min/avg/max are over MPI ranks; %varavg measures how unequal they are (load imbalance).</li>
    <li><b>Nghost</b>: ghost atoms each rank keeps; <b>Neighs</b>: pairs in the half list, each counted once.</li></ul></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Check every Performance line</button></div>
    <div class="prose pf-asum"></div>
  </section>

  <section class="block">
    <h2>2 · Modify the skin: ten runs, seven times each</h2>
    <div class="prose"><p>The list holds all pairs closer than r<sub>c</sub> + skin, and it stays valid until some atom has moved half a skin. A larger skin means more pairs to test on every step (Pair grows) but fewer rebuilds (Neigh shrinks). Ten skins, each run seven times, interleaved; the table shows the fastest of the seven:</p></div>
    <div class="tbl pf-scan"></div>
    <div class="prose"><p>On a shared machine the same run varies by up to 20 % (the spread column). Compare minima, repeat runs, and do not trust differences smaller than that. Even so the picture is clear. Below 0.15 σ the run spends more time rebuilding than computing forces. Above 0.4 σ it tests ever more pairs that are too far apart to interact: at skin 1.0, nearly two-thirds of the listed pairs are beyond r<sub>c</sub>. The optimum is broad and sits at the default, 0.3 σ.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Fit the cost model</button></div>
    <div class="prose pf-bsum"></div>
    <div class="tbl pf-btab"></div>
  </section>

  <section class="block">
    <h2>3 · Match: rebuilds, ghosts and neighbors by hand</h2>
    <div class="prose"><p><b>Rebuilds.</b> With <code>check yes</code> the list is rebuilt when some atom has moved more than skin/2 since the last build (page E.8). The fastest atom decides, so the number of builds should be about steps × 2·v<sub>max</sub>·dt / skin, with v<sub>max</sub> the largest speed among 4000 atoms at T ≈ ${NATIVE.T.toFixed(2)}.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Predict the builds from Maxwell</button></div>
    <div class="prose pf-csum"></div>
    <div class="prose"><p><b>Ghosts.</b> Each rank copies every atom within the communication cutoff r<sub>c</sub> + skin of its subdomain. For a box of sides L<sub>x</sub>, L<sub>y</sub>, L<sub>z</sub> at density ρ: N<sub>ghost</sub> ≈ ρ[(L<sub>x</sub>+2c)(L<sub>y</sub>+2c)(L<sub>z</sub>+2c) − L<sub>x</sub>L<sub>y</sub>L<sub>z</sub>].</p></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Estimate Nghost for all runs</button></div>
    <div class="prose pf-dsum"></div>
    <div class="prose"><p><b>Neighbors, exactly.</b> In a liquid no formula gives the pair count to the last pair, but the positions do. Your browser melts 864 atoms, dumps the positions at the setup of a <code>run 0</code> (where the list is rebuilt), and counts every pair closer than 2.8 σ:</p></div>
    <pre class="lmp-mini">${esc(EXACT)}</pre>
    <div class="row"><button class="btn sm primary" data-k="e">▶ Count the pairs</button><span class="pf-ste"></span></div>
    <div class="prose pf-esum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: system size and a second rank</h2>
    <div class="tbl pf-n"></div>
    <div class="prose"><p>Cost per atom per step is constant to within the noise from 864 to 23 328 atoms: with neighbor lists and binning, LAMMPS is O(N). Double the atoms, double the time.</p></div>
    <div class="tbl pf-p"></div>
    <div class="prose pf-psum"></div>
  </section>

  <section class="block">
    <h2>5 · Your turn in the browser</h2>
    <div class="prose"><p>The same skin scan, smaller (864 atoms, 600 timed steps), in WebAssembly on your device. The absolute times differ from the native ones; the shape should not.</p></div>
    <div class="row"><button class="btn sm primary" data-k="f">▶ Scan five skins</button><span class="pf-stf"></span></div>
    <div class="tbl pf-ftab"></div>
    <div class="prose pf-fsum"></div>
  </section>

  <section class="block">
    <h2>6 · Knobs that matter</h2>
    <div class="tbl pf-tips"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> At skin 0.05 there were exactly 1000 builds in 2000 steps. Why exactly every second step?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The fastest atom moves about v<sub>max</sub>·dt ≈ 3.4 × 0.005 = 0.017 σ per step. After one step that is below skin/2 = 0.025, after two steps (0.034) it is above. So the check triggers on every second step, and 2000/2 = 1000. Steps between builds come in whole numbers, so for small skins the formula of section 3 rounds up.</div></div>
    <div class="ex"><p><b>2.</b> The two MPI runs ended with 411 691 and 411 698 neighbors for the same input. Which one is right?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Both. With two ranks the forces are summed in a different order, the trajectories differ in the last bit and then diverge (page F.3). After 1500 steps they are two different liquids with slightly different pair counts. A parallel run is not bitwise the same as a serial one, and need not be.</div></div>
    <div class="ex"><p><b>3.</b> Your production system has 50 000 atoms. Should you change the skin from the default?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Probably not by much. The optimum depends on the speeds of the fastest atoms (temperature, mass, time step) and on the cost of a pair (the potential), not on N; the scan here varies by less than 10 % between 0.2 and 0.4. If you want to know, run 2–3 skins for a few thousand steps, repeat each, and compare the minima, as in section 2.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    const S = NATIVE.S.map(x => ({ s: x.s, T: x.r[0].tot, B: x.r[0].builds, G: x.r[0].ghost, P: min(x.r, 'pair'), Ne: min(x.r, 'neigh'), L: min(x.r, 'loop'), Lmax: Math.max(...x.r.map(r => r.loop)), O: Math.min(...x.r.map(r => r.loop - r.pair - r.neigh)) }));
    const ncut = s => RHO * 4 / 3 * Math.PI * (RC + s) ** 3 / 2;
    $('.pf-scan').innerHTML = `<table>${head(['skin', 'loop (s)', 'spread', 'Pair', 'Neigh', 'builds', 'neighbors/atom', 'uniform estimate'])}${S.map(r => row([r.s, r.L.toFixed(2), `+${(100 * (r.Lmax / r.L - 1)).toFixed(0)} %`, r.P.toFixed(2), r.Ne.toFixed(2), r.B, (r.T / 4000).toFixed(1), ncut(r.s).toFixed(1)])).join('')}</table>`;
    $('.pf-n').innerHTML = `<table>${head(['atoms', 'loop, 1000 steps (s)', 'µs per atom-step'])}${NATIVE.N.map(x => row([x.N, min(x.r, 'loop').toFixed(3), (1e6 * min(x.r, 'loop') / x.N / 1000).toFixed(3)])).join('')}</table>`;
    const p1 = NATIVE.P.np1[0], p2 = NATIVE.P.np2[0], t1 = min(NATIVE.P.np1, 'loop'), t2 = min(NATIVE.P.np2, 'loop');
    $('.pf-p').innerHTML = `<table>${head(['10 976 atoms, 1000 steps', '1 rank', '2 ranks'])}${[['loop (s)', t1.toFixed(2), t2.toFixed(2)], ['Pair %', (100 * p1.pair / p1.loop).toFixed(1), (100 * p2.pair / p2.loop).toFixed(1)], ['Neigh %', (100 * p1.neigh / p1.loop).toFixed(1), (100 * p2.neigh / p2.loop).toFixed(1)], ['Comm %', (100 * p1.comm / p1.loop).toFixed(1), (100 * p2.comm / p2.loop).toFixed(1)], ['Nlocal per rank', 10976, 5488], ['Nghost per rank', p1.ghost, p2.ghost]].map(row).join('')}</table>`;
    $('.pf-psum').innerHTML = `<p>Speed-up ${(t1 / t2).toFixed(2)} on 2 ranks, not 2. The processor grid is 1 × 1 × 2: each rank owns half the box but still needs a ghost layer on both cut faces, now real messages instead of local copies, so Comm rises from ${(100 * p1.comm / p1.loop).toFixed(1)} % to ${(100 * p2.comm / p2.loop).toFixed(1)} %, and the ghosts per rank (${p2.ghost}) outnumber the atoms it owns. The surface-to-volume ratio sets the limit; with millions of atoms per rank it is small, with a few thousand it dominates.</p>`;
    $('.pf-tips').innerHTML = `<table>${head(['knob', 'what to know'])}${TIPS.map(row).join('')}</table>`;

    $('[data-k=a]').onclick = () => {
      // finish.cpp:128-160: loop printed with {:.6g}; t_step = loop/steps; tau/day = 86400/t_step*dt; all with {:.3f}
      const all = [...NATIVE.S.flatMap(x => x.r.map(r => [r, 2000, 4000])), ...NATIVE.N.flatMap(x => x.r.map(r => [r, 1000, x.N])), ...NATIVE.P.np1.map(r => [r, 1000, 10976]), ...NATIVE.P.np2.map(r => [r, 1000, 10976])];
      let ok = 0, ex = null;
      for (const [r, st, N] of all) {
        const m = r.perf.match(/([\d.]+) tau\/day, ([\d.]+) timesteps\/s, ([\d.]+) (k|M)atom-step\/s/), u = m[4] === 'M' ? 1e6 : 1e3;
        const e = Math.floor(Math.log10(r.loop)) - 5, d = 0.5 * 10 ** e;          // half a unit in the 6th digit
        let lo = r.loop - d, hi = r.loop + d;
        const sp = Number(m[2]); lo = Math.max(lo, st / (sp + 5e-4)); hi = Math.min(hi, st / (sp - 5e-4));   // and from timesteps/s
        const inn = (v, f) => v >= f(hi) - 5e-4 && v <= f(lo) + 5e-4;
        const good = lo <= hi && inn(Number(m[1]), T => 86400 * st * DT / T) && inn(Number(m[3]), T => N * st / T / u);
        ok += good; if (!ex && r === NATIVE.S[4].r[0]) ex = [lo, hi];
      }
      const b = NATIVE.S[4].r[0];
      $('.pf-asum').innerHTML = `<p>All three numbers of a Performance line come from one wall time (${gh('src/finish.cpp', 142)}): t<sub>step</sub> = loop/steps, τ/day = 86 400 · dt / t<sub>step</sub>, atom-steps/s = N / t<sub>step</sub>. The loop time is printed with only 6 digits, so it cannot reproduce the 9-digit τ/day by itself (86 400 × 2000 × 0.005 / ${b.loop} = ${(86400 * 2000 * DT / b.loop).toFixed(3)}, printed ${b.perf.split(' ')[0]}). Combined with the printed timesteps/s it pins the true loop time to [${ex[0].toFixed(7)}, ${ex[1].toFixed(7)}] s, and every printed number must then fall inside the range that interval implies: ${ok} of ${all.length} Performance lines consistent ${ok === all.length ? OK : BAD}. The %total column is time / loop: Pair ${b.pair} / ${b.loop} = ${(100 * b.pair / b.loop).toFixed(2)} %.</p>`;
    };

    const fit = (X, Y) => { const n = X.length; let sx = 0, sy = 0, sxx = 0, sxy = 0; for (let i = 0; i < n; i++) { sx += X[i]; sy += Y[i]; sxx += X[i] * X[i]; sxy += X[i] * Y[i]; } const a = (n * sxy - sx * sy) / (n * sxx - sx * sx); return [a, (sy - a * sx) / n]; };
    $('[data-k=b]').onclick = () => {
      const [a] = fit(S.map(r => r.T), S.map(r => r.P / 2000)), aa = S.reduce((s, r) => s + r.P * r.T, 0) / S.reduce((s, r) => s + 2000 * r.T * r.T, 0);
      const [c, d] = fit(S.map(r => r.T), S.map(r => r.Ne / r.B));
      const mod = r => 2000 * aa * r.T + r.B * (c * r.T + d) + r.O;
      const best = S.reduce((p, q) => (q.L < p.L ? q : p)), bm = S.reduce((p, q) => (mod(q) < mod(p) ? q : p));
      $('.pf-bsum').innerHTML = `<p>Two costs, each proportional to the list size N<sub>list</sub>. Pair: every step tests every listed pair, ${(1e9 * aa).toFixed(1)} ns per pair per step (a fit through zero). Neigh: each build costs ${(1e9 * c).toFixed(1)} ns per listed pair plus ${(1e3 * d).toFixed(2)} ms fixed (binning, the ghost exchange that goes with it). So T ≈ steps · ${(1e9 * aa).toFixed(1)} ns · N<sub>list</sub> + builds · (${(1e9 * c).toFixed(1)} ns · N<sub>list</sub> + ${(1e3 * d).toFixed(2)} ms) + other. Fastest measured skin: ${best.s} σ (${best.L.toFixed(2)} s); fastest in the model: ${bm.s} σ ${best.s === bm.s ? OK : '(within the noise)'}.</p>`;
      $('.pf-btab').innerHTML = `<table>${head(['skin', 'measured (s)', 'model (s)', 'difference'])}${S.map(r => row([r.s, r.L.toFixed(2), mod(r).toFixed(2), `${(100 * (mod(r) / r.L - 1)).toFixed(0)} %`])).join('')}</table>`;
    };

    $('[data-k=c]').onclick = () => {
      // Maxwell speeds, 3 components of variance T: P(|v| > u) = erfc(x/√2) + √(2/π) x e^{-x²/2}, x = u/√T
      const erfc = x => { const t = 1 / (1 + 0.5 * Math.abs(x)), y = t * Math.exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))); return x >= 0 ? y : 2 - y; };
      const tail = x => erfc(x / Math.SQRT2) + Math.sqrt(2 / Math.PI) * x * Math.exp(-x * x / 2);
      let lo = 0, hi = 10; for (let k = 0; k < 100; k++) { const m = (lo + hi) / 2; if (tail(m) > 1 / 4000) lo = m; else hi = m; }
      const vmax = lo * Math.sqrt(NATIVE.T);
      const pr = s => { const per = s / 2 / (vmax * DT); return per < 3 ? 2000 / Math.ceil(per) : 2000 * 2 * vmax * DT / s; };
      let good = 0;
      const rows = S.map(r => { const p = pr(r.s), e = Math.abs(p / r.B - 1); good += e < 0.1; return [r.s, r.B, p.toFixed(0), `${(100 * (p / r.B - 1)).toFixed(0)} %`]; });
      $('.pf-csum').innerHTML = `<p>The fastest of 4000 atoms sits where P(|v| &gt; v<sub>max</sub>) = 1/4000: v<sub>max</sub> = ${lo.toFixed(3)}·√T = ${vmax.toFixed(2)} σ/τ. Predicted builds = 2000 × 2 × ${vmax.toFixed(2)} × 0.005 / skin = ${(2000 * 2 * vmax * DT).toFixed(1)} / skin (for the smallest skins, 2000 divided by the whole number of steps the fastest atom needs to cover skin/2). Measured: builds × skin = ${S.slice(1).map(r => (r.B * r.s).toFixed(0)).join(', ')}.</p><div class="tbl"><table>${head(['skin', 'builds', 'predicted', 'difference'])}${rows.map(row).join('')}</table></div><p>${good} of ${S.length} within 10 %. The prediction is an estimate: the fastest atom changes, and it does not move in a straight line.</p>`;
    };

    $('[data-k=d]').onclick = () => {
      const est = (Lx, Ly, Lz, c) => RHO * ((Lx + 2 * c) * (Ly + 2 * c) * (Lz + 2 * c) - Lx * Ly * Lz);
      const rows = S.map(r => { const L = 10 * A, e = est(L, L, L, RC + r.s); return [`skin ${r.s}, 1 rank`, r.G, e.toFixed(0), `${(100 * (e / r.G - 1)).toFixed(1)} %`]; });
      const L14 = 14 * A;
      rows.push(['10 976 atoms, 1 rank', NATIVE.P.np1[0].ghost, est(L14, L14, L14, 2.8).toFixed(0), `${(100 * (est(L14, L14, L14, 2.8) / NATIVE.P.np1[0].ghost - 1)).toFixed(1)} %`]);
      rows.push(['10 976 atoms, each of 2 ranks', NATIVE.P.np2[0].ghost, est(L14, L14, L14 / 2, 2.8).toFixed(0), `${(100 * (est(L14, L14, L14 / 2, 2.8) / NATIVE.P.np2[0].ghost - 1)).toFixed(1)} %`]);
      const w = Math.max(...rows.map(r => Math.abs(parseFloat(r[3]))));
      $('.pf-dsum').innerHTML = `<div class="tbl"><table>${head(['run', 'Nghost', 'ρ × shell volume', 'difference'])}${rows.map(row).join('')}</table></div><p>All within ${w.toFixed(1)} % ${w < 2 ? OK : BAD}; the rest is the liquid's density fluctuation in a shell only 2.5–3.5 σ thick. Note the 1 × 1 × 2 grid: each rank's slab is half as thick in z, yet keeps 70 % of the serial ghost count.</p>`;
    };

    $('[data-k=e]').onclick = async ev => {
      const b = ev.currentTarget, st = $('.pf-ste'); b.disabled = true; st.textContent = ' Running…';
      const R = await runLammps(EXACT, { read: ['e.dump'] }); if (!alive) return; b.disabled = false;
      if (R.error || !R.files['e.dump']) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no dump')}</span>`; return; } st.textContent = '';
      const tl = R.lines.filter(l => /Total # of neighbors/.test(l)), lam = Number(tl[tl.length - 1].split('=')[1]);
      const L = R.files['e.dump'].split('\n'), n = +L[3], lo = L[5].trim().split(/\s+/).map(Number), Lb = lo[1] - lo[0], X = new Float64Array(3 * n);
      for (let k = 0; k < n; k++) { const w = L[9 + k].trim().split(/\s+/); for (let c = 0; c < 3; c++) X[3 * k + c] = Number(w[1 + c]); }
      let np = 0, nin = 0, near = 1;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
        let dx = X[3 * i] - X[3 * j], dy = X[3 * i + 1] - X[3 * j + 1], dz = X[3 * i + 2] - X[3 * j + 2];
        dx -= Lb * Math.round(dx / Lb); dy -= Lb * Math.round(dy / Lb); dz -= Lb * Math.round(dz / Lb);
        const r2 = dx * dx + dy * dy + dz * dz; if (r2 < 2.8 * 2.8) np++; if (r2 < 6.25) nin++; near = Math.min(near, Math.abs(Math.sqrt(r2) - 2.8));
      }
      $('.pf-esum').innerHTML = `<p>"Total # of neighbors" printed by LAMMPS after run 0: ${lam}. Pairs closer than 2.8 σ counted here from the dumped positions: ${np} ${np === lam ? OK : BAD}. The closest pair to the list boundary is ${near.toExponential(1)} σ from it, so no rounding question arises. Of these pairs, ${nin} (${(100 * nin / np).toFixed(0)} %) are inside r<sub>c</sub> = 2.5 and contribute a force; the other ${np - nin} are tested every step and skipped. The uniform estimate ρ·(4/3)π·2.8³/2 per atom gives ${(n * ncut(0.3)).toFixed(0)}: ${(100 * (n * ncut(0.3) / np - 1)).toFixed(1)} % high, because a liquid has fewer neighbours than uniform at these distances (g(r) dips below 1 near 1.5 σ and 2.5 σ).</p>`;
    };

    $('[data-k=f]').onclick = async ev => {
      const b = ev.currentTarget, st = $('.pf-stf'); b.disabled = true; const out = [];
      for (const s of [0.05, 0.15, 0.3, 0.6, 1.0]) {
        st.textContent = ` skin ${s}…`;
        const R = await runLammps(SCAN(s, 6)); if (!alive) return;
        if (R.error) { st.innerHTML = `<span class="lmp-err">${esc(R.error)}</span>`; b.disabled = false; return; }
        const i = R.lines.findIndex(l => /for 600 steps/.test(l)), blk = R.lines.slice(i), g = k => { const l = blk.find(x => x.startsWith(k)); return l ? Number(l.split('|')[2]) : NaN; };
        const bl = blk.find(l => /Neighbor list builds/.test(l)), tn = blk.find(l => /Total # of neighbors/.test(l));
        out.push([s, Number(R.lines[i].match(/Loop time of (\S+)/)[1]), g('Pair'), g('Neigh'), Number(bl.split('=')[1]), Number(tn.split('=')[1])]);
        $('.pf-ftab').innerHTML = `<table>${head(['skin', 'loop (s)', 'Pair', 'Neigh', 'builds', 'neighbors'])}${out.map(r => row([r[0], r[1].toFixed(3), r[2].toFixed(3), r[3].toFixed(3), r[4], r[5]])).join('')}</table>`;
      }
      b.disabled = false; st.textContent = '';
      const best = out.reduce((p, q) => (q[1] < p[1] ? q : p));
      $('.pf-fsum').innerHTML = `<p>Fastest in your browser: skin ${best[0]} σ. Pair rises and Neigh falls with the skin, as natively; where they cross depends on the relative cost of a force and a rebuild on your device, but the default 0.3 should be close to the bottom of a flat curve.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
