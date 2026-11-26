// Research · F.8: multi-replica methods. -partition, world variables, parallel tempering (temper) and the nudged
// elastic band (neb), run natively with the REPLICA package. Every swap decision is reproduced from the replicas'
// energies and LAMMPS's Park–Miller generator; the NEB path's energies and reaction coordinates are recomputed here.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const nat = () => import('./multireplica-native.js').then(m => m.NATIVE);

// random_park.cpp: Park–Miller "minimal standard" with Schrage's trick, exactly as LAMMPS does it
class RanPark {
  constructor(seed) { this.seed = seed; }
  uniform() { const k = Math.trunc(this.seed / 127773); this.seed = 16807 * (this.seed - k * 127773) - 2836 * k; if (this.seed < 0) this.seed += 2147483647; return (1.0 / 2147483647) * this.seed; }
}
const TEMPS = [0.70, 0.76, 0.82, 0.88];
const RUN = `mpirun -np 4 lmp -partition 4x1 -in in.temper
# 4 partitions ("worlds") of 1 MPI rank each; each world runs the whole input
# output: log.lammps (universe: one status line per swap), log.lammps.0..3 and screen.0..3 (one per world)`;
const NRUN = `mpirun -np 4 lmp -partition 4x1 -in in.neb.hop1     # examples/neb, 4 images`;
const SP = k => `dimension 2
boundary p s p
atom_style atomic
atom_modify map array sort 0 0.0
read_data initial.hop1
pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
pair_modify shift yes
read_dump final.${k}.dump 124 x y box yes
thermo_style custom step pe
thermo_modify format float %.17g
run 0`;

export default {
  id: 'multi-replica',
  title: 'Multi-replica methods',
  minutes: 15,
  summary: 'Some methods need several copies of the system running side by side: parallel tempering swaps temperatures between replicas, the nudged elastic band strings images along a reaction path. LAMMPS runs them with -partition, one input for all, with "world" variables to make the copies different. Both run natively here with the REPLICA package. You then reproduce every one of the 40 temperature swaps from the replicas\' energies and LAMMPS\'s own random-number generator, and recompute the band\'s energies, reaction coordinates and barrier, including the factor of N that the printed barrier hides.',
  learn: [
    '-partition, worlds and the universe; world, universe and uloop variables',
    'temper: the swap rule and its random numbers, reproduced swap by swap',
    'neb: what the final line reports, and the per-atom barrier trap',
    'What the browser can and cannot run (one world only)',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: one input, several worlds</h2>
    <div class="prose"><p><code>-partition 4x1</code> splits the MPI ranks into 4 independent LAMMPS instances ("worlds"), each with its own atoms, its own log.lammps.N and screen.N. They all read the same input, so something has to make them differ. A <code>variable t world 0.70 0.76 0.82 0.88</code> takes the N-th value in world N. <code>universe</code> and <code>uloop</code> variables hand out values to whichever world asks next, for running many independent jobs. Commands from the REPLICA package (temper, neb, prd, tad, hyper) then make the worlds talk, through their root ranks. The browser runs one LAMMPS instance in one worker, so it has neither -partition nor REPLICA; the runs below are native, and the browser checks them.</p></div>
    <pre class="lmp-mini">${esc(RUN)}</pre>
  </section>

  <section class="block">
    <h2>2 · Parallel tempering</h2>
    <div class="prose"><p>Four 256-atom LJ liquids at T = 0.70, 0.76, 0.82, 0.88 under fix nvt. Every 100 steps neighbouring temperatures may swap. ${gh('src/REPLICA/temper.cpp', 240)}: with seed_swap = 0, even swaps pair temperatures (0,1),(2,3) and odd swaps pair (1,2). Of each pair, the world with the lower rank decides. It computes Δ = (U<sub>mine</sub> − U<sub>partner</sub>)(1/kT<sub>mine</sub> − 1/kT<sub>partner</sub>) from the two total potential energies, and swaps if Δ ≥ 0, otherwise if a uniform random number is below e<sup>Δ</sup>. The random numbers come from one Park–Miller generator per rank, seeded seed_boltz + rank and discarded 100 times first (line 149). On a swap the velocities are scaled by √(T<sub>new</sub>/T<sub>old</sub>) and the thermostat target is reset.</p></div>
    <div class="row"><button class="btn sm ghost" data-k="in">Show the input</button></div>
    <pre class="lmp-mini mr-in" hidden></pre>
    <div class="prose"><p>The universe log prints, after each swap, the temperature index of each world. The header "T0 T1 T2 T3" means worlds 0–3, not temperatures.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Reproduce all 40 swaps</button><span class="mr-sta"></span></div>
    <div class="prose mr-asum"></div>
    <div class="tbl mr-atab"></div>
  </section>

  <section class="block">
    <h2>3 · Nudged elastic band</h2>
    <div class="prose"><p>examples/neb/in.neb.hop1: a 2-d LJ surface (420 atoms) and an atom hopping out of it to become an adatom. Four images, initial and final fixed, the two in between relaxed by quickmin with springs along the path (fix neb on the 21 atoms near the hop); then a climbing-image stage pushes the highest image to the saddle.</p></div>
    <pre class="lmp-mini">${esc(NRUN)}</pre>
    <div class="row"><button class="btn sm ghost" data-k="nl">Show the final NEB output</button></div>
    <pre class="lmp-mini mr-nl" hidden></pre>
    <div class="prose"><p>${gh('src/REPLICA/neb.cpp', 640)}: PE<sub>i</sub> is each image's potential energy, divided by N when thermo normalizes (LJ units do by default). RD<sub>i</sub> is the cumulative distance along the path, ${gh('src/REPLICA/neb.cpp', 644)}: the sum of |R<sub>i</sub> − R<sub>i−1</sub>| over the fix neb atoms (minimum image), divided by the total length RDT. EBF and EBR are the climbing image's energy minus the first and last image's. Check all of it from the four final configurations:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Recompute the band in the browser</button><span class="mr-stb"></span></div>
    <div class="prose mr-bsum"></div>
    <div class="tbl mr-btab"></div>
  </section>

  <section class="block prose">
    <h2>4 · In practice</h2>
    <ul><li><b>Ranks per world.</b> <code>-partition 8x4</code> gives 8 replicas of 4 ranks each. temper and neb only exchange a few numbers per swap, so the cost is the replicas themselves.</li>
    <li><b>Temperature ladders.</b> Acceptance depends on how much neighbouring energy distributions overlap. Their width grows as √N and their spacing as N·ΔT, so larger systems need more, closer temperatures (roughly √N of them).</li>
    <li><b>Trajectories by temperature.</b> Each world's dump follows one replica wandering in temperature. tools/replica/reorder_remd_traj.py sorts frames back by temperature.</li>
    <li><b>Restarting temper.</b> Write a restart in each world; the optional 7th argument of temper (a world variable) tells each world which temperature index it had.</li>
    <li><b>NEB inputs.</b> The final (and optionally intermediate) images come as coordinate files; check the reaction coordinates are evenly spaced and the climbing image really is the maximum. A barrier is a total energy: multiply the printed per-atom EBF by N.</li></ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why does the world with the lower rank make the decision, and not both worlds?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Both must agree. If each drew its own random number, one could accept and the other reject, and two worlds would end up at the same temperature. So one decides and sends the result to the other (temper.cpp: the higher rank sends its energy, the lower replies with swap = 0 or 1). That is also why only the deciding ranks' generators advance, and why reproducing the run needs to know who drew when.</div></div>
    <div class="ex"><p><b>2.</b> The NEB barrier printed is EBF = 0.0071 for the hop. In ε, how high is it?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>EBF is per atom because thermo normalizes in LJ units: 0.0071061754 × 420 ≈ 2.98 ε. Use <code>thermo_modify norm no</code> before neb to print totals. In real or metal units (no normalization by default) the printed number is already the total.</div></div>
    <div class="ex"><p><b>3.</b> With seed_swap = 0, odd and even swaps alternate. What would a non-zero seed_swap change?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Every world then draws one number from its own copy of a common generator (same seed everywhere) and picks even or odd pairs with probability ½. The sequence is random but identical in every world, so they still agree on who is partnered with whom. Detailed balance is slightly cleaner that way; in practice alternation works fine.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    const tog = (k, sel, fill) => { const btn = $(`[data-k=${k}]`); btn.onclick = async () => { const p = $(sel); if (!p.textContent) p.textContent = await fill(); p.hidden = !p.hidden; btn.textContent = (p.hidden ? 'Show' : 'Hide') + btn.textContent.slice(4); }; };
    tog('in', '.mr-in', async () => (await nat()).tin);
    tog('nl', '.mr-nl', async () => (await nat()).neblog.slice(-6).join('\n'));

    $('[data-k=a]').onclick = async () => {
      const N = await nat(); if (!alive) return;
      const pe = [0, 1, 2, 3].map(w => N.pe[w].map(Number));   // total PE of world w at steps 1000, 1100, ...
      const rb = [0, 1, 2, 3].map(w => { const r = new RanPark(58728 + w); for (let i = 0; i < 100; i++) r.uniform(); return r; });
      let w2t = [0, 1, 2, 3]; const lines = ['1000 0 1 2 3'], att = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], draws = [0, 0, 0, 0]; let expected = 0;
      for (let isw = 0; isw < 40; isw++) {
        const step = 1000 + 100 * (isw + 1), which = isw % 2, t2w = []; w2t.forEach((t, w) => { t2w[t] = w; });
        const nw = w2t.slice(), pairs = which === 0 ? [[0, 1], [2, 3]] : [[1, 2]];
        for (const [a, b] of pairs) {
          const wa = t2w[a], wb = t2w[b], lo = Math.min(wa, wb), hi = Math.max(wa, wb), my = w2t[lo], pt = w2t[hi], k = isw + 1;
          const bf = (pe[lo][k] - pe[hi][k]) * (1.0 / (1.0 * TEMPS[my]) - 1.0 / (1.0 * TEMPS[pt]));
          let sw = bf >= 0; if (!sw) { draws[lo]++; sw = rb[lo].uniform() < Math.exp(bf); }
          expected += Math.min(1, Math.exp(bf)); att[a][0]++; if (sw) { att[a][1]++; nw[lo] = pt; nw[hi] = my; }
        }
        w2t = nw; lines.push(`${step} ${w2t.join(' ')}`);
      }
      const same = lines.filter((l, i) => l === N.status[i]).length, acc = att.reduce((s, a) => s + a[1], 0);
      const temps = [0, 1, 2, 3].map(w => N.status.map(l => +l.split(' ')[1 + w]));
      $('.mr-atab').innerHTML = `<table>${head(['pair of temperatures', 'attempts', 'accepted'])}${att.map((a, i) => row([`${TEMPS[i]} ↔ ${TEMPS[i + 1]}`, a[0], a[1]])).join('')}</table>`;
      $('.mr-asum').innerHTML = `<p>${same} of ${N.status.length} status lines reproduced exactly ${same === N.status.length ? OK : BAD}, from the four replicas' potential energies (thermo, 17 digits, norm no) and the generator ported from ${gh('src/random_park.cpp', 43)}. Worlds 0–3 drew ${draws.join(', ')} random numbers. They were needed only when Δ &lt; 0, and only by the deciding world, so each world's stream advanced at its own pace, and the order matters.</p>
      <p>${acc} of 60 attempts accepted. The sum of min(1, e<sup>Δ</sup>) over the attempts, the number expected from these energies, is ${expected.toFixed(1)}. Every world visited ${[0, 1, 2, 3].map(w => new Set(temps[w]).size).join(', ')} of the 4 temperatures. World 0, for example, climbed from the lowest to the highest temperature between steps 2400 and 2700, and was back at the lowest by step 4300. That walk is what lets a replica cross barriers at high temperature and come back to sample the low one.</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.mr-stb'); b.disabled = true;
      const N = await nat(); if (!alive) return;
      const hdr = N.neblog.filter(l => /^\s*Step/.test(l)).pop().trim().split(/\s+/), last = N.neblog[N.neblog.length - 1].trim().split(/\s+/), col = k => Number(last[hdr.indexOf(k)]);
      const D = N.dumps.map(t => { const L = t.split('\n'), n = +L[3], box = [5, 6, 7].map(k => L[k].trim().split(/\s+/).map(Number)), at = new Map(); for (let a = 0; a < n; a++) { const w = L[9 + a].trim().split(/\s+/).map(Number); at.set(w[0], w); } return { box, at, n }; });
      const Lx = D[0].box[0][1] - D[0].box[0][0], neb = [...D[0].at.values()].filter(w => w[1] === 3).map(w => w[0]);
      const dist = (A, B) => { let s = 0; for (const id of neb) { const p = A.at.get(id), q = B.at.get(id); let dx = p[2] - q[2]; dx -= Lx * Math.round(dx / Lx); const dy = p[3] - q[3], dz = p[4] - q[4]; s += dx * dx + dy * dy + dz * dz; } return Math.sqrt(s); };
      const seg = [1, 2, 3].map(i => dist(D[i], D[i - 1])), cum = [0, seg[0], seg[0] + seg[1], seg[0] + seg[1] + seg[2]], RDT = cum[3];
      const pes = [];
      for (let k = 1; k <= 4; k++) {
        st.textContent = ` single point ${k}/4…`;
        const R = await runLammps(SP(k), { files: { 'initial.hop1': N.data, [`final.${k}.dump`]: N.dumps[k - 1] } }); if (!alive) return;
        if (R.error) { st.innerHTML = `<span class="lmp-err">${esc(R.error)}</span>`; b.disabled = false; return; }
        pes.push(Number(R.thermo.rows[R.thermo.rows.length - 1][1]));
      }
      b.disabled = false; st.textContent = '';
      const g8 = x => Number(x.toPrecision(8)), climb = pes.indexOf(Math.max(...pes));
      const rows = [0, 1, 2, 3].map(i => { const rd = cum[i] / RDT, pp = col(`PE${i + 1}`), rr = col(`RD${i + 1}`); return [i + 1, s17(pes[i]), pp, g8(pes[i]) === pp ? OK : BAD, rd.toFixed(8), rr, Math.abs(rd - rr) < 5e-8 ? OK : BAD]; });
      const ebf = pes[climb] - pes[0], ebr = pes[climb] - pes[3];
      $('.mr-btab').innerHTML = `<table>${head(['image', 'PE, single point here', 'PE printed', '', 'RD by hand', 'RD printed', ''])}${rows.map(row).join('')}</table>`;
      $('.mr-bsum').innerHTML = `<p>Each image's final configuration was read into your browser's LAMMPS (read_data of the starting file, then read_dump of the positions) and evaluated with <code>run 0</code>. The energies agree with the NEB output to all 8 printed digits, and the reaction coordinates computed here from the ${neb.length} fix neb atoms agree to the printed 8 decimals. Total path length RDT: ${RDT.toFixed(7)} here, ${col('RDT')} printed ${Math.abs(RDT - col('RDT')) < 5e-7 ? OK : BAD}. The images are not evenly spaced (segments ${seg.map(x => x.toFixed(3)).join(', ')}): the climbing image moved up to the saddle and away from the middle.</p>
      <p>Climbing image ${climb + 1}: EBF = ${ebf.toPrecision(8)} (printed ${col('EBF')}), EBR = ${ebr.toPrecision(8)} (printed ${col('EBR')}). Per atom. The hop's barrier is EBF × 420 = ${(ebf * 420).toFixed(3)} ε, and the reverse ${(ebr * 420).toFixed(3)} ε. A barrier quoted without that factor is off by 420. With these energies and a prefactor you get a rate by transition-state theory; prd and tad, also in REPLICA, find such events by themselves.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
