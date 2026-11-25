// Research · D.4: communication between processors. CommBrick::setup and borders emulated swap by swap: ghost counts for any comm cutoff (live) and any processor grid (native mpirun), multi-hop swaps, and newton off by hand.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const grab = (lines, re) => { for (const l of lines) { const m = l.match(re); if (m) return m; } return null; };

const A = 1.6795961913825073, L = 10 * A, EPS = 1.0e-6, BIG = 1e20;
const X = []; for (let i = 0; i < 10; i++) for (const b of [0, 0.5]) X.push([(i + b) * A, 2 * i + (b ? 1 : 0)]);   // 20 planes: [coordinate, half-lattice index]

// one axis split over P procs: owned planes (create_atoms rule, page B.7), then the borders() swaps of comm_brick.cpp
function axis(P, c) {
  const lo = [], hi = [], own = [];
  for (let k = 0; k < P; k++) { lo.push(0 + L * (k * 1.0 / P)); hi.push(k < P - 1 ? 0 + L * ((k + 1) * 1.0 / P) : L); }
  for (let k = 0; k < P; k++) { const clo = lo[k] - (k === 0 ? EPS : 0), chi = hi[k] - (k === P - 1 ? 2 * EPS : 0); own.push(X.filter(([x]) => x >= clo && x < chi)); }
  const co = own.map(o => [...o]), m = Math.floor(c * P / L) + 1, nf = Array(P).fill(0), nl = Array(P).fill(0), sent = [];
  for (let ineed = 0; ineed < 2 * m; ineed++) {
    if (ineed % 2 === 0) for (let k = 0; k < P; k++) { nf[k] = nl[k]; nl[k] = co[k].length; }
    const recv = Array.from({ length: P }, () => []); let ns = 0;
    for (let k = 0; k < P; k++) {
      let a, b, sh, dest;
      if (ineed % 2 === 0) { a = ineed < 2 ? -BIG : 0.5 * (lo[k] + hi[k]); b = lo[k] + c; sh = k === 0 ? L : 0.0; dest = (k - 1 + P) % P; }
      else { a = hi[k] - c; b = ineed < 2 ? BIG : 0.5 * (lo[k] + hi[k]); sh = k === P - 1 ? -L : 0.0; dest = (k + 1) % P; }
      const s = co[k].slice(nf[k], nl[k]).filter(([x]) => x >= a && x <= b).map(([x, p]) => [x + sh, p]);
      ns += s.length; recv[dest].push(...s);
    }
    for (let k = 0; k < P; k++) co[k].push(...recv[k]);
    sent.push(ns / P);
  }
  return { own, co, m, sent };
}
const sites = (a, b, c) => { let n = 0; for (const [, i] of a) for (const [, j] of b) for (const [, k] of c) if ((i + j + k) % 2 === 0) n++; return n; };
function emulate(g, c = 2.8) {
  const ax = g.map(P => axis(P, c)), nl = [], ng = [];
  for (let a = 0; a < g[0]; a++) for (let b = 0; b < g[1]; b++) for (let d = 0; d < g[2]; d++) {
    const o = sites(ax[0].own[a], ax[1].own[b], ax[2].own[d]); nl.push(o); ng.push(sites(ax[0].co[a], ax[1].co[b], ax[2].co[d]) - o);
  }
  const st = v => [v.reduce((s, x) => s + x, 0) / v.length, Math.max(...v), Math.min(...v)];
  return { nl: st(nl), ng: st(ng), m: ax.map(x => x.m), planes: ax.map(x => x.co.map(q => q.length)), ax };
}
const f6 = x => String(+x.toPrecision(6));

// measured natively: melt, run 0, LAMMPS 22 Jul 2025, MPICH, 2-core Linux; 1×1×P grids forced with "processors 1 1 P"
const NATIVE = [
  [[1, 1, 1], [4000, 4000, 4000], [5841, 5841, 5841]], [[1, 1, 2], [2000, 2000, 2000], [4196, 4196, 4196]],
  [[1, 1, 3], [1333.33, 1400, 1200], [3648.33, 3903, 3339]], [[1, 2, 2], [1000, 1000, 1000], [2901, 2901, 2901]],
  [[1, 1, 4], [1000, 1200, 800], [3374, 3574, 3174]], [[1, 1, 5], [800, 1000, 600], [3209, 3409, 3009]],
  [[1, 2, 3], [666.667, 700, 600], [2470, 2613, 2284]], [[1, 1, 6], [666.667, 800, 600], [3100, 3410, 2845]],
  [[2, 2, 2], [500, 500, 500], [1956, 1956, 1956]], [[1, 1, 8], [500, 600, 400], [2780.5, 2881, 2680]],
  [[2, 2, 3], [333.333, 350, 300], [1641.67, 1723, 1529]], [[2, 2, 4], [250, 300, 200], [1484, 1534, 1434]],
  [[1, 1, 16], [250, 400, 200], [2666, 2716, 2516]], [[3, 3, 3], [148.148, 172, 108], [1128.19, 1264, 928]],
];

const MELT = `units           lj
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
neigh_modify    every 20 delay 0 check no
fix             1 all nve
`;
// pairs whose minimum image crosses the box: each is stored twice with newton off
function crossingPairs() {
  const offs = [];
  for (let i = -4; i <= 4; i++) for (let j = -4; j <= 4; j++) for (let k = -4; k <= 4; k++) if ((i || j || k) && ((i + j + k) % 2 + 2) % 2 === 0 && (i * i + j * j + k * k) * (A / 2) ** 2 < 2.8 * 2.8) offs.push([i, j, k]);
  let cross = 0;
  for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) for (let k = 0; k < 20; k++) { if ((i + j + k) % 2) continue; for (const [a, b, c] of offs) { const t = [i + a, j + b, k + c]; if (t.some(v => v < 0 || v > 19)) cross++; } }
  return { per: offs.length, contacts: cross, pairs: cross / 2 };
}

export default {
  id: 'communication',
  title: 'Communication between processors',
  minutes: 25,
  summary: 'Ghost atoms are made in 2 to 2·(n+1) swaps per dimension, each sending a slab of atoms to a neighbouring processor. You emulate comm_brick.cpp swap by swap and predict the exact number of ghosts for any communication cutoff (checked live, up to swaps of ghosts of ghosts) and for 14 processor grids (checked against real mpirun runs, including bricks thinner than the cutoff). Then newton off: by hand, 26 352 pairs cross the box and are stored twice.',
  learn: [
    'CommBrick::setup: maxneed, nswap and the slab of each swap',
    'borders(): x, then y, then z, so edges and corners come for free',
    'When a brick is thinner than the cutoff: two or three hops',
    'Each step: forward_comm sends ghost positions, reverse_comm returns ghost forces',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: swaps and slabs</h2>
    <div class="prose"><p>${gh('src/comm_brick.cpp', 193)} decides how many processors away ghosts must come from, per dimension:</p>
<pre class="lmp-mini">maxneed[0] = static_cast&lt;int&gt;(cutghost[0] * procgrid[0] / prd[0]) + 1;   // 310
nswap = 2 * (maxneed[0] + maxneed[1] + maxneed[2]);                        // 394</pre>
<p>Each swap sends a slab: to the west, atoms with x ≤ sublo + cutghost; to the east, atoms with x ≥ subhi − cutghost. Later swaps (ineed ≥ 2) send on what just arrived, cut at the brick's midpoint (lines 427–458). A slab that crosses the periodic boundary is shifted by ± the box length on the way. ${gh('src/comm_brick.cpp', 775)}, trimmed:</p>
<pre class="lmp-mini">for (dim = 0; dim &lt; 3; dim++) {
  for (ineed = 0; ineed &lt; 2*maxneed[dim]; ineed++) {
    if (ineed % 2 == 0) { nfirst = nlast; nlast = nlocal + nghost; }   // first pair: owned + ghosts so far
    for (i = nfirst; i &lt; nlast; i++)
      if (x[i][dim] &gt;= lo &amp;&amp; x[i][dim] &lt;= hi) sendlist[iswap][nsend++] = i;
    ... send, receive, append as ghosts ...
  }
}</pre>
<p>Because y is done after x, the y-slabs include the x-ghosts, so edge and corner ghosts arrive without extra messages. The sendlists are kept: every step, ${gh('src/comm_brick.cpp', 540)} sends the same atoms' new positions (forward), and with newton on ${gh('src/comm_brick.cpp', 600)} sends the forces on ghosts back (reverse). ${gh('src/comm_brick.cpp', 654)} (atoms that left the brick) runs only on rebuild steps.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: ghosts for any communication cutoff</h2>
    <div class="prose"><p>One processor still does all the swaps, sending to itself across the periodic boundary. <code>comm_modify cutoff c</code> sets cutghost = max(c, 2.8). Melt at step 0: atoms sit on 20 planes per axis. Predict the planes after the swaps, then count FCC sites:</p></div>
    <div class="row"><label class="lc-in">comm cutoff<select class="lc-u cm-c">${[2.8, 4, 6, 10, 17, 20, 30].map(c => `<option${c === 2.8 ? ' selected' : ''}>${c}</option>`).join('')}</select></label><button class="btn sm primary" data-k="a">▶ run 0</button><span class="lmp-status muted cm-sta"></span></div>
    <div class="calc cm-pred"></div>
    <div class="prose cm-outa"></div>
  </section>

  <section class="block">
    <h2>3 · Match: newton off</h2>
    <div class="prose"><p>With <code>newton off</code> nothing is sent back (no reverse_comm), so a pair between an owned atom and a ghost must be computed on both sides: it is stored twice. On one processor, those are exactly the pairs whose closest image crosses the box. Total # of neighbors = 156 000 + (crossing pairs).</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ run 0 with newton off</button><span class="lmp-status muted cm-stb"></span></div>
    <div class="calc cm-predb"></div>
    <div class="prose cm-outb"></div>
  </section>

  <section class="block">
    <h2>4 · Match: real processor grids</h2>
    <div class="prose"><p>The same emulation, now with P bricks per axis, against <code>mpirun -np P</code> on a 2-core computer (run 0, numbers from the log). The 1 × 1 × P grids were forced with <code>processors 1 1 P</code>: at P = 6 each brick is 2.7993 thick, just under the 2.8 cutoff, so two hops are needed; at P = 16, three.</p></div>
    <div class="tbl cm-outc"></div>
    <div class="prose cm-sumc"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Grid 1 × 1 × 6 for melt: maxneed and nswap?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>z: int(2.8 × 6 / 16.796) + 1 = int(1.0002) + 1 = 2; x, y: 1. nswap = 2 × (1 + 1 + 2) = <b>8</b>. A 0.02 % thinner cutoff (2.7993) would give 6.</div></div>
    <div class="ex"><p><b>2.</b> How many bytes does forward_comm send per step for P = 8 (2 × 2 × 2)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each processor sends the positions of the atoms that are ghosts elsewhere: 1956 per processor × 3 doubles × 8 bytes ≈ <b>47 kB</b> per processor per step; reverse_comm returns the same amount of forces.</div></div>
    <div class="ex"><p><b>3.</b> Why does newton off need no reverse communication but more pair work?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Each processor computes the owned–ghost pair force itself and keeps only its own half; the other processor does the same pair. 26 352 of the 182 352 pairs here are computed twice (17 % more), in exchange for skipping one message round per step.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const predA = () => {
      const c = +$('.cm-c').value, cg = Math.max(c, 2.8), E = emulate([1, 1, 1], cg), ax = E.ax[2];
      const lo = ax.co[0].filter(([x]) => x < 0).length, hiN = ax.co[0].filter(([x]) => x >= L).length;
      $('.cm-pred').textContent = `cutghost = ${cg};  maxneed = int(${cg} × 1 / ${L.toFixed(4)}) + 1 = ${E.m[2]};  nswap = 2 × 3 × ${E.m[2]} = ${6 * E.m[2]}\nper axis: 20 owned planes + ${lo} below 0 + ${hiN} above L = ${E.planes[2][0]} planes${E.m[2] > 1 ? '  (the second hop sends ghosts of ghosts)' : ''}\nFCC sites in ${E.planes[2][0]}³ planes − 4000 owned = ${E.ng[0]} ghosts`;
      $('.cm-outa').innerHTML = ''; return E;
    };
    $('.cm-c').onchange = predA; predA();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.cm-sta'), c = +$('.cm-c').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(`${MELT}${c !== 2.8 ? `comm_modify     cutoff ${c}\n` : ''}run             0`); b.disabled = false; if (res.error) return err(st, res); st.textContent = '';
      const E = predA(), m = grab(res.lines, /Nghost:\s+(\d+)/), gc = grab(res.lines, /ghost atom cutoff = ([\d.]+)/);
      $('.cm-outa').innerHTML = `<p>LAMMPS: ghost atom cutoff = ${gc ? gc[1] : '?'}, Nghost = <b>${m ? m[1] : '?'}</b> ${m && +m[1] === E.ng[0] ? OK : BAD}</p>`;
    };

    const C = crossingPairs();
    $('.cm-predb').textContent = `neighbours within 2.8 per atom: ${C.per}\ncontacts that leave the 20 × 20 × 20 planes (cross the box): ${C.contacts}\nso ${C.pairs} crossing pairs; newton on stores 156000 pairs, newton off ${156000 + C.pairs}`;
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.cm-stb'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(`${MELT.replace('atom_style      atomic', 'newton          off\natom_style      atomic')}run             0`); b.disabled = false; if (res.error) return err(st, res); st.textContent = '';
      const t = grab(res.lines, /Total # of neighbors = (\d+)/), pb = grab(res.lines, /pair build: (\S+)/);
      $('.cm-outb').innerHTML = `<p>LAMMPS: pair build ${pb ? `<code>${esc(pb[1])}</code>` : ''}, Total # of neighbors = <b>${t ? t[1] : '?'}</b> ${t && +t[1] === 156000 + C.pairs ? OK : BAD}</p>`;
    };

    let ok = 0;
    $('.cm-outc').innerHTML = `<table>${head(['grid', 'maxneed', 'nswap', 'Nlocal ave/max/min by hand', 'Nghost ave/max/min by hand', 'mpirun', ''])}${NATIVE.map(([g, nl, ng]) => {
      const E = emulate(g), hn = E.nl.map(f6), hg = E.ng.map(f6), same = hn.every((v, i) => +v === nl[i]) && hg.every((v, i) => +v === ng[i]); if (same) ok++;
      return row([g.join(' × '), E.m.join(' '), 2 * E.m.reduce((s, x) => s + x, 0), hn.join(' / '), hg.join(' / '), `${nl.join(' / ')}<br>${ng.join(' / ')}`, same ? OK : BAD]); }).join('')}</table>`;
    $('.cm-sumc').innerHTML = `<p>${ok} of ${NATIVE.length} grids: every average, maximum and minimum identical (6 significant digits, as LAMMPS prints them). The uneven splits come from lattice planes landing exactly on brick boundaries (page B.7).</p>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
