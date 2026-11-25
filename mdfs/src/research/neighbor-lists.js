// Research · D.3: neighbour lists. When to rebuild (decide, check_distance), how big the bins are, how the half list is filled; rebuild counts predicted by hand, and the box-corner rule that is not the "two largest" its comment promises.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const grab = (lines, re) => { for (const l of lines) { const m = l.match(re); if (m) return m; } return null; };

// A: one moving atom, one far away
const scriptA = p => `units           lj
atom_style      atomic
region          box block 0 20 0 20 0 20 units box
create_box      1 box
create_atoms    1 single 5.0 10.0 10.0 units box
create_atoms    1 single 15.0 10.0 10.0 units box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
set             atom 1 vx ${p.v}
neighbor        ${p.skin} bin
neigh_modify    every ${p.every} delay ${p.delay} check ${p.check ? 'yes' : 'no'}
fix             1 all nve
thermo          200
run             200`;
function handA(p) {
  let x = 5.0, hold = 5.0, ago = 0, builds = 0, danger = 0; const at = [];
  const trig = 0.5 * p.skin, trigsq = trig * trig;
  for (let s = 1; s <= 200; s++) {
    x += 0.005 * p.v; ago++;
    let rebuild = 0;
    if (ago >= p.delay && ago % p.every === 0) {
      if (!p.check) rebuild = 1;
      else { const d = x - hold, rsq = d * d; const flag = rsq > trigsq ? 1 : 0; if (flag && ago === Math.max(p.every, p.delay)) danger++; rebuild = flag; }
    }
    if (rebuild) { builds++; at.push(s); hold = x; ago = 0; }
  }
  return { builds, danger, at };
}

// B: triclinic box growing with lo fixed, atom remapped with it
const scriptB = p => `units           lj
atom_style      atomic
region          box prism 0 10 0 10 0 10 0 0 0 units box
create_box      1 box
create_atoms    1 single ${p.f.map(f => 10 * f).join(' ')} units box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        1.0 bin
neigh_modify    every 1 delay 0 check yes
fix             1 all nve
fix             2 all deform 1 x final 0 ${p.L[0]} y final 0 ${p.L[1]} z final 0 ${p.L[2]} remap x
thermo          1000
run             1000`;
function handB(p, rule) {
  const N = 1000, L0 = [10, 10, 10], box = s => [0, 1, 2].map(k => L0[k] + (p.L[k] - L0[k]) * s / N);
  const corners = L => { const c = []; for (const k of [0, 1]) for (const j of [0, 1]) for (const i of [0, 1]) c.push([i * L[0], j * L[1], k * L[2]]); return c; };
  let hc = corners(box(0)), hx = p.f.map((f, k) => f * L0[k]), builds = 0; const at = [];
  for (let s = 1; s <= N; s++) {
    const L = box(s), x = p.f.map((f, k) => f * L[k]), c = corners(L);
    const ds = c.map((q, i) => Math.hypot(q[0] - hc[i][0], q[1] - hc[i][1], q[2] - hc[i][2]));
    let d1 = 0, d2 = 0;
    if (rule === 'code') { for (const d of ds) { if (d > d1) d1 = d; else if (d > d2) d2 = d; } }
    else { const sd = [...ds].sort((a, b) => b - a); d1 = sd[0]; d2 = sd[1]; }
    let delta = 0.5 * (1.0 - (d1 + d2)); if (delta < 0) delta = 0;
    const rsq = x.reduce((m, q, k) => m + (q - hx[k]) ** 2, 0);
    if (rsq > delta * delta) { builds++; at.push(s); hc = c; hx = x; }
  }
  return { builds, at };
}

export default {
  id: 'neighbor-lists',
  title: 'Neighbour lists',
  minutes: 25,
  summary: 'Three questions the Neighbor class answers every step: rebuild now? how big are the bins? which pairs go in the list? You predict the number of rebuilds and "dangerous builds" for any every/delay/check/skin by hand and match LAMMPS; read the binning and the half-list rules; and find that for a changing triclinic box the "sum of the 2 largest corner moves" in the comment is not what the code computes, which is shown live by LAMMPS rebuilding fewer times than the stated rule would.',
  learn: [
    'decide() and check_distance(): every, delay, check, the half-skin trigger, dangerous builds',
    'Bins of half the list cutoff, and why the log\'s "bins =" is not the bin count used',
    'The half list with newton on: own bin forward, ghosts above-right, upper-half stencil',
    'A box-change check that is less conservative than documented',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: rebuild now?</h2>
    <div class="prose"><p>${gh('src/neighbor.cpp', 2363)} is called at stage 3 of every step (page C.2):</p>
<pre class="lmp-mini">ago++;
if (ago &gt;= delay &amp;&amp; ago % every == 0) {
  if (dist_check == 0) return 1;            // check no: rebuild
  return check_distance();
} else return 0;</pre>
<p>${gh('src/neighbor.cpp', 2393)}, for a box that does not change:</p>
<pre class="lmp-mini">rsq = |x[i] − xhold[i]|²;                  // xhold = positions at the last build
if (rsq &gt; triggersq) { flag = 1; break; }   // triggersq = (skin/2)²
if (flagall &amp;&amp; ago == MAX(every,delay)) ndanger++;   // 2443: "dangerous"</pre>
<p>Half the skin, because two atoms can each move skin/2 towards each other. A build is <b>dangerous</b> when it is needed on the very first step a check was allowed: the atoms may have crossed the trigger earlier, unseen.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: count the rebuilds</h2>
    <div class="prose"><p>Atom 1 moves at constant vx; atom 2 is 10 away, so there are no forces. |vx| ≤ 7 keeps atom 1 more than 2.5 from atom 2 (and its periodic image) for all 200 steps; faster, and the force would change its path. 200 steps, dt = 0.005. The build in setup is not counted.</p></div>
    <div class="row">
      <label class="lc-in">vx<input type="number" class="lmp-in sm nl-v" value="4" step="0.5" min="-7" max="7"></label>
      <label class="lc-in">skin<select class="lc-u nl-skin"><option>0.1</option><option selected>0.3</option><option>0.6</option><option>1.0</option></select></label>
      <label class="lc-in">every<select class="lc-u nl-every"><option selected>1</option><option>2</option><option>5</option><option>10</option></select></label>
      <label class="lc-in">delay<select class="lc-u nl-delay"><option selected>0</option><option>3</option><option>10</option></select></label>
      <label class="lc-in">check<select class="lc-u nl-check"><option value="1">yes</option><option value="0">no</option></select></label>
      <button class="btn sm primary" data-k="a">▶ Run 200 steps</button><span class="lmp-status muted nl-sta"></span></div>
    <div class="calc nl-pred"></div>
    <div class="prose nl-outa"></div>
  </section>

  <section class="block">
    <h2>3 · Understand: bins and the half list</h2>
    <div class="prose"><p><b>Bins.</b> ${gh('src/nbin_standard.cpp', 124)}: binsize = ½ (cutoff + skin), and</p>
<pre class="lmp-mini">nbinx = static_cast&lt;int&gt;(bbox[0]*binsizeinv);   // 139: rounds DOWN
binsizex = bbox[0]/nbinx;                         // 154: bins a bit larger than asked</pre>
<p>The log line comes from a different formula, ${gh('src/neighbor.cpp', 1807)}: <code>ceil(bbox/binsize)</code>, rounded up. For melt: 16.796 / 1.4 = 11.997, so the log says <b>bins = 12 12 12</b> while the code builds <b>11</b> bins of 1.527 per direction (plus ghost bins). The log line is a report, not the bin count.</p>
<p><b>The half list</b> with newton on (${gh('src/npair_bin.cpp', 112)}): for atom i, loop over a stencil of bins that covers the upper half of space. In i's own bin, take only the atoms after i in the bin's linked list, and a ghost only if it is "above and to the right" of i (higher z; same z and higher y; same z and y and higher x). In the other stencil bins, take everything; then keep the pair if r² &lt; (cutoff + skin)². Each pair is stored exactly once.</p></div>
    <div class="tbl nl-bins"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: a deforming triclinic box</h2>
    <div class="prose"><p>When the box changes, the trigger shrinks by how far the box moved (${gh('src/neighbor.cpp', 2381)}). The comment says: "reduce skin distance by sum of 2 largest of the 8 values". The code, ${gh('src/neighbor.cpp', 2414)}:</p>
<pre class="lmp-mini">for (int i = 0; i &lt; 8; i++) {
  delta = |corner[i] − corner_hold[i]|;
  if (delta &gt; delta1) delta1 = delta;        // the old delta1 is lost here
  else if (delta &gt; delta2) delta2 = delta;   // so delta2 is not always 2nd largest
}
delta = 0.5 * (skin - (delta1+delta2));</pre>
<p>Corners come in the order (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), … When the displacements grow along that order, each new maximum overwrites delta1 and the old one is not moved to delta2. delta2 ends up too small, the trigger too large, and LAMMPS rebuilds less often than the stated rule. Test it: a prism box 10³ grows to the sizes below with <code>fix deform … remap x</code>, skin 1.0, check every step, over 1000 steps. You predict the build count both ways:</p></div>
    <div class="row">
      <label class="lc-in">final box<select class="lc-u nl-L"><option value="10.5,11,12" selected>10.5 × 11 × 12</option><option value="12,11,10.5">12 × 11 × 10.5</option><option value="11,12,10.5">11 × 12 × 10.5</option></select></label>
      <label class="lc-in">atom at (fraction)<select class="lc-u nl-f"><option value="0.5,0.5,0.5" selected>0.5 0.5 0.5</option><option value="0.9,0.9,0.9">0.9 0.9 0.9</option><option value="0.3,0.6,0.9">0.3 0.6 0.9</option></select></label>
      <button class="btn sm primary" data-k="b">▶ Run 1000 steps</button><span class="lmp-status muted nl-stb"></span></div>
    <div class="tbl nl-outb"></div>
    <div class="prose nl-sumb"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> vx = 4, skin 0.3, every 1, delay 0, check yes. On which step is the first rebuild?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The atom moves 0.02 per step; it must move more than 0.15: after 7 steps it has moved 0.14, after 8 steps 0.16, so <b>step 8</b>. Then every 8 steps: 25 builds in 200 steps.</div></div>
    <div class="ex"><p><b>2.</b> Same, but every 10. Why are all builds dangerous?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The first check is at ago = 10, when the atom has moved 0.2 &gt; 0.15: it needed a rebuild at step 8 already. Every build happens on the first allowed check, so every build counts as dangerous: pairs could have been missed in steps 8–9.</div></div>
    <div class="ex"><p><b>3.</b> Corner moves arrive as 1, 2, 3, 0.5 (in that order). What do delta1 and delta2 end up as, in the code and by the comment?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Code: 1 → d1 = 1; 2 → d1 = 2; 3 → d1 = 3; 0.5 → d2 = 0.5. Sum 3.5. Comment ("two largest"): 3 + 2 = 5. The code reduces the trigger by 1.5 less than promised.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pA = () => { let v = +$('.nl-v').value; if (!isFinite(v)) v = 4; v = Math.max(-7, Math.min(7, v)); return { v, skin: +$('.nl-skin').value, every: +$('.nl-every').value, delay: +$('.nl-delay').value, check: $('.nl-check').value === '1' }; };
    const predA = () => { const p = pA(); if (p.delay % p.every !== 0) { $('.nl-pred').innerHTML = `delay ${p.delay} is not a multiple of every ${p.every}: LAMMPS will stop with "Neighbor delay must be 0 or multiple of every setting" (${gh('src/neighbor.cpp', 327)}). Run it to see.`; $('.nl-outa').innerHTML = ''; return null; } const H = handA(p); $('.nl-pred').textContent = `moves ${Math.abs(0.005 * p.v).toFixed(4)} per step; trigger ${p.check ? `|Δx| > skin/2 = ${p.skin / 2}` : 'none (check no)'}; first check when ago ≥ ${p.delay} and ago % ${p.every} = 0\nyour prediction: ${H.builds} builds, ${p.check ? H.danger + ' dangerous' : 'dangerous not checked'}\nbuild steps: ${H.at.slice(0, 12).join(', ')}${H.at.length > 12 ? ', …' : ''}`; $('.nl-outa').innerHTML = ''; return H; };
    root.querySelectorAll('.nl-v,.nl-skin,.nl-every,.nl-delay,.nl-check').forEach(e => e.oninput = e.onchange = predA); predA();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.nl-sta'), p = pA(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptA(p)); b.disabled = false; const H = predA(); if (res.error) { err(st, res); if (!H) $('.nl-outa').innerHTML = `<p>As predicted: LAMMPS refused the input. ${OK}</p>`; return; } st.textContent = '';
      const nb = grab(res.lines, /Neighbor list builds = (\d+)/), dg = grab(res.lines, /Dangerous builds = (\d+)/), bi = grab(res.lines, /binsize = ([\d.]+), bins = (\d+) (\d+) (\d+)/);
      const okb = nb && +nb[1] === H.builds, okd = p.check ? (dg && +dg[1] === H.danger) : !dg;
      $('.nl-outa').innerHTML = `<p>LAMMPS: Neighbor list builds = <b>${nb ? nb[1] : '?'}</b> ${okb ? OK : BAD}; ${dg ? `Dangerous builds = <b>${dg[1]}</b>` : 'Dangerous builds not checked'} ${okd ? OK : BAD}.${bi ? ` Log: binsize = ${bi[1]}, bins = ${bi[2]} ${bi[3]} ${bi[4]}.` : ''}</p>`;
    };

    const bins = () => {
      const L = 10 * 1.6795961913825073;
      $('.nl-bins').innerHTML = `<table>${head(['melt, skin', 'binsize = (2.5 + skin)/2', 'L / binsize', 'log says (ceil)', 'code builds (int)', 'actual bin width'])}${[0.1, 0.3, 0.5, 1.0].map(s => { const bs = 0.5 * (2.5 + s), q = L * (1 / bs), n = Math.max(1, Math.floor(q)); return row([s, bs, q.toFixed(4), Math.ceil(L / bs), n, (L / n).toFixed(4)]); }).join('')}</table><p class="muted">The "code builds" column counts bins across the box itself; the code then adds bins for the ghost region on each side.</p>`;
    };
    bins();

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.nl-stb'), p = { L: $('.nl-L').value.split(',').map(Number), f: $('.nl-f').value.split(',').map(Number) }; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptB(p)); b.disabled = false; if (res.error) return err(st, res); st.textContent = '';
      const nb = grab(res.lines, /Neighbor list builds = (\d+)/), C = handB(p, 'code'), D = handB(p, 'doc'), n = nb ? +nb[1] : NaN;
      $('.nl-outb').innerHTML = `<table>${head(['rule', 'builds by hand', 'build steps', 'LAMMPS'])}${row(['as the code is written', C.builds, C.at.join(', '), n === C.builds ? `${n} ${OK}` : n])}${row(['as the comment says (2 largest)', D.builds, D.at.join(', '), n === D.builds ? `${n} ${OK}` : n])}</table>`;
      $('.nl-sumb').innerHTML = C.builds === D.builds ? '<p>For this box the two rules agree. Try 10.5 × 11 × 12: z grows most, and the z corners come last.</p>' : `<p>LAMMPS rebuilt ${n} times: the code's rule. The documented rule asks for ${D.builds}. Missed rebuilds mean the list can lack pairs that came within the cutoff. It only matters for triclinic boxes that change shape or size (npt, deform). The same lines are in the current development version; as far as I found, it has not been reported.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
