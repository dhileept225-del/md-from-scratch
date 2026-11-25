// Research · D.2: domain and periodic boundaries. Domain::pbc, image flags in 10 bits, unwrapped coordinates, triclinic h and h_inv with the lamda round trip; all reproduced bit for bit, including an image-flag overflow.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const IMGMASK = 1023, IMGMAX = 512, IMGBITS = 10;

const parseDump = t => t.split('ITEM: TIMESTEP').slice(1).map(f => { const ls = f.trim().split('\n'); return { step: +ls[0], w: ls[ls.length - 1].trim().split(/\s+/).map(Number) }; });

// orthogonal box 0..5, one free atom
const scriptO = (vx, vy, ev, n, cols = 'id x y ix iy xu yu', de = 1) => `units           lj
atom_style      atomic
region          box block 0 5 0 5 0 5 units box
create_box      1 box
create_atoms    1 single 4.9 2.5 2.5 units box
mass            1 1.0
velocity        all set ${vx} ${vy} 0.0 units box
comm_modify     cutoff 1.0
neigh_modify    every ${ev} delay 0 check no
fix             1 all nve
dump            d all custom ${de} o.dump ${cols}
dump_modify     d format float %.17g
run             ${n}`;
// the image int is an unsigned 10-bit field per dimension; wrap like the C++ does
const imgInc = (im, d, s) => { const sh = d * IMGBITS, idim = (im >>> sh) & IMGMASK, other = (im ^ (idim << sh)) >>> 0; return (other | ((((idim + s) & IMGMASK) << sh) >>> 0)) >>> 0; };
const ix = (im, d) => ((im >>> (d * IMGBITS)) & IMGMASK) - IMGMAX;
const IMG0 = ((IMGMAX << 20) | (IMGMAX << 10) | IMGMAX) >>> 0;
function handO(vx, vy, ev, n, L = 5.0) {
  const x = [4.9, 2.5, 2.5], v = [vx, vy, 0.0]; let im = IMG0; const out = [];
  const rec = s => out.push({ s, x: x[0], y: x[1], ix: ix(im, 0), iy: ix(im, 1), xu: x[0] + ix(im, 0) * L, yu: x[1] + ix(im, 1) * L, im });
  rec(0);
  for (let s = 1; s <= n; s++) {
    for (let k = 0; k < 3; k++) x[k] += 0.005 * v[k];
    if (s % ev === 0) for (let k = 0; k < 3; k++) {
      if (x[k] < 0.0) { x[k] += L; im = imgInc(im, k, -1); }
      if (x[k] >= L) { x[k] -= L; x[k] = Math.max(x[k], 0.0); im = imgInc(im, k, +1); }
    }
    rec(s);
  }
  return out;
}

// triclinic prism 0..5 with tilt xy
const scriptT = (xy, vy) => `units           lj
atom_style      atomic
region          box prism 0 5 0 5 0 5 ${xy} 0.0 0.0 units box
create_box      1 box
create_atoms    1 single 3.3 4.6 2.5 units box
mass            1 1.0
velocity        all set 3.0 ${vy} 0.0 units box
comm_modify     cutoff 1.0
neigh_modify    every 1 delay 0 check no
fix             1 all nve
dump            d all custom 1 t.dump id x y xs ys ix iy
dump_modify     d format float %.17g
run             6`;
function handT(xy, vy) {
  const h = [5.0, 5.0, 5.0, 0.0, 0.0, xy];
  const hi = [1.0 / h[0], 1.0 / h[1], 1.0 / h[2], -h[3] / (h[1] * h[2]), (h[3] * h[5] - h[1] * h[4]) / (h[0] * h[1] * h[2]), -h[5] / (h[0] * h[1])];
  let x = [3.3, 4.6, 2.5]; const v = [3.0, vy, 0.0], im = [0, 0, 0], out = [];
  const xs = () => [hi[0] * (x[0] - 0) + hi[5] * (x[1] - 0) + hi[4] * (x[2] - 0), hi[1] * (x[1] - 0) + hi[3] * (x[2] - 0)];
  const wrap = () => { // Verlet::setup and every rebuild step: x2lamda, pbc in lamda, lamda2x
    const d = [x[0] - 0.0, x[1] - 0.0, x[2] - 0.0], l = [hi[0] * d[0] + hi[5] * d[1] + hi[4] * d[2], hi[1] * d[1] + hi[3] * d[2], hi[2] * d[2]];
    for (let k = 0; k < 3; k++) { if (l[k] < 0.0) { l[k] += 1.0; im[k]--; } if (l[k] >= 1.0) { l[k] -= 1.0; l[k] = Math.max(l[k], 0.0); im[k]++; } }
    x = [h[0] * l[0] + h[5] * l[1] + h[4] * l[2] + 0.0, h[1] * l[1] + h[3] * l[2] + 0.0, h[2] * l[2] + 0.0];
  };
  let plain = x[0]; wrap();
  out.push({ s: 0, x: x[0], y: x[1], xs: xs(), ix: im[0], iy: im[1], plain });
  for (let s = 1; s <= 6; s++) {
    for (let k = 0; k < 3; k++) x[k] += 0.005 * v[k];
    plain = x[0]; wrap();
    out.push({ s, x: x[0], y: x[1], xs: xs(), ix: im[0], iy: im[1], plain });
  }
  return { out, hi };
}

export default {
  id: 'domain-pbc',
  title: 'Domain and periodic boundaries',
  minutes: 20,
  summary: 'Domain::pbc is twenty lines: add or subtract one box length and count it in an image flag packed into 10 bits. You follow one fast atom through an orthogonal box and a tilted one, reproduce every wrapped and unwrapped coordinate bit for bit, see a triclinic coordinate change in its last digit just from the trip into box fractions and back, and make an image flag overflow live: after 512 crossings the unwrapped position jumps by 1024 box lengths.',
  learn: [
    'Domain::pbc: one shift per dimension, only on rebuild steps',
    'Image flags: three 10-bit fields in one int, and unwrapped xu = x + ix · Lx',
    'Triclinic: h, h_inv, and the x → λ → x round trip every rebuild',
    'The default build can count only −512 … 511 crossings',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: wrap, count, unwrap</h2>
    <div class="prose"><p>${gh('src/domain.cpp', 769)}, the x part (y and z are the same with shifted bits):</p>
<pre class="lmp-mini">if (x[i][0] &lt; lo[0]) {                 // 808
  x[i][0] += period[0];
  idim = image[i] &amp; IMGMASK;            // the low 10 bits
  otherdims = image[i] ^ idim;
  idim--;  idim &amp;= IMGMASK;             // count down, modulo 1024
  image[i] = otherdims | idim;
}
if (x[i][0] &gt;= hi[0]) {                // 817
  x[i][0] -= period[0];
  x[i][0] = MAX(x[i][0],lo[0]);
  ...  idim++;  idim &amp;= IMGMASK; ...
}</pre>
<p>It runs only on neighbour-rebuild steps (${gh('src/verlet.cpp', 277)}). The image int (${gh('src/lmptype.h', 122)}): IMGMASK = 1023, IMGMAX = 512, IMGBITS = 10, so</p>
<pre class="lmp-mini">image = (iz + 512) &lt;&lt; 20  |  (iy + 512) &lt;&lt; 10  |  (ix + 512)
xu    = x + ((image &amp; IMGMASK) - IMGMAX) * xprd        // dump_custom.cpp:2665</pre>
<p>For a tilted (triclinic) box, ${gh('src/domain.cpp', 257)} builds h = (Lx, Ly, Lz, yz, xz, xy) and its inverse, and pbc works in box fractions λ = h⁻¹(x − lo), period 1:</p>
<pre class="lmp-mini">h_inv[5] = -h[5] / (h[0]*h[1]);                         // 280
x2lamda: λx = h_inv[0]·Δx + h_inv[5]·Δy + h_inv[4]·Δz   // 2241
lamda2x: x  = h[0]·λx + h[5]·λy + h[4]·λz + boxlo[0]     // 2212</pre></div>
  </section>

  <section class="block">
    <h2>2 · Match: one atom through an orthogonal box</h2>
    <div class="prose"><p>Box 0…5, the atom starts at (4.9, 2.5), no forces, so it moves dt·v per step; pbc acts every <i>N</i> steps.</p></div>
    <div class="row">
      <label class="lc-in">vx<input type="number" class="lmp-in sm dp-vx" value="37" step="1"></label>
      <label class="lc-in">vy<input type="number" class="lmp-in sm dp-vy" value="-11" step="1"></label>
      <label class="lc-in">rebuild every<select class="lc-u dp-ev"><option>1</option><option>3</option><option>5</option></select></label>
      <button class="btn sm primary" data-k="o">▶ Run 10 steps</button><span class="lmp-status muted dp-sto"></span></div>
    <div class="tbl dp-outo"></div>
    <div class="prose dp-sumo"></div>
  </section>

  <section class="block">
    <h2>3 · Match: a tilted box</h2>
    <div class="prose"><p><code>region box prism 0 5 0 5 0 5 xy 0 0</code>: the top face (y = 5) is shifted by xy in x. An atom leaving through the top comes back at the bottom with x reduced by xy. Every step here is a rebuild step, so every step does x → λ → x, and so does the setup before step 0 (${gh('src/verlet.cpp', 117)}).</p></div>
    <div class="row">
      <label class="lc-in">xy<select class="lc-u dp-xy"><option>0.0</option><option>1.0</option><option selected>2.0</option><option>2.5</option></select></label>
      <label class="lc-in">vy<input type="number" class="lmp-in sm dp-tvy" value="41" step="1"></label>
      <button class="btn sm primary" data-k="t">▶ Run 6 steps</button><span class="lmp-status muted dp-stt"></span></div>
    <div class="tbl dp-outt"></div>
    <div class="prose dp-sumt"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: the image flag overflows</h2>
    <div class="prose"><p>vx = 800 moves the atom 4 per step through a box of 5: it crosses 4 boundaries every 5 steps. After 640 steps that is 512 crossings, one more than a 10-bit field can count. Predict ix every 50 steps, then run 700.</p></div>
    <div class="row"><button class="btn sm primary" data-k="v">▶ Run 700 steps</button><span class="lmp-status muted dp-stv"></span></div>
    <div class="tbl dp-outv"></div>
    <div class="prose dp-sumv"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> An atom has image = 537395713. What are ix, iy, iz?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>537395713 = 512·2²⁰ + 512·2¹⁰ + 513, so iz = 0, iy = 0, <b>ix = 1</b>.</div></div>
    <div class="ex"><p><b>2.</b> Box 0…5 with xy = 2. An atom at (3.0, 4.999) moves to y = 5.001. Where is it after pbc?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>λy goes from 1.0002 to 0.0002 and λx is unchanged, so x drops by xy: <b>(1.0, 0.001)</b>, iy = 1.</div></div>
    <div class="ex"><p><b>3.</b> Water at 300 K, D = 2.3 × 10⁻⁹ m²/s, box 30 Å, a 1 µs run. Is the image flag at risk?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>RMS displacement √(6Dt) ≈ 1200 Å per molecule ≈ 40 box lengths, and the farthest molecules a few times more: safe. A light gas over many µs, or a small box, can reach 512. Then build with <code>-DLAMMPS_BIGBIG</code> (21 bits, ±1 048 576), or dump unwrapped coordinates often enough to stitch them yourself.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const num = (s, d) => { const v = +$(s).value; return isFinite(v) && Math.abs(v) <= 190 ? v : d; };

    $('[data-k=o]').onclick = async e => {
      const b = e.currentTarget, st = $('.dp-sto'), vx = num('.dp-vx', 37), vy = num('.dp-vy', -11), ev = +$('.dp-ev').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptO(vx, vy, ev, 10), { read: ['o.dump'] }); b.disabled = false;
      const t = res.files?.['o.dump']; if (res.error || !t) return err(st, res); st.textContent = '';
      const H = handO(vx, vy, ev, 10), D = parseDump(t); let n = 0, ok = 0;
      $('.dp-outo').innerHTML = `<table>${head(['step', 'x', 'y', 'ix', 'iy', 'xu', 'image int', 'LAMMPS'])}${D.map(({ step, w }) => { const h = H[step], hv = [h.x, h.y, h.ix, h.iy, h.xu, h.yu], same = hv.map((q, i) => q === w[i + 1]); same.forEach(s => { n++; if (s) ok++; }); return row([step, `<span class="mono">${s17(h.x)}</span>`, `<span class="mono">${s17(h.y)}</span>`, h.ix, h.iy, `<span class="mono">${s17(h.xu)}</span>`, h.im, same.every(Boolean) ? `${OK} 6 of 6` : BAD]); }).join('')}</table>`;
      let xc = 4.9; for (let s = 0; s < 10; s++) xc += 0.005 * vx;
      $('.dp-sumo').innerHTML = `<p>${ok} of ${n} numbers identical. Atoms outside 0…5 between rebuilds are allowed (page C.2): try "every 5". The same 10 steps in a box too big to wrap would end at x = ${s17(xc)}; the unwrapped xu here is ${s17(H[10].xu)}${xc === H[10].xu ? ', the same.' : ': not the same double. Each wrap adds or subtracts 5 and rounds, so xu is the continuous path only to round-off.'}</p>`;
    };

    $('[data-k=t]').onclick = async e => {
      const b = e.currentTarget, st = $('.dp-stt'), xy = +$('.dp-xy').value, vy = num('.dp-tvy', 41); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptT(xy, vy), { read: ['t.dump'] }); b.disabled = false;
      const t = res.files?.['t.dump']; if (res.error || !t) return err(st, res); st.textContent = '';
      const { out: H, hi } = handT(xy, vy), D = parseDump(t); let n = 0, ok = 0, rt = 0;
      $('.dp-outt').innerHTML = `<table>${head(['step', 'x by hand', 'x before x→λ→x', 'y', 'xs', 'ys', 'ix iy', 'LAMMPS'])}${D.map(({ step, w }) => { const h = H[step], hv = [h.x, h.y, h.xs[0], h.xs[1], h.ix, h.iy], same = hv.map((q, i) => q === w[i + 1]); same.forEach(s => { n++; if (s) ok++; }); const moved = (step === 0 || h.iy === H[step - 1].iy) && h.plain !== h.x; if (moved) rt++; return row([step, `<span class="mono">${s17(h.x)}</span>`, `<span class="mono">${s17(h.plain)}</span>${moved ? ' ≠' : ''}`, `<span class="mono">${s17(h.y)}</span>`, `<span class="mono">${s17(h.xs[0])}</span>`, `<span class="mono">${s17(h.xs[1])}</span>`, `${h.ix} ${h.iy}`, same.every(Boolean) ? `${OK} 6 of 6` : BAD]); }).join('')}</table>`;
      $('.dp-sumt').innerHTML = `<p>${ok} of ${n} identical. h_inv[5] = −xy/(Lx·Ly) = ${s17(hi[5])}. On ${rt} step${rt === 1 ? '' : 's'} without any crossing, the trip into box fractions and back still changed x in its last digits (marked ≠)${xy === 0 ? '; with xy = 0 there is no tilt term, and only the division and multiplication by 5 can round' : ''}.</p>`;
    };

    $('[data-k=v]').onclick = async e => {
      const b = e.currentTarget, st = $('.dp-stv'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptO(800, 0, 1, 700, 'id x ix xu', 50), { read: ['o.dump'] }); b.disabled = false;
      const t = res.files?.['o.dump']; if (res.error || !t) return err(st, res); st.textContent = '';
      const H = handO(800, 0, 1, 700), D = parseDump(t); let ok = 0;
      let cross = 0; const cr = [0]; { let x = 4.9; for (let s = 1; s <= 700; s++) { x += 0.005 * 800; if (x >= 5) { x -= 5; cross++; } cr.push(cross); } }
      $('.dp-outv').innerHTML = `<table>${head(['step', 'crossings', 'ix by hand', 'xu by hand', 'LAMMPS ix', 'LAMMPS xu', ''])}${D.map(({ step, w }) => { const h = H[step], s = h.ix === w[2] && h.xu === w[3]; if (s) ok++; return row([step, cr[step], h.ix, `<span class="mono">${s17(h.xu)}</span>`, w[2], `<span class="mono">${s17(w[3])}</span>`, s ? OK : BAD]); }).join('')}</table>`;
      $('.dp-sumv').innerHTML = `<p>${ok} of ${D.length} rows identical. At crossing 512 the 10-bit field rolls over from 511 to −512, with no warning: xu drops by 1024 × 5 = 5120. MSD computed from xu would show a sudden jump. compute msd unwraps with the same image flags, so it jumps too.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
