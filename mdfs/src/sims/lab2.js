// Lab 2 widget: 9 atoms in a 2D periodic box (L = 3.6), exactly the system of Lab 2.7.
// Start from the lattice or from the random positions, prepare velocities step by step, look at periodic images,
// switch the minimum image off, choose the cutoff and skin, and watch the neighbour list being rebuilt as atoms move.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const L = 3.6, N = 9, DT = 0.005, T0 = 0.72;
const LAT = [0, 1, 2].flatMap(j => [0, 1, 2].map(i => [(i + 0.5) * 1.2, (j + 0.5) * 1.2]));
const RND = [[3.40, 1.84], [3.51, 0.29], [2.19, 1.36], [2.89, 0.63], [3.14, 1.96], [3.25, 1.72], [1.55, 2.84], [3.54, 1.33], [3.49, 3.34]];
const GAU = [[-0.33, -0.88], [-0.66, -0.67], [0.38, -0.11], [1.48, -1.83], [-0.00, -0.89], [0.78, -2.12], [-0.34, 0.21], [-1.48, 0.99], [0.18, 1.01]];
const U = r => 4 * (r ** -12 - r ** -6);

export function mountLab2(root) {
  const st = { start: 'lat', vel: 'scaled', images: true, mic: true, rc: 1.8, skin: 0.3, sel: 0, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Tap an atom to see its neighbours</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="lat" class="on">Lattice start</button><button data-s="rnd">Random start</button></div>
    <div class="seg"><button data-v="raw">Raw Gaussian</button><button data-v="nodrift">Drift removed</button><button data-v="scaled" class="on">Scaled to T</button></div>
    <div class="seg"><button data-m="1" class="on">Minimum image on</button><button data-m="0">off</button></div>
    <div class="seg"><button data-i="1" class="on">Show images</button><button data-i="0">Hide</button></div>
    <label><span>Cutoff r<sub>c</sub> (≤ L/2 = 1.8)</span><output data-o="rc"></output><input type="range" data-k="rc" min="1.0" max="1.8" step="0.05" value="1.8"></label>
    <label><span>Neighbour-list skin</span><output data-o="sk"></output><input type="range" data-k="sk" min="0" max="0.6" step="0.05" value="0.3"></label>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span>Potential energy (per atom)</span><b data-o="E"></b></div>
      <div><span><i class="nc">T</i> with N<sub>f</sub> = 16 (2N = 18)</span><b data-o="T"></b></div>
      <div><span>Drift velocity</span><b data-o="dr"></b></div>
      <div><span>List rebuilds / steps · pairs</span><b data-o="nl"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let x, v, f, list, ref, steps, rebuilds;
  const mi = d => st.mic ? d - L * Math.round(d / L) : d;
  function velocities() {
    let vv = GAU.map(p => [...p]);
    if (st.vel !== 'raw') { const m = [0, 1].map(k => vv.reduce((s, p) => s + p[k], 0) / N); vv = vv.map(p => [p[0] - m[0], p[1] - m[1]]); }
    if (st.vel === 'scaled') { const s2 = vv.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2, 0), lam = Math.sqrt(T0 * 16 / s2); vv = vv.map(p => [p[0] * lam, p[1] * lam]); }
    return vv;
  }
  function build() { list = []; const rl = st.rc + st.skin; for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { const r = Math.hypot(mi(x[j][0] - x[i][0]), mi(x[j][1] - x[i][1])); if (r < rl) list.push([i, j]); } ref = x.map(p => [...p]); rebuilds++; }
  function forces() {
    f = x.map(() => [0, 0]); let e = 0;
    for (const [i, j] of list) { const dx = mi(x[j][0] - x[i][0]), dy = mi(x[j][1] - x[i][1]), r = Math.hypot(dx, dy); if (r >= st.rc) continue;
      e += U(r); const ff = 24 * (2 * r ** -14 - r ** -8); f[i][0] -= ff * dx; f[i][1] -= ff * dy; f[j][0] += ff * dx; f[j][1] += ff * dy; }
    return e;
  }
  let Epot = 0;
  function reset() { x = (st.start === 'lat' ? LAT : RND).map(p => [...p]); v = velocities(); steps = 0; rebuilds = 0; build(); Epot = forces(); }
  function step() {
    for (let i = 0; i < N; i++) for (let k = 0; k < 2; k++) { v[i][k] += 0.5 * DT * f[i][k]; x[i][k] += DT * v[i][k]; if (st.mic) x[i][k] -= L * Math.floor(x[i][k] / L); }
    const moved = Math.max(...x.map((p, i) => Math.hypot(mi(p[0] - ref[i][0]), mi(p[1] - ref[i][1]))));
    if (moved > st.skin / 2) build();
    Epot = forces(); for (let i = 0; i < N; i++) for (let k = 0; k < 2; k++) v[i][k] += 0.5 * DT * f[i][k];
    steps++;
  }
  const seg = (a, fn) => root.querySelectorAll(`[data-${a}]`).forEach(b => b.onclick = () => { root.querySelectorAll(`[data-${a}]`).forEach(z => z.classList.toggle('on', z === b)); fn(b.dataset[a]); draw(); });
  seg('s', val => { st.start = val; reset(); }); seg('v', val => { st.vel = val; reset(); });
  seg('m', val => { st.mic = val === '1'; reset(); }); seg('i', val => { st.images = val === '1'; });
  $('[data-k=rc]').oninput = e => { st.rc = +e.target.value; reset(); draw(); };
  $('[data-k=sk]').oninput = e => { st.skin = +e.target.value; reset(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600, S = 100, ox = 0, oy = 0;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  cv.addEventListener('pointerdown', e => { const r = cv.getBoundingClientRect(), px = (e.clientX - r.left - ox) / S, py = L - (e.clientY - r.top - oy) / S;
    let best = 0, bd = 1e9; x.forEach((p, i) => { const d = Math.hypot(mi(p[0] - px), mi(p[1] - py)); if (d < bd) { bd = d; best = i; } }); st.sel = best; draw(); });
  const toS = (px, py) => [ox + px * S, oy + (L - py) * S];
  function draw() {
    if (!W || !x) return;
    const span = st.images ? 3 * L : 1.3 * L; S = Math.min(W - 20, H - 20) / span; ox = (W - L * S) / 2; oy = (H - L * S) / 2;
    cx.clearRect(0, 0, W, H);
    if (st.images) for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) { if (!a && !b) continue;
      cx.strokeStyle = 'rgba(143,163,199,.18)'; cx.strokeRect(ox + a * L * S, oy - b * L * S, L * S, L * S);
      x.forEach(p => { const [X, Y] = toS(p[0] + a * L, p[1] + b * L); cx.fillStyle = 'rgba(56,225,255,.25)'; cx.beginPath(); cx.arc(X, Y, 0.3 * S, 0, 7); cx.fill(); }); }
    cx.strokeStyle = 'rgba(56,225,255,.7)'; cx.lineWidth = 2; cx.strokeRect(ox, oy, L * S, L * S);
    // selected atom: cutoff and list circles, neighbour lines
    const p0 = x[st.sel], [SX, SY] = toS(p0[0], p0[1]);
    cx.setLineDash([5, 5]); cx.strokeStyle = 'rgba(255,194,75,.8)'; cx.beginPath(); cx.arc(SX, SY, st.rc * S, 0, 7); cx.stroke();
    cx.strokeStyle = 'rgba(255,194,75,.35)'; cx.beginPath(); cx.arc(SX, SY, (st.rc + st.skin) * S, 0, 7); cx.stroke(); cx.setLineDash([]);
    for (let j = 0; j < N; j++) { if (j === st.sel) continue; const dx = mi(x[j][0] - p0[0]), dy = mi(x[j][1] - p0[1]), r = Math.hypot(dx, dy); if (r >= st.rc + st.skin) continue;
      const [X, Y] = toS(p0[0] + dx, p0[1] + dy); cx.strokeStyle = r < st.rc ? COL.amber : 'rgba(255,194,75,.35)'; cx.lineWidth = r < st.rc ? 2 : 1; cx.beginPath(); cx.moveTo(SX, SY); cx.lineTo(X, Y); cx.stroke(); }
    x.forEach((p, i) => { const [X, Y] = toS(p[0], p[1]); cx.fillStyle = i === st.sel ? COL.amber : COL.cyan; cx.beginPath(); cx.arc(X, Y, 0.3 * S, 0, 7); cx.fill();
      cx.fillStyle = '#060C1A'; cx.font = `700 ${Math.max(10, 0.3 * S)}px "Space Grotesk", sans-serif`; cx.textAlign = 'center'; cx.fillText(i + 1, X, Y + 0.11 * S);
      const vs = 0.25 * S; cx.strokeStyle = COL.coral; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(X, Y); cx.lineTo(X + v[i][0] * vs, Y - v[i][1] * vs); cx.stroke(); });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('red arrows: velocities · solid lines: pairs within r_c · faint: in the list (r_c + skin)', 8, H - 8);
    // readouts
    const s2 = v.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2, 0), m = [0, 1].map(k => v.reduce((s, p) => s + p[k], 0) / N);
    o.rc.textContent = st.rc.toFixed(2); o.sk.textContent = st.skin.toFixed(2);
    o.E.textContent = Math.abs(Epot) > 1e5 ? `${Epot.toExponential(2)} (overlap!)` : `${Epot.toFixed(4)} (${(Epot / N).toFixed(4)})`;
    o.T.textContent = `${(s2 / 16).toFixed(4)} (${(s2 / 18).toFixed(4)})`;
    o.dr.textContent = `(${m[0].toFixed(3)}, ${m[1].toFixed(3)})`;
    o.nl.textContent = `${rebuilds} / ${steps} · ${list.length}`;
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let q = 0; q < 3; q++) step(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
