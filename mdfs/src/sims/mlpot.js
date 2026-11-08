// MLIP lab on a 2D energy surface (Mueller-Brown, scaled): the 'reference' is the analytic surface (think DFT),
// the model is Gaussian-process regression on the points computed so far. A Langevin walker runs on the MODEL.
// Where the model has no data it guesses badly; its own uncertainty tells you where to compute next (active learning).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const MB = { A: [-200, -100, -170, 15], a: [-1, -1, -6.5, 0.7], b: [0, 0, 11, 0.6], c: [-10, -10, -6.5, 0.7], x0: [1, 0, -0.5, -1], y0: [0, 0.5, 1.5, 1] };
const X0 = -1.6, X1 = 1.2, Y0 = -0.5, Y1 = 2.1, GN = 90;
function ref(x, y) { let e = 0; for (let k = 0; k < 4; k++) { const dx = x - MB.x0[k], dy = y - MB.y0[k]; e += MB.A[k] * Math.exp(MB.a[k] * dx * dx + MB.b[k] * dx * dy + MB.c[k] * dy * dy); } return e / 100; }
function refGrad(x, y) { let gx = 0, gy = 0; for (let k = 0; k < 4; k++) { const dx = x - MB.x0[k], dy = y - MB.y0[k], ex = MB.A[k] * Math.exp(MB.a[k] * dx * dx + MB.b[k] * dx * dy + MB.c[k] * dy * dy) / 100; gx += ex * (2 * MB.a[k] * dx + MB.b[k] * dy); gy += ex * (MB.b[k] * dx + 2 * MB.c[k] * dy); } return [gx, gy]; }

export function mountMlpot(root) {
  const st = { kT: 0.25, ell: 0.25, view: 'model', running: true, walkOn: 'model' };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Tap the map to compute a reference point there</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-v="model" class="on">Model</button><button data-v="ref">Reference</button><button data-v="err">Error</button><button data-v="sd">Uncertainty</button></div>
    <div class="seg"><button data-w="model" class="on">Walker on model</button><button data-w="ref">Walker on reference</button></div>
    <label><span>Temperature <i class="nc">kT</i></span><output data-o="kT"></output><input type="range" data-k="kT" min="0.05" max="0.6" step="0.01" value="0.25"></label>
    <label><span>Kernel length scale</span><output data-o="ell"></output><input type="range" data-k="ell" min="0.1" max="0.6" step="0.01" value="0.25"></label>
    <div class="row"><button class="btn sm primary" data-k="al">Active-learning step</button><button class="btn sm ghost" data-k="al5">×5</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span>Reference calculations</span><b data-o="n"></b></div>
      <div><span>RMS error (low-energy region)</span><b data-o="rms"></b></div>
      <div><span>Walker: model / true energy</span><b data-o="wk"></b></div>
      <div><span>Largest uncertainty on path</span><b data-o="sd"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 5; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let P = [], alpha = [], Lch = null, mean0 = 0, walker = [-0.55, 1.45], path = [], img = null, dirty = true;
  const SF2 = 1.0, NOISE = 1e-4;
  const k = (a, b) => SF2 * Math.exp(-((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) / (2 * st.ell * st.ell));
  function train() {   // GP: alpha = K^-1 (y - m), Cholesky kept for the variance
    const n = P.length; mean0 = -0.5;   // prior mean: the model's guess far from any data
    const K = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => k(P[i], P[j]) + (i === j ? NOISE : 0)));
    Lch = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { let s = K[i][j]; for (let q = 0; q < j; q++) s -= Lch[i][q] * Lch[j][q]; Lch[i][j] = i === j ? Math.sqrt(Math.max(s, 1e-12)) : s / Lch[j][j]; }
    const y = P.map(p => p[2] - mean0), z = new Float64Array(n); alpha = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = y[i]; for (let q = 0; q < i; q++) s -= Lch[i][q] * z[q]; z[i] = s / Lch[i][i]; }
    for (let i = n - 1; i >= 0; i--) { let s = z[i]; for (let q = i + 1; q < n; q++) s -= Lch[q][i] * alpha[q]; alpha[i] = s / Lch[i][i]; }
    dirty = true;
  }
  function predict(x, y, wantVar) {
    const n = P.length, ks = new Float64Array(n); let mu = mean0, gx = 0, gy = 0;
    for (let i = 0; i < n; i++) { ks[i] = k([x, y], P[i]); mu += ks[i] * alpha[i]; const f = -ks[i] * alpha[i] / (st.ell * st.ell); gx += f * (x - P[i][0]); gy += f * (y - P[i][1]); }
    let sd = 0; if (wantVar) { const v = new Float64Array(n); let s2 = 0; for (let i = 0; i < n; i++) { let s = ks[i]; for (let q = 0; q < i; q++) s -= Lch[i][q] * v[q]; v[i] = s / Lch[i][i]; s2 += v[i] * v[i]; } sd = Math.sqrt(Math.max(0, SF2 - s2)); }
    return { mu, gx, gy, sd };
  }
  function addPoint(x, y) { P.push([x, y, ref(x, y)]); train(); }
  function reset() { P = []; [[-0.56, 1.44], [-0.45, 1.5], [-0.65, 1.35], [-0.55, 1.6], [-0.7, 1.5], [-0.4, 1.35]].forEach(([x, y]) => P.push([x, y, ref(x, y)])); train(); walker = [-0.55, 1.45]; path = []; }
  function activeStep() {   // compute the reference where the model is least sure, among recently visited points
    let best = null, bs = -1; for (const p of path) { const s = predict(p[0], p[1], true).sd; if (s > bs) { bs = s; best = p; } }
    if (best) addPoint(best[0], best[1]);
  }
  function walk() {   // overdamped Langevin (Brownian) dynamics on the model or the reference: dx = -grad E * s^2/(2kT) + s * noise
    const sw = 0.035;
    for (let q = 0; q < 4; q++) {
      const g = st.walkOn === 'model' ? predict(walker[0], walker[1], false) : (([gx, gy]) => ({ gx, gy }))(refGrad(walker[0], walker[1]));
      const cap = v => Math.max(-40, Math.min(40, v));
      walker[0] += -cap(g.gx) * sw * sw / (2 * st.kT) + sw * gauss();
      walker[1] += -cap(g.gy) * sw * sw / (2 * st.kT) + sw * gauss();
      if (walker[0] < X0) walker[0] = 2 * X0 - walker[0]; if (walker[0] > X1) walker[0] = 2 * X1 - walker[0];
      if (walker[1] < Y0) walker[1] = 2 * Y0 - walker[1]; if (walker[1] > Y1) walker[1] = 2 * Y1 - walker[1];
    }
    path.push([walker[0], walker[1]]); if (path.length > 300) path.shift();
  }
  // controls
  root.querySelectorAll('[data-v]').forEach(b => b.onclick = () => { st.view = b.dataset.v; root.querySelectorAll('[data-v]').forEach(z => z.classList.toggle('on', z === b)); dirty = true; });
  root.querySelectorAll('[data-w]').forEach(b => b.onclick = () => { st.walkOn = b.dataset.w; root.querySelectorAll('[data-w]').forEach(z => z.classList.toggle('on', z === b)); });
  $('[data-k=kT]').oninput = e => { st.kT = +e.target.value; };
  $('[data-k=ell]').oninput = e => { st.ell = +e.target.value; train(); };
  $('[data-k=al]').onclick = () => activeStep(); $('[data-k=al5]').onclick = () => { for (let q = 0; q < 5; q++) { for (let s = 0; s < 40; s++) walk(); activeStep(); } };
  $('[data-k=reset]').onclick = () => reset();

  let W = 600, H = 600, map = { x: 0, y: 0, w: 1, h: 1 };
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); dirty = true; };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const toS = (x, y) => [map.x + (x - X0) / (X1 - X0) * map.w, map.y + (Y1 - y) / (Y1 - Y0) * map.h];
  cv.addEventListener('pointerdown', e => { const r = cv.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
    const x = X0 + (px - map.x) / map.w * (X1 - X0), y = Y1 - (py - map.y) / map.h * (Y1 - Y0); if (x > X0 && x < X1 && y > Y0 && y < Y1) addPoint(x, y); });

  function colour(v, lo, hi, diverge) {   // simple perceptual-ish ramps
    const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo)));
    if (diverge) return [Math.round(10 + 245 * t), Math.round(20 + 90 * t), Math.round(40 + 70 * t)];
    const a = [[8, 16, 40], [20, 80, 140], [56, 225, 255], [255, 194, 75], [255, 107, 107]], s = t * 4, i = Math.min(3, Math.floor(s)), f = s - i;
    return a[i].map((c, j) => Math.round(c + (a[i + 1][j] - c) * f));
  }
  let rmsCache = 0;
  function rebuild() {
    const side = Math.min(W - 20, H - 20); map = { x: (W - side * (X1 - X0) / Math.max(X1 - X0, Y1 - Y0)) / 2, y: 10, w: side * (X1 - X0) / Math.max(X1 - X0, Y1 - Y0), h: side * (Y1 - Y0) / Math.max(X1 - X0, Y1 - Y0) };
    const off = document.createElement('canvas'); off.width = GN; off.height = GN; const oc = off.getContext('2d'), id = oc.createImageData(GN, GN); let se = 0, ne = 0;
    for (let j = 0; j < GN; j++) for (let i = 0; i < GN; i++) {
      const x = X0 + (i + 0.5) / GN * (X1 - X0), y = Y1 - (j + 0.5) / GN * (Y1 - Y0), r = ref(x, y), p = predict(x, y, st.view === 'sd');
      if (r < 0.5) { se += (p.mu - r) ** 2; ne++; }
      const c = st.view === 'model' ? colour(p.mu, -1.5, 1.0) : st.view === 'ref' ? colour(r, -1.5, 1.0) : st.view === 'err' ? (r > 1 ? [28, 34, 52] : colour(Math.abs(p.mu - r), 0, 0.8, true)) : colour(p.sd, 0, 1, true);
      const q = 4 * (j * GN + i); id.data[q] = c[0]; id.data[q + 1] = c[1]; id.data[q + 2] = c[2]; id.data[q + 3] = 255;
    }
    oc.putImageData(id, 0, 0); img = off; rmsCache = Math.sqrt(se / ne); dirty = false;
  }
  function draw() {
    if (!W) return;
    if (dirty || !img) rebuild();
    cx.clearRect(0, 0, W, H); cx.imageSmoothingEnabled = true; cx.drawImage(img, map.x, map.y, map.w, map.h);
    cx.strokeStyle = 'rgba(242,245,250,.25)'; cx.strokeRect(map.x, map.y, map.w, map.h);
    // path, training points, walker
    cx.strokeStyle = 'rgba(242,245,250,.55)'; cx.lineWidth = 1.2; cx.beginPath(); path.forEach((p, q) => { const [X, Y] = toS(p[0], p[1]); q ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke();
    P.forEach(p => { const [X, Y] = toS(p[0], p[1]); cx.fillStyle = '#060C1A'; cx.beginPath(); cx.arc(X, Y, 4.5, 0, 7); cx.fill(); cx.strokeStyle = COL.ink; cx.lineWidth = 1.5; cx.stroke(); });
    const [WX, WY] = toS(walker[0], walker[1]); cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(WX, WY, 7, 0, 7); cx.fill();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.fillText({ model: 'model energy (GP)', ref: 'reference energy ("DFT")', err: '|model − reference| (grey: reference > 1, irrelevant)', sd: 'model uncertainty (GP std)' }[st.view], map.x + 8, map.y + map.h - 10);
    const p = predict(walker[0], walker[1], true), r = ref(walker[0], walker[1]);
    o.kT.textContent = st.kT.toFixed(2); o.ell.textContent = st.ell.toFixed(2); o.n.textContent = P.length;
    o.rms.textContent = rmsCache.toFixed(3); o.wk.textContent = `${p.mu.toFixed(2)} / ${r.toFixed(2)}${Math.abs(p.mu - r) > 0.3 ? '  ⚠' : ''}`;
    let ms = 0; for (let q = 0; q < path.length; q += 10) ms = Math.max(ms, predict(path[q][0], path[q][1], true).sd); o.sd.textContent = ms.toFixed(2);
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) walk(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
