// Lab 7.4 widget, two modes.
// ML fit: fit a Morse curve from a few training points with a rigid LJ-form model, a parabola, or K Gaussian basis functions;
//         tap the plot to add a training point (active learning).
// Metadynamics: Langevin walker in U(s) = 5 (s^2 - 1)^2 (kT = 1) with standard or well-tempered hills.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const TRUE = r => (1 - Math.exp(-2 * (r - 1.12))) ** 2 - 1;
const R0 = 0.88, R1 = 3.0;
function solve(A, b) {   // Gaussian elimination with partial pivoting
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) { let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r; [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; } }
  const x = new Array(n).fill(0); for (let r = n - 1; r >= 0; r--) { let s = M[r][n]; for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k]; x[r] = s / M[r][r]; } return x;
}

export function mountLab7(root) {
  const st = { mode: 'ml', model: 'lj', K: 8, pts: [1.0, 1.2, 1.5], w: 0.5, wt: 0, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view bx-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="ml" class="on">ML fit</button><button data-m="md">Metadynamics</button></div>
    <div data-p="ml">
      <div class="seg"><button data-f="lj" class="on">A/r¹² − B/r⁶</button><button data-f="par">Parabola</button><button data-f="rbf">K Gaussians</button></div>
      <label data-kk hidden><span>Basis functions K</span><output data-o="K"></output><input type="range" data-k="K" min="2" max="20" step="1" value="8"></label>
      <div class="row"><button class="btn sm ghost" data-k="three">3 points</button><button class="btn sm ghost" data-k="twenty">20 points</button><span class="note" style="font-size:12px;color:var(--muted)">tap the plot to add a point</span></div>
    </div>
    <div data-p="md" hidden>
      <div class="seg"><button data-t="0" class="on">Standard</button><button data-t="4">Well-tempered ΔT = 4</button></div>
      <label><span>Hill height w</span><output data-o="w"></output><input type="range" data-k="w" min="0.05" max="1" step="0.05" value="0.5"></label>
      <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    </div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.bx-view'), $ = s => root.querySelector(s);
  let rng = 9; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  // ---------- ML fit ----------
  let model = null;
  function fitModel() {
    const P = st.pts, E = P.map(TRUE); let phi;
    if (st.model === 'lj') phi = r => [r ** -12, -(r ** -6)];
    else if (st.model === 'par') phi = r => [1, r, r * r];
    else { const K = st.K, mu = Array.from({ length: K }, (_, k) => 0.9 + 1.3 * k / Math.max(1, K - 1)), wd = 1.3 / K; phi = r => mu.map(m => Math.exp(-((r - m) ** 2) / (2 * wd * wd))); }
    const n = phi(1).length, A = Array.from({ length: n }, () => new Array(n).fill(0)), b = new Array(n).fill(0);
    P.forEach((r, i) => { const f = phi(r); for (let a = 0; a < n; a++) { b[a] += f[a] * E[i]; for (let c = 0; c < n; c++) A[a][c] += f[a] * f[c]; } });
    for (let a = 0; a < n; a++) A[a][a] += 1e-8;
    const c = solve(A, b); model = r => phi(r).reduce((s, f, k) => s + f * c[k], 0);
  }
  // ---------- metadynamics ----------
  const G0 = -2.5, NG = 1001, DG = 5 / (NG - 1), SG = 0.2, DT = 0.001, GAM = 5;
  let s, v, Vg, dVg, nh, cross, side, t, trace;
  function resetMD() { s = -1; v = 0; Vg = new Float64Array(NG); dVg = new Float64Array(NG); nh = 0; cross = 0; side = -1; t = 0; trace = []; }
  const interp = (arr, x) => { const u = (x - G0) / DG, i = Math.min(Math.max(Math.floor(u), 0), NG - 2); return arr[i] + (arr[i + 1] - arr[i]) * (u - i); };
  const dU = x => 20 * x * (x * x - 1), U = x => 5 * (x * x - 1) ** 2;
  function mdStep() {
    const c = Math.exp(-GAM * DT), sn = Math.sqrt(1 - c * c);
    v += -0.5 * DT * (dU(s) + interp(dVg, s)); s += DT * v; v += -0.5 * DT * (dU(s) + interp(dVg, s)); v = c * v + sn * gauss(); t++;
    if (s * side < 0 && Math.abs(s) > 0.8) { cross++; side = -side; }
    if (t % 200 === 0) { const ww = st.wt ? st.w * Math.exp(-interp(Vg, s) / st.wt) : st.w; nh++;
      for (let k = 0; k < NG; k++) { const x = G0 + k * DG, g = Math.exp(-((x - s) ** 2) / (2 * SG * SG)); Vg[k] += ww * g; dVg[k] += -ww * (x - s) / (SG * SG) * g; } }
    if (t % 10000 === 0) { trace.push(barrier()); if (trace.length > 400) trace.shift(); }
  }
  function barrier() { const f = st.wt ? 1 + 1 / st.wt : 1; let mn = Infinity; for (let k = 0; k < NG; k++) { const x = G0 + k * DG; if (x > -1.5 && x < -0.5) mn = Math.min(mn, -f * Vg[k]); } return -f * interp(Vg, 0) - mn; }

  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); draw(); });
  root.querySelectorAll('[data-f]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-f]').forEach(z => z.classList.toggle('on', z === b)); st.model = b.dataset.f; $('[data-kk]').hidden = st.model !== 'rbf'; fitModel(); draw(); });
  root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-t]').forEach(z => z.classList.toggle('on', z === b)); st.wt = +b.dataset.t; resetMD(); draw(); });
  $('[data-k=K]').oninput = e => { st.K = +e.target.value; fitModel(); draw(); };
  $('[data-k=w]').oninput = e => { st.w = +e.target.value; resetMD(); draw(); };
  $('[data-k=three]').onclick = () => { st.pts = [1.0, 1.2, 1.5]; fitModel(); draw(); };
  $('[data-k=twenty]').onclick = () => { st.pts = Array.from({ length: 20 }, (_, k) => 0.95 + 1.05 * k / 19); fitModel(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { resetMD(); draw(); };
  cv.addEventListener('click', e => { if (st.mode !== 'ml') return; const b = cv.getBoundingClientRect(), r = R0 + (e.clientX - b.left - 50) / (W - 62) * (R1 - R0);
    if (r > R0 && r < R1 && !st.pts.some(p => Math.abs(p - r) < 0.01)) { st.pts.push(r); st.pts.sort((a, c) => a - c); fitModel(); draw(); } });

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (tx, x, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(tx, x, y); };
  const curve = (fn, X, Y, a, b, c, w = 2, n = 300) => { cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); for (let k = 0; k <= n; k++) { const x = a + (b - a) * k / n; k ? cx.lineTo(X(x), Y(fn(x))) : cx.moveTo(X(x), Y(fn(x))); } cx.stroke(); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const Rr = $('[data-ro]');
    if (st.mode === 'ml') {
      const x0 = 50, x1 = W - 12, y0 = 24, y1 = H - 30, lo = -1.4, hi = 1.0, X = r => x0 + (r - R0) / (R1 - R0) * (x1 - x0), Y = e => y1 - (Math.min(hi, Math.max(lo, e)) - lo) / (hi - lo) * (y1 - y0);
      axes(x0, y0, x1, y1); label('energy vs r: truth (amber), model (cyan), training points', x0 + 6, y0 - 8, COL.ink, 'left', W < 500 ? '600 10px' : '600 12px');
      const pmin = Math.min(...st.pts), pmax = Math.max(...st.pts); cx.fillStyle = 'rgba(56,225,255,.07)'; cx.fillRect(X(pmin), y0, X(pmax) - X(pmin), y1 - y0);
      cx.strokeStyle = 'rgba(143,163,199,.35)'; cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
      curve(TRUE, X, Y, R0, R1, COL.amber, 2.5); curve(model, X, Y, R0, R1, COL.cyan, 2);
      st.pts.forEach(r => { cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(X(r), Y(TRUE(r)), 5, 0, 7); cx.fill(); });
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 1.5, 2, 2.5, 3].forEach(r => cx.fillText(r, X(r), y1 + 14));
      cx.textAlign = 'right'; [-1, 0, 0.5].forEach(e => cx.fillText(e, x0 - 4, Y(e) + 4));
      let tr = 0, te = 0, out = 0; st.pts.forEach(r => tr = Math.max(tr, Math.abs(model(r) - TRUE(r))));
      for (let k = 0; k <= 200; k++) { const r = 0.95 + 1.05 * k / 200; te = Math.max(te, Math.abs(model(r) - TRUE(r))); }
      for (let k = 0; k <= 200; k++) { const r = R0 + (R1 - R0) * k / 200; if (r < pmin || r > pmax) out = Math.max(out, Math.abs(model(r) - TRUE(r))); }
      Rr.innerHTML = `<div><span>Training points</span><b>${st.pts.length}</b></div><div><span>Max training error</span><b>${tr.toFixed(4)}</b></div>
        <div><span>Max test error, r 0.95–2.0</span><b>${te.toFixed(4)}</b></div><div><span>Max error outside the data</span><b>${out ? out.toFixed(3) : '—'}</b></div>`;
      $('[data-o=K]').textContent = st.K;
    } else {
      const half = H * 0.55, x0 = 46, x1 = W - 12, y0 = 24, y1 = half, X = x => x0 + (x + 2) / 4 * (x1 - x0), Emax = 9, Y = e => y1 - Math.min(Emax, Math.max(-0.5, e)) / Emax * (y1 - y0);
      axes(x0, y0, x1, y1); label(W < 500 ? 'U (white), bias (grey), U + bias (cyan)' : 'U(s) (white), bias V (grey), U + V (cyan), free energy from the bias (amber)', x0 + 6, y0 - 8);
      const f = st.wt ? 1 + 1 / st.wt : 1, Vat = x => interp(Vg, x), Fm = (() => { let m = Infinity; for (let x = -1.5; x < -0.5; x += 0.01) m = Math.min(m, -f * Vat(x)); return m; })();
      curve(U, X, Y, -2, 2, 'rgba(242,245,250,.8)', 1.5); curve(Vat, X, Y, -2, 2, 'rgba(143,163,199,.7)', 1.5); curve(x => U(x) + Vat(x), X, Y, -2, 2, COL.cyan, 2);
      if (nh) { cx.setLineDash([5, 4]); curve(x => -f * Vat(x) - Fm, X, Y, -1.6, 1.6, COL.amber, 2); cx.setLineDash([]); }
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X(s), Y(U(s) + Vat(s)), 6, 0, 7); cx.fill();
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [-1, 0, 1].forEach(x => cx.fillText(x, X(x), y1 + 14));
      // barrier estimate trace
      const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 24, BY = b => by1 - Math.min(10, Math.max(0, b)) / 10 * (by1 - by0), BX = k => bx0 + k / 399 * (bx1 - bx0);
      axes(bx0, by0, bx1, by1); label('barrier estimate vs time (dashed: exact 5)', bx0 + 6, by0 - 10);
      cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(bx0, BY(5)); cx.lineTo(bx1, BY(5)); cx.stroke(); cx.setLineDash([]);
      cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); trace.forEach((b, k) => k ? cx.lineTo(BX(k), BY(b)) : cx.moveTo(BX(k), BY(b))); cx.stroke();
      cx.fillStyle = COL.muted; cx.textAlign = 'right'; [0, 5, 10].forEach(b => cx.fillText(b, bx0 - 4, BY(b) + 4));
      Rr.innerHTML = `<div><span>Hills</span><b>${nh}</b></div><div><span>Barrier from the bias (exact 5)</span><b>${nh ? barrier().toFixed(2) : '…'}</b></div>
        <div><span>Crossings</span><b>${cross}</b></div><div><span>Steps</span><b>${t.toLocaleString()}</b></div>`;
      $('[data-o=w]').textContent = st.w.toFixed(2);
    }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running || st.mode !== 'md') return; for (let q = 0; q < 2000; q++) mdStep(); draw(); };
  fitModel(); resetMD(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
