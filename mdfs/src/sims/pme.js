// PME lab: smooth particle-mesh Ewald (reciprocal part) for 216 ions in a periodic cube, in the browser.
// Choose the grid size and B-spline order; compare with the exact Ewald Fourier sum (accuracy and time).
// Pictures: one slice of the spread charge grid, and the same slice of the smooth reciprocal-space potential.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 216, RHO = 0.5, L = Math.cbrt(N / RHO), ALPHA = 0.8, KREF = 10;

function fft1(re, im, inv) {                       // in-place radix-2 complex FFT
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; } }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (inv ? 2 : -2) * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) { let cr = 1, ci = 0; for (let k = 0; k < len / 2; k++) {
      const a = i + k, b = a + len / 2, tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
      re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti; const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t; } }
  }
}
function fft3(re, im, K, inv) {
  const r = new Float64Array(K), m = new Float64Array(K), id = (a, b, c) => (a * K + b) * K + c;
  for (let ax = 0; ax < 3; ax++) for (let p = 0; p < K; p++) for (let q = 0; q < K; q++) {
    for (let t = 0; t < K; t++) { const k = ax === 0 ? id(t, p, q) : ax === 1 ? id(p, t, q) : id(p, q, t); r[t] = re[k]; m[t] = im[k]; }
    fft1(r, m, inv);
    for (let t = 0; t < K; t++) { const k = ax === 0 ? id(t, p, q) : ax === 1 ? id(p, t, q) : id(p, q, t); re[k] = r[t]; im[k] = m[t]; }
  }
}
function bspline(w, n) {                           // M_n(w + j), j = 0..n-1
  let M = [w, 1 - w];
  for (let p = 3; p <= n; p++) { const Mn = []; for (let j = 0; j < p; j++) { const x = w + j, left = j < p - 1 ? M[j] : 0, right = j >= 1 ? M[j - 1] : 0; Mn.push((x * left + (p - x) * right) / (p - 1)); } M = Mn; }
  return M;
}

export function mountPme(root) {
  const st = { K: 16, n: 4 };
  let X = [], Q = [], Eref = 0, tRef = 0;
  function makeSystem(seed) { let s = seed; const R = () => (s = (s * 16807) % 2147483647) / 2147483647; X = []; Q = []; for (let i = 0; i < N; i++) { X.push([R() * L, R() * L, R() * L]); Q.push(i % 2 ? -1 : 1); } }
  function direct() {
    const t0 = performance.now(); let E = 0;
    for (let a = -KREF; a <= KREF; a++) for (let b = -KREF; b <= KREF; b++) for (let c = 0; c <= KREF; c++) {
      if (c === 0 && (b < 0 || (b === 0 && a <= 0))) continue;
      const kx = 2 * Math.PI * a / L, ky = 2 * Math.PI * b / L, kz = 2 * Math.PI * c / L, k2 = kx * kx + ky * ky + kz * kz;
      let Sr = 0, Si = 0; for (let i = 0; i < N; i++) { const ph = kx * X[i][0] + ky * X[i][1] + kz * X[i][2]; Sr += Q[i] * Math.cos(ph); Si += Q[i] * Math.sin(ph); }
      E += 2 * Math.exp(-k2 / (4 * ALPHA * ALPHA)) / k2 * (Sr * Sr + Si * Si);
    }
    tRef = performance.now() - t0; return 2 * Math.PI / L ** 3 * E;
  }
  function pme(K, n) {
    const t0 = performance.now(), G = K * K * K, re = new Float64Array(G), im = new Float64Array(G), id = (a, b, c) => (((a % K) + K) % K * K + (((b % K) + K) % K)) * K + (((c % K) + K) % K);
    for (let i = 0; i < N; i++) {
      const u = X[i].map(x => x / L * K), f = u.map(Math.floor), W = u.map((x, d) => bspline(x - f[d], n));
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) { const wab = Q[i] * W[0][a] * W[1][b]; for (let c = 0; c < n; c++) re[id(f[0] - a, f[1] - b, f[2] - c)] += wab * W[2][c]; }
    }
    const charge = Float64Array.from(re);
    fft3(re, im, K, false);
    const Mint = bspline(0, n).slice(1), b2 = new Float64Array(K);
    for (let m = 0; m < K; m++) { const mm = m < K / 2 ? m : m - K; let dr = 0, di = 0; for (let k = 0; k < n - 1; k++) { dr += Mint[k] * Math.cos(2 * Math.PI * mm * k / K); di += Mint[k] * Math.sin(2 * Math.PI * mm * k / K); } b2[m] = 1 / (dr * dr + di * di); }
    let E = 0;
    for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) for (let c = 0; c < K; c++) {
      const k = (a * K + b) * K + c; if (!k) { re[k] = im[k] = 0; continue; }
      const ma = a < K / 2 ? a : a - K, mb = b < K / 2 ? b : b - K, mc = c < K / 2 ? c : c - K, m2 = (ma * ma + mb * mb + mc * mc) / (L * L);
      const g = Math.exp(-Math.PI * Math.PI * m2 / (ALPHA * ALPHA)) / m2 * b2[a] * b2[b] * b2[c];
      E += g * (re[k] * re[k] + im[k] * im[k]); re[k] *= g; im[k] *= g;
    }
    fft3(re, im, K, true);                          // convolution: the smooth potential on the grid (arbitrary scale)
    return { E: E / (2 * Math.PI * L ** 3), t: performance.now() - t0, charge, pot: re };
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg" data-g="K"><button data-grid="8">8³ grid</button><button data-grid="16" class="on">16³</button><button data-grid="32">32³</button></div>
    <div class="seg" data-g="n"><button data-n="2">order 2</button><button data-n="4" class="on">order 4</button><button data-n="6">order 6</button><button data-n="8">order 8</button></div>
    <div class="row"><button class="btn sm ghost" data-k="new">New configuration</button></div>
    <div class="readouts">
      <div><span>PME Fourier energy</span><b data-o="E"></b></div>
      <div><span>Exact Ewald value</span><b data-o="ref"></b></div>
      <div><span>Relative error</span><b data-o="err"></b></div>
      <div><span>Grid spacing</span><b data-o="h"></b></div>
      <div class="wide"><span>Time: PME / direct Ewald (this browser)</span><b data-o="t"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let res = null;
  const run = () => { res = pme(st.K, st.n); draw(); };
  root.querySelectorAll('[data-grid]').forEach(b => b.onclick = () => { st.K = +b.dataset.grid; root.querySelectorAll('[data-grid]').forEach(z => z.classList.toggle('on', z === b)); run(); });
  root.querySelectorAll('[data-n]').forEach(b => b.onclick = () => { st.n = +b.dataset.n; root.querySelectorAll('[data-n]').forEach(z => z.classList.toggle('on', z === b)); run(); });
  let seed = 11; const init = () => { makeSystem(seed); Eref = direct(); run(); };
  $('[data-k=new]').onclick = () => { seed += 97; init(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function heat(grid, K, x0, y0, s, title) {
    const z = 0; let mx = 0; for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) mx = Math.max(mx, Math.abs(grid[(a * K + b) * K + z]));
    const c = s / K;
    for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) { const v = grid[(a * K + b) * K + z] / (mx || 1); cx.fillStyle = v > 0 ? `rgba(255,107,107,${Math.abs(v)})` : `rgba(56,225,255,${Math.abs(v)})`; cx.fillRect(x0 + a * c, y0 + b * c, c + 0.5, c + 0.5); }
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.strokeRect(x0, y0, s, s);
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(title, x0 + s / 2, y0 - 8);
  }
  function draw() {
    if (!W || !res) return;
    cx.clearRect(0, 0, W, H);
    const wide = W > 520, s = wide ? Math.min((W - 60) / 2, H - 70) : Math.min(W - 40, (H - 90) / 2);
    const yy = wide ? Math.max(40, (H - s) / 2) : Math.max(30, (H - 2 * s - 30) / 2), p1 = wide ? [20, yy] : [(W - s) / 2, yy], p2 = wide ? [40 + s, yy] : [(W - s) / 2, yy + s + 30];
    heat(res.charge, st.K, p1[0], p1[1], s, 'charges spread onto the grid (one slice)');
    // ions near this slice
    const sh = st.n / 2 - 0.5;                       // B-spline weights are centred n/2 - 1/2 cells before the atom
    X.forEach((p, i) => { const uz = p[2] / L * st.K; if (uz > st.n) return; const ux = ((p[0] / L * st.K - sh) % st.K + st.K) % st.K, uy = ((p[1] / L * st.K - sh) % st.K + st.K) % st.K; cx.fillStyle = Q[i] > 0 ? COL.coral : COL.cyan; cx.strokeStyle = COL.ink; cx.lineWidth = 1; cx.beginPath(); cx.arc(p1[0] + (ux + 0.5) / st.K * s, p1[1] + (uy + 0.5) / st.K * s, 3.5, 0, 7); cx.fill(); cx.stroke(); });
    heat(res.pot, st.K, p2[0], p2[1], s, 'smooth potential after FFT, filter, inverse FFT');
    const err = Math.abs(res.E - Eref) / Math.abs(Eref);
    o.E.textContent = res.E.toFixed(6); o.ref.textContent = Eref.toFixed(6);
    o.err.innerHTML = `<span style="color:${err < 1e-4 ? COL.cyan : err < 1e-2 ? COL.amber : COL.coral}">${err.toExponential(1)}</span>`;
    o.h.textContent = (L / st.K).toFixed(3) + ` σ  (${st.K}³ = ${(st.K ** 3).toLocaleString('en-US')} points)`;
    o.t.textContent = `${res.t.toFixed(1)} ms / ${tRef.toFixed(1)} ms`;
  }
  init();
  return () => ro.disconnect();
}
