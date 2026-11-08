// Structure-factor lab widget: 2D soft disks (WCA, kT = 1) evolved by Monte Carlo. S(k) directly from positions at the
// allowed k vectors vs the 2D transform of g(r), S(k) = 1 + 2 pi rho int r (g - 1) J0(kr) dr, with or without a window.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 120, KMAX = 10, RC2 = 2 ** (1 / 3);
function J0(x) {   // Abramowitz & Stegun 9.4.1 / 9.4.3
  const ax = Math.abs(x);
  if (ax <= 3) { const y = (x / 3) ** 2; return 1 - 2.2499997 * y + 1.2656208 * y ** 2 - 0.3163866 * y ** 3 + 0.0444479 * y ** 4 - 0.0039444 * y ** 5 + 0.00021 * y ** 6; }
  const y = 3 / ax, f = 0.79788456 - 0.00000077 * y - 0.0055274 * y ** 2 - 0.00009512 * y ** 3 + 0.00137237 * y ** 4 - 0.00072805 * y ** 5 + 0.00014476 * y ** 6;
  const t = ax - 0.78539816 - 0.04166397 * y - 0.00003954 * y ** 2 + 0.00262573 * y ** 3 - 0.00054125 * y ** 4 - 0.00029333 * y ** 5 + 0.00013558 * y ** 6;
  return f * Math.cos(t) / Math.sqrt(ax);
}

export function mountLab13(root) {
  const st = { rho: 0.7, win: false, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view sk-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Density ρ</span><output data-o="rho"></output><input type="range" data-k="rho" min="0.1" max="0.85" step="0.05" value="0.7"></label>
    <div class="seg"><button data-w="0" class="on">Raw transform</button><button data-w="1">With window</button></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span>Smallest allowed k</span><b data-o="kmin"></b></div>
      <div><span>Direct S at the smallest k</span><b data-o="s0"></b></div>
      <div><span>Transform at k → 0</span><b data-o="t0"></b></div>
      <div><span>Main peak (k, S)</span><b data-o="pk"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.sk-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let seed = 31; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  let L, x, kv, kmag, Sacc, nS, gh, ng, sweeps, step = 0.15;
  const NBR = 100;
  function reset() {
    L = Math.sqrt(N / st.rho); const n = Math.ceil(Math.sqrt(N)), a = L / n; x = [];
    for (let i = 0; i < N; i++) x.push([(i % n + 0.5) * a, (Math.floor(i / n) + 0.5) * a]);
    const nmax = Math.floor(KMAX * L / (2 * Math.PI)); kv = []; kmag = [];
    for (let p = 0; p <= nmax; p++) for (let q = -nmax; q <= nmax; q++) { if (p === 0 && q <= 0) continue; const kx = 2 * Math.PI * p / L, ky = 2 * Math.PI * q / L, k = Math.hypot(kx, ky); if (k <= KMAX) { kv.push([kx, ky]); kmag.push(k); } }
    Sacc = new Float64Array(kv.length); nS = 0; gh = new Float64Array(NBR); ng = 0; sweeps = 0;
  }
  const pairU = r2 => { if (r2 >= RC2) return 0; const s6 = 1 / (r2 * r2 * r2); return 4 * (s6 * s6 - s6) + 1; };
  function energyOf(i, p) { let u = 0; for (let j = 0; j < N; j++) { if (j === i) continue; let dx = x[j][0] - p[0], dy = x[j][1] - p[1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); u += pairU(dx * dx + dy * dy); } return u; }
  function sweep() {
    let acc = 0;
    for (let m = 0; m < N; m++) { const i = Math.floor(rand() * N), old = x[i], nw = [(old[0] + step * (rand() - 0.5) + L) % L, (old[1] + step * (rand() - 0.5) + L) % L];
      const dU = energyOf(i, nw) - energyOf(i, old); if (dU <= 0 || rand() < Math.exp(-dU)) { x[i] = nw; acc++; } }
    step *= acc / N > 0.4 ? 1.05 : 0.95; step = Math.min(step, L / 4); sweeps++;
    if (sweeps > 50 && sweeps % 2 === 0) measure();
  }
  function measure() {
    for (let v = 0; v < kv.length; v++) { let c = 0, s = 0; const [kx, ky] = kv[v]; for (let i = 0; i < N; i++) { const ph = kx * x[i][0] + ky * x[i][1]; c += Math.cos(ph); s += Math.sin(ph); } Sacc[v] += (c * c + s * s) / N; }
    nS++; const dr = L / 2 / NBR;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { let dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); const b = Math.floor(Math.hypot(dx, dy) / dr); if (b < NBR) gh[b]++; }
    ng++;
  }
  function gofr() { const dr = L / 2 / NBR; return Array.from({ length: NBR }, (_, b) => { const r0 = b * dr, r1 = r0 + dr; return ng ? gh[b] / ng / (st.rho * Math.PI * (r1 * r1 - r0 * r0) * N / 2) : 0; }); }
  function transform(g, k) { const dr = L / 2 / NBR, R = L / 2; let s = 0; for (let b = 0; b < NBR; b++) { const r = (b + 0.5) * dr, w = st.win ? Math.sin(Math.PI * r / R) / (Math.PI * r / R) : 1; s += r * (g[b] - 1) * J0(k * r) * w * dr; } return 1 + 2 * Math.PI * st.rho * s; }
  $('[data-k=rho]').oninput = e => { st.rho = +e.target.value; reset(); draw(); };
  root.querySelectorAll('[data-w]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-w]').forEach(z => z.classList.toggle('on', z === b)); st.win = b.dataset.w === '1'; draw(); });
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !x) return;
    cx.clearRect(0, 0, W, H);
    const half = H * 0.42, bw = Math.min(W * 0.4, half - 16), S = bw / L;
    cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.strokeRect(8, 8, bw, bw);
    x.forEach(p => { cx.fillStyle = 'rgba(56,225,255,.75)'; cx.beginPath(); cx.arc(8 + p[0] * S, 8 + p[1] * S, Math.max(1.5, 0.5 * S), 0, 7); cx.fill(); });
    // g(r)
    const g = gofr(), gx0 = bw + 50, gx1 = W - 12, gy0 = 22, gy1 = half, gmax = Math.max(2, ...g) * 1.05, GX = r => gx0 + r / (L / 2) * (gx1 - gx0), GY = v => gy1 - v / gmax * (gy1 - gy0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(gx0, gy0); cx.lineTo(gx0, gy1); cx.lineTo(gx1, gy1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`g(r) up to L/2 = ${(L / 2).toFixed(1)}`, gx0 + 6, gy0 - 6);
    cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.beginPath(); cx.moveTo(gx0, GY(1)); cx.lineTo(gx1, GY(1)); cx.stroke();
    if (ng) { cx.strokeStyle = COL.cyan; cx.lineWidth = 1.8; cx.beginPath(); g.forEach((v, b) => { const X = GX((b + 0.5) * L / 2 / NBR); b ? cx.lineTo(X, GY(v)) : cx.moveTo(X, GY(v)); }); cx.stroke(); }
    // S(k)
    const sx0 = 46, sx1 = W - 12, sy0 = half + 46, sy1 = H - 28, nb = 40, bins = new Float64Array(nb), cnt = new Float64Array(nb);
    if (nS) kmag.forEach((k, v) => { const b = Math.min(nb - 1, Math.floor(k / KMAX * nb)); bins[b] += Sacc[v] / nS; cnt[b]++; });
    const dir = Array.from({ length: nb }, (_, b) => cnt[b] ? bins[b] / cnt[b] : NaN), tr = Array.from({ length: 101 }, (_, i) => transform(g, 1e-3 + KMAX * i / 100));
    const smax = Math.max(2, ...dir.filter(isFinite), ...(ng ? tr : [0])) * 1.08, smin = Math.min(0, ...(ng ? tr : [0])), SX = k => sx0 + k / KMAX * (sx1 - sx0), SY = v => sy1 - (v - smin) / (smax - smin) * (sy1 - sy0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(sx0, sy0); cx.lineTo(sx0, sy1); cx.lineTo(sx1, sy1); cx.stroke();
    cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.beginPath(); cx.moveTo(sx0, SY(0)); cx.lineTo(sx1, SY(0)); cx.stroke();
    cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.fillText(W < 500 ? 'S(k): direct (dots), from g(r) (line)' : 'S(k): directly from positions (dots) and transform of g(r) (line)', sx0 + 6, sy0 - 10);
    const kmin = 2 * Math.PI / L; cx.fillStyle = 'rgba(255,107,107,.12)'; cx.fillRect(sx0, sy0, SX(kmin) - sx0, sy1 - sy0);
    if (ng) { cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); tr.forEach((v, i) => { const X = SX(KMAX * i / 100); i ? cx.lineTo(X, SY(v)) : cx.moveTo(X, SY(v)); }); cx.stroke(); }
    dir.forEach((v, b) => { if (isFinite(v)) { cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(SX((b + 0.5) * KMAX / nb), SY(v), 3.5, 0, 7); cx.fill(); } });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0, 2, 4, 6, 8, 10].forEach(k => cx.fillText(k, SX(k), sy1 + 14));
    cx.textAlign = 'right'; [0, 1, 2].forEach(v => cx.fillText(v, sx0 - 4, SY(v) + 4));
    const s0 = nS ? kmag.reduce((a, k, v) => Math.abs(k - kmin) < 1e-6 ? a.concat(Sacc[v] / nS) : a, []) : [];
    let pk = 0; dir.forEach((v, b) => { if (isFinite(v) && v > (dir[pk] || 0)) pk = b; });
    o.rho.textContent = st.rho.toFixed(2); o.kmin.textContent = kmin.toFixed(3);
    o.s0.textContent = s0.length ? (s0.reduce((a, b) => a + b, 0) / s0.length).toFixed(3) : '…'; o.t0.textContent = ng ? tr[0].toFixed(3) : '…';
    o.pk.textContent = nS ? `${((pk + 0.5) * KMAX / nb).toFixed(2)}, ${dir[pk].toFixed(2)}` : '…';
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let k = 0; k < 6; k++) sweep(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
