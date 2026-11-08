// Ewald lab: 32 ions (16+, 16−) in a periodic cubic box of side L = 1 (units: q²/L, 1/(4πε0) = 1).
// Choose the splitting α, the real-space cutoff and the number of k vectors; see the three parts, the error
// against a converged reference, and the error as a function of α (the classic U-shaped curve).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 32, L = 1;
function erfc(z) {                    // accurate to ~1e-14: series for small z, continued fraction for large z
  if (z < 2.5) { let term = z, sum = z; for (let n = 1; n < 80; n++) { term *= -z * z / n; const add = term / (2 * n + 1); sum += add; if (Math.abs(add) < 1e-17) break; } return 1 - 2 / Math.sqrt(Math.PI) * sum; }
  let f = 0; for (let n = 60; n >= 1; n--) f = n / 2 / (z + f); return Math.exp(-z * z) / Math.sqrt(Math.PI) / (z + f);
}

export function mountEwald(root) {
  const st = { alpha: 7, rc: 0.45, kmax: 7 };
  let X = [], Q = [], ref = 0;
  function makeSystem(seed) {
    let s = seed; const R = () => (s = (s * 16807) % 2147483647) / 2147483647;
    X = []; Q = [];
    while (X.length < N) { const p = [R(), R(), R()]; if (X.every(q => { let d = 0; for (let k = 0; k < 3; k++) { let x = p[k] - q[k]; x -= Math.round(x); d += x * x; } return d > 0.14 ** 2; })) { X.push(p); Q.push(X.length % 2 ? 1 : -1); } }
  }
  function ewald(alpha, rc, kmax) {
    let Er = 0, npairs = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let d2 = 0; for (let k = 0; k < 3; k++) { let x = X[i][k] - X[j][k]; x -= L * Math.round(x / L); d2 += x * x; }
      const r = Math.sqrt(d2); if (r < rc) { Er += Q[i] * Q[j] * erfc(alpha * r) / r; npairs++; }
    }
    let Ek = 0, nk = 0; const V = L ** 3;
    for (let a = -kmax; a <= kmax; a++) for (let b = -kmax; b <= kmax; b++) for (let c = 0; c <= kmax; c++) {
      if (c === 0 && (b < 0 || (b === 0 && a <= 0))) continue;           // use k and -k symmetry
      const kx = 2 * Math.PI * a / L, ky = 2 * Math.PI * b / L, kz = 2 * Math.PI * c / L, k2 = kx * kx + ky * ky + kz * kz;
      let Sr = 0, Si = 0; for (let i = 0; i < N; i++) { const ph = kx * X[i][0] + ky * X[i][1] + kz * X[i][2]; Sr += Q[i] * Math.cos(ph); Si += Q[i] * Math.sin(ph); }
      Ek += 2 * (2 * Math.PI / V) * Math.exp(-k2 / (4 * alpha * alpha)) / k2 * (Sr * Sr + Si * Si); nk += 2;
    }
    const Es = -alpha / Math.sqrt(Math.PI) * N;
    return { Er, Ek, Es, E: Er + Ek + Es, npairs, nk };
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Splitting <i class="nc">α</i> (1/L)</span><output data-o="al"></output><input type="range" data-k="al" min="2" max="14" step="0.25" value="7"></label>
    <label><span>Real-space cutoff r<sub>c</sub> (L)</span><output data-o="rc"></output><input type="range" data-k="rc" min="0.2" max="0.5" step="0.01" value="0.45"></label>
    <label><span>k vectors per direction k<sub>max</sub></span><output data-o="km"></output><input type="range" data-k="km" min="1" max="10" step="1" value="7"></label>
    <div class="row"><button class="btn sm ghost" data-k="new">New configuration</button></div>
    <div class="readouts">
      <div><span>Real-space part</span><b data-o="er"></b></div>
      <div><span>Reciprocal part</span><b data-o="ek"></b></div>
      <div><span>Self part</span><b data-o="es"></b></div>
      <div><span>Total</span><b data-o="E"></b></div>
      <div class="wide"><span>Error vs converged reference</span><b data-o="err"></b></div>
      <div class="wide"><span>Work: real pairs / k vectors</span><b data-o="cost"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let curve = [];
  const recurve = () => { curve = []; for (let a = 2; a <= 14.001; a += 0.4) curve.push([a, Math.abs(ewald(a, st.rc, st.kmax).E - ref)]); };
  const init = seed => { makeSystem(seed); ref = ewald(8, 0.5, 14).E; recurve(); draw(); };
  $('[data-k=al]').oninput = e => { st.alpha = +e.target.value; draw(); };
  $('[data-k=rc]').oninput = e => { st.rc = +e.target.value; recurve(); draw(); };
  $('[data-k=km]').oninput = e => { st.kmax = +e.target.value; recurve(); draw(); };
  let seed = 7; $('[data-k=new]').onclick = () => init(seed += 101);

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function draw() {
    if (!W || !X.length) return;
    cx.clearRect(0, 0, W, H);
    const res = ewald(st.alpha, st.rc, st.kmax), top = H * 0.46;
    // left: box projection with cutoff around ion 0; right: bars
    const bs = Math.min(W * 0.42, top - 30), bx = 14, by = 24;
    cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.lineWidth = 1.5; cx.strokeRect(bx, by, bs, bs);
    X.forEach((p, i) => { cx.fillStyle = Q[i] > 0 ? COL.coral : COL.cyan; cx.beginPath(); cx.arc(bx + p[0] * bs, by + p[1] * bs, i === 0 ? 6 : 4, 0, 7); cx.fill(); });
    cx.save(); cx.beginPath(); cx.rect(bx, by, bs, bs); cx.clip(); cx.strokeStyle = COL.amber; cx.setLineDash([5, 4]);
    for (const sx of [-1, 0, 1]) for (const sy of [-1, 0, 1]) { cx.beginPath(); cx.arc(bx + (X[0][0] + sx) * bs, by + (X[0][1] + sy) * bs, st.rc * bs, 0, 7); cx.stroke(); }
    cx.restore(); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('box seen from above; r_c around one ion', bx, by - 8);
    const parts = [['real', res.Er, COL.cyan], ['recip', res.Ek, COL.amber], ['self', res.Es, COL.coral], ['total', res.E, COL.ink]];
    const gx0 = bx + bs + 30, gx1 = W - 10, mid = by + bs / 2, sc = (bs / 2 - 14) / Math.max(...parts.map(p => Math.abs(p[1])), 1), bw = (gx1 - gx0) / 4;
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(gx0, mid); cx.lineTo(gx1, mid); cx.stroke();
    parts.forEach(([nm, v, c], k) => { const x = gx0 + k * bw + bw * 0.15; cx.fillStyle = c; cx.fillRect(x, mid, bw * 0.7, -v * sc); cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '11px "Space Grotesk", sans-serif'; cx.fillText(nm, x + bw * 0.35, v < 0 ? mid - 6 : mid + 14); });
    // error vs alpha
    const x0 = 50, x1 = W - 12, y0 = top + 30, y1 = H - 26, lo = -10, hi = 1;
    const Xa = a => x0 + (a - 2) / 12 * (x1 - x0), Ye = e => y1 - (Math.max(lo, Math.min(hi, Math.log10(Math.max(e, 1e-12)))) - lo) / (hi - lo) * (y1 - y0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [0, -3, -6, -9].forEach(p => cx.fillText('1e' + p, x0 - 5, Ye(10 ** p) + 4));
    cx.textAlign = 'center'; [2, 4, 6, 8, 10, 12, 14].forEach(a => cx.fillText(a, Xa(a), y1 + 14));
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.cyan; cx.fillText('|error| vs α for the current r_c and k_max (log scale)', x0 + 6, y0 - 8);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.fillText('← real-space error', x0 + 6, y0 + 12); cx.textAlign = 'right'; cx.fillText('reciprocal error →', x1 - 4, y0 + 12);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); curve.forEach(([a, e], k) => (k ? cx.lineTo(Xa(a), Ye(e)) : cx.moveTo(Xa(a), Ye(e)))); cx.stroke();
    const err = Math.abs(res.E - ref); cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(Xa(st.alpha), Ye(err), 6, 0, 7); cx.fill();
    // readouts
    o.al.textContent = st.alpha.toFixed(2); o.rc.textContent = st.rc.toFixed(2) + `  (erfc(α r_c) = ${erfc(st.alpha * st.rc).toExponential(1)})`; o.km.textContent = st.kmax;
    o.er.textContent = res.Er.toFixed(5); o.ek.textContent = res.Ek.toFixed(5); o.es.textContent = res.Es.toFixed(5); o.E.textContent = res.E.toFixed(6);
    o.err.innerHTML = `<span style="color:${err < 1e-4 ? COL.cyan : COL.coral}">${err.toExponential(1)}  (reference ${ref.toFixed(6)})</span>`;
    o.cost.textContent = `${res.npairs} pairs / ${res.nk} k vectors`;
  }
  init(seed);
  return () => ro.disconnect();
}
