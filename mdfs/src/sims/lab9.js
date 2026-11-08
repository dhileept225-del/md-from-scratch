// Free-energy lab widget, two modes.
// FEP: spring k0 = 1 -> k1 (kT = 1); forward and backward Zwanzig estimates over 30 repeats vs exact (1/2) ln(k1/k0).
// Umbrella: double well U = 5 (s^2 - 1)^2, n windows with spring k; Langevin windows, WHAM, F(s) vs exact.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const U = s => 5 * (s * s - 1) ** 2, dU = s => 20 * s * (s * s - 1);
const NB = 40, B0 = -1.6, BW = 3.2 / NB, BM = Array.from({ length: NB }, (_, i) => B0 + (i + 0.5) * BW);

export function mountLab9(root) {
  const st = { mode: 'fep', k1: 100, n: 1000, nw: 9, kw: 30, steps: 50000 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view fe-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="fep" class="on">FEP: overlap</button><button data-m="us">Umbrella + WHAM</button></div>
    <div data-p="fep">
      <label><span>Final stiffness k₁ (start k₀ = 1)</span><output data-o="k1"></output><input type="range" data-k="k1" min="0" max="3" step="0.05" value="2"></label>
      <label><span>Samples per estimate</span><output data-o="n"></output><input type="range" data-k="n" min="1" max="5" step="1" value="3"></label>
    </div>
    <div data-p="us" hidden>
      <label><span>Number of windows</span><output data-o="nw"></output><input type="range" data-k="nw" min="3" max="21" step="1" value="9"></label>
      <label><span>Window spring k</span><output data-o="kw"></output><input type="range" data-k="kw" min="5" max="100" step="5" value="30"></label>
      <label><span>Steps per window</span><output data-o="steps"></output><input type="range" data-k="steps" min="5000" max="200000" step="5000" value="50000"></label>
    </div>
    <div class="row"><button class="btn sm primary" data-k="go">Run</button></div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.fe-view'), $ = s => root.querySelector(s);
  let seed = 0x9E3779B9; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };   // mulberry32
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let fw = [], bw = [], us = null;
  function runFep() {
    fw = []; bw = []; const k1 = st.k1, n = st.n;
    for (let r = 0; r < 30; r++) {
      let a = 0, b = 0; for (let i = 0; i < n; i++) { const x0 = gauss(), x1 = gauss() / Math.sqrt(k1); a += Math.exp(-0.5 * (k1 - 1) * x0 * x0); b += Math.exp(0.5 * (k1 - 1) * x1 * x1); }
      fw.push(-Math.log(a / n)); bw.push(Math.log(b / n));
    }
  }
  function runUs() {
    const cen = Array.from({ length: st.nw }, (_, i) => -1.2 + 2.4 * i / (st.nw - 1)), H = [], Nw = [], dt = 0.001, c = Math.exp(-5 * dt), sn = Math.sqrt(1 - c * c);
    cen.forEach(s0 => { const h = new Float64Array(NB); let s = s0, v = 0; const F = x => -dU(x) - st.kw * (x - s0);
      for (let k = 0; k < st.steps + 2000; k++) { v += 0.5 * dt * F(s); s += dt * v; v += 0.5 * dt * F(s); v = c * v + sn * gauss();
        if (k >= 2000) { const b = Math.floor((s - B0) / BW); if (b >= 0 && b < NB) h[b]++; } }
      H.push(h); Nw.push(h.reduce((a, q) => a + q, 0)); });
    const bias = cen.map(s0 => BM.map(x => 0.5 * st.kw * (x - s0) ** 2)); let f = new Array(cen.length).fill(0), P = new Array(NB).fill(0), it = 0;
    for (it = 0; it < 5000; it++) {
      for (let j = 0; j < NB; j++) { let num = 0, den = 0; for (let i = 0; i < cen.length; i++) { num += H[i][j]; den += Nw[i] * Math.exp(-(bias[i][j] - f[i])); } P[j] = num / den; }
      const fn = cen.map((_, i) => -Math.log(P.reduce((a, p, j) => a + p * Math.exp(-bias[i][j]), 0))), d0 = fn[0]; let diff = 0;
      fn.forEach((q, i) => { fn[i] = q - d0; diff = Math.max(diff, Math.abs(fn[i] - f[i])); }); f = fn; if (diff < 1e-6) break;
    }
    const F = P.map(p => p > 0 ? -Math.log(p) : NaN); let mnL = Infinity; BM.forEach((x, j) => { if (x > -1.3 && x < -0.7 && isFinite(F[j])) mnL = Math.min(mnL, F[j]); });
    const top = F[BM.reduce((bi, x, j) => Math.abs(x) < Math.abs(BM[bi]) ? j : bi, 0)];
    us = { cen, H, F: F.map(q => q - mnL), barrier: top - mnL, it, empty: BM.filter((x, j) => Math.abs(x) < 1.2 && H.every(h => h[j] === 0)).length };
  }
  const fresh = () => { st.mode === 'fep' ? runFep() : runUs(); draw(); };
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); if (st.mode === 'us' && !us) runUs(); draw(); });
  $('[data-k=k1]').oninput = e => { st.k1 = 10 ** +e.target.value; runFep(); draw(); };
  $('[data-k=n]').oninput = e => { st.n = 10 ** +e.target.value; runFep(); draw(); };
  ['nw', 'kw', 'steps'].forEach(k => $(`[data-k=${k}]`).oninput = e => { st[k] = +e.target.value; $(`[data-o=${k}]`).textContent = st[k].toLocaleString(); });
  $('[data-k=go]').onclick = fresh;

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x, y); };
  const mean = a => a.reduce((s, q) => s + q, 0) / a.length, sd = a => { const m = mean(a); return Math.sqrt(a.reduce((s, q) => s + (q - m) ** 2, 0) / (a.length - 1)); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const R = $('[data-ro]'), half = H * 0.48;
    if (st.mode === 'fep') {
      const k1 = st.k1, ex = 0.5 * Math.log(k1);
      // distributions of x in both states
      const x0 = 46, x1 = W - 12, y0 = 22, y1 = half, X = x => x0 + (x + 3.5) / 7 * (x1 - x0), pmax = Math.sqrt(Math.max(1, k1) / (2 * Math.PI)), Y = p => y1 - Math.min(p, pmax) / pmax * (y1 - y0);
      axes(x0, y0, x1, y1); label('where each state lives: k₀ = 1 (cyan), k₁ (coral)', x0 + 6, y0 - 6);
      [[1, COL.cyan], [k1, COL.coral]].forEach(([k, col]) => { cx.strokeStyle = col; cx.lineWidth = 2; cx.beginPath(); for (let i = 0; i <= 300; i++) { const x = -3.5 + 7 * i / 300, p = Math.sqrt(k / (2 * Math.PI)) * Math.exp(-0.5 * k * x * x); i ? cx.lineTo(X(x), Y(p)) : cx.moveTo(X(x), Y(p)); } cx.stroke(); });
      // estimates
      const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 28, all = [...fw, ...bw, ex], lo = Math.min(...all) - 0.2, hi = Math.max(...all) + 0.2, EX = v => bx0 + (v - lo) / (hi - lo) * (bx1 - bx0);
      axes(bx0, by0, bx1, by1); label(`30 estimates of ΔF with ${st.n.toLocaleString()} samples each`, bx0 + 6, by0 - 10);
      cx.setLineDash([5, 4]); cx.strokeStyle = COL.ink; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(EX(ex), by0); cx.lineTo(EX(ex), by1); cx.stroke(); cx.setLineDash([]);
      const yF = by0 + (by1 - by0) * 0.3, yB = by0 + (by1 - by0) * 0.7;
      fw.forEach(v => { cx.fillStyle = 'rgba(56,225,255,.7)'; cx.beginPath(); cx.arc(EX(v), yF + (rand() - 0.5) * 18, 4, 0, 7); cx.fill(); });
      bw.forEach(v => { cx.fillStyle = 'rgba(255,107,107,.7)'; cx.beginPath(); cx.arc(EX(v), yB + (rand() - 0.5) * 18, 4, 0, 7); cx.fill(); });
      label('forward (sample k₀)', bx0 + 6, yF - 14, COL.cyan, 'left', '11px'); label('backward (sample k₁)', bx0 + 6, yB - 14, COL.coral, 'left', '11px');
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(`exact ${ex.toFixed(3)}`, EX(ex), by1 + 14);
      R.innerHTML = `<div><span>Exact free-energy change</span><b>${ex.toFixed(3)}</b></div><div><span>Forward</span><b>${mean(fw).toFixed(3)} ± ${sd(fw).toFixed(3)}</b></div>
        <div><span>Backward</span><b>${mean(bw).toFixed(3)} ± ${sd(bw).toFixed(3)}</b></div><div><span>Stiffness ratio</span><b>${k1.toFixed(1)}</b></div>`;
      $('[data-o=k1]').textContent = k1.toFixed(1); $('[data-o=n]').textContent = st.n.toLocaleString();
    } else {
      const x0 = 46, x1 = W - 12, y0 = 22, y1 = half, X = s => x0 + (s + 1.6) / 3.2 * (x1 - x0);
      let hmax = 1; us.H.forEach(h => h.forEach(q => hmax = Math.max(hmax, q)));
      axes(x0, y0, x1, y1); label(`histograms of ${us.cen.length} windows`, x0 + 6, y0 - 6);
      us.H.forEach((h, i) => { cx.strokeStyle = `hsla(${190 + 150 * i / us.H.length},90%,65%,.85)`; cx.lineWidth = 1.5; cx.beginPath(); h.forEach((q, j) => { const yy = y1 - q / hmax * (y1 - y0); j ? cx.lineTo(X(BM[j]), yy) : cx.moveTo(X(BM[j]), yy); }); cx.stroke(); });
      const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 28, FY = f => by1 - Math.min(8, Math.max(-0.5, f)) / 8 * (by1 - by0);
      axes(bx0, by0, bx1, by1); label('free energy: exact (white) vs WHAM (amber)', bx0 + 6, by0 - 10);
      cx.strokeStyle = 'rgba(242,245,250,.8)'; cx.lineWidth = 1.5; cx.beginPath(); for (let i = 0; i <= 200; i++) { const s = -1.6 + 3.2 * i / 200; i ? cx.lineTo(X(s), FY(U(s))) : cx.moveTo(X(s), FY(U(s))); } cx.stroke();
      cx.fillStyle = COL.amber; BM.forEach((x, j) => { if (isFinite(us.F[j])) { cx.beginPath(); cx.arc(X(x), FY(us.F[j]), 3.5, 0, 7); cx.fill(); } });
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [-1, 0, 1].forEach(s => cx.fillText(s, X(s), by1 + 14));
      R.innerHTML = `<div><span>Barrier from WHAM (exact 5)</span><b>${us.barrier.toFixed(2)}</b></div><div><span>WHAM iterations</span><b>${us.it}</b></div>
        <div><span>Empty bins between windows</span><b>${us.empty}</b></div><div><span>Total steps</span><b>${(us.cen.length * st.steps).toLocaleString()}</b></div>`;
    }
  }
  ['nw', 'kw', 'steps'].forEach(k => $(`[data-o=${k}]`).textContent = st[k].toLocaleString());
  runFep(); draw();
  return () => ro.disconnect();
}
