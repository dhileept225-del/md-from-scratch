// Metadynamics lab: Brownian dynamics (kT = 1) on a 2D surface with two basins, U = H (x^2 - 1)^2 + 6 (y - 0.5 x)^2.
// Gaussian hills are added along the chosen collective variable (x, or the poorer y). The bottom panel compares the
// exact free energy along that CV, F(cv) = -ln of the integral of exp(-U) over the other coordinate, with the estimate.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const XL = -2, XH = 2, YL = -1.6, YH = 1.6, NB = 200, DT = 1e-3;

export function mountMetad(root) {
  const st = { method: 'wt', cv: 'x', H: 6, w: 0.4, s: 0.12, stride: 200, gamma: 8, running: true };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="plain">Plain MD</button><button data-m="meta">Metadynamics</button><button data-m="wt" class="on">Well-tempered</button></div>
    <div class="seg"><button data-c="x" class="on">CV = x (good)</button><button data-c="y">CV = y (poor)</button></div>
    <label><span>Barrier height (kT)</span><output data-o="H"></output><input type="range" data-k="H" min="2" max="12" step="0.5" value="6"></label>
    <label><span>Hill height (kT)</span><output data-o="w"></output><input type="range" data-k="w" min="0.05" max="1.5" step="0.05" value="0.4"></label>
    <label><span>Hill width</span><output data-o="s"></output><input type="range" data-k="s" min="0.04" max="0.4" step="0.01" value="0.12"></label>
    <label><span>Bias factor <i class="nc">γ</i> (well-tempered)</span><output data-o="g"></output><input type="range" data-k="g" min="2" max="30" step="1" value="8"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span>Hills / time</span><b data-o="nh"></b></div>
      <div><span>Barrier crossings</span><b data-o="cr"></b></div>
      <div><span>Barrier along CV: estimate / exact</span><b data-o="bar"></b></div>
      <div><span>RMS error of F (basins)</span><b data-o="err"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 17; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  const U = (x, y) => st.H * (x * x - 1) ** 2 + 6 * (y - 0.5 * x) ** 2;
  const dU = (x, y) => [4 * st.H * x * (x * x - 1) - 3 * (y - 0.5 * x), 12 * (y - 0.5 * x)];
  let p, V, dV, hills, steps, side, crossings, trail, Fex, mapImg;
  const lo = () => st.cv === 'x' ? XL : YL, hi = () => st.cv === 'x' ? XH : YH, dq = () => (hi() - lo()) / (NB - 1), qOf = (k) => lo() + k * dq();
  function exactF() {   // F(cv) = -ln integral exp(-U) d(other)
    const F = new Float64Array(NB);
    for (let k = 0; k < NB; k++) { const q = qOf(k); let z = 0; for (let j = 0; j < 400; j++) { const t = (st.cv === 'x' ? YL : XL) + j / 399 * ((st.cv === 'x' ? YH : XH) - (st.cv === 'x' ? YL : XL)); z += Math.exp(-(st.cv === 'x' ? U(q, t) : U(t, q))); } F[k] = -Math.log(z + 1e-300); }
    const m = Math.min(...F); return F.map(f => f - m);
  }
  function reset() { p = [-1, -0.5]; V = new Float64Array(NB); dV = new Float64Array(NB); hills = 0; steps = 0; side = -1; crossings = 0; trail = []; Fex = exactF(); mapImg = null; }
  function step() {
    const q = st.cv === 'x' ? p[0] : p[1], kf = (q - lo()) / dq(), k = Math.max(0, Math.min(NB - 2, Math.floor(kf))), f = Math.max(0, Math.min(1, kf - k)), inside = q >= lo() && q <= hi();
    const g = dU(p[0], p[1]), bf = inside && st.method !== 'plain' ? dV[k] * (1 - f) + dV[k + 1] * f : 0;
    if (st.cv === 'x') g[0] += bf; else g[1] += bf;
    p[0] += -g[0] * DT + Math.sqrt(2 * DT) * gauss(); p[1] += -g[1] * DT + Math.sqrt(2 * DT) * gauss();
    steps++;
    if (st.method !== 'plain' && inside && steps % st.stride === 0) {
      const Vq = V[k] * (1 - f) + V[k + 1] * f, w = st.method === 'wt' ? st.w * Math.exp(-Vq / (st.gamma - 1)) : st.w;   // dT = (gamma - 1) kT
      for (let j = 0; j < NB; j++) { const d = qOf(j) - q, gg = w * Math.exp(-d * d / (2 * st.s * st.s)); V[j] += gg; dV[j] += -gg * d / (st.s * st.s); }
      hills++;
    }
    if (p[0] > 0.8 && side < 0) { side = 1; crossings++; } else if (p[0] < -0.8 && side > 0) { side = -1; crossings++; }
    if (steps % 20 === 0) { trail.push([p[0], p[1]]); if (trail.length > 400) trail.shift(); }
  }
  const Fest = () => { const fac = st.method === 'wt' ? st.gamma / (st.gamma - 1) : 1, F = Array.from(V, v => -fac * v); const m = Math.min(...F.filter((_, k) => Fex[k] < 12)); return F.map(f => f - m); };

  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.method = b.dataset.m; root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); reset(); });
  root.querySelectorAll('[data-c]').forEach(b => b.onclick = () => { st.cv = b.dataset.c; root.querySelectorAll('[data-c]').forEach(z => z.classList.toggle('on', z === b)); reset(); });
  const bind = (k, key) => { $(`[data-k=${k}]`).oninput = e => { st[key] = +e.target.value; reset(); }; };
  bind('H', 'H'); bind('w', 'w'); bind('s', 's'); bind('g', 'gamma');
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => reset();

  let W = 600, Hc = 600, M = { x: 0, y: 0, w: 1, h: 1 };
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; Hc = view.clientHeight; cv.width = W * d; cv.height = Hc * d; cx.setTransform(d, 0, 0, d, 0, 0); mapImg = null; draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function buildMap() {
    const topH = Hc * 0.5 - 10, wmap = Math.min(W - 20, topH * (XH - XL) / (YH - YL)); M = { x: (W - wmap) / 2, y: 6, w: wmap, h: wmap * (YH - YL) / (XH - XL) };
    const n = 100, m = Math.round(n * (YH - YL) / (XH - XL)), off = document.createElement('canvas'); off.width = n; off.height = m; const oc = off.getContext('2d'), id = oc.createImageData(n, m);
    for (let j = 0; j < m; j++) for (let i = 0; i < n; i++) { const x = XL + (i + 0.5) / n * (XH - XL), y = YH - (j + 0.5) / m * (YH - YL), t = Math.min(1, U(x, y) / (st.H + 6)), q = 4 * (j * n + i);
      const lo_ = [40, 150, 230], hi_ = [10, 18, 40]; for (let c = 0; c < 3; c++) id.data[q + c] = lo_[c] * (1 - t) + hi_[c] * t; id.data[q + 3] = 255; }
    oc.putImageData(id, 0, 0); mapImg = off;
  }
  const toS = (x, y) => [M.x + (x - XL) / (XH - XL) * M.w, M.y + (YH - y) / (YH - YL) * M.h];
  function draw() {
    if (!W) return;
    if (!mapImg) buildMap();
    cx.clearRect(0, 0, W, Hc); cx.imageSmoothingEnabled = true; cx.drawImage(mapImg, M.x, M.y, M.w, M.h);
    cx.strokeStyle = 'rgba(242,245,250,.6)'; cx.lineWidth = 1.2; cx.beginPath(); trail.forEach((t, q) => { const [X, Y] = toS(t[0], t[1]); q ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke();
    const [PX, PY] = toS(p[0], p[1]); cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(PX, PY, 6, 0, 7); cx.fill();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`U(x, y): two basins; hills are added along ${st.cv}`, M.x + 6, M.y + M.h - 8);
    // free-energy panel
    const x0 = 46, x1 = W - 14, yt = M.y + M.h + 40, yb = Hc - 26, Fm = 14, X = q => x0 + (q - lo()) / (hi() - lo()) * (x1 - x0), Y = f => yb - Math.min(f, Fm) / Fm * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillText(`free energy along ${st.cv} (kT)`, x0 + 6, yt - 10);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [5, 10].forEach(f => cx.fillText(f, x0 - 5, Y(f) + 4)); cx.textAlign = 'center'; [-1, 0, 1].forEach(q => cx.fillText(q, X(q), yb + 14));
    const line = (arr, c, w, dash) => { cx.setLineDash(dash || []); cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); arr.forEach((f, k) => { k ? cx.lineTo(X(qOf(k)), Y(f)) : cx.moveTo(X(qOf(k)), Y(f)); }); cx.stroke(); cx.setLineDash([]); };
    line(Fex, COL.cyan, 2.5);
    let est = null; if (st.method !== 'plain' && hills > 3) { est = Fest(); line(est, COL.amber, 2.2, [6, 4]); line(Array.from(V), 'rgba(255,107,107,.6)', 1.5); }
    cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [[COL.cyan, 'exact'], ...(est ? [[COL.amber, 'estimate from the bias'], ['rgba(255,107,107,.9)', 'bias V']] : [])].forEach(([c, t], j) => { cx.fillStyle = c; cx.fillText(t, x1, yt + 4 + 14 * j); });
    // readouts
    o.H.textContent = st.H; o.w.textContent = st.w.toFixed(2); o.s.textContent = st.s.toFixed(2); o.g.textContent = st.gamma;
    o.nh.textContent = `${hills} / ${(steps * DT).toFixed(0)}`; o.cr.textContent = crossings;
    const kmid = Math.round((0 - lo()) / dq()), basins = Fex.map((f, k) => f < 4 ? k : -1).filter(k => k >= 0);
    const exB = Fex[kmid] - 0;
    if (est) { const rms = Math.sqrt(basins.reduce((s, k) => s + (est[k] - Fex[k]) ** 2, 0) / basins.length); o.bar.textContent = `${(est[kmid]).toFixed(1)} / ${exB.toFixed(1)} kT`; o.err.textContent = rms.toFixed(2) + ' kT'; }
    else { o.bar.textContent = `— / ${exB.toFixed(1)} kT`; o.err.textContent = st.method === 'plain' ? 'no bias' : '…'; }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 60; k++) step(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
