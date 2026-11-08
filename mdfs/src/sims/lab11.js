// Viscosity lab widget: running Green-Kubo integral of model shear-stress signals.
// Each stress component is an Ornstein-Uhlenbeck signal with the size (C(0) = 0.0722) and decay time of the LJ liquid
// of the lab (eta = V/kT C(0) tau = 2.42, V = 320). Independent pieces show the error bar; 1 or 5 components.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const V = 320, C0 = 0.07218, ETA = 2.42, TAU = ETA / (V * C0), DT = 0.02;
const TS = Array.from({ length: 36 }, (_, k) => 0.02 * Math.pow(400, k / 35));   // upper limits 0.02 ... 8

export function mountLab11(root) {
  const st = { len: 80, comp: 5, pieces: 8 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view vi-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Length of each run (time units)</span><output data-o="len"></output><input type="range" data-k="len" min="10" max="400" step="10" value="80"></label>
    <div class="seg"><button data-c="1">1 component</button><button data-c="5" class="on">5 components</button></div>
    <div class="row"><button class="btn sm primary" data-k="go">New runs</button></div>
    <div class="readouts">
      <div><span>Viscosity at t = 1.5, mean of runs</span><b data-o="m"></b></div>
      <div><span>Spread between runs</span><b data-o="s"></b></div>
      <div><span>True value</span><b>2.42</b></div>
      <div><span>Steps per run (dt = 0.004)</span><b data-o="n"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.vi-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let seed = 4242; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };
  const gauss = () => Math.sqrt(-2 * Math.log(rand())) * Math.cos(2 * Math.PI * rand());
  let curves = [];
  function runningIntegral(P) {   // I(T) = V/kT * mean_n P_n * sum_{k=0..K} P_{n+k} dt  (rectangle rule)
    const n = P.length, S = new Float64Array(n + 1); for (let i = 0; i < n; i++) S[i + 1] = S[i] + P[i];
    return TS.map(T => { const K = Math.max(1, Math.round(T / DT)); let s = 0; const m = n - K; for (let i = 0; i < m; i++) s += P[i] * (S[i + K] - S[i] - 0.5 * P[i]); return V * s / m * DT; });
  }
  function simulate() {
    const n = Math.round(st.len / DT), a = Math.exp(-DT / TAU), sn = Math.sqrt(C0 * (1 - a * a)); curves = [];
    for (let p = 0; p < st.pieces; p++) {
      const acc = new Array(TS.length).fill(0);
      for (let c = 0; c < st.comp; c++) { const P = new Float64Array(n); P[0] = Math.sqrt(C0) * gauss(); for (let i = 1; i < n; i++) P[i] = a * P[i - 1] + sn * gauss();
        runningIntegral(P).forEach((v, k) => acc[k] += v / st.comp); }
      curves.push(acc);
    }
  }
  $('[data-k=len]').oninput = e => { st.len = +e.target.value; o.len.textContent = st.len; };
  root.querySelectorAll('[data-c]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-c]').forEach(z => z.classList.toggle('on', z === b)); st.comp = +b.dataset.c; simulate(); draw(); });
  $('[data-k=go]').onclick = () => { simulate(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !curves.length) return;
    cx.clearRect(0, 0, W, H);
    const x0 = 46, x1 = W - 12, y0 = 24, y1 = H * 0.62, X = t => x0 + Math.log(t / 0.02) / Math.log(400) * (x1 - x0), ymax = 5, Y = v => y1 - Math.max(0, Math.min(ymax, v)) / ymax * (y1 - y0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('running Green–Kubo integral η(t), log time axis', x0 + 6, y0 - 8);
    cx.setLineDash([5, 4]); cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(x0, Y(ETA)); cx.lineTo(x1, Y(ETA)); cx.stroke();
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.beginPath(); cx.moveTo(X(1.5), y0); cx.lineTo(X(1.5), y1); cx.stroke(); cx.setLineDash([]);
    curves.forEach(c => { cx.strokeStyle = 'rgba(143,163,199,.55)'; cx.lineWidth = 1.2; cx.beginPath(); c.forEach((v, k) => k ? cx.lineTo(X(TS[k]), Y(v)) : cx.moveTo(X(TS[k]), Y(v))); cx.stroke(); });
    const mean = TS.map((_, k) => curves.reduce((s, c) => s + c[k], 0) / curves.length);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2.5; cx.beginPath(); mean.forEach((v, k) => k ? cx.lineTo(X(TS[k]), Y(v)) : cx.moveTo(X(TS[k]), Y(v))); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0.02, 0.1, 0.5, 1.5, 8].forEach(t => cx.fillText(t, X(t), y1 + 14));
    cx.textAlign = 'right'; [0, 2.42, 5].forEach(v => cx.fillText(v, x0 - 4, Y(v) + 4));
    // plateau values
    const k15 = TS.reduce((bi, t, k) => Math.abs(t - 1.5) < Math.abs(TS[bi] - 1.5) ? k : bi, 0), vals = curves.map(c => c[k15]);
    const m = vals.reduce((a, b) => a + b, 0) / vals.length, s = Math.sqrt(vals.reduce((a, b) => a + (b - m) ** 2, 0) / (vals.length - 1));
    const bx0 = 46, bx1 = W - 12, by = H * 0.62 + 60, BX = v => bx0 + Math.max(0, Math.min(6, v)) / 6 * (bx1 - bx0);
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`η at t = 1.5 from ${vals.length} independent runs`, bx0, by - 22);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(bx0, by); cx.lineTo(bx1, by); cx.stroke();
    cx.setLineDash([5, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(BX(ETA), by - 14); cx.lineTo(BX(ETA), by + 14); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = 'rgba(56,225,255,.25)'; cx.fillRect(BX(m - s), by - 8, BX(m + s) - BX(m - s), 16);
    vals.forEach(v => { cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(BX(v), by, 5, 0, 7); cx.fill(); });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0, 2, 4, 6].forEach(v => cx.fillText(v, BX(v), by + 26));
    o.m.textContent = m.toFixed(2); o.s.textContent = `± ${s.toFixed(2)} (${(100 * s / ETA).toFixed(0)} %)`; o.n.textContent = Math.round(st.len / 0.004).toLocaleString(); o.len.textContent = st.len;
  }
  simulate(); draw();
  return () => ro.disconnect();
}
