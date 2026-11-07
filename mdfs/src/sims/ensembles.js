// Ensemble lab: one 2D Lennard-Jones fluid run in NVE, NVT (Andersen thermostat) or NPT (plus Berendsen barostat).
// Live traces of energy, temperature and density show which quantity is fixed and which fluctuates.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RC = 2.5, DT = 0.005, SUB = 4, HIST = 400, NU = 2.0, KP = 0.004;
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountEnsembles(root) {
  const st = { ens: 'NVT', T0: 1.0, P0: 1.0, N: 100, running: true };
  let N, L, x, v, f, U = 0, Wv = 0, hist = [], acc;
  const Nf = () => 2 * N - 2;
  function forces() {
    f.fill(0); U = 0; Wv = 0; const rc2 = RC * RC, Ush = 4 * (Math.pow(RC, -12) - Math.pow(RC, -6));
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[2 * i] - x[2 * j], dy = x[2 * i + 1] - x[2 * j + 1];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r2 = dx * dx + dy * dy; if (r2 >= rc2) continue;
      const i6 = 1 / (r2 * r2 * r2), rF = 24 * (2 * i6 * i6 - i6), ff = rF / r2;
      U += 4 * (i6 * i6 - i6) - Ush; Wv += rF;
      f[2 * i] += ff * dx; f[2 * i + 1] += ff * dy; f[2 * j] -= ff * dx; f[2 * j + 1] -= ff * dy;
    }
  }
  const kin = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return 0.5 * s; };
  const resetAcc = () => { acc = { n: 0, s: [0, 0, 0, 0], q: [0, 0, 0, 0] }; };
  function restart() {
    N = st.N; L = Math.sqrt(N / 0.6); x = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = Math.ceil(Math.sqrt(N)), a = L / n;
    for (let i = 0; i < N; i++) { x[2 * i] = (i % n + 0.5) * a; x[2 * i + 1] = (Math.floor(i / n) + 0.5) * a; v[2 * i] = gauss() * Math.sqrt(st.T0); v[2 * i + 1] = gauss() * Math.sqrt(st.T0); }
    let px = 0, py = 0; for (let i = 0; i < N; i++) { px += v[2 * i]; py += v[2 * i + 1]; }
    for (let i = 0; i < N; i++) { v[2 * i] -= px / N; v[2 * i + 1] -= py / N; }
    forces(); hist = []; resetAcc();
  }
  function step() {
    let P = 0, T = 0;
    for (let s = 0; s < SUB; s++) {
      for (let i = 0; i < 2 * N; i++) { v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];
      if (st.ens !== 'NVE') for (let i = 0; i < N; i++) if (Math.random() < NU * DT) { const sd = Math.sqrt(st.T0); v[2 * i] = gauss() * sd; v[2 * i + 1] = gauss() * sd; }
      T = 2 * kin() / Nf(); P = (N * T + Wv / 2) / (L * L);   // 2D virial pressure
      if (st.ens === 'NPT') { const mu = Math.sqrt(Math.max(0.98, Math.min(1.02, 1 - KP * (st.P0 - P)))); L *= mu; for (let i = 0; i < 2 * N; i++) x[i] *= mu; }
    }
    const E = (kin() + U) / N, rho = N / (L * L), q = [E, T, P, rho];
    hist.push(q); if (hist.length > HIST) hist.shift();
    acc.n++; q.forEach((y, k) => { acc.s[k] += y; acc.q[k] += y * y; });
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view en-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-e="NVE">NVE</button><button data-e="NVT" class="on">NVT</button><button data-e="NPT">NPT</button></div>
    <label data-t>Thermostat T <output data-o="T0"></output><input type="range" data-k="T0" min="0.5" max="2" step="0.1" value="1"></label>
    <label data-p hidden>Barostat P <output data-o="P0"></output><input type="range" data-k="P0" min="0.2" max="3" step="0.1" value="1"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button><button class="btn sm ghost" data-k="stats">Reset averages</button></div>
    <div class="readouts">
      <div><span>Energy per atom</span><b data-o="E"></b></div>
      <div><span>Temperature</span><b data-o="T"></b></div>
      <div><span>Pressure</span><b data-o="P"></b></div>
      <div><span>Density</span><b data-o="rho"></b></div>
      <div class="wide"><span>Fixed in this ensemble</span><b data-o="fixed"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.en-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-e]').forEach(b => b.onclick = () => {
    st.ens = b.dataset.e; root.querySelectorAll('[data-e]').forEach(z => z.classList.toggle('on', z === b));
    $('[data-t]').hidden = st.ens === 'NVE'; $('[data-p]').hidden = st.ens !== 'NPT'; resetAcc();
  });
  $('[data-k=T0]').oninput = e => { st.T0 = +e.target.value; resetAcc(); };
  $('[data-k=P0]').oninput = e => { st.P0 = +e.target.value; resetAcc(); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = restart; $('[data-k=stats]').onclick = resetAcc;

  let W = 600, H = 640;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  const series = [['energy per atom', 0, COL.cyan, 'NVE'], ['temperature', 1, COL.coral, 'NVT'], ['density', 3, COL.amber, 'NPT']];
  function draw() {
    cx.clearRect(0, 0, W, H);
    // box (fixed scale, so the NPT box visibly grows and shrinks)
    const boxH = H * 0.42, S = boxH / Math.sqrt(st.N / 0.3), bs = L * S, ox = W / 2 - bs / 2, oy = 10 + (boxH - bs) / 2;
    cx.strokeStyle = st.ens === 'NPT' ? COL.amber : 'rgba(56,225,255,.6)'; cx.lineWidth = st.ens === 'NPT' ? 3 : 2; cx.strokeRect(ox, oy, bs, bs);
    if (st.ens !== 'NVE') { cx.strokeStyle = 'rgba(255,194,75,.25)'; cx.lineWidth = 6; cx.strokeRect(ox - 7, oy - 7, bs + 14, bs + 14); }
    const r = Math.max(2, 0.42 * S), T0 = Math.max(0.3, st.ens === 'NVE' ? 1 : st.T0);
    for (let i = 0; i < N; i++) {
      const e = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / (4 * T0);
      cx.fillStyle = e > 1 ? COL.coral : e > 0.45 ? COL.amber : COL.cyan;
      cx.beginPath(); cx.arc(ox + x[2 * i] * S, oy + x[2 * i + 1] * S, r, 0, 7); cx.fill();
    }
    // three strip charts
    const top = boxH + 30, h = (H - top - 10) / 3, x0 = 50, x1 = W - 12;
    series.forEach(([name, k, c, fixedIn], s) => {
      const y0 = top + s * h + 18, y1 = top + (s + 1) * h - 8, vals = hist.map(q => q[k]);
      cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
      const fixed = (fixedIn === 'NVE' && st.ens === 'NVE') || (k === 1 && st.ens !== 'NVE') || (k === 3 && st.ens !== 'NPT');
      cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = c;
      cx.fillText(name + (fixed ? (k === 1 ? '  (fixed on average)' : '  (fixed)') : '  (free to fluctuate)'), x0 + 6, y0 - 4);
      if (vals.length < 2) return;
      let lo = Math.min(...vals), hi = Math.max(...vals); const pad = Math.max((hi - lo) * 0.15, k === 3 ? 0.02 : 0.05); lo -= pad; hi += pad;
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right';
      cx.fillText(hi.toFixed(2), x0 - 5, y0 + 9); cx.fillText(lo.toFixed(2), x0 - 5, y1);
      cx.strokeStyle = c; cx.lineWidth = 1.6; cx.beginPath();
      vals.forEach((y, j) => { const X = x0 + (x1 - x0) * j / (HIST - 1), Y = y1 - (y - lo) / (hi - lo) * (y1 - y0); j ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke();
    });
  }
  function readouts() {
    o.T0.textContent = st.T0.toFixed(1); o.P0.textContent = st.P0.toFixed(1);
    const fmt = (k, d) => { if (acc.n < 30) return 'measuring…'; const m = acc.s[k] / acc.n, sd = Math.sqrt(Math.max(0, acc.q[k] / acc.n - m * m)); return `${m.toFixed(d)} ± ${sd.toFixed(d)}`; };
    o.E.textContent = fmt(0, 3); o.T.textContent = fmt(1, 3); o.P.textContent = fmt(2, 2); o.rho.textContent = fmt(3, 3);
    o.fixed.textContent = { NVE: 'N, V, E: temperature and pressure fluctuate', NVT: 'N, V, T: energy and pressure fluctuate', NPT: 'N, P, T: energy and density fluctuate' }[st.ens];
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) step(); draw(); if (frame++ % 8 === 0) readouts(); };
  restart(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
