// Barostat lab: 2D Lennard-Jones fluid with a Langevin thermostat and a choice of barostat:
// none (NVT), Berendsen, C-rescale (stochastic cell rescaling), or a Parrinello–Rahman-style piston with mass.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RC = 2.5, DT = 0.005, SUB = 4, HIST = 450, N = 100, T0 = 1.0, GAMMA = 1.0, BETA = 0.3;
const TAUS = [0.1, 0.2, 0.5, 1, 2, 5];
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountBarostats(root) {
  const st = { b: 'crescale', P0: 1.0, tau: 0.5, running: true };
  let L, x, v, f, Wv = 0, eta = 0, t = 0, hist = [], statFrom = 0;
  function forces() {
    f.fill(0); Wv = 0; const rc2 = RC * RC;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[2 * i] - x[2 * j], dy = x[2 * i + 1] - x[2 * j + 1];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r2 = dx * dx + dy * dy; if (r2 >= rc2) continue;
      const i6 = 1 / (r2 * r2 * r2), rF = 24 * (2 * i6 * i6 - i6), ff = rF / r2; Wv += rF;
      f[2 * i] += ff * dx; f[2 * i + 1] += ff * dy; f[2 * j] -= ff * dx; f[2 * j + 1] -= ff * dy;
    }
  }
  const kin = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return 0.5 * s; };
  const pressure = () => (kin() + Wv / 2) / (L * L);          // 2D virial pressure (K = N kT for Nf = 2N)
  function restart() {
    L = Math.sqrt(N / 0.6); x = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = 10, a = L / n;
    for (let i = 0; i < N; i++) { x[2 * i] = (i % n + 0.5) * a; x[2 * i + 1] = (Math.floor(i / n) + 0.5) * a; v[2 * i] = gauss(); v[2 * i + 1] = gauss(); }
    forces(); eta = 0; t = 0; hist = []; statFrom = 3;
  }
  function scaleBox(de) {                                       // de = change of ln(area)
    de = Math.max(-0.02, Math.min(0.02, de)); const mu = Math.exp(de / 2);
    L *= mu; for (let i = 0; i < 2 * N; i++) x[i] *= mu;
    if (st.b !== 'berendsen') for (let i = 0; i < 2 * N; i++) v[i] /= mu;
  }
  function step() {
    const c1 = Math.exp(-GAMMA * DT), c2 = Math.sqrt((1 - c1 * c1) * T0);
    for (let s = 0; s < SUB; s++) {
      for (let i = 0; i < 2 * N; i++) { v[i] += 0.5 * DT * f[i]; x[i] += 0.5 * DT * v[i]; }
      for (let i = 0; i < 2 * N; i++) v[i] = c1 * v[i] + c2 * gauss();
      for (let i = 0; i < 2 * N; i++) { x[i] += 0.5 * DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];
      const P = pressure(), A = L * L;
      if (st.b === 'berendsen') scaleBox(-BETA / st.tau * (st.P0 - P) * DT);
      else if (st.b === 'crescale') scaleBox(-BETA / st.tau * (st.P0 - P - T0 / A) * DT + Math.sqrt(2 * T0 * BETA / (A * st.tau) * DT) * gauss());
      else if (st.b === 'piston') { const W = A * st.tau * st.tau / (4 * Math.PI * Math.PI * BETA); eta += DT * A * (P - st.P0) / W; scaleBox(eta * DT); }
      t += DT;
    }
    hist.push({ t, P: pressure(), A: L * L }); if (hist.length > HIST) hist.shift();
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view en-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-b="none">None</button><button data-b="berendsen">Berendsen</button><button data-b="crescale" class="on">C-rescale</button><button data-b="piston">Piston</button></div>
    <label data-bo>Target pressure <output data-o="P0"></output><input type="range" data-k="P0" min="0.2" max="3" step="0.1" value="1"></label>
    <label data-bo><span>Coupling time <i class="nc">τ</i><sub>P</sub></span><output data-o="tau"></output><input type="range" data-k="tau" min="0" max="5" step="1" value="2"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button></div>
    <div class="readouts">
      <div><span>Average pressure</span><b data-o="P"></b></div>
      <div><span>Spread of instant P</span><b data-o="sP"></b></div>
      <div><span>Average density</span><b data-o="rho"></b></div>
      <div><span>Spread of box area</span><b data-o="sA"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.en-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const resetStats = w => { statFrom = t + w; };
  root.querySelectorAll('[data-b]').forEach(bt => bt.onclick = () => {
    st.b = bt.dataset.b; root.querySelectorAll('[data-b]').forEach(z => z.classList.toggle('on', z === bt));
    root.querySelectorAll('[data-bo]').forEach(e => e.hidden = st.b === 'none'); eta = 0; resetStats(2 + 3 * st.tau);
  });
  $('[data-k=P0]').oninput = e => { st.P0 = +e.target.value; resetStats(2 + 3 * st.tau); };
  $('[data-k=tau]').oninput = e => { st.tau = TAUS[+e.target.value]; eta = 0; resetStats(2 + 3 * st.tau); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = restart;

  let W = 600, H = 660;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function panel(y0, y1, key, color, name, target) {
    const x0 = 50, x1 = W - 12, vals = hist.map(q => q[key]);
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = color; cx.fillText(name, x0 + 6, y0 - 4);
    if (vals.length < 2) return;
    let lo = Math.min(...vals, target ?? Infinity), hi = Math.max(...vals, target ?? -Infinity); const pad = (hi - lo) * 0.1 + 1e-3; lo -= pad; hi += pad;
    const Y = q => y1 - (q - lo) / (hi - lo) * (y1 - y0), X = k => x0 + (x1 - x0) * k / (HIST - 1);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText(hi.toFixed(1), x0 - 5, y0 + 9); cx.fillText(lo.toFixed(1), x0 - 5, y1);
    if (target != null) { cx.strokeStyle = COL.ink; cx.setLineDash([5, 5]); cx.beginPath(); cx.moveTo(x0, Y(target)); cx.lineTo(x1, Y(target)); cx.stroke(); cx.setLineDash([]); }
    cx.strokeStyle = color; cx.lineWidth = 1.5; cx.beginPath(); vals.forEach((q, k) => (k ? cx.lineTo(X(k), Y(q)) : cx.moveTo(X(k), Y(q)))); cx.stroke();
  }
  function draw() {
    cx.clearRect(0, 0, W, H);
    const boxH = H * 0.44, S = boxH / Math.sqrt(N / 0.25), bs = L * S, ox = W / 2 - bs / 2, oy = 8 + (boxH - bs) / 2;
    cx.strokeStyle = st.b === 'none' ? 'rgba(143,163,199,.7)' : COL.amber; cx.lineWidth = st.b === 'none' ? 1.5 : 3; cx.strokeRect(ox, oy, bs, bs);
    for (let i = 0; i < N; i++) { const e = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / (4 * T0); cx.fillStyle = e > 1 ? COL.coral : e > 0.45 ? COL.amber : COL.cyan; cx.beginPath(); cx.arc(ox + x[2 * i] * S, oy + x[2 * i + 1] * S, Math.max(2, 0.42 * S), 0, 7); cx.fill(); }
    const top = boxH + 34, h = (H - top - 8) / 2;
    panel(top, top + h - 22, 'P', COL.amber, 'instantaneous pressure' + (st.b === 'none' ? '' : '  (dashed: target)'), st.b === 'none' ? null : st.P0);
    panel(top + h + 14, H - 10, 'A', COL.cyan, 'box area', null);
  }
  function readouts() {
    o.P0.textContent = st.P0.toFixed(1); o.tau.textContent = st.tau;
    const s = hist.filter(q => q.t >= statFrom), ok = s.length > 80;
    if (!ok) { o.P.textContent = o.sP.textContent = o.rho.textContent = o.sA.textContent = 'settling…'; return; }
    const m = k => s.reduce((a, q) => a + q[k], 0) / s.length, sd = k => { const mm = m(k); return Math.sqrt(s.reduce((a, q) => a + (q[k] - mm) ** 2, 0) / s.length); };
    o.P.textContent = m('P').toFixed(3); o.sP.textContent = '± ' + sd('P').toFixed(2);
    o.rho.textContent = (s.reduce((a, q) => a + N / q.A, 0) / s.length).toFixed(3); o.sA.textContent = st.b === 'none' ? '0 (fixed box)' : '± ' + sd('A').toFixed(2);
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) step(); draw(); if (frame++ % 8 === 0) readouts(); };
  restart(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
