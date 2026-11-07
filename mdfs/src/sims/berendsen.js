// Thermostat lab: 2D Lennard-Jones fluid with no thermostat, exact rescaling, Berendsen, or Bussi stochastic rescaling.
// Shows how fast T reaches the target and whether the fluctuations have the canonical size.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RHO = 0.7, RC = 2.5, DT = 0.005, SUB = 4, HIST = 500, N0 = 100;
const TAUS = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5];
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());
function heat(x) {
  x = Math.min(1, Math.max(0, x));
  const S = [[56, 225, 255], [255, 194, 75], [255, 107, 107]], k = x < 0.5 ? 0 : 1, f = x < 0.5 ? x * 2 : x * 2 - 1;
  return `rgb(${S[k].map((c, i) => Math.round(c + (S[k + 1][i] - c) * f)).join(',')})`;
}

export function mountBerendsen(root) {
  const st = { th: 'berendsen', tau: 0.5, T0: 1.0, running: true };
  const N = N0, Nf = 2 * N - 2;
  let Lx, Ly, x, v, f, hist = [], t = 0, change = null, statFrom = 0;
  function forces() {
    f.fill(0); const rc2 = RC * RC;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[2 * i] - x[2 * j], dy = x[2 * i + 1] - x[2 * j + 1];
      dx -= Lx * Math.round(dx / Lx); dy -= Ly * Math.round(dy / Ly);
      const r2 = dx * dx + dy * dy; if (r2 >= rc2) continue;
      const i6 = 1 / (r2 * r2 * r2), ff = 24 * (2 * i6 * i6 - i6) / r2;
      f[2 * i] += ff * dx; f[2 * i + 1] += ff * dy; f[2 * j] -= ff * dx; f[2 * j + 1] -= ff * dy;
    }
  }
  const kin = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return 0.5 * s; };
  function lambda(K) {
    const K0 = 0.5 * Nf * st.T0;
    if (st.th === 'rescale') return Math.sqrt(K0 / K);
    if (st.th === 'berendsen') return Math.sqrt(Math.max(0, 1 + DT / st.tau * (K0 / K - 1)));
    if (st.th === 'bussi') {
      const c = Math.exp(-DT / st.tau), R1 = gauss(); let S = 0; for (let k = 0; k < Nf - 1; k++) { const g = gauss(); S += g * g; }
      const a2 = c + (1 - c) * K0 / (Nf * K) * (R1 * R1 + S) + 2 * R1 * Math.sqrt(c * (1 - c) * K0 / (Nf * K));
      return Math.sign(R1 + Math.sqrt(c * Nf * K / ((1 - c) * K0))) * Math.sqrt(a2);
    }
    return 1;
  }
  function restart() {
    x = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = 10, a = Math.sqrt(2 / (Math.sqrt(3) * RHO)), h = a * Math.sqrt(3) / 2; Lx = n * a; Ly = n * h;
    for (let i = 0; i < N; i++) { const r = Math.floor(i / n); x[2 * i] = (i % n + 0.25 + 0.5 * (r % 2)) * a; x[2 * i + 1] = (r + 0.5) * h; v[2 * i] = gauss(); v[2 * i + 1] = gauss(); }
    let px = 0, py = 0; for (let i = 0; i < N; i++) { px += v[2 * i]; py += v[2 * i + 1]; }
    for (let i = 0; i < N; i++) { v[2 * i] -= px / N; v[2 * i + 1] -= py / N; }
    const s = Math.sqrt(st.T0 / (2 * kin() / Nf)); for (let i = 0; i < 2 * N; i++) v[i] *= s;
    forces(); hist = []; t = 0; change = null; statFrom = 0;
  }
  function step() {
    for (let s = 0; s < SUB; s++) {
      for (let i = 0; i < 2 * N; i++) { const L = i % 2 ? Ly : Lx; v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];
      const lam = lambda(kin()); if (lam !== 1) for (let i = 0; i < 2 * N; i++) v[i] *= lam;
      t += DT;
    }
    const T = 2 * kin() / Nf; hist.push({ t, T, T0: st.T0 }); if (hist.length > HIST) hist.shift();
    if (change && change.t63 == null && hist.length >= 5) {
      const sm = hist.slice(-5).reduce((s2, q) => s2 + q.T, 0) / 5;
      if (Math.abs(sm - change.Tstart) >= 0.63 * Math.abs(st.T0 - change.Tstart)) change.t63 = t - change.t;
    }
  }
  const settle = () => { statFrom = t + (st.th === 'none' ? 1 : Math.min(3 * st.tau, 15) + 1); };

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-t="none">None</button><button data-t="rescale">Rescale</button><button data-t="berendsen" class="on">Berendsen</button><button data-t="bussi">Bussi</button></div>
    <label>Target temperature <output data-o="T0"></output><input type="range" data-k="T0" min="0.4" max="2" step="0.1" value="1"></label>
    <label data-tau><span>Coupling time <i class="nc">τ</i></span><output data-o="tau"></output><input type="range" data-k="tau" min="0" max="7" step="1" value="4"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button></div>
    <div class="readouts">
      <div><span>Average T</span><b data-o="avg"></b></div>
      <div><span>63% time after a jump</span><b data-o="t63"></b></div>
      <div><span>Spread of T (measured)</span><b data-o="sd"></b></div>
      <div><span>Spread if canonical</span><b data-o="sdc"></b></div>
      <div class="wide"><span>Measured ÷ canonical spread</span><b data-o="ratio"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => {
    st.th = b.dataset.t; root.querySelectorAll('[data-t]').forEach(z => z.classList.toggle('on', z === b));
    $('[data-tau]').hidden = !(st.th === 'berendsen' || st.th === 'bussi'); change = null; settle();
  });
  $('[data-k=T0]').oninput = e => {
    st.T0 = +e.target.value; const T = hist.length ? hist[hist.length - 1].T : st.T0;
    change = st.th === 'none' ? null : { t, Tstart: T, t63: null }; settle();
  };
  $('[data-k=tau]').oninput = e => { st.tau = TAUS[+e.target.value]; change = null; settle(); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = () => { restart(); settle(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function draw() {
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(160, H * 0.38), S = Math.min((W - 24) / Lx, (H - plotH - 40) / Ly), ox = (W - Lx * S) / 2, oy = 12;
    cx.strokeStyle = st.th === 'none' ? 'rgba(143,163,199,.7)' : COL.amber; cx.lineWidth = st.th === 'none' ? 1.5 : 3; cx.strokeRect(ox, oy, Lx * S, Ly * S);
    const r = Math.max(2.5, 0.42 * S);
    for (let i = 0; i < N; i++) { const e = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / 2; cx.fillStyle = heat(e / (2.5 * Math.max(0.4, st.T0))); cx.beginPath(); cx.arc(ox + x[2 * i] * S, oy + x[2 * i + 1] * S, r, 0, 7); cx.fill(); }
    const px0 = 46, px1 = W - 14, py1 = H - 22, py0 = oy + Ly * S + 44;
    if (hist.length < 2) return;
    const hi = Math.max(2.2, ...hist.map(q => q.T)) * 1.05, Y = T => py1 - T / hi * (py1 - py0), X = k => px0 + (px1 - px0) * k / (HIST - 1);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(px0, py0); cx.lineTo(px0, py1); cx.lineTo(px1, py1); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right';
    [0, 1, 2].forEach(T => cx.fillText(T.toFixed(1), px0 - 6, Y(T) + 4));
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.coral; cx.fillText('instant T', px0 + 8, py0 - 6); cx.fillStyle = COL.muted; cx.fillText('target', px0 + 80, py0 - 6);
    if (st.th !== 'none') { cx.strokeStyle = COL.muted; cx.setLineDash([6, 5]); cx.lineWidth = 1.5; cx.beginPath(); hist.forEach((q, k) => { k ? cx.lineTo(X(k), Y(q.T0)) : cx.moveTo(X(k), Y(q.T0)); }); cx.stroke(); cx.setLineDash([]); }
    cx.strokeStyle = COL.coral; cx.lineWidth = 1.6; cx.beginPath(); hist.forEach((q, k) => k ? cx.lineTo(X(k), Y(q.T)) : cx.moveTo(X(k), Y(q.T))); cx.stroke();
  }
  function readouts() {
    o.T0.textContent = st.T0.toFixed(1); o.tau.textContent = st.tau;
    const s = hist.filter(q => q.t >= statFrom), ok = s.length > 60;
    const m = ok ? s.reduce((a, q) => a + q.T, 0) / s.length : 0, sd = ok ? Math.sqrt(s.reduce((a, q) => a + (q.T - m) ** 2, 0) / s.length) : 0;
    const sdc = Math.sqrt(2 / Nf) * (ok ? m : st.T0);
    o.avg.textContent = ok ? m.toFixed(3) : 'settling…'; o.sd.textContent = ok ? '± ' + sd.toFixed(3) : 'settling…'; o.sdc.textContent = '± ' + sdc.toFixed(3);
    o.ratio.innerHTML = !ok ? 'settling…' : `<span style="color:${sd / sdc < 0.75 ? COL.coral : COL.cyan}">${(sd / sdc).toFixed(2)}${sd / sdc < 0.75 ? ' (too small)' : ' (about right)'}</span>`;
    o.t63.textContent = st.th === 'none' ? '— (no thermostat)' : !change ? 'move the target' : change.t63 == null ? 'waiting…' : change.t63.toFixed(2);
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) step(); draw(); if (frame++ % 8 === 0) readouts(); };
  restart(); settle(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
