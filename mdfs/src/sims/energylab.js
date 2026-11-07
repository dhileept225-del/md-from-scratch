// Energy-check lab: 2D Lennard-Jones liquid in NVE, with switches that break energy conservation.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const NS = 8, N = NS * NS, L = 9.6, T0 = 0.6, HIST = 300;

export function mountEnergyLab(root) {
  const st = { dt: 0.005, rc: 2.5, shift: true, bug: false, thermo: false, running: true, steps: 0, hist: [], E0: null };
  const x = new Float64Array(N * 2), v = new Float64Array(N * 2), f = new Float64Array(N * 2);
  let U = 0;

  function potential(pos, wantForce) {
    const rc2 = st.rc * st.rc, vc = st.shift ? 4 * (st.rc ** -12 - st.rc ** -6) : 0; let u = 0;
    if (wantForce) f.fill(0);
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = pos[2 * i] - pos[2 * j], dy = pos[2 * i + 1] - pos[2 * j + 1];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r2 = dx * dx + dy * dy; if (r2 > rc2) continue;
      const i2 = 1 / r2, i6 = i2 * i2 * i2;
      u += 4 * (i6 * i6 - i6) - vc;
      if (wantForce) {
        const ff = st.bug ? 24 * i2 * (i6 * i6 - i6) : 24 * i2 * (2 * i6 * i6 - i6);   // bug: a missing factor of 2
        f[2 * i] += ff * dx; f[2 * i + 1] += ff * dy; f[2 * j] -= ff * dx; f[2 * j + 1] -= ff * dy;
      }
    }
    return u;
  }
  const K = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return 0.5 * s; };
  function reset() {
    const a = L / NS; let k = 0;
    for (let i = 0; i < NS; i++) for (let j = 0; j < NS; j++, k++) { x[2 * k] = (i + 0.5) * a; x[2 * k + 1] = (j + 0.5) * a; }
    for (let i = 0; i < 2 * N; i++) v[i] = (Math.random() * 2 - 1);
    for (let d = 0; d < 2; d++) { let m = 0; for (let i = 0; i < N; i++) m += v[2 * i + d]; m /= N; for (let i = 0; i < N; i++) v[2 * i + d] -= m; }
    const s = Math.sqrt(T0 * (2 * N - 2) / (2 * K())); for (let i = 0; i < 2 * N; i++) v[i] *= s;
    U = potential(x, true); st.steps = 0; st.hist = []; st.E0 = null; msg.hidden = true; blown = false;
  }
  let blown = false;
  function step() {
    const dt = st.dt, h = 0.5 * dt;
    for (let i = 0; i < 2 * N; i++) { v[i] += h * f[i]; x[i] += dt * v[i]; x[i] -= L * Math.floor(x[i] / L); }
    U = potential(x, true);
    for (let i = 0; i < 2 * N; i++) v[i] += h * f[i];
    if (st.thermo) { const T = 2 * K() / (2 * N - 2), lam = Math.sqrt(1 + (dt / 0.2) * (T0 / Math.max(T, 1e-6) - 1)); for (let i = 0; i < 2 * N; i++) v[i] *= lam; }
    st.steps++;
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view lab-view"><canvas></canvas><div class="sim-msg" hidden></div></div>
  <div class="sim-ctrl">
    <label>Time step Δt <output data-o="dt"></output><input type="range" data-k="dt" min="0.001" max="0.02" step="0.001" value="0.005"></label>
    <label>Cutoff rc <output data-o="rc"></output><input type="range" data-k="rc" min="1.3" max="3" step="0.05" value="2.5"></label>
    <label class="check"><input type="checkbox" data-k="shift" checked> Shift potential so V(rc) = 0</label>
    <label class="check"><input type="checkbox" data-k="bug"> <span style="color:${COL.coral}">Bug:</span> force missing a factor of 2</label>
    <label class="check"><input type="checkbox" data-k="thermo"> Thermostat on (not NVE)</label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">↺ Reset</button><button class="btn sm ghost" data-k="fc">Force check</button></div>
    <div class="readouts">
      <div><span>Drift / atom / time</span><b data-o="drift">–</b></div>
      <div><span>Fluctuation ratio</span><b data-o="ratio">–</b></div>
      <div class="wide"><span>Verdict</span><b data-o="verdict">–</b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.lab-view'), msg = root.querySelector('.sim-msg');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const restartStats = () => { st.hist = []; st.E0 = null; U = potential(x, true); if (!blown) msg.hidden = true; };
  $('[data-k=dt]').oninput = e => { st.dt = +e.target.value; restartStats(); };
  $('[data-k=rc]').oninput = e => { st.rc = +e.target.value; restartStats(); };
  $('[data-k=shift]').onchange = e => { st.shift = e.target.checked; restartStats(); };
  $('[data-k=bug]').onchange = e => { st.bug = e.target.checked; restartStats(); };
  $('[data-k=thermo]').onchange = e => { st.thermo = e.target.checked; restartStats(); };
  const runBtn = $('[data-k=run]');
  runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = reset;
  $('[data-k=fc]').onclick = () => {   // finite-difference check on atom 0, x direction
    const hh = 1e-5, p = Float64Array.from(x);
    potential(x, true); const Fa = f[0];
    p[0] = x[0] + hh; const Up = potential(p, false); p[0] = x[0] - hh; const Um = potential(p, false);
    const Fn = -(Up - Um) / (2 * hh), rel = Math.abs(Fa - Fn) / Math.max(1e-8, Math.abs(Fn)), ok = rel < 1e-4;
    U = potential(x, true);
    msg.hidden = false;
    msg.innerHTML = `<b>Force check on atom 1 (x):</b> analytic F = ${Fa.toFixed(5)}, numeric −dU/dx = ${Fn.toFixed(5)}<br>Relative error ${rel.toExponential(1)}: <b style="color:${ok ? COL.cyan : COL.coral}">${ok ? 'PASS, forces match the energy' : 'FAIL, forces do not match the energy'}</b>`;
  };

  let W = 600, H = 460;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  // legend sits ABOVE the plot box, so it never overlaps the data
  function plot(x0, y0, w, h, series, legend) {
    cx.font = '500 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.textBaseline = 'alphabetic';
    let lx = x0;
    legend.forEach(([txt, c]) => { const sw = c !== COL.muted; if (sw) { cx.fillStyle = c; cx.fillRect(lx, y0 - 11, 14, 3); } cx.fillStyle = c; cx.fillText(txt, lx + (sw ? 18 : 0), y0 - 6); lx += (sw ? 26 : 8) + cx.measureText(txt).width; });
    cx.strokeStyle = 'rgba(143,163,199,.45)'; cx.lineWidth = 1; cx.strokeRect(x0, y0, w, h);
    let lo = Infinity, hi = -Infinity; series.forEach(sr => sr.data.forEach(q => { lo = Math.min(lo, q); hi = Math.max(hi, q); }));
    if (!isFinite(lo)) return; const pad = Math.max((hi - lo) * 0.15, series.minSpan || 1e-3); lo -= pad; hi += pad;
    cx.save(); cx.beginPath(); cx.rect(x0 + 1, y0 + 1, w - 2, h - 2); cx.clip();
    series.forEach(sr => {
      cx.strokeStyle = sr.color; cx.lineWidth = sr.w || 2; cx.beginPath();
      sr.data.forEach((q, i) => { const X = x0 + 4 + i / (HIST - 1) * (w - 8), Yv = y0 + 4 + (h - 8) - (q - lo) / (hi - lo) * (h - 8); i ? cx.lineTo(X, Yv) : cx.moveTo(X, Yv); }); cx.stroke();
    });
    cx.restore();
  }
  function draw() {
    cx.clearRect(0, 0, W, H);
    const wide = W > 520, bs = wide ? Math.min(H - 20, W * 0.48) : Math.min(W - 20, H * 0.5), bx = 10, by = 10, s = bs / L;
    cx.strokeStyle = 'rgba(56,225,255,.4)'; cx.lineWidth = 1.5; cx.strokeRect(bx, by, bs, bs);
    cx.save(); cx.beginPath(); cx.rect(bx, by, bs, bs); cx.clip();
    for (let i = 0; i < N; i++) {
      const sp = Math.hypot(v[2 * i], v[2 * i + 1]) / 2.2;
      cx.fillStyle = sp < 0.5 ? COL.cyan : sp < 1 ? COL.amber : COL.coral;
      cx.beginPath(); cx.arc(bx + x[2 * i] * s, by + x[2 * i + 1] * s, Math.max(3, 0.48 * s), 0, Math.PI * 2); cx.fill();
    }
    cx.restore();
    const px = wide ? bx + bs + 16 : 10, py = wide ? 26 : by + bs + 30, pw = wide ? W - px - 10 : W - 20, ph = wide ? (H - 26 - 36 - 10) / 2 : (H - py - 36 - 10) / 2;
    const hK = st.hist.map(q => q[0]), hU = st.hist.map(q => q[1]), hE = st.hist.map(q => q[0] + q[1]);
    plot(px, py, pw, ph, [{ data: hK, color: COL.cyan }, { data: hU, color: COL.coral }, { data: hE, color: COL.amber, w: 3 }], [['kinetic', COL.cyan], ['potential', COL.coral], ['total', COL.amber], ['(per atom)', COL.muted]]);
    const dev = Object.assign([{ data: hE.map(e => e - (st.E0 ?? e)), color: COL.amber, w: 2.5 }], { minSpan: 2e-3 });
    plot(px, py + ph + 36, pw, ph, dev, [['total energy change, zoomed in', COL.amber]]);
  }
  function stats() {
    o.dt.textContent = st.dt.toFixed(3); o.rc.textContent = st.rc.toFixed(2) + ' σ';
    const n = st.hist.length; if (n < 30) { o.drift.textContent = o.ratio.textContent = '…'; o.verdict.textContent = 'collecting data…'; return; }
    const E = st.hist.map(q => q[0] + q[1]), Kh = st.hist.map(q => q[0]), tt = st.hist.map((_, i) => i * 10 * st.dt);
    const mean = a => a.reduce((s, q) => s + q, 0) / a.length, sd = a => { const m = mean(a); return Math.sqrt(mean(a.map(q => (q - m) ** 2))); };
    const mt = mean(tt), mE = mean(E); let num = 0, den = 0; tt.forEach((t, i) => { num += (t - mt) * (E[i] - mE); den += (t - mt) ** 2; });
    const drift = num / den, ratio = sd(E) / Math.max(1e-9, sd(Kh));
    o.drift.textContent = drift.toExponential(1); o.ratio.textContent = ratio.toFixed(3);
    o.verdict.innerHTML = st.thermo ? `<span style="color:${COL.amber}">Thermostat on: energy is not meant to be conserved</span>`
      : ratio < 0.02 && Math.abs(drift) < 2e-3 ? `<span style="color:${COL.cyan}">Healthy: energy conserved</span>`
      : `<span style="color:${COL.coral}">Not conserved: check Δt, cutoff, forces</span>`;
  }

  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  function loop() {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running && !blown) {
      for (let k = 0; k < 10; k++) step();
      const k = K() / N, u = U / N;
      if (!isFinite(k + u) || k > 50) { blown = true; msg.hidden = false; msg.innerHTML = '<b>The simulation blew up.</b> Lower Δt and press Reset.'; }
      else { if (st.E0 === null) st.E0 = k + u; st.hist.push([k, u]); if (st.hist.length > HIST) st.hist.shift(); }
    }
    draw(); if (frame++ % 10 === 0) stats();
  }
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
