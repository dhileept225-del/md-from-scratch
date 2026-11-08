// Constraints lab: a small cluster of SPC water molecules in the plane (nm, ps, amu, kJ/mol), run NVE with
// flexible bonds or with rigid SHAKE/RATTLE constraints. Change the time step and the SHAKE tolerance and watch
// the energy, the number of SHAKE sweeps and the bond-length error.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const COUL = 138.935, QO = -0.82, QH = 0.41, SIG = 0.3166, EPS = 0.650, R0 = 0.1, TH0 = 109.47 * Math.PI / 180, KB = 345000, KT = 383;
const DHH = 2 * R0 * Math.sin(TH0 / 2), NW = 4, N = 3 * NW, RW = 0.32, KW = 4000;
const MASS = Array.from({ length: N }, (_, i) => (i % 3 ? 1.008 : 15.999)), Q = Array.from({ length: N }, (_, i) => (i % 3 ? QH : QO));
const DTS = [0.5, 1, 2, 3, 4], TOLS = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], FS_PER_FRAME = 1.2, HIST = 400;
const CONS = []; for (let w = 0; w < NW; w++) { const o = 3 * w; CONS.push([o, o + 1, R0], [o, o + 2, R0], [o + 1, o + 2, DHH]); }
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountShake(root) {
  const st = { rigid: true, dt: 2, tol: 1e-8, running: true };
  let x, v, F, U = 0, E0 = null, hist = [], acc = 0, tfs = 0, sweeps = 0, lenErr = 0, blown = false;
  function forces() {
    F = Array.from({ length: N }, () => [0, 0]); U = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      if (Math.floor(i / 3) === Math.floor(j / 3)) continue;                   // intramolecular non-bonded: excluded
      const dx = x[i][0] - x[j][0], dy = x[i][1] - x[j][1], r2 = dx * dx + dy * dy, r = Math.sqrt(r2);
      let fr = COUL * Q[i] * Q[j] / (r2 * r); U += COUL * Q[i] * Q[j] / r;
      if (i % 3 === 0 && j % 3 === 0) { const s6 = (SIG * SIG / r2) ** 3; U += 4 * EPS * (s6 * s6 - s6); fr += 24 * EPS * (2 * s6 * s6 - s6) / r2; }
      F[i][0] += fr * dx; F[i][1] += fr * dy; F[j][0] -= fr * dx; F[j][1] -= fr * dy;
    }
    for (let w = 0; w < NW; w++) {                                            // soft circular wall keeps the cluster together
      const o = 3 * w, r = Math.hypot(x[o][0], x[o][1]);
      if (r > RW) { U += 0.5 * KW * (r - RW) ** 2; const f = -KW * (r - RW) / r; F[o][0] += f * x[o][0]; F[o][1] += f * x[o][1]; }
    }
    if (!st.rigid) for (let w = 0; w < NW; w++) {
      const o = 3 * w;
      for (const h of [o + 1, o + 2]) { const dx = x[h][0] - x[o][0], dy = x[h][1] - x[o][1], r = Math.hypot(dx, dy), f = -KB * (r - R0) / r; U += 0.5 * KB * (r - R0) ** 2; F[h][0] += f * dx; F[h][1] += f * dy; F[o][0] -= f * dx; F[o][1] -= f * dy; }
      const ax = x[o + 1][0] - x[o][0], ay = x[o + 1][1] - x[o][1], bx = x[o + 2][0] - x[o][0], by = x[o + 2][1] - x[o][1], la = Math.hypot(ax, ay), lb = Math.hypot(bx, by);
      const c = Math.max(-0.9999, Math.min(0.9999, (ax * bx + ay * by) / (la * lb))), th = Math.acos(c), g = KT * (th - TH0) / Math.sin(th);
      const fax = g * (bx / (la * lb) - c * ax / (la * la)), fay = g * (by / (la * lb) - c * ay / (la * la)), fbx = g * (ax / (la * lb) - c * bx / (lb * lb)), fby = g * (ay / (la * lb) - c * by / (lb * lb));
      U += 0.5 * KT * (th - TH0) ** 2; F[o + 1][0] += fax; F[o + 1][1] += fay; F[o + 2][0] += fbx; F[o + 2][1] += fby; F[o][0] -= fax + fbx; F[o][1] -= fay + fby;
    }
  }
  function shake(xn, xo) {
    for (let it = 1; it <= 1000; it++) {
      let worst = 0;
      for (const [i, j, d] of CONS) {
        const sx = xn[i][0] - xn[j][0], sy = xn[i][1] - xn[j][1], diff = sx * sx + sy * sy - d * d; worst = Math.max(worst, Math.abs(diff) / (d * d));
        const rx = xo[i][0] - xo[j][0], ry = xo[i][1] - xo[j][1], g = diff / (2 * (1 / MASS[i] + 1 / MASS[j]) * (sx * rx + sy * ry));
        xn[i][0] -= g * rx / MASS[i]; xn[i][1] -= g * ry / MASS[i]; xn[j][0] += g * rx / MASS[j]; xn[j][1] += g * ry / MASS[j];
      }
      if (worst < st.tol) return it;
    }
    return 1000;
  }
  function rattle() {
    for (let it = 0; it < 1000; it++) {
      let worst = 0;
      for (const [i, j, d] of CONS) {
        const rx = x[i][0] - x[j][0], ry = x[i][1] - x[j][1], dot = rx * (v[i][0] - v[j][0]) + ry * (v[i][1] - v[j][1]), k = dot / ((1 / MASS[i] + 1 / MASS[j]) * d * d);
        worst = Math.max(worst, Math.abs(dot)); v[i][0] -= k * rx / MASS[i]; v[i][1] -= k * ry / MASS[i]; v[j][0] += k * rx / MASS[j]; v[j][1] += k * ry / MASS[j];
      }
      if (worst < 1e-10) return;
    }
  }
  const kinetic = () => v.reduce((s, q, i) => s + 0.5 * MASS[i] * (q[0] ** 2 + q[1] ** 2), 0);
  function restart() {
    x = []; const centres = [[-0.15, -0.14], [0.15, -0.14], [-0.15, 0.14], [0.15, 0.14]], angs = [0.3, 2.4, -1.2, 4.0];
    centres.forEach((c, w) => { const a = angs[w]; x.push([c[0], c[1]], [c[0] + R0 * Math.cos(a - TH0 / 2), c[1] + R0 * Math.sin(a - TH0 / 2)], [c[0] + R0 * Math.cos(a + TH0 / 2), c[1] + R0 * Math.sin(a + TH0 / 2)]); });
    v = x.map((_, i) => { const s = Math.sqrt(2.494 / MASS[i]); return [gauss() * s, gauss() * s]; });
    const P = [0, 1].map(d => v.reduce((s, q, i) => s + MASS[i] * q[d], 0) / MASS.reduce((s, q) => s + q, 0)); v = v.map(q => [q[0] - P[0], q[1] - P[1]]);
    if (st.rigid) rattle(); forces(); E0 = null; hist = []; tfs = 0; blown = false;
  }
  function step() {
    const dt = st.dt / 1000;
    for (let i = 0; i < N; i++) { v[i][0] += 0.5 * dt * F[i][0] / MASS[i]; v[i][1] += 0.5 * dt * F[i][1] / MASS[i]; }
    const xn = x.map((q, i) => [q[0] + dt * v[i][0], q[1] + dt * v[i][1]]);
    if (st.rigid) { sweeps = shake(xn, x); v = xn.map((q, i) => [(q[0] - x[i][0]) / dt, (q[1] - x[i][1]) / dt]); }
    x = xn; forces();
    for (let i = 0; i < N; i++) { v[i][0] += 0.5 * dt * F[i][0] / MASS[i]; v[i][1] += 0.5 * dt * F[i][1] / MASS[i]; }
    if (st.rigid) rattle();
    lenErr = CONS.reduce((m, [i, j, d]) => Math.max(m, Math.abs(Math.hypot(x[i][0] - x[j][0], x[i][1] - x[j][1]) - d) / d), 0);
    tfs += st.dt; const E = U + kinetic(); if (E0 === null) E0 = E;
    if (!isFinite(E) || Math.abs(E - E0) > 500) { blown = true; return; }
    hist.push({ t: tfs, dE: E - E0 }); if (hist.length > HIST) hist.shift();
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="flex">Flexible bonds</button><button data-m="rigid" class="on">Rigid (SHAKE/RATTLE)</button></div>
    <label><span>Time step</span><output data-o="dt"></output><input type="range" data-k="dt" min="0" max="4" step="1" value="2"></label>
    <label data-tol><span>SHAKE tolerance</span><output data-o="tol"></output><input type="range" data-k="tol" min="0" max="4" step="1" value="3"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button></div>
    <div class="readouts">
      <div><span>Energy spread</span><b data-o="sp"></b></div>
      <div><span>Energy drift</span><b data-o="dr"></b></div>
      <div><span>SHAKE sweeps per step</span><b data-o="sw"></b></div>
      <div><span>Largest bond-length error</span><b data-o="le"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.rigid = b.dataset.m === 'rigid'; root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); $('[data-tol]').hidden = !st.rigid; restart(); });
  $('[data-k=dt]').oninput = e => { st.dt = DTS[+e.target.value]; restart(); };
  $('[data-k=tol]').oninput = e => { st.tol = TOLS[+e.target.value]; restart(); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = restart;

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  function draw() {
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(150, H * 0.32), areaH = H - plotH - 20;
    const cmx = x.reduce((s, q, i) => s + MASS[i] * q[0], 0) / MASS.reduce((s, q) => s + q, 0), cmy = x.reduce((s, q, i) => s + MASS[i] * q[1], 0) / MASS.reduce((s, q) => s + q, 0);
    const S = Math.min(W, areaH) / 0.8, P = q => [W / 2 + q[0] * S, areaH / 2 + 8 - q[1] * S];
    cx.strokeStyle = 'rgba(143,163,199,.3)'; cx.setLineDash([6, 6]); cx.lineWidth = 1.5; cx.beginPath(); cx.arc(W / 2, areaH / 2 + 8, RW * S, 0, 7); cx.stroke(); cx.setLineDash([]);
    for (let w = 0; w < NW; w++) {
      const Pm = [0, 1, 2].map(k => P(x[3 * w + k]));
      cx.strokeStyle = st.rigid ? COL.ink : COL.muted; cx.lineWidth = st.rigid ? 6 : 4; cx.lineCap = 'round';
      if (st.rigid) { cx.beginPath(); cx.moveTo(Pm[1][0], Pm[1][1]); cx.lineTo(Pm[0][0], Pm[0][1]); cx.lineTo(Pm[2][0], Pm[2][1]); cx.stroke(); cx.strokeStyle = 'rgba(255,194,75,.6)'; cx.lineWidth = 2; cx.setLineDash([5, 4]); cx.beginPath(); cx.moveTo(Pm[1][0], Pm[1][1]); cx.lineTo(Pm[2][0], Pm[2][1]); cx.stroke(); cx.setLineDash([]); }
      else for (const k of [1, 2]) { const ax = Pm[0][0], ay = Pm[0][1], bx = Pm[k][0], by = Pm[k][1], L = Math.hypot(bx - ax, by - ay), nx = -(by - ay) / L, ny = (bx - ax) / L; cx.beginPath(); cx.moveTo(ax, ay); for (let s = 1; s < 10; s++) { const f = s / 10, a_ = (s % 2 ? 1 : -1) * 4; cx.lineTo(ax + (bx - ax) * f + nx * a_, ay + (by - ay) * f + ny * a_); } cx.lineTo(bx, by); cx.stroke(); }
      cx.fillStyle = COL.ink; [1, 2].forEach(k => { cx.beginPath(); cx.arc(Pm[k][0], Pm[k][1], 0.016 * S, 0, 7); cx.fill(); });
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(Pm[0][0], Pm[0][1], 0.026 * S, 0, 7); cx.fill();
    }
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = COL.muted; cx.fillText(`t = ${(tfs / 1000).toFixed(2)} ps`, 12, 22);
    if (blown) { cx.fillStyle = COL.coral; cx.font = '700 18px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText('Unstable: the energy blew up. Smaller time step, or constraints?', W / 2, areaH / 2); }
    const x0 = 46, x1 = W - 12, top = areaH + 18, y1 = H - 12;
    if (hist.length < 2) return;
    const lim = Math.max(0.01, ...hist.map(q => Math.abs(q.dE))) * 1.1, Y = e => (top + y1) / 2 - e / lim * (y1 - top) / 2, X = k => x0 + (x1 - x0) * k / (HIST - 1);
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(x0, top); cx.lineTo(x0, y1); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
    cx.fillStyle = COL.amber; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`total energy change (±${lim.toPrecision(2)} kJ/mol)`, x0 + 6, top - 4);
    cx.strokeStyle = COL.amber; cx.lineWidth = 1.6; cx.beginPath(); hist.forEach((q, k) => (k ? cx.lineTo(X(k), Y(q.dE)) : cx.moveTo(X(k), Y(q.dE)))); cx.stroke();
  }
  function readouts() {
    o.dt.textContent = st.dt + ' fs  (' + (10.4 / st.dt).toFixed(1) + ' steps per O–H vibration)'; o.tol.textContent = st.tol.toExponential(0);
    if (hist.length < 5) { o.sp.textContent = o.dr.textContent = '–'; }
    else { const e = hist.map(q => q.dE), m = e.reduce((s, q) => s + q, 0) / e.length, sd = Math.sqrt(e.reduce((s, q) => s + (q - m) ** 2, 0) / e.length);
      o.sp.innerHTML = `<span style="color:${sd > 0.2 ? COL.coral : COL.cyan}">± ${sd.toPrecision(2)} kJ/mol</span>`; o.dr.textContent = (e[e.length - 1]).toPrecision(2) + ' kJ/mol'; }
    if (blown) o.sp.innerHTML = `<span style="color:${COL.coral}">blew up</span>`;
    o.sw.textContent = st.rigid ? sweeps : '— (no constraints)'; o.le.textContent = lenErr.toExponential(1);
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running && !blown) { acc += FS_PER_FRAME; let n = 0; while (acc >= st.dt && n++ < 40) { step(); acc -= st.dt; if (blown) break; } }
    draw(); if (frame++ % 8 === 0) readouts();
  };
  restart(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
