// Bonds and angles lab: a flexible water molecule in 2D (real units: nm, ps, amu, kJ/mol).
// Drag an atom to stretch or bend, release to watch it vibrate. Bond and angle forces are drawn separately;
// the time step and the bond potential (harmonic / Morse) can be changed.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const M = [15.999, 1.008, 1.008], KB = 345000, R0 = 0.1, KT = 383, T0 = 109.47 * Math.PI / 180, DE = 458;
const AM = Math.sqrt(KB / (2 * DE)), DTS = [0.25, 0.5, 1, 2, 3], FS_PER_FRAME = 0.12, HIST = 400;

export function mountBonds(root) {
  const st = { dt: 0.5, morse: false, drag: -1, running: true };
  let x, v, F, Fb, Fa, U = 0, Ub = 0, Ua = 0, E0 = null, hist = [], acc = 0, tfs = 0;
  function bondTerm(i, j) {
    const dx = x[i][0] - x[j][0], dy = x[i][1] - x[j][1], r = Math.hypot(dx, dy);
    let u, dU;
    if (st.morse) { const e = Math.exp(-AM * (r - R0)); u = DE * (1 - e) ** 2; dU = 2 * DE * AM * e * (1 - e); }
    else { u = 0.5 * KB * (r - R0) ** 2; dU = KB * (r - R0); }
    const fx = -dU * dx / r, fy = -dU * dy / r;
    Fb[i][0] += fx; Fb[i][1] += fy; Fb[j][0] -= fx; Fb[j][1] -= fy; return u;
  }
  function angleTerm(i, j, k) {
    const ax = x[i][0] - x[j][0], ay = x[i][1] - x[j][1], bx = x[k][0] - x[j][0], by = x[k][1] - x[j][1];
    const la = Math.hypot(ax, ay), lb = Math.hypot(bx, by), c = Math.max(-0.9999, Math.min(0.9999, (ax * bx + ay * by) / (la * lb)));
    const th = Math.acos(c), dU = KT * (th - T0), s = Math.sin(th);
    const fix = dU / s * (bx / (la * lb) - c * ax / (la * la)), fiy = dU / s * (by / (la * lb) - c * ay / (la * la));
    const fkx = dU / s * (ax / (la * lb) - c * bx / (lb * lb)), fky = dU / s * (ay / (la * lb) - c * by / (lb * lb));
    Fa[i][0] += fix; Fa[i][1] += fiy; Fa[k][0] += fkx; Fa[k][1] += fky; Fa[j][0] -= fix + fkx; Fa[j][1] -= fiy + fky;
    return 0.5 * KT * (th - T0) ** 2;
  }
  function forces() {
    Fb = [[0, 0], [0, 0], [0, 0]]; Fa = [[0, 0], [0, 0], [0, 0]];
    Ub = bondTerm(1, 0) + bondTerm(2, 0); Ua = angleTerm(1, 0, 2); U = Ub + Ua;
    F = Fb.map((f, i) => [f[0] + Fa[i][0], f[1] + Fa[i][1]]);
  }
  const kinetic = () => v.reduce((s, q, i) => s + 0.5 * M[i] * (q[0] ** 2 + q[1] ** 2), 0);
  const geom = () => { const r1 = Math.hypot(x[1][0] - x[0][0], x[1][1] - x[0][1]), r2 = Math.hypot(x[2][0] - x[0][0], x[2][1] - x[0][1]);
    const c = ((x[1][0] - x[0][0]) * (x[2][0] - x[0][0]) + (x[1][1] - x[0][1]) * (x[2][1] - x[0][1])) / (r1 * r2); return { r1, r2, th: Math.acos(Math.max(-1, Math.min(1, c))) }; };
  function setup(r1, r2, th) {
    x = [[0, 0], [r1 * Math.cos(-Math.PI / 2 - th / 2), r1 * Math.sin(-Math.PI / 2 - th / 2)], [r2 * Math.cos(-Math.PI / 2 + th / 2), r2 * Math.sin(-Math.PI / 2 + th / 2)]];
    const cm = [0, 1].map(d => x.reduce((s, q, i) => s + M[i] * q[d], 0) / 18.015); x = x.map(q => [q[0] - cm[0], q[1] - cm[1]]);
    v = [[0, 0], [0, 0], [0, 0]]; forces(); E0 = null; hist = []; tfs = 0;
  }
  function step() {
    const dt = st.dt / 1000;
    for (let i = 0; i < 3; i++) { v[i][0] += 0.5 * dt * F[i][0] / M[i]; v[i][1] += 0.5 * dt * F[i][1] / M[i]; x[i][0] += dt * v[i][0]; x[i][1] += dt * v[i][1]; }
    forces();
    for (let i = 0; i < 3; i++) { v[i][0] += 0.5 * dt * F[i][0] / M[i]; v[i][1] += 0.5 * dt * F[i][1] / M[i]; }
    tfs += st.dt; const E = U + kinetic(); if (E0 === null) E0 = E;
    hist.push({ t: tfs, r1: geom().r1, E }); if (hist.length > HIST) hist.shift();
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Drag an atom, then let go</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-p="h" class="on">Harmonic bonds</button><button data-p="m">Morse bonds</button></div>
    <label><span>Time step</span><output data-o="dt"></output><input type="range" data-k="dt" min="0" max="4" step="1" value="1"></label>
    <div class="row"><button class="btn sm primary" data-k="stretch">Stretch a bond</button><button class="btn sm ghost" data-k="bend">Bend</button><button class="btn sm ghost" data-k="reset">Rest</button></div>
    <div class="readouts">
      <div><span>O–H bond lengths</span><b data-o="r"></b></div>
      <div><span>H–O–H angle</span><b data-o="th"></b></div>
      <div><span>Bond energy / angle energy</span><b data-o="u"></b></div>
      <div><span>Total energy change</span><b data-o="dE"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-p]').forEach(b => b.onclick = () => { st.morse = b.dataset.p === 'm'; root.querySelectorAll('[data-p]').forEach(z => z.classList.toggle('on', z === b)); setup(0.115, R0, T0); });
  $('[data-k=dt]').oninput = e => { st.dt = DTS[+e.target.value]; E0 = null; hist = []; };
  $('[data-k=stretch]').onclick = () => setup(0.115, R0, T0);
  $('[data-k=bend]').onclick = () => setup(R0, R0, T0 - 0.35);
  $('[data-k=reset]').onclick = () => setup(R0, R0, T0);

  let W = 600, H = 600, S = 1000, ox = 300, oy = 200;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const toWorld = (px, py) => [(px - ox) / S, (py - oy) / S];
  cv.addEventListener('pointerdown', e => {
    const r = cv.getBoundingClientRect(), p = toWorld(e.clientX - r.left, e.clientY - r.top);
    let best = -1, bd = 0.04; x.forEach((q, i) => { const d = Math.hypot(q[0] - p[0], q[1] - p[1]); if (d < bd) { bd = d; best = i; } });
    if (best >= 0) { st.drag = best; cv.setPointerCapture(e.pointerId); }
  });
  cv.addEventListener('pointermove', e => {
    if (st.drag < 0) return; const r = cv.getBoundingClientRect(), p = toWorld(e.clientX - r.left, e.clientY - r.top);
    x[st.drag] = p; v = [[0, 0], [0, 0], [0, 0]]; forces(); E0 = null; hist = [];
  });
  const up = () => { if (st.drag >= 0) { st.drag = -1; E0 = null; hist = []; tfs = 0; } };
  cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);

  function arrow(x0, y0, dx, dy, c) {
    const L = Math.hypot(dx, dy); if (L < 3) return; const ux = dx / L, uy = dy / L;
    cx.strokeStyle = c; cx.fillStyle = c; cx.lineWidth = 3; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0 + dx - ux * 8, y0 + dy - uy * 8); cx.stroke();
    cx.beginPath(); cx.moveTo(x0 + dx, y0 + dy); cx.lineTo(x0 + dx - ux * 12 - uy * 6, y0 + dy - uy * 12 + ux * 6); cx.lineTo(x0 + dx - ux * 12 + uy * 6, y0 + dy - uy * 12 - ux * 6); cx.closePath(); cx.fill();
  }
  function draw() {
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(150, H * 0.32), areaH = H - plotH - 20;
    S = Math.min(W, areaH) / 0.34; ox = W / 2; oy = areaH * 0.68;
    const P = x.map(q => [ox + q[0] * S, oy + q[1] * S]);
    cx.strokeStyle = 'rgba(143,163,199,.7)'; cx.lineWidth = 6; cx.beginPath(); cx.moveTo(P[1][0], P[1][1]); cx.lineTo(P[0][0], P[0][1]); cx.lineTo(P[2][0], P[2][1]); cx.stroke();
    const g = geom(); cx.strokeStyle = COL.amber; cx.lineWidth = 2.5; const a1 = Math.atan2(x[1][1] - x[0][1], x[1][0] - x[0][0]), a2 = Math.atan2(x[2][1] - x[0][1], x[2][0] - x[0][0]);
    cx.beginPath(); cx.arc(P[0][0], P[0][1], 0.04 * S, Math.min(a1, a2), Math.max(a1, a2)); cx.stroke();
    [[P[1], 0.016, COL.ink], [P[2], 0.016, COL.ink], [P[0], 0.026, COL.coral]].forEach(([p, r, c]) => { cx.fillStyle = c; cx.beginPath(); cx.arc(p[0], p[1], r * S, 0, 7); cx.fill(); });
    const fs = 0.06 * S / 2000;           // 2000 kJ/mol/nm ~ length of 0.06 nm on screen
    for (let i = 0; i < 3; i++) { arrow(P[i][0], P[i][1], Fb[i][0] * fs, Fb[i][1] * fs, COL.cyan); arrow(P[i][0], P[i][1], Fa[i][0] * fs, Fa[i][1] * fs, COL.amber); }
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = COL.cyan; cx.fillText('bond forces', 12, 22); cx.fillStyle = COL.amber; cx.fillText('angle forces', 12, 40);
    cx.fillStyle = COL.muted; cx.fillText(`t = ${tfs.toFixed(1)} fs`, 12, 58);
    // trace: bond length and energy change
    const x0 = 46, x1 = W - 12, top = areaH + 16, mid = top + plotH * 0.5, y1 = H - 10;
    if (hist.length < 2) { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.fillText(st.drag >= 0 ? 'release to start' : '', W / 2, top + plotH / 2); return; }
    const X = k => x0 + (x1 - x0) * k / (HIST - 1);
    const rv = hist.map(q => q.r1), rlo = Math.min(...rv), rhi = Math.max(...rv) + 1e-6;
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(x0, top); cx.lineTo(x0, mid - 8); cx.lineTo(x1, mid - 8); cx.moveTo(x0, mid + 14); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.cyan; cx.textAlign = 'left'; cx.fillText('O–H₁ length', x0 + 6, top + 2);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 1.6; cx.beginPath(); hist.forEach((q, k) => { const Y = mid - 8 - (q.r1 - rlo) / (rhi - rlo) * (mid - 18 - top); k ? cx.lineTo(X(k), Y) : cx.moveTo(X(k), Y); }); cx.stroke();
    const ev = hist.map(q => q.E - E0), elim = Math.max(0.05, ...ev.map(Math.abs)) * 1.1, YE = e => (mid + 14 + y1) / 2 - e / elim * (y1 - mid - 18) / 2;
    cx.fillStyle = COL.amber; cx.fillText(`total energy change (±${elim.toFixed(2)} kJ/mol)`, x0 + 6, mid + 14);
    cx.strokeStyle = COL.amber; cx.beginPath(); hist.forEach((q, k) => (k ? cx.lineTo(X(k), YE(q.E - E0)) : cx.moveTo(X(k), YE(q.E - E0)))); cx.stroke();
  }
  function readouts() {
    const g = geom(); o.dt.textContent = st.dt + ' fs';
    o.r.textContent = `${g.r1.toFixed(4)} / ${g.r2.toFixed(4)} nm`; o.th.textContent = (g.th * 180 / Math.PI).toFixed(1) + '°';
    o.u.textContent = `${Ub.toFixed(1)} / ${Ua.toFixed(1)} kJ/mol`;
    if (!hist.length || E0 === null) { o.dE.textContent = '–'; return; }
    const ev = hist.map(q => q.E - E0), sp = Math.max(...ev) - Math.min(...ev), Evib = E0 - 0;
    const pct = 50 * sp / Math.max(Evib, 1e-9);
    o.dE.innerHTML = `<span style="color:${pct > 3 ? COL.coral : COL.cyan}">± ${(sp / 2).toFixed(3)} kJ/mol (${(50 * sp / Math.max(Evib, 1e-9)).toFixed(2)} % of ${Evib.toFixed(1)})</span>`;
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running && st.drag < 0) { acc += FS_PER_FRAME; while (acc >= st.dt) { step(); acc -= st.dt; } }
    draw(); if (frame++ % 6 === 0) readouts();
  };
  setup(0.115, R0, T0); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
