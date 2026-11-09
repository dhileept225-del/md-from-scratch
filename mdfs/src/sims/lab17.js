// Integrator extras lab widget.
// Shadow: harmonic oscillator phase space, Verlet (points stay on the shadow ellipse) vs RK4 (spirals inward).
// Precision: 7-atom 2D LJ cluster with positions stored in float32 vs float64, cluster offset from the origin.
// RESPA: x'' = -(wf^2 + ws^2) x, wf = 10, ws = 1; outer step H, inner step ~0.005; resonance at H = k pi / omega.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const WF = 10, WS = 1, WT = Math.sqrt(WF * WF + WS * WS);

export function mountLab17(root) {
  const st = { mode: 'sh', h: 0.5, rk: false, off: 2, H: 0.1, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view ix-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="sh" class="on">Shadow</button><button data-m="pr">Precision</button><button data-m="re">RESPA</button></div>
    <div data-p="sh"><label><span>Time step h (×1/ω)</span><output data-o="h"></output><input type="range" data-k="h" min="0.05" max="2.1" step="0.05" value="0.5"></label>
      <div class="seg"><button data-r="0" class="on">Verlet</button><button data-r="1">RK4</button></div></div>
    <div data-p="pr" hidden><label><span>Distance of the cluster from the origin</span><output data-o="off"></output><input type="range" data-k="off" min="0" max="4" step="1" value="2"></label></div>
    <div data-p="re" hidden><label><span>Outer step H (slow force)</span><output data-o="H"></output><input type="range" data-k="H" min="0.02" max="0.7" step="0.0025" value="0.1"></label></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.ix-view'), $ = s => root.querySelector(s);
  const f32 = Math.fround;
  // ---- shadow ----
  let sx, sp, sn, pts;
  const vv = (x, p, h) => { p -= 0.5 * h * x; x += h * p; p -= 0.5 * h * x; return [x, p]; };
  const rk4 = (x, p, h) => { const k1 = [p, -x], k2 = [p + h / 2 * k1[1], -(x + h / 2 * k1[0])], k3 = [p + h / 2 * k2[1], -(x + h / 2 * k2[0])], k4 = [p + h * k3[1], -(x + h * k3[0])];
    return [x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), p + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])]; };
  function resetSh() { sx = 1; sp = 0; sn = 0; pts = [[1, 0]]; }
  // ---- precision: 7-atom LJ cluster (2D) ----
  const OFFS = [0, 10, 100, 1000, 10000], H7 = 0.005;
  let c32, c64, tP, Etr;
  const base = [[0, 0], ...Array.from({ length: 6 }, (_, k) => [1.12 * Math.cos(k * Math.PI / 3), 1.12 * Math.sin(k * Math.PI / 3)])];
  const v0 = [[0.31, -0.22], [-0.4, 0.15], [0.12, 0.45], [-0.28, -0.36], [0.44, 0.08], [-0.05, -0.3], [-0.14, 0.2]];
  function ljF(x) { const n = x.length, F = x.map(() => [0, 0]); let U = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const dx = x[i][0] - x[j][0], dy = x[i][1] - x[j][1], r2 = dx * dx + dy * dy, s6 = 1 / (r2 * r2 * r2), f = 24 * (2 * s6 * s6 - s6) / r2;
      U += 4 * (s6 * s6 - s6); F[i][0] += f * dx; F[i][1] += f * dy; F[j][0] -= f * dx; F[j][1] -= f * dy; } return [F, U]; }
  function mk(single) { const o = OFFS[st.off], r = single ? f32 : (q => q); const x = base.map(p => [r(p[0] + o), r(p[1] + o)]), p = v0.map(q => [r(q[0]), r(q[1])]); const s = { x, p, r, F: ljF(x)[0] }; s.E0 = en(s); return s; }
  function en(s) { const [, U] = ljF(s.x); return U + 0.5 * s.p.reduce((a, q) => a + q[0] * q[0] + q[1] * q[1], 0); }
  function stepC(s) { const r = s.r, h = H7;
    s.p = s.p.map((q, i) => [r(q[0] + r(0.5 * h * s.F[i][0])), r(q[1] + r(0.5 * h * s.F[i][1]))]);
    s.x = s.x.map((q, i) => [r(q[0] + r(h * s.p[i][0])), r(q[1] + r(h * s.p[i][1]))]); s.F = ljF(s.x)[0];
    s.p = s.p.map((q, i) => [r(q[0] + r(0.5 * h * s.F[i][0])), r(q[1] + r(0.5 * h * s.F[i][1]))]); }
  function resetPr() { c32 = mk(true); c64 = mk(false); tP = 0; Etr = []; }
  // ---- RESPA ----
  let rx, rp, vx, vp, tR, Rtr, rmax, vmax, vdead;
  const Er = (x, p) => 0.5 * p * p + 0.5 * WT * WT * x * x;
  function respa(x, p, H) { const n = Math.ceil(H / 0.005), h = H / n; p -= 0.5 * H * WS * WS * x; for (let k = 0; k < n; k++) { p -= 0.5 * h * WF * WF * x; x += h * p; p -= 0.5 * h * WF * WF * x; } p -= 0.5 * H * WS * WS * x; return [x, p]; }
  function resetRe() { rx = 1; rp = 0; vx = 1; vp = 0; tR = 0; Rtr = []; rmax = 0; vmax = 0; vdead = false; }

  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); draw(); });
  root.querySelectorAll('[data-r]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-r]').forEach(z => z.classList.toggle('on', z === b)); st.rk = b.dataset.r === '1'; resetSh(); draw(); });
  $('[data-k=h]').oninput = e => { st.h = +e.target.value; resetSh(); draw(); };
  $('[data-k=off]').oninput = e => { st.off = +e.target.value; resetPr(); draw(); };
  $('[data-k=H]').oninput = e => { st.H = +e.target.value; resetRe(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { ({ sh: resetSh, pr: resetPr, re: resetRe })[st.mode](); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const label = (t, x_, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x_, y); };
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const fmt = q => !isFinite(q) ? '∞' : Math.abs(q) >= 1e4 || (Math.abs(q) < 1e-3 && q !== 0) ? q.toExponential(2) : q.toFixed(4);
  function series(tr, x0, y0, x1, y1, keys, cols, sym) {
    if (tr.length < 2) return; let m = 1e-12; tr.forEach(r => keys.forEach(k => { if (isFinite(r[k])) m = Math.max(m, Math.abs(r[k])); })); if (!sym) m *= 1.1;
    const ta = tr[0][0], tb = tr[tr.length - 1][0] || 1, X = t => x0 + (t - ta) / (tb - ta || 1) * (x1 - x0), Y = q => (y0 + y1) / 2 - Math.max(-m, Math.min(m, q)) / m * (y1 - y0) / 2;
    cx.strokeStyle = 'rgba(143,163,199,.35)'; cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
    keys.forEach((k, j) => { cx.strokeStyle = cols[j]; cx.lineWidth = 1.6; cx.beginPath(); tr.forEach((r, i) => i ? cx.lineTo(X(r[0]), Y(r[k])) : cx.moveTo(X(r[0]), Y(r[k]))); cx.stroke(); });
    label(`±${m.toExponential(1)}`, x0 + 4, y0 + 12, COL.muted, 'left', '10px');
  }
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const R = $('[data-ro]');
    if (st.mode === 'sh') {
      const sz = Math.min(W, H) * 0.42, ox = W / 2, oy = H / 2, S = sz / 1.3;
      cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.beginPath(); cx.moveTo(ox - sz, oy); cx.lineTo(ox + sz, oy); cx.moveTo(ox, oy - sz); cx.lineTo(ox, oy + sz); cx.stroke();
      label('x', ox + sz - 8, oy - 6, COL.muted); label('p', ox + 6, oy - sz + 10, COL.muted);
      cx.strokeStyle = 'rgba(242,245,250,.8)'; cx.lineWidth = 1.5; cx.beginPath(); cx.arc(ox, oy, S, 0, 7); cx.stroke();
      const q = 1 - st.h * st.h / 4;
      if (!st.rk && q > 0) { const Hs = 0.5 * q; cx.strokeStyle = COL.amber; cx.setLineDash([5, 4]); cx.beginPath(); cx.ellipse(ox, oy, Math.sqrt(2 * Hs / q) * S, Math.sqrt(2 * Hs) * S, 0, 0, 7); cx.stroke(); cx.setLineDash([]); }
      cx.fillStyle = st.rk ? COL.coral : COL.cyan; pts.forEach(p => { if (Math.abs(p[0]) < 3 && Math.abs(p[1]) < 3) { cx.beginPath(); cx.arc(ox + p[0] * S, oy - p[1] * S, 2, 0, 7); cx.fill(); } });
      label(st.rk ? 'white: exact orbit · coral: RK4' : 'white: exact · amber dashed: shadow · cyan: Verlet', 10, H - 10, COL.muted, 'left', '11px');
      const E = 0.5 * (sx * sx + sp * sp), es = pts.map(p => 0.5 * (p[0] * p[0] + p[1] * p[1]));
      R.innerHTML = `<div><span>Steps</span><b>${sn.toLocaleString()}</b></div><div><span>Energy now (start 0.5)</span><b>${fmt(E)}</b></div>
        <div><span>${st.rk ? 'Lowest energy so far' : 'Shadow energy (constant)'}</span><b>${st.rk ? fmt(Math.min(...es)) : q > 0 ? fmt(0.5 * sp * sp + 0.5 * sx * sx * q) : 'none: unstable'}</b></div>
        <div><span>${st.rk ? 'Energy lost' : 'Shadow frequency / true'}</span><b>${st.rk ? (100 * (1 - E / 0.5)).toFixed(3) + ' %' : st.h < 2 ? (2 / st.h * Math.asin(st.h / 2)).toFixed(4) : 'unstable (h ≥ 2)'}</b></div>`;
      $('[data-o=h]').textContent = st.h.toFixed(2);
    } else if (st.mode === 'pr') {
      const x0 = 50, x1 = W - 12, y0 = 26, y1 = H - 26; axes(x0, y0, x1, y1); label('E − E0 vs time: cyan float64, coral float32', x0 + 6, y0 - 8);
      series(Etr, x0, y0, x1, y1, [1, 2], [COL.cyan, COL.coral], true);
      const o = OFFS[st.off], sp32 = o ? 2 ** (Math.floor(Math.log2(o)) - 23) : 2 ** -23 * 1.12;
      const last = Etr[Etr.length - 1];
      R.innerHTML = `<div><span>Time · steps</span><b>${tP.toFixed(1)} · ${Math.round(tP / H7).toLocaleString()}</b></div><div><span>Float32 spacing at the cluster</span><b>${sp32.toExponential(2)}</b></div>
        <div><span>Distance moved per step (about)</span><b>${(0.4 * H7).toExponential(1)}</b></div><div><span>E − E0: float64 · float32</span><b>${last ? fmt(last[1]) + ' · ' + fmt(last[2]) : '…'}</b></div>`;
      $('[data-o=off]').textContent = OFFS[st.off].toLocaleString();
    } else {
      const x0 = 50, x1 = W - 12, y0 = 26, y1 = H - 26; axes(x0, y0, x1, y1); label('E/E0 − 1 vs time: amber RESPA, coral Verlet with step H', x0 + 6, y0 - 8);
      series(Rtr, x0, y0, x1, y1, vdead ? [1] : [1, 2], [COL.amber, COL.coral], true);
      const n = Math.ceil(st.H / 0.005);
      R.innerHTML = `<div><span>H × omega / pi (resonance at 1, 2, …)</span><b>${(st.H * WT / Math.PI).toFixed(3)}</b></div><div><span>RESPA: max |E/E0 − 1|</span><b>${fmt(rmax)}</b></div>
        <div><span>Verlet at step H</span><b>${st.H * WT >= 2 ? 'unstable (H × omega ≥ 2)' : fmt(vmax)}</b></div><div><span>Slow-force evaluations per unit time · fast</span><b>${(1 / st.H).toFixed(0)} · ${(n / st.H).toFixed(0)}</b></div>`;
      $('[data-o=H]').textContent = st.H.toFixed(4);
    }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return;
    if (st.mode === 'sh') { for (let k = 0; k < (st.rk ? 20 : 3); k++) { [sx, sp] = (st.rk ? rk4 : vv)(sx, sp, st.h); sn++; if (!st.rk || sn % 2 === 0) pts.push([sx, sp]); } if (pts.length > 1500) pts.splice(0, pts.length - 1500); }
    else if (st.mode === 'pr') { for (let k = 0; k < 60; k++) { stepC(c32); stepC(c64); tP += H7; } Etr.push([tP, en(c64) - c64.E0, en(c32) - c32.E0]); if (Etr.length > 4000) Etr.shift(); }
    else { const E0 = Er(1, 0); for (let k = 0; k < Math.max(1, Math.round(0.2 / st.H)); k++) { [rx, rp] = respa(rx, rp, st.H); if (!vdead) { [vx, vp] = vv(vx, vp, st.H * WT); } tR += st.H;
        const er = Er(rx, rp) / E0 - 1, ev = 0.5 * (vx * vx + vp * vp) / 0.5 - 1; rmax = Math.max(rmax, Math.abs(er)); if (!vdead) vmax = Math.max(vmax, Math.abs(ev)); if (!(Math.abs(ev) < 1e6)) vdead = true;
        if (!(Math.abs(er) < 1e9)) { st.running = false; $('[data-k=run]').textContent = '▶ Run'; break; }
        Rtr.push([tR, er, vdead ? 0 : ev]); } if (Rtr.length > 3000) Rtr.splice(0, Rtr.length - 3000); }
    draw(); };
  resetSh(); resetPr(); resetRe(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
