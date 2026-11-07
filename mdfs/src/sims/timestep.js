// Time-step explorer: one bond vibrating, integrated with velocity Verlet at the chosen Δt (in fs).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const BONDS = { oh: { name: 'O–H', wn: 3700 }, ch: { name: 'C–H', wn: 3000 }, cc: { name: 'C–C', wn: 1000 } };
const WINDOW = 60;   // fs shown

export function mountTimestep(root) {
  const st = { bond: 'ch', dt: 1.0, born: performance.now(), res: null };
  const period = () => 1 / (BONDS[st.bond].wn * 2.998e-5);   // fs

  function simulate() {
    const T = period(), w = 2 * Math.PI / T, dt = st.dt, x0 = 1;
    let x = x0, v = 0, a = -w * w * x; const E0 = 0.5 * w * w * x0 * x0;
    const pts = [[0, x, 0]]; let worst = 0, blown = false;
    for (let t = dt; t <= WINDOW + 1e-9; t += dt) {
      v += 0.5 * a * dt; x += v * dt; a = -w * w * x; v += 0.5 * a * dt;
      const err = (0.5 * v * v + 0.5 * w * w * x * x - E0) / E0;
      worst = Math.max(worst, Math.abs(err)); pts.push([t, x, err]);
      if (Math.abs(x) > 6) { blown = true; break; }
    }
    return { T, w, pts, worst, blown };
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view ts-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg">${Object.entries(BONDS).map(([k, b]) => `<button data-b="${k}" class="${k === 'ch' ? 'on' : ''}">${b.name}</button>`).join('')}</div>
    <label>Time step Δt <output data-o="dt"></output><input type="range" data-k="dt" min="0.1" max="6" step="0.05" value="1"></label>
    <div class="readouts">
      <div><span>Vibration period</span><b data-o="T"></b></div>
      <div><span>Steps per vibration</span><b data-o="spv"></b></div>
      <div><span>Energy error</span><b data-o="err"></b></div>
      <div><span>Steps for 1 ns</span><b data-o="cost"></b></div>
      <div class="wide"><span>Verdict</span><b data-o="verdict"></b></div>
    </div>
    <p class="muted small">Top: bond length over time. <b style="color:${COL.amber}">Exact</b> curve vs <b style="color:${COL.cyan}">computed</b> steps. Bottom: this bond's energy error. In a full system the total energy fluctuates far less, because fast bonds hold only a small part of it.</p>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.ts-view');
  const o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const recompute = () => { st.res = simulate(); st.born = performance.now(); readouts(); };
  root.querySelectorAll('[data-b]').forEach(b => b.onclick = () => { st.bond = b.dataset.b; root.querySelectorAll('[data-b]').forEach(x => x.classList.toggle('on', x === b)); recompute(); });
  root.querySelector('[data-k=dt]').oninput = e => { st.dt = +e.target.value; recompute(); };

  function readouts() {
    const r = st.res, spv = r.T / st.dt;
    o.dt.textContent = st.dt.toFixed(2) + ' fs'; o.T.textContent = r.T.toFixed(1) + ' fs'; o.spv.textContent = spv.toFixed(1);
    o.err.textContent = r.blown ? '∞' : (100 * r.worst).toFixed(1) + ' %';
    o.cost.textContent = Math.round(1e6 / st.dt).toLocaleString('en-US');
    o.verdict.innerHTML = r.blown ? `<span style="color:${COL.coral}">Explodes: Δt is longer than T/π</span>`
      : r.worst < 0.10 ? `<span style="color:${COL.cyan}">Good: accurate and stable</span>`
      : r.worst < 0.35 ? `<span style="color:${COL.amber}">Rough: stable but inaccurate</span>`
      : `<span style="color:${COL.coral}">Bad: big energy errors</span>`;
  }

  let W = 600, H = 420;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function panel(y0, y1, label) {
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1.5; const yc = (y0 + y1) / 2;
    cx.beginPath(); cx.moveTo(50, y0); cx.lineTo(50, y1); cx.moveTo(50, yc); cx.lineTo(W - 18, yc); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '500 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(label, 58, y0 + 12);
    return yc;
  }
  const X = t => 50 + t / WINDOW * (W - 70);
  function draw(now) {
    cx.clearRect(0, 0, W, H);
    const r = st.res, grow = Math.min(1, (now - st.born) / 900), tEnd = grow * WINDOW;
    // top: bond length
    const y0 = 20, y1 = H * 0.58, yc = panel(y0, y1, 'bond stretch'), A = (y1 - y0) / 2 * 0.78;
    cx.save(); cx.beginPath(); cx.rect(50, y0, W - 60, y1 - y0); cx.clip();
    cx.strokeStyle = COL.amber; cx.lineWidth = 2.5; cx.globalAlpha = 0.75; cx.beginPath();
    for (let t = 0; t <= tEnd; t += 0.1) { const xx = X(t), yy = yc - A * Math.cos(r.w * t); t ? cx.lineTo(xx, yy) : cx.moveTo(xx, yy); } cx.stroke(); cx.globalAlpha = 1;
    const c = r.blown ? COL.coral : r.worst > 0.35 ? COL.coral : COL.cyan;
    const vis = r.pts.filter(p => p[0] <= tEnd);
    cx.strokeStyle = c; cx.fillStyle = c; cx.lineWidth = 2; cx.beginPath();
    vis.forEach(([t, x], i) => { const xx = X(t), yy = yc - A * x; i ? cx.lineTo(xx, yy) : cx.moveTo(xx, yy); }); cx.stroke();
    if (vis.length < 260) vis.forEach(([t, x]) => { cx.beginPath(); cx.arc(X(t), yc - A * x, 3.2, 0, Math.PI * 2); cx.fill(); });
    cx.restore();
    // bottom: energy error
    const e0 = H * 0.66, e1 = H - 30, ec = panel(e0, e1, 'energy error'), sc = (e1 - e0) / 2 / 0.5;   // ±50 % full scale
    cx.save(); cx.beginPath(); cx.rect(50, e0, W - 60, e1 - e0); cx.clip();
    cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath();
    vis.forEach(([t, , err], i) => { const xx = X(t), yy = ec - sc * err; i ? cx.lineTo(xx, yy) : cx.moveTo(xx, yy); }); cx.stroke();
    cx.restore();
    cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '500 12px "Space Grotesk", sans-serif';
    for (let t = 0; t <= WINDOW; t += 10) cx.fillText(t + (t === WINDOW ? ' fs' : ''), X(t), H - 10);
    cx.textAlign = 'right'; cx.fillText('+50%', 46, e0 + 10); cx.fillText('−50%', 46, e1);
  }
  let raf = 0;
  const loop = now => { raf = requestAnimationFrame(loop); draw(now); };
  recompute(); raf = requestAnimationFrame(loop);
  return () => { cancelAnimationFrame(raf); ro.disconnect(); };
}
