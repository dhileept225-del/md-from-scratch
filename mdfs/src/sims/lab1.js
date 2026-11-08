// Lab 1 widget: the hand calculation of Lab 1.7, one step at a time. Two atoms on a line, LJ in reduced units.
// Euler and velocity Verlet run side by side from the same start; the table shows every number of the chosen recipe.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const U = r => 4 * (r ** -12 - r ** -6), F = r => 24 * (2 * r ** -13 - r ** -7);

export function mountLab1(root) {
  const st = { dt: 0.01, show: 'verlet', running: false, r0: 1.05 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-r="euler">Table: Euler</button><button data-r="verlet" class="on">Table: Verlet</button></div>
    <div class="seg">${[0.005, 0.01, 0.03, 0.06, 0.12].map(d => `<button data-dt="${d}" class="${d === 0.01 ? 'on' : ''}">${d}</button>`).join('')}</div>
    <label><span>Starting distance r₀</span><output data-o="r0"></output><input type="range" data-k="r0" min="0.95" max="1.5" step="0.01" value="1.05"></label>
    <div class="row"><button class="btn sm primary" data-k="step">Step</button><button class="btn sm ghost" data-k="ten">×10</button><button class="btn sm ghost" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="tbl" style="font-size:12px"><table data-o="tab"></table></div>
    <div class="readouts">
      <div><span>Euler energy error</span><b data-o="ee"></b></div>
      <div><span>Verlet energy error</span><b data-o="ev"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let S, E0, n, hist, last;
  const energy = s => 0.5 * s.v1 ** 2 + 0.5 * s.v2 ** 2 + U(s.x2 - s.x1);
  function reset() { const s = { x1: 0, x2: st.r0, v1: 0, v2: 0 }; S = { euler: { ...s }, verlet: { ...s } }; E0 = energy(s); n = 0; hist = { euler: [[0, E0, st.r0]], verlet: [[0, E0, st.r0]] }; last = { euler: null, verlet: null }; rows = { euler: [rowOf('euler', null)], verlet: [rowOf('verlet', null)] }; }
  let rows = {};
  function rowOf(k, extra) { const s = S[k]; return { n, ...s, r: s.x2 - s.x1, F: F(s.x2 - s.x1), E: energy(s), extra }; }
  function step() {
    const dt = st.dt;
    { const s = S.euler, f = F(s.x2 - s.x1);   // Euler: old velocity moves, old force kicks
      const x1 = s.x1 + s.v1 * dt, x2 = s.x2 + s.v2 * dt; s.v1 -= f * dt; s.v2 += f * dt; s.x1 = x1; s.x2 = x2; last.euler = `used F(old) = ${f.toFixed(4)}`; }
    { const s = S.verlet, f = F(s.x2 - s.x1);   // Verlet: half kick, move, new force, half kick
      const v2h = s.v2 + 0.5 * f * dt, v1h = s.v1 - 0.5 * f * dt; s.x1 += v1h * dt; s.x2 += v2h * dt; const fn = F(s.x2 - s.x1); s.v1 = v1h - 0.5 * fn * dt; s.v2 = v2h + 0.5 * fn * dt;
      last.verlet = `v2_half = ${v2h.toFixed(5)}, F(new) = ${fn.toFixed(4)}`; }
    n++;
    for (const k of ['euler', 'verlet']) { const s = S[k]; hist[k].push([n * dt, energy(s), s.x2 - s.x1]); if (hist[k].length > 4000) hist[k].shift(); rows[k].push(rowOf(k, last[k])); if (rows[k].length > 6) rows[k].shift(); }
  }
  root.querySelectorAll('[data-r]').forEach(b => b.onclick = () => { st.show = b.dataset.r; root.querySelectorAll('[data-r]').forEach(z => z.classList.toggle('on', z === b)); draw(); });
  root.querySelectorAll('[data-dt]').forEach(b => b.onclick = () => { st.dt = +b.dataset.dt; root.querySelectorAll('[data-dt]').forEach(z => z.classList.toggle('on', z === b)); reset(); draw(); });
  $('[data-k=r0]').oninput = e => { st.r0 = +e.target.value; reset(); draw(); };
  $('[data-k=step]').onclick = () => { step(); draw(); };
  $('[data-k=ten]').onclick = () => { for (let q = 0; q < 10; q++) step(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !S) return;
    cx.clearRect(0, 0, W, H);
    // atoms: two rows (Euler, Verlet)
    const sc = Math.min(140, W / 4.5), cxm = W / 2;
    [['euler', COL.coral, 40], ['verlet', COL.cyan, 110]].forEach(([k, c, y]) => {
      const s = S[k], r = s.x2 - s.x1, mid = (s.x1 + s.x2) / 2, X1 = cxm + (s.x1 - mid) * sc, X2 = cxm + (s.x2 - mid) * sc, broke = r > 3;
      cx.fillStyle = COL.muted; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(k === 'euler' ? 'Euler' : 'Verlet', 10, y + 4);
      cx.strokeStyle = 'rgba(242,245,250,.25)'; cx.lineWidth = 2; if (!broke) { cx.beginPath(); cx.moveTo(X1, y); cx.lineTo(X2, y); cx.stroke(); }
      [X1, X2].forEach(X => { cx.fillStyle = c; cx.beginPath(); cx.arc(Math.max(-30, Math.min(W + 30, X)), y, 0.45 * sc * 0.5, 0, 7); cx.fill(); });
      cx.fillStyle = broke ? COL.coral : COL.ink; cx.textAlign = 'right'; cx.fillText(broke ? 'molecule broke!' : `r = ${r.toFixed(4)}`, W - 10, y + 4);
    });
    // energy plot
    const x0 = 50, x1 = W - 14, yt = 175, yb = H - 26, tmax = Math.max(10 * st.dt, ...hist.verlet.map(h => h[0]));
    const all = [...hist.euler, ...hist.verlet].map(h => h[1] - E0), lo = Math.min(-0.002, ...all), hi = Math.max(0.002, ...all);
    const X = t => x0 + t / tmax * (x1 - x0), Y = e => yb - (e - lo) / (hi - lo) * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('energy error E − E₀', x0 + 6, yt - 10);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [lo, 0, hi].forEach(e => cx.fillText(e.toPrecision(2), x0 - 4, Y(e) + 4)); cx.fillText(`time ${tmax.toFixed(2)}`, x1, yb + 14);
    [['euler', COL.coral], ['verlet', COL.cyan]].forEach(([k, c]) => { cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath(); hist[k].forEach((h, q) => q ? cx.lineTo(X(h[0]), Y(h[1] - E0)) : cx.moveTo(X(h[0]), Y(h[1] - E0))); cx.stroke(); });
    cx.textAlign = 'left'; cx.fillStyle = COL.coral; cx.fillText('Euler', x0 + 10, yt + 14); cx.fillStyle = COL.cyan; cx.fillText('Verlet', x0 + 60, yt + 14);
    // table
    const k = st.show, f = v => (Math.abs(v) < 1e-12 ? 0 : v).toFixed(5);
    o.tab.innerHTML = `<tr><th>n</th><th>x1</th><th>x2</th><th>v2</th><th>r</th><th>F</th><th>E</th></tr>` + rows[k].map(r => `<tr><td>${r.n}</td><td>${f(r.x1)}</td><td>${f(r.x2)}</td><td>${f(r.v2)}</td><td>${f(r.r)}</td><td>${r.F.toFixed(4)}</td><td>${r.E.toFixed(6)}</td></tr>`).join('') +
      (last[k] ? `<tr><td colspan="7" style="text-align:left;color:#8FA3C7">last step: ${last[k]}</td></tr>` : '');
    o.r0.textContent = st.r0.toFixed(2);
    o.ee.textContent = (energy(S.euler) - E0).toExponential(2); o.ev.textContent = (energy(S.verlet) - E0).toExponential(2);
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let q = 0; q < 2; q++) step(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
