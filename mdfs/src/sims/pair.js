// Two-atom Lennard-Jones explorer: distance slider, drag, energy/force curves, and "release" (1D velocity Verlet).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RMIN_PLOT = 0.8, RMAX_PLOT = 3.2;

export function mountPair(root) {
  const st = { r: 1.6, eps: 1, sig: 1, showF: false, running: false, v: 0, msg: '' };
  const V = r => 4 * st.eps * ((st.sig / r) ** 12 - (st.sig / r) ** 6);
  const F = r => 24 * st.eps * (2 * (st.sig / r) ** 12 - (st.sig / r) ** 6) / r;   // + = push

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view pair-view"><canvas></canvas><div class="sim-msg" hidden></div>
    <div class="sim-tip">Drag the <b style="color:${COL.coral}">red atom</b> or use the slider</div></div>
  <div class="sim-ctrl">
    <label>Distance r <output data-o="r"></output><input type="range" data-k="r" min="0.85" max="3.0" step="0.005"></label>
    <label>Well depth ε <output data-o="eps"></output><input type="range" data-k="eps" min="0.3" max="2" step="0.05" value="1"></label>
    <label>Atom size σ <output data-o="sig"></output><input type="range" data-k="sig" min="0.8" max="1.3" step="0.01" value="1"></label>
    <label class="check"><input type="checkbox" data-k="showF"> Show force curve F(r)</label>
    <div class="row"><button class="btn sm primary" data-k="go">▶ Release atoms</button><button class="btn sm ghost" data-k="min">Jump to minimum</button></div>
    <div class="readouts">
      <div><span>Energy V</span><b data-o="V">–</b></div>
      <div><span>Force F</span><b data-o="F">–</b></div>
      <div class="wide"><span>What the atoms feel</span><b data-o="state">–</b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.pair-view'), msg = root.querySelector('.sim-msg');
  const $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const slider = $('[data-k=r]');
  const stop = () => { st.running = false; $('[data-k=go]').textContent = '▶ Release atoms'; };
  slider.oninput = e => { stop(); msg.hidden = true; st.r = +e.target.value; };
  $('[data-k=eps]').oninput = e => { st.eps = +e.target.value; };
  $('[data-k=sig]').oninput = e => { st.sig = +e.target.value; };
  $('[data-k=showF]').onchange = e => { st.showF = e.target.checked; };
  $('[data-k=min]').onclick = () => { stop(); msg.hidden = true; st.r = 2 ** (1 / 6) * st.sig; };
  $('[data-k=go]').onclick = () => {
    if (st.running) return stop();
    msg.hidden = true; st.v = 0; st.running = true; $('[data-k=go]').textContent = '❚❚ Stop';
  };

  // layout (CSS pixels)
  let W = 600, H = 400;
  const resize = () => {
    const dpr = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight;
    cv.width = W * dpr; cv.height = H * dpr; cx.setTransform(dpr, 0, 0, dpr, 0, 0);
  };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const plot = () => ({ x0: 56, x1: W - 24, y0: H * 0.44, y1: H - 34 });
  const X = r => { const p = plot(); return p.x0 + (r - RMIN_PLOT) / (RMAX_PLOT - RMIN_PLOT) * (p.x1 - p.x0); };
  const Rx = x => { const p = plot(); return RMIN_PLOT + (x - p.x0) / (p.x1 - p.x0) * (RMAX_PLOT - RMIN_PLOT); };
  const vlo = -2.3, vhi = 2.5;
  const Y = v => { const p = plot(); return p.y1 - (v - vlo) / (vhi - vlo) * (p.y1 - p.y0); };

  // drag the red atom
  let drag = false;
  const atomY = () => H * 0.21;
  cv.addEventListener('pointerdown', e => {
    const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
    const bx = atomX(st.r).b;
    if (Math.hypot(x - bx, y - atomY()) < atomRad() * 1.6 || y > plot().y0) { drag = true; cv.setPointerCapture(e.pointerId); stop(); msg.hidden = true; moveTo(x, y); }
  });
  cv.addEventListener('pointermove', e => { if (drag) { const b = cv.getBoundingClientRect(); moveTo(e.clientX - b.left, e.clientY - b.top); } });
  cv.addEventListener('pointerup', () => { drag = false; });
  function moveTo(x, y) {
    if (y > plot().y0) st.r = Rx(x);                 // on the plot: x maps to r
    else st.r = 2 * (x - W / 2) / k();                // on the atoms: the pair stays centred
    st.r = Math.max(0.85, Math.min(3.0, st.r));
  }
  const k = () => Math.min(H * 0.16, (W - 80) / 3.6);          // pixels per sigma in the atom row
  const atomRad = () => 0.5 * st.sig * k();
  const atomX = r => ({ a: W / 2 - r * k() / 2, b: W / 2 + r * k() / 2 });

  function arrow(x, y, len, color) {
    if (Math.abs(len) < 3) return;
    const d = Math.sign(len); cx.strokeStyle = cx.fillStyle = color; cx.lineWidth = 5; cx.lineCap = 'round';
    cx.beginPath(); cx.moveTo(x, y); cx.lineTo(x + len - d * 10, y); cx.stroke();
    cx.beginPath(); cx.moveTo(x + len, y); cx.lineTo(x + len - d * 16, y - 9); cx.lineTo(x + len - d * 16, y + 9); cx.closePath(); cx.fill();
  }
  function ball(x, y, rad, c) {
    const g = cx.createRadialGradient(x - rad * 0.35, y - rad * 0.4, rad * 0.1, x, y, rad);
    g.addColorStop(0, '#ffffff'); g.addColorStop(0.18, c); g.addColorStop(1, shade(c));
    cx.fillStyle = g; cx.beginPath(); cx.arc(x, y, rad, 0, Math.PI * 2); cx.fill();
  }
  const shade = c => c === COL.cyan ? '#137A96' : '#9C2F3A';

  function curve(fn, color, width, dash) {
    const p = plot(); cx.save(); cx.beginPath(); cx.rect(p.x0, p.y0 - 6, p.x1 - p.x0, p.y1 - p.y0 + 12); cx.clip();
    cx.strokeStyle = color; cx.lineWidth = width; cx.setLineDash(dash || []); cx.beginPath();
    for (let x = p.x0; x <= p.x1; x += 2) { const r = Rx(x), y = Y(Math.max(vlo - 1, Math.min(vhi + 1, fn(r)))); x === p.x0 ? cx.moveTo(x, y) : cx.lineTo(x, y); }
    cx.stroke(); cx.restore();
  }

  function draw() {
    cx.clearRect(0, 0, W, H);
    const p = plot(), r = st.r, f = F(r), v = V(r);
    // atoms
    const { a, b } = atomX(r), y = atomY(), rad = atomRad();
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 2;
    cx.beginPath(); cx.moveTo(a, y + rad + 16); cx.lineTo(b, y + rad + 16); cx.moveTo(a, y + rad + 8); cx.lineTo(a, y + rad + 24); cx.moveTo(b, y + rad + 8); cx.lineTo(b, y + rad + 24); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '500 14px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(`r = ${r.toFixed(2)}`, (a + b) / 2, y + rad + 42);
    ball(a, y, rad, COL.cyan); ball(b, y, rad, COL.coral);
    const L = Math.sign(f) * Math.min(90, 34 * Math.log1p(Math.abs(f) * 2));
    arrow(a, y - rad - 16, -L, COL.amber); arrow(b, y - rad - 16, L, COL.amber);
    // axes
    cx.strokeStyle = 'rgba(143,163,199,.55)'; cx.lineWidth = 1.5;
    cx.beginPath(); cx.moveTo(p.x0, p.y0 - 6); cx.lineTo(p.x0, p.y1); cx.moveTo(p.x0, Y(0)); cx.lineTo(p.x1, Y(0)); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '500 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.textAlign = 'right'; cx.fillText('energy', p.x0 - 6, p.y0 + 10); cx.fillText('r (units of σ₀)', p.x1, Y(0) - 8);
    for (const t of [1, 1.5, 2, 2.5, 3]) { cx.textAlign = 'center'; cx.fillText(t, X(t), p.y1 + 18); cx.beginPath(); cx.moveTo(X(t), p.y1 - 3); cx.lineTo(X(t), p.y1 + 3); cx.stroke(); }
    // curves
    if (st.showF) curve(F, COL.cyan, 2.5, [8, 6]);
    curve(V, COL.amber, 3.5);
    // total energy line while running
    if (st.running) { const E = 0.25 * st.v * st.v + v; cx.strokeStyle = 'rgba(242,245,250,.35)'; cx.setLineDash([4, 6]); cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(p.x0, Y(E)); cx.lineTo(p.x1, Y(E)); cx.stroke(); cx.setLineDash([]); }
    // marker
    cx.fillStyle = COL.ink; cx.shadowColor = 'rgba(255,255,255,.8)'; cx.shadowBlur = 12;
    cx.beginPath(); cx.arc(X(r), Y(Math.max(vlo, Math.min(vhi, v))), 7, 0, Math.PI * 2); cx.fill(); cx.shadowBlur = 0;
    // legend
    cx.textAlign = 'right'; cx.font = '500 13px "Space Grotesk", sans-serif';
    cx.fillStyle = COL.amber; cx.fillText('V(r)', p.x1, p.y0 + 4);
    if (st.showF) { cx.fillStyle = COL.cyan; cx.fillText('F(r)', p.x1 - 44, p.y0 + 4); }
    // readouts
    o.r.textContent = r.toFixed(3) + ' σ₀'; o.eps.textContent = st.eps.toFixed(2); o.sig.textContent = st.sig.toFixed(2);
    if (!drag && document.activeElement !== slider) slider.value = r;
    o.V.textContent = v.toFixed(3) + ' ε₀'; o.F.textContent = f.toFixed(3);
    const rm = 2 ** (1 / 6) * st.sig, near = Math.abs(r - rm) < 0.015 * st.sig;
    o.state.innerHTML = near ? `<span style="color:${COL.amber}">Balanced: no net force</span>` : f > 0 ? `<span style="color:${COL.coral}">Pushing apart (too close)</span>` : `<span style="color:${COL.cyan}">Pulling together</span>`;
  }

  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  function loop() {
    raf = requestAnimationFrame(loop);
    if (!visible) return;
    if (st.running) {
      const mu = 0.5, dt = 0.0015;      // reduced mass of two unit-mass atoms
      for (let k = 0; k < 6; k++) {
        st.v += 0.5 * dt * F(st.r) / mu; st.r += dt * st.v; st.v += 0.5 * dt * F(st.r) / mu;
      }
      if (st.r > 3.0) {
        st.r = 3.0; stop(); msg.hidden = false;
        msg.innerHTML = '<b>They flew apart!</b> The atoms started with more energy than the well depth ε, so the pull could not hold them. Start them between r = 1.05 and 1.5 and release again.';
      }
    }
    draw();
  }
  loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
