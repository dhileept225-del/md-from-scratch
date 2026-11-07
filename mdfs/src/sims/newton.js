// Newton's law, one step at a time: F -> a = F/m -> v += a*dt -> x += v*dt, with every number shown.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const START = { spring: { x: [4, 0], v: [0, 2.2] }, push: { x: [-5, 3.2], v: [1.6, 0.8] } };

export function mountNewton(root) {
  const st = { mode: 'spring', m: 1, dt: 0.15, x: [0, 0], v: [0, 0], F: [0, 0], a: [0, 0], n: 0, trail: [], running: false, acc: 0 };
  const force = x => st.mode === 'spring' ? [-x[0], -x[1]] : [0, -1];   // spring k = 1, or constant push of size 1
  function reset() {
    const s = START[st.mode]; st.x = [...s.x]; st.v = [...s.v]; st.n = 0; st.trail = [[...st.x]];
    st.F = force(st.x); st.a = [st.F[0] / st.m, st.F[1] / st.m]; msg.hidden = true;
  }
  function step() {
    st.F = force(st.x);                                   // 1. force
    st.a = [st.F[0] / st.m, st.F[1] / st.m];              // 2. a = F / m
    st.v = [st.v[0] + st.a[0] * st.dt, st.v[1] + st.a[1] * st.dt];   // 3. update velocity
    st.x = [st.x[0] + st.v[0] * st.dt, st.x[1] + st.v[1] * st.dt];   // 4. move
    st.n++; st.trail.push([...st.x]); if (st.trail.length > 400) st.trail.shift();
    const far = Math.hypot(...st.x) > 9 || st.x[1] < -6;
    if (far) {
      pause(); msg.hidden = false;
      msg.innerHTML = st.mode === 'spring'
        ? '<b>The atom escaped!</b> A light atom vibrates fast, and this time step is too big to follow it: every jump overshoots. Lower Δt or raise the mass, then press Reset. This is why light hydrogen atoms limit the time step in real MD (Lesson 1.4).'
        : '<b>The atom left the box.</b> Press Reset to throw it again. Try a heavier mass: the same push bends its path less.';
    }
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view newton-view"><canvas></canvas><div class="sim-msg" hidden></div>
    <div class="sim-tip">Dots = one dot per time step · <b style="color:${COL.amber}">force</b> · <b style="color:${COL.ink}">velocity</b></div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-mode="spring" class="on">Spring (a bond)</button><button data-mode="push">Constant push</button></div>
    <label>Mass m <output data-o="m"></output><input type="range" data-k="m" min="0.25" max="4" step="0.05" value="1"></label>
    <label>Time step Δt <output data-o="dt"></output><input type="range" data-k="dt" min="0.02" max="1.2" step="0.01" value="0.15"></label>
    <div class="row"><button class="btn sm ghost" data-k="step">Step ▸</button><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">↺ Reset</button></div>
    <table class="mem"><caption>In memory after step <b data-o="n">0</b></caption>
      <tr><th></th><th>x</th><th>y</th></tr>
      <tr><td style="color:${COL.amber}">force F</td><td data-o="Fx"></td><td data-o="Fy"></td></tr>
      <tr><td>accel a = F/m</td><td data-o="ax"></td><td data-o="ay"></td></tr>
      <tr><td>velocity v</td><td data-o="vx"></td><td data-o="vy"></td></tr>
      <tr><td style="color:${COL.cyan}">position</td><td data-o="px"></td><td data-o="py"></td></tr>
    </table>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.newton-view'), msg = root.querySelector('.sim-msg');
  const $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const runBtn = $('[data-k=run]');
  function pause() { st.running = false; runBtn.textContent = '▶ Run'; }
  root.querySelectorAll('[data-mode]').forEach(b => b.onclick = () => {
    st.mode = b.dataset.mode; root.querySelectorAll('[data-mode]').forEach(x => x.classList.toggle('on', x === b)); pause(); reset();
  });
  $('[data-k=m]').oninput = e => { st.m = +e.target.value; };
  $('[data-k=dt]').oninput = e => { st.dt = +e.target.value; };
  $('[data-k=step]').onclick = () => { pause(); if (msg.hidden) step(); };
  runBtn.onclick = () => { if (st.running) return pause(); if (!msg.hidden) reset(); st.running = true; runBtn.textContent = '❚❚ Pause'; };
  $('[data-k=reset]').onclick = () => { pause(); reset(); };

  let W = 600, H = 420;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const S = () => Math.min(W, H) / 13;
  const P = x => [W / 2 + x[0] * S(), H / 2 - x[1] * S()];

  function arrow(from, vec, color, k) {
    const [x0, y0] = P(from), dx = vec[0] * k * S(), dy = -vec[1] * k * S(), L = Math.hypot(dx, dy);
    if (L < 4) return;
    const ux = dx / L, uy = dy / L, x1 = x0 + dx, y1 = y0 + dy;
    cx.strokeStyle = cx.fillStyle = color; cx.lineWidth = 4; cx.lineCap = 'round';
    cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x1 - ux * 10, y1 - uy * 10); cx.stroke();
    cx.beginPath(); cx.moveTo(x1, y1); cx.lineTo(x1 - ux * 15 - uy * 8, y1 - uy * 15 + ux * 8); cx.lineTo(x1 - ux * 15 + uy * 8, y1 - uy * 15 - ux * 8); cx.closePath(); cx.fill();
  }
  function draw() {
    cx.clearRect(0, 0, W, H);
    // background grid
    cx.strokeStyle = 'rgba(143,163,199,.08)'; cx.lineWidth = 1;
    for (let g = -12; g <= 12; g++) { const [gx] = P([g, 0]), [, gy] = P([0, g]); cx.beginPath(); cx.moveTo(gx, 0); cx.lineTo(gx, H); cx.moveTo(0, gy); cx.lineTo(W, gy); cx.stroke(); }
    if (st.mode === 'spring') {
      const [ax, ay] = P([0, 0]), [bx, by] = P(st.x), n = 14, L = Math.hypot(bx - ax, by - ay) || 1, ux = (bx - ax) / L, uy = (by - ay) / L;
      cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(ax, ay);
      for (let i = 1; i < n; i++) { const s = i / n, w = (i % 2 ? 7 : -7); cx.lineTo(ax + ux * L * s - uy * w, ay + uy * L * s + ux * w); }
      cx.lineTo(bx, by); cx.stroke();
      cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(ax, ay, 6, 0, Math.PI * 2); cx.fill();
    } else {
      cx.fillStyle = 'rgba(255,194,75,.18)';
      for (let gx = -6; gx <= 6; gx += 2) for (let gy = -4; gy <= 4; gy += 2) { const [px, py] = P([gx, gy]); cx.beginPath(); cx.moveTo(px - 5, py - 6); cx.lineTo(px + 5, py - 6); cx.lineTo(px, py + 6); cx.closePath(); cx.fill(); }
    }
    // trail: one dot per step
    cx.strokeStyle = 'rgba(56,225,255,.25)'; cx.lineWidth = 1.5; cx.beginPath();
    st.trail.forEach((p, i) => { const [x, y] = P(p); i ? cx.lineTo(x, y) : cx.moveTo(x, y); }); cx.stroke();
    cx.fillStyle = COL.cyan; st.trail.forEach(p => { const [x, y] = P(p); cx.beginPath(); cx.arc(x, y, 3, 0, Math.PI * 2); cx.fill(); });
    // atom + arrows
    const F = force(st.x);
    arrow(st.x, F, COL.amber, 0.9); arrow(st.x, st.v, COL.ink, 0.9);
    const [x, y] = P(st.x), r = 10 + 5 * Math.cbrt(st.m);
    const g = cx.createRadialGradient(x - r * 0.35, y - r * 0.4, r * 0.1, x, y, r);
    g.addColorStop(0, '#fff'); g.addColorStop(0.2, COL.cyan); g.addColorStop(1, '#137A96');
    cx.fillStyle = g; cx.beginPath(); cx.arc(x, y, r, 0, Math.PI * 2); cx.fill();
    // numbers
    const f = v => (Math.abs(v) < 5e-3 ? 0 : v).toFixed(2);
    o.m.textContent = st.m.toFixed(2); o.dt.textContent = st.dt.toFixed(2); o.n.textContent = st.n;
    o.Fx.textContent = f(st.F[0]); o.Fy.textContent = f(st.F[1]); o.ax.textContent = f(st.a[0]); o.ay.textContent = f(st.a[1]);
    o.vx.textContent = f(st.v[0]); o.vy.textContent = f(st.v[1]); o.px.textContent = f(st.x[0]); o.py.textContent = f(st.x[1]);
  }

  let raf = 0, last = performance.now(), visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  function loop(now) {
    raf = requestAnimationFrame(loop);
    const dtReal = Math.min(0.1, (now - last) / 1000); last = now;
    if (!visible) return;
    if (st.running) { st.acc += dtReal * 12; while (st.acc >= 1 && st.running) { st.acc -= 1; step(); } }   // 12 steps per second
    draw();
  }
  reset(); raf = requestAnimationFrame(loop);
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
