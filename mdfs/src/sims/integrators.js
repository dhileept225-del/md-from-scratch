// Integrator race: explicit Euler vs Euler–Cromer (Lesson 1.3) vs velocity Verlet vs leapfrog on the same 2D spring orbit.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', green: '#7CE38B', ink: '#F2F5FA', muted: '#8FA3C7' };
const METHODS = [
  { id: 'eu', name: 'Euler', color: COL.coral, step(s, dt) { const a = [-s.x[0], -s.x[1]]; s.x = [s.x[0] + s.v[0] * dt, s.x[1] + s.v[1] * dt]; s.v = [s.v[0] + a[0] * dt, s.v[1] + a[1] * dt]; } },
  { id: 'ec', name: 'Euler–Cromer', color: COL.amber, step(s, dt) { s.v = [s.v[0] - s.x[0] * dt, s.v[1] - s.x[1] * dt]; s.x = [s.x[0] + s.v[0] * dt, s.x[1] + s.v[1] * dt]; } },
  { id: 'vv', name: 'Velocity Verlet', color: COL.cyan, step(s, dt) { s.v = [s.v[0] - 0.5 * s.x[0] * dt, s.v[1] - 0.5 * s.x[1] * dt]; s.x = [s.x[0] + s.v[0] * dt, s.x[1] + s.v[1] * dt]; s.v = [s.v[0] - 0.5 * s.x[0] * dt, s.v[1] - 0.5 * s.x[1] * dt]; } },
  // leapfrog: positions at whole steps, velocities at half steps. vm = v(t − dt/2), v = v(t + dt/2)
  { id: 'lf', name: 'Leapfrog', color: COL.green, step(s, dt) { s.x = [s.x[0] + s.v[0] * dt, s.x[1] + s.v[1] * dt]; s.vm = s.v; s.v = [s.v[0] - s.x[0] * dt, s.v[1] - s.x[1] * dt]; } },
];
const X0 = [3, 0], V0 = [0, 2.2];
const E = s => 0.5 * (s.v[0] ** 2 + s.v[1] ** 2 + s.x[0] ** 2 + s.x[1] ** 2);
// leapfrog has no velocity at time t: estimate the kinetic energy from the two half-step velocities
const ELF = (s, mode) => {
  const U = 0.5 * (s.x[0] ** 2 + s.x[1] ** 2);
  if (mode === 'vavg') { const a = (s.vm[0] + s.v[0]) / 2, b = (s.vm[1] + s.v[1]) / 2; return 0.5 * (a * a + b * b) + U; }
  return 0.25 * (s.vm[0] ** 2 + s.vm[1] ** 2 + s.v[0] ** 2 + s.v[1] ** 2) + U;   // GROMACS integrator = md
};

export function mountIntegrators(root) {
  const st = { dt: 0.1, running: true, on: { eu: true, ec: true, vv: true, lf: true }, ke: 'keavg', sys: {}, n: 0, reversed: false, lastRev: 0, stopAt: 0 };
  const E0 = 0.5 * (V0[0] ** 2 + V0[1] ** 2 + X0[0] ** 2 + X0[1] ** 2);
  function reset() {
    METHODS.forEach(m => { st.sys[m.id] = { x: [...X0], v: [...V0], trail: [[...X0]], hist: [1], gone: false }; });
    // leapfrog starts from the same x(0), v(0): v(−dt/2) and v(+dt/2) from half a kick either way
    const l = st.sys.lf, a0 = [-X0[0], -X0[1]]; l.vm = [V0[0] - 0.5 * st.dt * a0[0], V0[1] - 0.5 * st.dt * a0[1]]; l.v = [V0[0] + 0.5 * st.dt * a0[0], V0[1] + 0.5 * st.dt * a0[1]];
    st.n = 0; st.reversed = false; st.lastRev = 0; st.stopAt = 0; msg.hidden = true; revBtn.textContent = '⇄ Reverse velocities';
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view int-view"><canvas></canvas><div class="sim-msg" hidden></div><div class="sim-tip">Same spring, same Δt, four update rules · ◎ = start · leapfrog (green) sits under Verlet (blue)</div></div>
  <div class="sim-ctrl">
    ${METHODS.map(m => `<label class="check"><input type="checkbox" data-m="${m.id}" checked> <b style="color:${m.color}">${m.name}</b></label>`).join('')}
    <label>Leapfrog kinetic energy <select data-k="ke"><option value="keavg">average of the two half-step KEs (GROMACS)</option><option value="vavg">KE of the average half-step velocity</option></select></label>
    <label>Time step Δt <output data-o="dt"></output><input type="range" data-k="dt" min="0.02" max="0.4" step="0.01" value="0.1"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">↺ Reset</button></div>
    <button class="btn sm ghost" data-k="rev">⇄ Reverse velocities</button>
    <div class="readouts">${METHODS.map(m => `<div><span>${m.name} energy</span><b data-o="${m.id}"></b></div>`).join('')}<div><span>Steps</span><b data-o="n"></b></div></div>
    <p class="muted small">Run a while, then press Reverse. Only Verlet and leapfrog find their way back to ◎.</p>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.int-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const runBtn = $('[data-k=run]'), revBtn = $('[data-k=rev]'), msg = root.querySelector('.sim-msg');
  root.querySelectorAll('[data-m]').forEach(c => c.onchange = () => { st.on[c.dataset.m] = c.checked; });
  $('[data-k=dt]').oninput = e => { st.dt = +e.target.value; reset(); };
  $('[data-k=ke]').onchange = e => { st.ke = e.target.value; reset(); };
  runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = reset;
  revBtn.onclick = () => {
    METHODS.forEach(m => { const s = st.sys[m.id]; if (m.id === 'lf') { const vm = s.vm; s.vm = [-s.v[0], -s.v[1]]; s.v = [-vm[0], -vm[1]]; } else s.v = [-s.v[0], -s.v[1]]; });
    st.reversed = !st.reversed; revBtn.textContent = st.reversed ? '⇄ Reverse again' : '⇄ Reverse velocities';
    st.stopAt = st.n + (st.n - st.lastRev); st.lastRev = st.n; msg.hidden = true;
    if (!st.running) runBtn.click();
  };

  let W = 600, H = 460;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function draw() {
    cx.clearRect(0, 0, W, H);
    const oh = H * 0.66, S = Math.min(W, oh) / 16, cxp = W / 2, cyp = oh / 2 + 6, P = p => [cxp + p[0] * S, cyp - p[1] * S];
    // exact orbit (ellipse) and start marker
    cx.strokeStyle = 'rgba(242,245,250,.25)'; cx.setLineDash([6, 6]); cx.lineWidth = 1.5; cx.beginPath();
    for (let k = 0; k <= 200; k++) { const th = k / 200 * Math.PI * 2, p = P([3 * Math.cos(th), 2.2 * Math.sin(th)]); k ? cx.lineTo(...p) : cx.moveTo(...p); } cx.stroke(); cx.setLineDash([]);
    const s0 = P(X0); cx.strokeStyle = COL.ink; cx.lineWidth = 2; cx.beginPath(); cx.arc(s0[0], s0[1], 10, 0, Math.PI * 2); cx.stroke(); cx.beginPath(); cx.arc(s0[0], s0[1], 3, 0, Math.PI * 2); cx.fillStyle = COL.ink; cx.fill();
    cx.save(); cx.beginPath(); cx.rect(0, 0, W, oh); cx.clip();
    METHODS.forEach(m => {
      if (!st.on[m.id]) return; const s = st.sys[m.id];
      cx.strokeStyle = m.color; cx.globalAlpha = 0.55; cx.lineWidth = 2; cx.beginPath();
      s.trail.forEach((p, i) => { const q = P(p); i ? cx.lineTo(...q) : cx.moveTo(...q); }); cx.stroke(); cx.globalAlpha = 1;
      const q = P(s.x); cx.fillStyle = m.color; cx.beginPath(); cx.arc(q[0], q[1], 8, 0, Math.PI * 2); cx.fill();
    });
    cx.restore();
    // energy panel
    const y0 = oh + 14, y1 = H - 22, x0 = 46, x1 = W - 16, emax = 3, Y = e => y1 - Math.min(e, emax) / emax * (y1 - y0);
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.setLineDash([4, 6]); cx.beginPath(); cx.moveTo(x0, Y(1)); cx.lineTo(x1, Y(1)); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '500 12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText('E/E₀', x0 - 6, y0 + 10); cx.fillText('1', x0 - 6, Y(1) + 4); cx.fillText('time →', x1, y1 + 15);
    METHODS.forEach(m => {
      if (!st.on[m.id]) return; const h = st.sys[m.id].hist;
      cx.strokeStyle = m.color; cx.lineWidth = 2; cx.beginPath();
      h.forEach((e, i) => { const X = x0 + i / 599 * (x1 - x0); i ? cx.lineTo(X, Y(e)) : cx.moveTo(X, Y(e)); }); cx.stroke();
    });
  }
  function readouts() {
    o.dt.textContent = st.dt.toFixed(2); o.n.textContent = st.n.toLocaleString('en-US');
    METHODS.forEach(m => { const s = st.sys[m.id], d = (en(m, s) / E0 - 1) * 100; o[m.id].textContent = s.gone ? 'escaped' : `${d >= 0 ? '+' : ''}${Math.abs(d) < 1e4 ? d.toFixed(2) : d.toExponential(1)} %`; });
  }

  const en = (m, s) => m.id === 'lf' ? ELF(s, st.ke) : E(s);
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  function loop() {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running) for (let k = 0; k < 2; k++) {
      st.n++;
      METHODS.forEach(m => {
        const s = st.sys[m.id]; if (s.gone) return;
        m.step(s, st.dt); s.trail.push([...s.x]); if (s.trail.length > 900) s.trail.shift();
        if (st.n % 2 === 0) { s.hist.push(en(m, s) / E0); if (s.hist.length > 600) s.hist.shift(); }
        if (Math.hypot(...s.x) > 40) s.gone = true;
      });
      if (st.stopAt && st.n >= st.stopAt) {   // ran the same number of steps backwards: who made it home?
        st.stopAt = 0; st.running = false; runBtn.textContent = '▶ Run';
        const d = m => { const s = st.sys[m.id]; if (s.gone) return 'escaped, never coming back'; const r = Math.hypot(s.x[0] - X0[0], s.x[1] - X0[1]); return r < 1e-6 ? `back home (off by ${r.toExponential(0)})` : `off by ${r.toFixed(3)}`; };
        msg.hidden = false;
        msg.innerHTML = '<b>Same number of steps backwards. Who is back at ◎?</b><br>' + METHODS.map(m => `<span style="color:${m.color}">${m.name}</span>: ${d(m)}`).join('<br>');
        break;
      }
    }
    draw(); readouts();
  }
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
