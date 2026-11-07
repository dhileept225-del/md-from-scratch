// Langevin lab: 2D Lennard-Jones fluid with the BAOAB Langevin integrator. Friction γ sets how strongly each atom is
// coupled to the bath. Shows T and its canonical spread, a tagged atom's path, and the diffusion coefficient vs γ.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RHO = 0.7, RC = 2.5, DT = 0.005, SUB = 4, HIST = 400, N = 100, TRAIL = 900;
const GAMMAS = [0, 0.1, 0.5, 1, 2, 5, 10, 20, 50];
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountLangevin(root) {
  const st = { g: 1, T0: 1.0, running: true };
  let Lx, Ly, x, xu, v, f, t = 0, hist = [], trail = [], ref = null, refT = 0, statFrom = 0, results = [];
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
  const Nf = () => (st.g > 0 ? 2 * N : 2 * N - 2);
  const temp = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return s / Nf(); };
  function restart() {
    x = new Float64Array(2 * N); xu = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = 10, a = Math.sqrt(2 / (Math.sqrt(3) * RHO)), h = a * Math.sqrt(3) / 2; Lx = n * a; Ly = n * h;
    for (let i = 0; i < N; i++) { const r = Math.floor(i / n); x[2 * i] = (i % n + 0.25 + 0.5 * (r % 2)) * a; x[2 * i + 1] = (r + 0.5) * h; v[2 * i] = gauss() * Math.sqrt(st.T0); v[2 * i + 1] = gauss() * Math.sqrt(st.T0); }
    xu.set(x); forces(); t = 0; hist = []; trail = []; results = []; markChange(2);
  }
  function markChange(wait) { statFrom = t + wait; ref = null; trail = []; }
  function step() {
    const c1 = Math.exp(-st.g * DT), c2 = Math.sqrt((1 - c1 * c1) * st.T0);
    for (let s = 0; s < SUB; s++) {
      for (let i = 0; i < 2 * N; i++) { v[i] += 0.5 * DT * f[i]; xu[i] += 0.5 * DT * v[i]; }          // B, A
      if (st.g > 0) for (let i = 0; i < 2 * N; i++) v[i] = c1 * v[i] + c2 * gauss();                     // O
      for (let i = 0; i < 2 * N; i++) { xu[i] += 0.5 * DT * v[i]; const L = i % 2 ? Ly : Lx; x[i] = xu[i] - L * Math.floor(xu[i] / L); } // A
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];                                // B
      t += DT;
    }
    hist.push({ t, T: temp() }); if (hist.length > HIST) hist.shift();
    trail.push([xu[0], xu[1]]); if (trail.length > TRAIL) trail.shift();
    if (t >= statFrom && !ref) { ref = Float64Array.from(xu); refT = t; }
  }
  function diffusion() {
    if (!ref || t - refT < 2) return null;
    let s = 0; for (let i = 0; i < 2 * N; i++) s += (xu[i] - ref[i]) ** 2;
    return s / N / (4 * (t - refT));          // 2D: MSD = 4 D t
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Friction <i class="nc">γ</i></span><output data-o="g"></output><input type="range" data-k="g" min="0" max="8" step="1" value="3"></label>
    <label>Bath temperature <output data-o="T0"></output><input type="range" data-k="T0" min="0.5" max="2" step="0.1" value="1"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button><button class="btn sm ghost" data-k="save">Add D to chart</button></div>
    <div class="readouts">
      <div><span>Average T</span><b data-o="avg"></b></div>
      <div><span>Spread of T ÷ canonical</span><b data-o="ratio"></b></div>
      <div class="wide"><span>Diffusion coefficient D (from MSD)</span><b data-o="D"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const save = () => { const D = diffusion(); if (D == null) return; results = results.filter(r => r.g !== st.g); results.push({ g: st.g, D }); };
  $('[data-k=g]').oninput = e => { save(); st.g = GAMMAS[+e.target.value]; markChange(st.g > 0 ? Math.min(3, 1 + 1 / st.g) : 1); };
  $('[data-k=T0]').oninput = e => { st.T0 = +e.target.value; results = []; markChange(st.g > 0 ? Math.min(5, 1 + 2 / st.g) : 1); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = restart; $('[data-k=save]').onclick = save;

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function draw() {
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(170, H * 0.36), S = Math.min((W - 24) / Lx, (H - plotH - 30) / Ly), bx = (W - Lx * S) / 2, by = 10;
    cx.strokeStyle = st.g > 0 ? COL.amber : 'rgba(143,163,199,.7)'; cx.lineWidth = 2; cx.strokeRect(bx, by, Lx * S, Ly * S);
    for (let i = 1; i < N; i++) { const e = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / (4 * st.T0); cx.fillStyle = e > 1 ? 'rgba(255,107,107,.55)' : e > 0.45 ? 'rgba(255,194,75,.55)' : 'rgba(56,225,255,.45)'; cx.beginPath(); cx.arc(bx + x[2 * i] * S, by + x[2 * i + 1] * S, Math.max(2, 0.4 * S), 0, 7); cx.fill(); }
    // tagged atom trail (unwrapped, folded into the box for display)
    cx.save(); cx.beginPath(); cx.rect(bx, by, Lx * S, Ly * S); cx.clip();
    cx.strokeStyle = COL.ink; cx.lineWidth = 1.6; cx.beginPath();
    const sx = Lx * Math.floor(xu[0] / Lx), sy = Ly * Math.floor(xu[1] / Ly);
    trail.forEach((p, k) => { const X = bx + (p[0] - sx) * S, Y = by + (p[1] - sy) * S; k ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke(); cx.restore();
    cx.fillStyle = COL.ink; cx.beginPath(); cx.arc(bx + x[0] * S, by + x[1] * S, Math.max(3, 0.5 * S), 0, 7); cx.fill();
    // T trace (left) and D vs gamma (right)
    const top = by + Ly * S + 34, y1 = H - 26, split = W * 0.58, x0 = 40, x1 = split - 16;
    if (hist.length > 1) {
      const hi = Math.max(2.2, ...hist.map(q => q.T)) * 1.05, Y = T => y1 - T / hi * (y1 - top), X = k => x0 + (x1 - x0) * k / (HIST - 1), sc = Math.sqrt(2 / Nf()) * st.T0;
      cx.fillStyle = 'rgba(56,225,255,.12)'; cx.fillRect(x0, Y(st.T0 + sc), x1 - x0, Y(st.T0 - sc) - Y(st.T0 + sc));
      cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, top); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
      cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [0, 1, 2].forEach(T => cx.fillText(T, x0 - 6, Y(T) + 4));
      cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.coral; cx.fillText('instant T', x0 + 6, top - 6); cx.fillStyle = COL.cyan; cx.fillText('canonical ±σ', x0 + 74, top - 6);
      cx.strokeStyle = COL.coral; cx.lineWidth = 1.5; cx.beginPath(); hist.forEach((q, k) => (k ? cx.lineTo(X(k), Y(q.T)) : cx.moveTo(X(k), Y(q.T)))); cx.stroke();
    }
    const dx0 = split + 34, dx1 = W - 12, cur = diffusion(), pts = [...results, ...(cur != null ? [{ g: st.g, D: cur, live: true }] : [])];
    const Dmax = Math.max(0.15, ...pts.map(p => p.D)) * 1.1, gx = g => dx0 + (dx1 - dx0) * (g === 0 ? 0 : (Math.log10(g) + 1.5) / 3.5), gy = D => y1 - D / Dmax * (y1 - top);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(dx0, top); cx.lineTo(dx0, y1); cx.lineTo(dx1, y1); cx.stroke();
    cx.fillStyle = COL.cyan; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('D vs γ', dx0 + 6, top - 6);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [[0, '0'], [0.1, '0.1'], [1, '1'], [10, '10']].forEach(([g, s]) => cx.fillText(s, gx(g), y1 + 14));
    pts.sort((p, q) => p.g - q.g).forEach(p => { cx.fillStyle = p.live ? COL.amber : COL.cyan; cx.beginPath(); cx.arc(gx(p.g), gy(p.D), p.live ? 6 : 5, 0, 7); cx.fill(); });
  }
  function readouts() {
    o.g.textContent = st.g === 0 ? '0 (plain NVE)' : st.g; o.T0.textContent = st.T0.toFixed(1);
    const s = hist.filter(q => q.t >= statFrom), ok = s.length > 80, m = ok ? s.reduce((a, q) => a + q.T, 0) / s.length : 0;
    const sd = ok ? Math.sqrt(s.reduce((a, q) => a + (q.T - m) ** 2, 0) / s.length) : 0, sdc = Math.sqrt(2 / Nf()) * st.T0;
    o.avg.textContent = ok ? m.toFixed(3) : 'settling…';
    o.ratio.innerHTML = !ok ? 'settling…' : `<span style="color:${Math.abs(sd / sdc - 1) < 0.25 ? COL.cyan : COL.coral}">${(sd / sdc).toFixed(2)}</span>`;
    const D = diffusion(); o.D.textContent = D == null ? 'measuring…' : D.toFixed(3) + (st.g > 0 ? `   (free particle: kT/γ = ${(st.T0 / st.g).toFixed(3)})` : '');
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) step(); draw(); if (frame++ % 8 === 0) readouts(); };
  restart(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
