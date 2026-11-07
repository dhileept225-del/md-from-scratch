// Nosé–Hoover lab. Liquid: 2D Lennard-Jones fluid with a Nosé–Hoover chain (length 1 or 3): T and ξ traces, a kick
// to show ringing, the extended-energy check. Oscillator: one harmonic oscillator, phase-space density for NH vs NHC.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RHO = 0.7, RC = 2.5, DT = 0.005, SUB = 4, HIST = 500, N = 100, WARM = 1.5;
const TAUS = [0.05, 0.1, 0.2, 0.5, 1, 2];
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

// Nosé–Hoover chain: integrate thermostat variables for time h and return the velocity scale factor.
// G0 uses the system kinetic energy (2K vs Nf kT); higher links thermostat the one below.
function chainHalf(c, twoK, h) {
  const M = c.xi.length, G = j => (j === 0 ? twoK - c.Nf * c.kT : c.Q[j - 1] * c.xi[j - 1] ** 2 - c.kT) / c.Q[j];
  for (let j = M - 1; j >= 0; j--) { if (j < M - 1) c.xi[j] *= Math.exp(-c.xi[j + 1] * h / 4); c.xi[j] += G(j) * h / 2; if (j < M - 1) c.xi[j] *= Math.exp(-c.xi[j + 1] * h / 4); }
  const s = Math.exp(-c.xi[0] * h); twoK *= s * s;
  for (let j = 0; j < M; j++) c.eta[j] += c.xi[j] * h;
  for (let j = 0; j < M; j++) { if (j < M - 1) c.xi[j] *= Math.exp(-c.xi[j + 1] * h / 4); c.xi[j] += G(j) * h / 2; if (j < M - 1) c.xi[j] *= Math.exp(-c.xi[j + 1] * h / 4); }
  return s;
}
const makeChain = (M, Nf, kT, tau) => ({ Nf, kT, xi: new Float64Array(M), eta: new Float64Array(M), Q: Array.from({ length: M }, (_, j) => (j === 0 ? Nf : 1) * kT * tau * tau) });
const chainEnergy = c => c.Q.reduce((s, q, j) => s + 0.5 * q * c.xi[j] ** 2, 0) + c.Nf * c.kT * c.eta[0] + c.kT * Array.from(c.eta).slice(1).reduce((s, e) => s + e, 0);

export function mountNoseHoover(root) {
  const st = { sys: 'liquid', M: 1, tau: 0.2, T0: 1.0, running: true };
  // ---- liquid ----
  const Nf = 2 * N - 2; let Lx, Ly, x, v, f, U = 0, t = 0, ch, hist = [], H0 = null, statFrom = 0;
  function forces() {
    f.fill(0); U = 0; const rc2 = RC * RC, ush = 4 * (Math.pow(RC, -12) - Math.pow(RC, -6));
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[2 * i] - x[2 * j], dy = x[2 * i + 1] - x[2 * j + 1];
      dx -= Lx * Math.round(dx / Lx); dy -= Ly * Math.round(dy / Ly);
      const r2 = dx * dx + dy * dy; if (r2 >= rc2) continue;
      const i6 = 1 / (r2 * r2 * r2), ff = 24 * (2 * i6 * i6 - i6) / r2; U += 4 * (i6 * i6 - i6) - ush;
      f[2 * i] += ff * dx; f[2 * i + 1] += ff * dy; f[2 * j] -= ff * dx; f[2 * j + 1] -= ff * dy;
    }
  }
  const twoK = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return s; };
  const newChain = () => { ch = makeChain(st.M, Nf, st.T0, st.tau); H0 = null; };
  function restartLiquid() {
    x = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = 10, a = Math.sqrt(2 / (Math.sqrt(3) * RHO)), h = a * Math.sqrt(3) / 2; Lx = n * a; Ly = n * h;
    for (let i = 0; i < N; i++) { const r = Math.floor(i / n); x[2 * i] = (i % n + 0.25 + 0.5 * (r % 2)) * a; x[2 * i + 1] = (r + 0.5) * h; v[2 * i] = gauss() * 0.5; v[2 * i + 1] = gauss() * 0.5; }
    let px = 0, py = 0; for (let i = 0; i < N; i++) { px += v[2 * i]; py += v[2 * i + 1]; }
    for (let i = 0; i < N; i++) { v[2 * i] -= px / N; v[2 * i + 1] -= py / N; }
    forces(); t = 0; hist = []; newChain(); statFrom = WARM + 3;
  }
  function stepLiquid() {
    for (let s = 0; s < SUB; s++) {
      const warm = t < WARM;
      if (!warm) { const sc = chainHalf(ch, twoK(), DT / 2); for (let i = 0; i < 2 * N; i++) v[i] *= sc; }
      for (let i = 0; i < 2 * N; i++) { const L = i % 2 ? Ly : Lx; v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];
      if (warm) { const lam = Math.sqrt(1 + DT / 0.05 * (Nf * st.T0 / twoK() - 1)); for (let i = 0; i < 2 * N; i++) v[i] *= lam; }
      else { const sc = chainHalf(ch, twoK(), DT / 2); for (let i = 0; i < 2 * N; i++) v[i] *= sc; }
      t += DT;
    }
    const K2 = twoK(), Hext = 0.5 * K2 + U + chainEnergy(ch);
    if (t >= WARM && H0 === null) H0 = Hext;
    hist.push({ t, T: K2 / Nf, xi: t < WARM ? 0 : ch.xi[0], warm: t < WARM, H: Hext }); if (hist.length > HIST) hist.shift();
  }
  // ---- oscillator ----
  const GRID = 120, RANGE = 3.6; let ox, op, och, dens, nS, m2, m4;
  function restartOsc() { ox = 1; op = 0; och = makeChain(st.M, 1, 1, st.tau); dens = new Float32Array(GRID * GRID); nS = 0; m2 = 0; m4 = 0; }
  function stepOsc() {
    const dt = 0.02;
    for (let s = 0; s < 400; s++) {
      op *= chainHalf(och, op * op, dt / 2); op -= 0.5 * dt * ox; ox += dt * op; op -= 0.5 * dt * ox; op *= chainHalf(och, op * op, dt / 2);
      const gx = Math.floor((ox / RANGE + 1) / 2 * GRID), gp = Math.floor((1 - op / RANGE) / 2 * GRID);
      if (gx >= 0 && gx < GRID && gp >= 0 && gp < GRID) dens[gp * GRID + gx]++;
      nS++; m2 += ox * ox; m4 += ox ** 4;
    }
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="liquid" class="on">Liquid (100 atoms)</button><button data-s="osc">One oscillator</button></div>
    <div class="seg"><button data-m="1" class="on">Nosé–Hoover</button><button data-m="3">NH chain (3)</button></div>
    <label><span>Thermostat period <i class="nc">τ</i></span><output data-o="tau"></output><input type="range" data-k="tau" min="0" max="5" step="1" value="2"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="restart">Restart</button><button class="btn sm ghost" data-k="kick" data-liq>Kick: v × 1.5</button></div>
    <div class="readouts" data-liq>
      <div><span>Average T (target 1.0)</span><b data-o="avg"></b></div>
      <div><span>Spread of T ÷ canonical</span><b data-o="ratio"></b></div>
      <div class="wide"><span>Extended energy drift (per atom)</span><b data-o="drift"></b></div>
    </div>
    <div class="readouts" data-osc hidden>
      <div><span><i class="nc">⟨x²⟩</i> (canonical: 1)</span><b data-o="x2"></b></div>
      <div><span><i class="nc">⟨x⁴⟩/3⟨x²⟩²</i> (canonical: 1)</span><b data-o="kurt"></b></div>
      <div class="wide"><span>Verdict</span><b data-o="verdict"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const restart = () => (st.sys === 'liquid' ? restartLiquid() : restartOsc());
  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => {
    st.sys = b.dataset.s; root.querySelectorAll('[data-s]').forEach(z => z.classList.toggle('on', z === b));
    root.querySelectorAll('[data-liq]').forEach(e => e.hidden = st.sys !== 'liquid'); root.querySelectorAll('[data-osc]').forEach(e => e.hidden = st.sys === 'liquid'); restart();
  });
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => {
    st.M = +b.dataset.m; root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b));
    if (st.sys === 'liquid') { newChain(); statFrom = t + 3; } else restartOsc();
  });
  $('[data-k=tau]').oninput = e => { st.tau = TAUS[+e.target.value]; if (st.sys === 'liquid') { newChain(); statFrom = t + 3; } else restartOsc(); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=restart]').onclick = restart;
  $('[data-k=kick]').onclick = () => { for (let i = 0; i < 2 * N; i++) v[i] *= 1.5; H0 = null; newChain(); statFrom = t + 1e9; hist.push({ ...hist[hist.length - 1], kick: true }); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const off = document.createElement('canvas'); off.width = off.height = GRID; const octx = off.getContext('2d'), img = octx.createImageData(GRID, GRID);

  function axisBox(x0, y0, x1, y1) { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); }
  function drawLiquid() {
    const plotH = Math.max(200, H * 0.45), S = Math.min((W - 24) / Lx, (H - plotH - 30) / Ly), bx = (W - Lx * S) / 2, by = 10;
    cx.strokeStyle = t < WARM ? COL.muted : COL.amber; cx.lineWidth = 2.5; cx.strokeRect(bx, by, Lx * S, Ly * S);
    for (let i = 0; i < N; i++) { const e = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / (4 * st.T0); cx.fillStyle = e > 1 ? COL.coral : e > 0.45 ? COL.amber : COL.cyan; cx.beginPath(); cx.arc(bx + x[2 * i] * S, by + x[2 * i + 1] * S, Math.max(2, 0.4 * S), 0, 7); cx.fill(); }
    if (t < WARM) { cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText('warm-up with Berendsen…', W / 2, by + Ly * S / 2); }
    if (hist.length < 2) return;
    const x0 = 46, x1 = W - 12, top = by + Ly * S + 34, mid = top + (H - top) * 0.62, X = k => x0 + (x1 - x0) * k / (HIST - 1);
    // T panel
    const ty1 = mid - 16, hiT = Math.max(2.2, ...hist.map(q => q.T)) * 1.05, YT = T => ty1 - T / hiT * (ty1 - top);
    const sc = Math.sqrt(2 / Nf) * st.T0; cx.fillStyle = 'rgba(56,225,255,.12)'; cx.fillRect(x0, YT(st.T0 + sc), x1 - x0, YT(st.T0 - sc) - YT(st.T0 + sc));
    axisBox(x0, top, x1, ty1); cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [0, 1, 2].forEach(T => cx.fillText(T.toFixed(0), x0 - 6, YT(T) + 4));
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.coral; cx.fillText('instant T', x0 + 8, top - 6); cx.fillStyle = COL.cyan; cx.fillText('canonical ±σ', x0 + 80, top - 6);
    hist.forEach((q, k) => { if (q.warm) { cx.fillStyle = 'rgba(143,163,199,.08)'; cx.fillRect(X(k), top, (x1 - x0) / HIST + 1, ty1 - top); } if (q.kick) { cx.strokeStyle = COL.coral; cx.beginPath(); cx.moveTo(X(k), top); cx.lineTo(X(k), ty1); cx.stroke(); } });
    cx.strokeStyle = COL.coral; cx.lineWidth = 1.5; cx.beginPath(); hist.forEach((q, k) => (k ? cx.lineTo(X(k), YT(q.T)) : cx.moveTo(X(k), YT(q.T)))); cx.stroke();
    // xi panel
    const zy0 = mid + 14, zy1 = H - 14, lim = Math.max(0.5, ...hist.map(q => Math.abs(q.xi))) * 1.1, YZ = z => (zy0 + zy1) / 2 - z / lim * (zy1 - zy0) / 2;
    axisBox(x0, zy0, x1, zy1); cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.beginPath(); cx.moveTo(x0, YZ(0)); cx.lineTo(x1, YZ(0)); cx.stroke();
    cx.fillStyle = COL.amber; cx.textAlign = 'left'; cx.fillText('friction ξ (above 0: braking, below: pushing)', x0 + 8, zy0 - 2);
    cx.strokeStyle = COL.amber; cx.lineWidth = 1.5; cx.beginPath(); hist.forEach((q, k) => (k ? cx.lineTo(X(k), YZ(q.xi)) : cx.moveTo(X(k), YZ(q.xi)))); cx.stroke();
  }
  function drawOsc() {
    const side = Math.min(W - 40, H * 0.58), sx = (W - side) / 2, sy = 24;
    let mx = 0; for (let k = 0; k < dens.length; k++) mx = Math.max(mx, dens[k]);
    for (let k = 0; k < dens.length; k++) { const q = mx ? Math.sqrt(dens[k] / mx) : 0, c = st.M === 1 ? [255, 107, 107] : [56, 225, 255]; img.data[4 * k] = c[0]; img.data[4 * k + 1] = c[1]; img.data[4 * k + 2] = c[2]; img.data[4 * k + 3] = Math.round(255 * q); }
    octx.putImageData(img, 0, 0); cx.imageSmoothingEnabled = true; cx.drawImage(off, sx, sy, side, side);
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(sx, sy + side / 2); cx.lineTo(sx + side, sy + side / 2); cx.moveTo(sx + side / 2, sy); cx.lineTo(sx + side / 2, sy + side); cx.stroke();
    cx.setLineDash([4, 4]); [1, 2].forEach(r => { cx.beginPath(); cx.arc(sx + side / 2, sy + side / 2, r / RANGE * side / 2, 0, 7); cx.stroke(); }); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText('position x', sx + side, sy + side / 2 - 6); cx.textAlign = 'left'; cx.fillText('momentum p', sx + side / 2 + 6, sy + 12);
    cx.fillText('dashed: canonical 1σ and 2σ', sx, sy + side + 16);
    // histogram of x vs Gaussian
    const x0 = 30, x1 = W - 30, y1 = H - 16, y0 = sy + side + 40, nb = 40, h = new Float64Array(nb);
    for (let gx = 0; gx < GRID; gx++) { let s = 0; for (let gp = 0; gp < GRID; gp++) s += dens[gp * GRID + gx]; h[Math.floor(gx / GRID * nb)] += s; }
    const tot = h.reduce((s, q) => s + q, 0) || 1, bw = 2 * RANGE / nb, gmax = 1 / Math.sqrt(2 * Math.PI), sc = (y1 - y0) / (gmax * 1.6);
    axisBox(x0, y0, x1, y1); const wpx = (x1 - x0) / nb;
    cx.fillStyle = st.M === 1 ? 'rgba(255,107,107,.7)' : 'rgba(56,225,255,.7)';
    h.forEach((q, k) => { const hh = q / tot / bw * sc; cx.fillRect(x0 + k * wpx + 1, y1 - hh, wpx - 2, hh); });
    cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.setLineDash([6, 4]); cx.beginPath();
    for (let s = 0; s <= 1; s += 0.01) { const xx = -RANGE + 2 * RANGE * s, Y = y1 - Math.exp(-xx * xx / 2) * gmax * sc; s ? cx.lineTo(x0 + s * (x1 - x0), Y) : cx.moveTo(x0, Y); } cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.amber; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('canonical (Gaussian)', x0 + 8, y0 + 4);
  }
  function readouts() {
    o.tau.textContent = st.tau;
    if (st.sys === 'liquid') {
      const s = hist.filter(q => q.t >= statFrom), ok = s.length > 80, m = ok ? s.reduce((a, q) => a + q.T, 0) / s.length : 0;
      const sd = ok ? Math.sqrt(s.reduce((a, q) => a + (q.T - m) ** 2, 0) / s.length) : 0, sdc = Math.sqrt(2 / Nf) * st.T0;
      o.avg.textContent = statFrom > 1e8 ? 'after a kick: watch it ring' : ok ? m.toFixed(3) : 'settling…';
      o.ratio.innerHTML = statFrom > 1e8 || !ok ? '–' : `<span style="color:${Math.abs(sd / sdc - 1) < 0.3 ? COL.cyan : COL.coral}">${(sd / sdc).toFixed(2)}</span>`;
      o.drift.textContent = H0 === null ? '–' : ((hist[hist.length - 1].H - H0) / N).toExponential(1);
    } else {
      if (nS < 2000) { o.x2.textContent = o.kurt.textContent = o.verdict.textContent = 'collecting…'; return; }
      const x2 = m2 / nS, k = m4 / nS / (3 * x2 * x2), good = Math.abs(x2 - 1) < 0.1 && Math.abs(k - 1) < 0.1;
      o.x2.textContent = x2.toFixed(3); o.kurt.textContent = k.toFixed(3);
      o.verdict.innerHTML = `<span style="color:${good ? COL.cyan : COL.coral}">${good ? 'canonical: the thermostat explores all of phase space' : 'not canonical: a single NH stays on a narrow band'}</span>`;
    }
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running) st.sys === 'liquid' ? stepLiquid() : stepOsc();
    cx.clearRect(0, 0, W, H); st.sys === 'liquid' ? drawLiquid() : drawOsc(); if (frame++ % 8 === 0) readouts();
  };
  restartLiquid(); restartOsc(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
