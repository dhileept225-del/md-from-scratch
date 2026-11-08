// Van Hove lab. Two systems: a live Lennard-Jones liquid (256 atoms, 3D) and a hopping model
// (particles rattle around sites and jump a distance a = 1 after random waiting times).
// Shows the self part 4 pi r^2 G_s(r, t) with the Gaussian of the same MSD, the non-Gaussian parameter,
// and (liquid) the distinct part G_d(r, t) / rho, all for a lag time t you choose.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 256, DT = 0.005, RC2 = 6.25, EVERY = 10, H = EVERY * DT, K = 200, RHO = 0.8442, TEMP = 0.722;
const M = 600, HOPDT = 0.05, TAGS = [3, 77, 150];

export function mountVanhove(root) {
  const st = { sys: 'lj', lag: 1.0, wait: 1.0, rattle: 0.1, running: true };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="lj" class="on">LJ liquid</button><button data-s="hop">Hopping model</button></div>
    <label><span>Lag time <i class="nc">t</i></span><output data-o="lag"></output><input type="range" data-k="lag" min="0.1" max="5" step="0.05" value="1.0"></label>
    <label data-hop hidden><span>Mean waiting time between hops</span><output data-o="wait"></output><input type="range" data-k="wait" min="0.3" max="3" step="0.1" value="1.0"></label>
    <label data-hop hidden><span>Rattling amplitude</span><output data-o="rat"></output><input type="range" data-k="rat" min="0.03" max="0.3" step="0.01" value="0.1"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Restart</button></div>
    <div class="readouts">
      <div><span>MSD at <i class="nc">t</i></span><b data-o="msd"></b></div>
      <div><span>Non-Gaussian <i class="nc">α</i>₂</span><b data-o="a2"></b></div>
      <div><span>Moved less than 0.5</span><b data-o="stay"></b></div>
      <div><span data-o="lab4">Shell peak now</span><b data-o="gd"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 99; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  // ---- LJ liquid ----
  let L = 1, x = new Float64Array(3 * N), u = new Float64Array(3 * N), v = new Float64Array(3 * N), f = new Float64Array(3 * N), steps = 0;
  // ---- hopping model ----
  let site = new Float64Array(3 * M), hp = new Float64Array(3 * M);
  let buf = [], tracks = [];
  const nP = () => st.sys === 'lj' ? N : M;
  const frameDt = () => st.sys === 'lj' ? H : HOPDT;
  const lagK = () => Math.max(1, Math.min(K / 2, Math.round(st.lag / frameDt())));   // keep at least K/2 time origins

  function build() {
    buf = []; tracks = TAGS.map(() => []); steps = 0;
    if (st.sys === 'lj') {
      const a = Math.cbrt(4 / RHO), B = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]; L = 4 * a; let k = 0;
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let l = 0; l < 4; l++) for (const b of B) { x[3 * k] = (i + b[0]) * a; x[3 * k + 1] = (j + b[1]) * a; x[3 * k + 2] = (l + b[2]) * a; k++; }
      const m = [0, 0, 0]; for (let i = 0; i < 3 * N; i++) { v[i] = gauss() * Math.sqrt(TEMP); m[i % 3] += v[i] / N; } for (let i = 0; i < 3 * N; i++) v[i] -= m[i % 3];
      forces(); u.set(x);
    } else { site.fill(0); for (let i = 0; i < 3 * M; i++) hp[i] = st.rattle * gauss(); }
  }
  function forces() {
    f.fill(0);
    for (let i = 0; i < N - 1; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
      const r2 = dx * dx + dy * dy + dz * dz; if (r2 > RC2) continue;
      const ir6 = 1 / (r2 * r2 * r2), ff = 24 * (2 * ir6 * ir6 - ir6) / r2;
      f[3 * i] += ff * dx; f[3 * i + 1] += ff * dy; f[3 * i + 2] += ff * dz; f[3 * j] -= ff * dx; f[3 * j + 1] -= ff * dy; f[3 * j + 2] -= ff * dz;
    }
  }
  function ljStep() {
    for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; u[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
    forces(); let ke = 0; for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; ke += v[i] * v[i]; }
    const Tn = ke / (3 * N - 3), lam = Math.sqrt(1 + DT / (steps < 600 ? 0.1 : 2.0) * (TEMP / Tn - 1)); for (let i = 0; i < 3 * N; i++) v[i] *= lam;
    steps++; if (steps === 600) { buf = []; u.set(x); tracks = TAGS.map(() => []); }
    if (steps > 600 && steps % EVERY === 0) push(u, x);
  }
  function hopStep() {   // one frame of the hopping model
    for (let i = 0; i < M; i++) {
      if (rand() < HOPDT / st.wait) { let a = gauss(), b = gauss(), c = gauss(); const n = Math.hypot(a, b, c); site[3 * i] += a / n; site[3 * i + 1] += b / n; site[3 * i + 2] += c / n; }
      for (let d = 0; d < 3; d++) hp[3 * i + d] = site[3 * i + d] + st.rattle * gauss();
    }
    push(hp, null);
  }
  function push(U, X) {
    buf.push({ u: Float64Array.from(U), x: X ? Float64Array.from(X) : null }); if (buf.length > K + 1) buf.shift();
    TAGS.forEach((t, j) => { tracks[j].push([U[3 * t], U[3 * t + 1]]); if (tracks[j].length > 400) tracks[j].shift(); });
  }
  function stats() {   // self part from all origins in the buffer, distinct part (liquid) from up to 12 origins
    const k = lagK(), n = nP(), nb = 80, rmax = st.sys === 'lj' ? Math.max(1.5, Math.min(3, 4 * Math.sqrt(st.lag * 0.18 + 0.02))) : 3.2, dr = rmax / nb;
    const hs = new Float64Array(nb); let c = 0, s2 = 0, s4 = 0, stay = 0;
    if (buf.length <= k) return null;
    for (let o = 0; o + k < buf.length; o++) { const A = buf[o].u, B = buf[o + k].u;
      for (let i = 0; i < n; i++) { const dx = B[3 * i] - A[3 * i], dy = B[3 * i + 1] - A[3 * i + 1], dz = B[3 * i + 2] - A[3 * i + 2], r2 = dx * dx + dy * dy + dz * dz, r = Math.sqrt(r2);
        s2 += r2; s4 += r2 * r2; c++; if (r < 0.5) stay++; const b = Math.floor(r / dr); if (b < nb) hs[b]++; } }
    const msd = s2 / c, a2 = 3 * (s4 / c) / (5 * msd * msd) - 1; for (let b = 0; b < nb; b++) hs[b] /= c * dr;
    let gd = null, gr = 3.3;
    if (st.sys === 'lj') {
      const nbd = 66, drd = gr / nbd, hd = new Float64Array(nbd), origins = []; for (let o = 0; o + k < buf.length && origins.length < 12; o += Math.max(1, Math.floor((buf.length - k) / 12))) origins.push(o);
      for (const o of origins) { const A = buf[o].x, B = buf[o + k].u, A0 = buf[o].u;
        for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { if (i === j) continue;
          // position of j at t0 + t (unwrapped -> wrapped consistently) relative to i at t0
          let dx = (A[3 * j] + B[3 * j] - A0[3 * j]) - A[3 * i], dy = (A[3 * j + 1] + B[3 * j + 1] - A0[3 * j + 1]) - A[3 * i + 1], dz = (A[3 * j + 2] + B[3 * j + 2] - A0[3 * j + 2]) - A[3 * i + 2];
          dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
          const r = Math.sqrt(dx * dx + dy * dy + dz * dz), b = Math.floor(r / drd); if (b < nbd) hd[b]++; } }
      gd = Array.from(hd, (h, b) => h / (origins.length * N * RHO * 4 / 3 * Math.PI * (((b + 1) * drd) ** 3 - (b * drd) ** 3))); gd.dr = drd;
    }
    return { hs, dr, rmax, msd, a2, stay: stay / c, gd, origins: buf.length - k };
  }

  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { st.sys = b.dataset.s; root.querySelectorAll('[data-s]').forEach(z => z.classList.toggle('on', z === b)); root.querySelectorAll('[data-hop]').forEach(e => e.hidden = st.sys !== 'hop'); build(); });
  $('[data-k=lag]').oninput = e => { st.lag = +e.target.value; draw(); };
  $('[data-k=wait]').oninput = e => { st.wait = +e.target.value; buf = []; };
  $('[data-k=rat]').oninput = e => { st.rattle = +e.target.value; buf = []; };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => build();

  let W = 600, Hc = 600, cache = null, frame = 0;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; Hc = view.clientHeight; cv.width = W * d; cv.height = Hc * d; cx.setTransform(d, 0, 0, d, 0, 0); cache = null; draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function axesP(yt, yb, xmax, ymax, title, xlab, ticks) {
    const x0 = 44, x1 = W - 14, X = r => x0 + r / xmax * (x1 - x0), Y = q => yb - Math.min(q, ymax * 1.04) / ymax * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(title, x0 + 6, yt - 8);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText(xlab, x1, yb + 14);
    cx.textAlign = 'center'; ticks.filter(t => t < xmax * 0.8).forEach(t => cx.fillText(t, X(t), yb + 14));
    return { X, Y };
  }
  function draw() {
    if (!W) return;
    if (!cache || frame % 4 === 0) cache = stats();
    const S = cache; cx.clearRect(0, 0, W, Hc);
    // tracks
    const topH = Hc * 0.26, tcol = [COL.cyan, COL.amber, COL.coral];
    let span = 0.8; tracks.forEach(tr => tr.forEach(p => { span = Math.max(span, Math.abs(p[0] - tr[0][0]) * 1.1, Math.abs(p[1] - tr[0][1]) * 1.1); }));
    const sc = Math.min(W / 3, topH) / 2 / span;
    tracks.forEach((tr, j) => { if (tr.length < 2) return; const ox = W * (j + 0.5) / 3, oy = topH / 2 + 8;
      cx.strokeStyle = tcol[j]; cx.lineWidth = 1.3; cx.beginPath(); tr.forEach((p, q) => { const X = ox + (p[0] - tr[0][0]) * sc, Y = oy - (p[1] - tr[0][1]) * sc; q ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke(); });
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(st.sys === 'lj' ? 'three atoms (x–y), smooth wandering' : 'three particles (x–y): rattle, then jump', 10, 16);
    const two = st.sys === 'lj', yb1 = two ? Hc * 0.62 : Hc - 26;
    if (!S) { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '13px "Space Grotesk", sans-serif'; cx.fillText(st.sys === 'lj' ? 'equilibrating / filling the time window…' : 'filling the time window…', W / 2, (topH + yb1) / 2); }
    else {
      const ymax = Math.max(...S.hs) * 1.15, P = axesP(topH + 30, yb1, S.rmax, ymax, `4πr² G_s(r, t) at t = ${st.lag.toFixed(2)} τ`, 'distance moved r', S.rmax > 2 ? [1, 2, 3] : [0.5, 1, 1.5]);
      cx.setLineDash([6, 4]); cx.strokeStyle = COL.amber; cx.lineWidth = 1.5; cx.beginPath();
      for (let b = 0; b <= 120; b++) { const r = b / 120 * S.rmax, g = 4 * Math.PI * r * r * (3 / (2 * Math.PI * S.msd)) ** 1.5 * Math.exp(-3 * r * r / (2 * S.msd)); b ? cx.lineTo(P.X(r), P.Y(g)) : cx.moveTo(P.X(r), P.Y(g)); } cx.stroke(); cx.setLineDash([]);
      cx.strokeStyle = COL.cyan; cx.lineWidth = 2.5; cx.beginPath(); S.hs.forEach((h, b) => { const r = (b + 0.5) * S.dr; b ? cx.lineTo(P.X(r), P.Y(h)) : cx.moveTo(P.X(r), P.Y(h)); }); cx.stroke();
      cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillStyle = COL.cyan; cx.fillText('simulation', W - 16, topH + 36); cx.fillStyle = COL.amber; cx.fillText('dashed: Gaussian, same MSD', W - 16, topH + 51);
      if (two && S.gd) {
        const P2 = axesP(yb1 + 44, Hc - 26, 3.3, 3.4, `G_d(r, t) / ρ  (t = 0 would be g(r))`, 'r (σ)', [1, 2, 3]);
        cx.setLineDash([4, 4]); cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.beginPath(); cx.moveTo(P2.X(0), P2.Y(1)); cx.lineTo(P2.X(3.3), P2.Y(1)); cx.stroke(); cx.setLineDash([]);
        cx.strokeStyle = COL.coral; cx.lineWidth = 2.2; cx.beginPath(); S.gd.forEach((g, b) => { const r = (b + 0.5) * S.gd.dr; b ? cx.lineTo(P2.X(r), P2.Y(g)) : cx.moveTo(P2.X(r), P2.Y(g)); }); cx.stroke();
      }
    }
    o.lag.textContent = st.lag.toFixed(2) + ' τ'; o.wait.textContent = st.wait.toFixed(1) + ' τ'; o.rat.textContent = st.rattle.toFixed(2);
    o.lab4.textContent = st.sys === 'lj' ? 'Shell peak now (3.0 at t = 0)' : 'Expected not hopped, e^(−t/wait)';
    if (S) {
      o.msd.textContent = S.msd.toFixed(3); o.a2.textContent = S.a2.toFixed(3) + (S.a2 > 0.3 ? ' (far from Gaussian)' : ''); o.stay.textContent = (100 * S.stay).toFixed(0) + ' %';
      o.gd.textContent = st.sys === 'lj' ? (S.gd ? Math.max(...S.gd.slice(Math.floor(0.8 / S.gd.dr), Math.floor(1.5 / S.gd.dr))).toFixed(2) : '…') : (100 * Math.exp(-st.lag / st.wait)).toFixed(0) + ' %';
    } else o.msd.textContent = o.a2.textContent = o.stay.textContent = o.gd.textContent = '…';
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) { if (st.sys === 'lj') for (let k = 0; k < 6; k++) ljStep(); else hopStep(); } frame++; draw(); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
