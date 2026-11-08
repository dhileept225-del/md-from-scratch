// MSD lab: 256 Lennard-Jones atoms in 3D (velocity Verlet, gentle Berendsen thermostat).
// Mean squared displacement averaged over atoms and all time origins, built up live from unwrapped
// coordinates (optionally also from wrapped ones, to show the classic mistake). D = slope / 6.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 256, DT = 0.005, RC2 = 6.25, EVERY = 10, K = 100, TAGS = [3, 77, 150];   // frame every 0.05 tau, lags up to 5 tau
const STATES = { gas: { rho: 0.05, T: 2.0 }, liquid: { rho: 0.8442, T: 0.722 }, solid: { rho: 1.1, T: 0.3 } };

export function mountMsd(root) {
  const st = { state: 'liquid', T: 0.722, running: true, wrapped: false, log: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="gas">Gas</button><button data-s="liquid" class="on">Liquid</button><button data-s="solid">Solid</button></div>
    <label><span>Temperature <i class="nc">T</i>*</span><output data-o="T0"></output><input type="range" data-k="T" min="0.2" max="2.5" step="0.05" value="0.722"></label>
    <div class="seg"><button data-ax="lin" class="on">Linear axes</button><button data-ax="log">Log–log</button></div>
    <div class="seg"><button data-w="0" class="on">Unwrapped</button><button data-w="1">+ wrapped (wrong)</button></div>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset average</button></div>
    <div class="readouts">
      <div><span>Diffusion coefficient <i class="nc">D</i>*</span><b data-o="D"></b></div>
      <div><span>Log–log slope (end)</span><b data-o="sl"></b></div>
      <div><span>Simulated time / origins</span><b data-o="nf"></b></div>
      <div><span>Measured <i class="nc">T</i>*</span><b data-o="T"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let L = 1, x = new Float64Array(3 * N), u = new Float64Array(3 * N), v = new Float64Array(3 * N), f = new Float64Array(3 * N), Tnow = 0, steps = 0;
  let bufU = [], bufW = [], sumU, sumW, cnt, nfr = 0, tracks = [];
  let rng = 4242; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());

  function build() {
    const s = STATES[st.state], a = Math.cbrt(4 / s.rho), B = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]; L = 4 * a; let k = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let l = 0; l < 4; l++) for (const b of B) { x[3 * k] = (i + b[0]) * a; x[3 * k + 1] = (j + b[1]) * a; x[3 * k + 2] = (l + b[2]) * a; k++; }
    const m = [0, 0, 0]; for (let i = 0; i < 3 * N; i++) { v[i] = gauss() * Math.sqrt(st.T); m[i % 3] += v[i] / N; } for (let i = 0; i < 3 * N; i++) v[i] -= m[i % 3];
    forces(); steps = 0; resetAvg();
  }
  function resetAvg() { bufU = []; bufW = []; sumU = new Float64Array(K + 1); sumW = new Float64Array(K + 1); cnt = new Float64Array(K + 1); nfr = 0; u.set(x); tracks = TAGS.map(() => []); }
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
  function step() {
    for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; u[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }   // u: unwrapped copy
    forces(); let ke = 0;
    for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; ke += v[i] * v[i]; }
    Tnow = ke / (3 * N - 3); const lam = Math.sqrt(1 + DT / (steps < 300 ? 0.1 : 2.0) * (st.T / Tnow - 1)); for (let i = 0; i < 3 * N; i++) v[i] *= lam;   // strong coupling while equilibrating, then weak (tau = 2) so the dynamics are barely perturbed
    steps++;
    if (steps === 300) resetAvg();                                    // equilibration done
    if (steps > 300 && steps % EVERY === 0) sample();
  }
  function sample() {   // add the new frame as the END point of every lag 1..K (all time origins)
    bufU.push(Float64Array.from(u)); bufW.push(Float64Array.from(x)); if (bufU.length > K + 1) { bufU.shift(); bufW.shift(); }
    const n = bufU.length - 1, cu = bufU[n], cw = bufW[n];
    for (let k = 1; k <= n; k++) {
      const pu = bufU[n - k], pw = bufW[n - k]; let a = 0, b = 0;
      for (let i = 0; i < 3 * N; i++) { const d1 = cu[i] - pu[i], d2 = cw[i] - pw[i]; a += d1 * d1; b += d2 * d2; }
      sumU[k] += a / N; sumW[k] += b / N; cnt[k]++;
    }
    nfr++; TAGS.forEach((t, j) => { tracks[j].push([u[3 * t], u[3 * t + 1]]); if (tracks[j].length > 600) tracks[j].shift(); });
  }
  const msd = (s) => Array.from({ length: K + 1 }, (_, k) => k && cnt[k] ? s[k] / cnt[k] : (k ? null : 0));
  function fitD(m) {   // least squares over 1 <= t <= 5 tau (lags 20..100), through the data (not the origin)
    let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let k = 20; k <= K; k++) { if (m[k] == null) continue; const t = k * EVERY * DT; n++; sx += t; sy += m[k]; sxx += t * t; sxy += t * m[k]; }
    return n > 10 ? (n * sxy - sx * sy) / (n * sxx - sx * sx) / 6 : null;
  }

  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { st.state = b.dataset.s; st.T = STATES[st.state].T; $('[data-k=T]').value = st.T; root.querySelectorAll('[data-s]').forEach(z => z.classList.toggle('on', z === b)); build(); });
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; steps = 0; };
  root.querySelectorAll('[data-ax]').forEach(b => b.onclick = () => { st.log = b.dataset.ax === 'log'; root.querySelectorAll('[data-ax]').forEach(z => z.classList.toggle('on', z === b)); draw(); });
  root.querySelectorAll('[data-w]').forEach(b => b.onclick = () => { st.wrapped = b.dataset.w === '1'; root.querySelectorAll('[data-w]').forEach(z => z.classList.toggle('on', z === b)); draw(); });
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => resetAvg();

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function draw() {
    if (!W) return;
    const mU = msd(sumU), mW = msd(sumW), D = fitD(mU);
    cx.clearRect(0, 0, W, H);
    // tracks of three tagged atoms (x-y projection of unwrapped positions)
    const topH = H * 0.42, tcol = [COL.cyan, COL.amber, COL.coral];
    let span = 1.5; tracks.forEach(tr => { if (!tr.length) return; tr.forEach(p => { span = Math.max(span, Math.abs(p[0] - tr[0][0]) * 1.1, Math.abs(p[1] - tr[0][1]) * 1.1); }); });
    const sc = Math.min(W / 3, topH) / 2 / span;
    tracks.forEach((tr, j) => {
      if (tr.length < 2) return; const ox = W * (j + 0.5) / 3, oy = topH / 2 + 8;
      cx.strokeStyle = tcol[j]; cx.lineWidth = 1.5; cx.beginPath(); tr.forEach((p, q) => { const X = ox + (p[0] - tr[0][0]) * sc, Y = oy - (p[1] - tr[0][1]) * sc; q ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }); cx.stroke();
      const e = tr[tr.length - 1]; cx.fillStyle = tcol[j]; cx.beginPath(); cx.arc(ox + (e[0] - tr[0][0]) * sc, oy - (e[1] - tr[0][1]) * sc, 4, 0, 7); cx.fill();
      cx.fillStyle = COL.ink; cx.beginPath(); cx.arc(ox, oy, 2.5, 0, 7); cx.fill();
    });
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`three atoms, last ${Math.min(30, nfr * EVERY * DT).toFixed(0)} τ (unwrapped, x–y)`, 10, 16);
    cx.strokeStyle = COL.muted; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(W - 10 - sc, topH - 4); cx.lineTo(W - 10, topH - 4); cx.stroke(); cx.textAlign = 'right'; cx.fillText('1 σ', W - 10, topH - 10);
    // MSD plot
    const x0 = 46, x1 = W - 14, yt = topH + 34, yb = H - 26, tmax = K * EVERY * DT;
    const shown = [mU, ...(st.wrapped ? [mW] : [])].flat().filter(q => q != null && q > 0), mmax = Math.max(0.05, ...shown, D ? 6 * D * tmax : 0) * 1.08;
    const lo = Math.log10(0.05), hi = Math.log10(tmax), ylo = Math.log10(Math.max(1e-4, Math.min(...(shown.length ? shown : [1e-3])) * 0.8)), yhi = Math.log10(mmax);
    const X = t => st.log ? x0 + (Math.log10(t) - lo) / (hi - lo) * (x1 - x0) : x0 + t / tmax * (x1 - x0);
    const Y = q => st.log ? yb - (Math.log10(q) - ylo) / (yhi - ylo) * (yb - yt) : yb - q / mmax * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center';
    (st.log ? [0.1, 1] : [1, 2, 3, 4]).forEach(t => cx.fillText(t, X(t), yb + 14)); cx.textAlign = 'right'; cx.fillText('t (τ)', x1, yb + 14);
    if (st.log) { for (let e = Math.ceil(ylo); e <= Math.floor(yhi); e++) cx.fillText('10' + ['⁻⁴', '⁻³', '⁻²', '⁻¹', '⁰', '¹', '²'][e + 4], x0 - 5, Y(10 ** e) + 4); }
    else { const raw = mmax / 4, p10 = 10 ** Math.floor(Math.log10(raw)), stp = [1, 2, 5, 10].map(m => m * p10).find(q => q >= raw); for (let q = stp; q < mmax; q += stp) cx.fillText(+q.toPrecision(3), x0 - 5, Y(q) + 4); }
    cx.textAlign = 'left'; cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillText('MSD (σ²)', x0 + 6, yt - 10);
    const curve = (m, c, w) => { cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); let on = false; for (let k = 1; k <= K; k++) { if (m[k] == null || m[k] <= 0) continue; const px = X(k * EVERY * DT), py = Y(m[k]); on ? cx.lineTo(px, py) : cx.moveTo(px, py); on = true; } cx.stroke(); };
    if (nfr > 1) {
      if (D != null && D > 0) { cx.setLineDash([6, 4]); cx.strokeStyle = COL.amber; cx.lineWidth = 1.5; cx.beginPath(); let on = false; for (let k = 2; k <= K; k++) { const t = k * EVERY * DT, q = 6 * D * t; if (q <= 0) continue; on ? cx.lineTo(X(t), Y(q)) : cx.moveTo(X(t), Y(q)); on = true; } cx.stroke(); cx.setLineDash([]); }
      if (st.log) { cx.setLineDash([3, 4]); cx.strokeStyle = 'rgba(255,107,107,.6)'; cx.lineWidth = 1.2; cx.beginPath(); let on = false; for (let k = 1; k <= 40; k++) { const t = k * EVERY * DT, q = 3 * st.T * t * t; if (q > mmax) break; on ? cx.lineTo(X(t), Y(q)) : cx.moveTo(X(t), Y(q)); on = true; } cx.stroke(); cx.setLineDash([]); }
      if (st.wrapped) curve(mW, COL.coral, 2);
      curve(mU, COL.cyan, 2.5);
      cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; let ly = yt + 4;
      [[COL.cyan, 'MSD (unwrapped)'], [COL.amber, 'dashed: 6D t (fit, 1–5 τ)'], ...(st.log ? [[COL.coral, 'dotted: ballistic 3T t²']] : []), ...(st.wrapped ? [[COL.coral, 'wrapped coordinates (wrong)']] : [])].forEach(([c, t]) => { cx.fillStyle = c; cx.fillText(t, x1, ly); ly += 15; });
    } else { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '13px "Space Grotesk", sans-serif'; cx.fillText('equilibrating…', (x0 + x1) / 2, (yt + yb) / 2); }
    // readouts
    o.T0.textContent = st.T.toFixed(2); o.T.textContent = Tnow.toFixed(2);
    const k2 = cnt[K] ? K : 0, sl = k2 && mU[K] > 0 && mU[K / 2] > 0 ? Math.log(mU[K] / mU[K / 2]) / Math.log(2) : null;
    o.D.textContent = D == null ? 'collecting…' : D < 2e-4 ? '≈ 0 (atoms only vibrate)' : sl != null && sl > 1.25 ? `${D.toFixed(3)}? not diffusive yet (slope > 1)` : D.toFixed(4);
    o.sl.textContent = sl == null ? '…' : `${sl.toFixed(2)} (1 = diffusive)`;
    o.nf.textContent = `${(nfr * EVERY * DT).toFixed(1)} τ / ${cnt[1] | 0}`;
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 6; k++) step(); draw(); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
