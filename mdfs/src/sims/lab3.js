// Lab 3 widget: thermostats on the 9-atom 2D LJ system of Labs 2 and 3, or on one harmonic oscillator.
// Step shows the thermostat quantities of each step; Run collects the statistics that decide whether a thermostat
// produces the correct (canonical) ensemble: average T, kinetic-energy fluctuation, energy in collective sliding modes.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const L = 3.6, N = 9, RC = 1.8, T0 = 0.72, NF = 16;
const LAT = [0, 1, 2].flatMap(j => [0, 1, 2].map(i => [(i + 0.5) * 1.2, (j + 0.5) * 1.2]));
const GAU = [[-0.33, -0.88], [-0.66, -0.67], [0.38, -0.11], [1.48, -1.83], [-0.00, -0.89], [0.78, -2.12], [-0.34, 0.21], [-1.48, 0.99], [0.18, 1.01]];
const NAMES = { nve: 'NVE', rescale: 'Rescale', ber: 'Berendsen', nh: 'Nosé–Hoover', nhc: 'NH chain', bussi: 'Bussi', lang: 'Langevin' };

export function mountLab3(root) {
  const st = { sys: 'atoms', th: 'ber', tau: 0.1, gamma: 1, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="atoms" class="on">9 atoms</button><button data-s="osc">1 oscillator</button></div>
    <div class="seg" style="flex-wrap:wrap">${Object.entries(NAMES).map(([k, n]) => `<button data-t="${k}" class="${k === 'ber' ? 'on' : ''}">${n}</button>`).join('')}</div>
    <label data-tau><span>Coupling time τ</span><output data-o="tau"></output><input type="range" data-k="tau" min="0.02" max="2" step="0.01" value="0.1"></label>
    <label data-gam hidden><span>Friction γ</span><output data-o="gam"></output><input type="range" data-k="gam" min="0.1" max="20" step="0.1" value="1"></label>
    <div class="row"><button class="btn sm primary" data-k="step">Step</button><button class="btn sm ghost" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="tbl" style="font-size:12px"><table data-o="tab"></table></div>
    <div class="readouts">
      <div><span>Average <i class="nc">T</i> (target)</span><b data-o="T"></b></div>
      <div><span>KE fluctuation (correct)</span><b data-o="fl"></b></div>
      <div><span data-o="l3">Sliding-mode KE (correct 0.25)</span><b data-o="sl"></b></div>
      <div><span>Time sampled</span><b data-o="ts"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 11; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let x, v, f, xi, xi2, t, rows, Tt, stats, hist, phase;
  const isOsc = () => st.sys === 'osc', DT = () => isOsc() ? 0.05 : 0.005, Tg = () => isOsc() ? 1 : T0, nf = () => isOsc() ? 1 : (st.th === 'lang' ? 2 * N : NF);
  function forces() {
    if (isOsc()) { f = [[-x[0][0], 0]]; return; }
    f = x.map(() => [0, 0]);
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { let dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r2 = dx * dx + dy * dy; if (r2 > RC * RC) continue; const ir6 = 1 / (r2 * r2 * r2), ff = 24 * (2 * ir6 * ir6 - ir6) / r2; f[i][0] -= ff * dx; f[i][1] -= ff * dy; f[j][0] += ff * dx; f[j][1] += ff * dy; }
  }
  const sumv2 = () => isOsc() ? v[0][0] ** 2 : v.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2, 0);
  const scaleV = l => v.forEach(p => { p[0] *= l; p[1] *= l; });
  function reset() {
    if (isOsc()) { x = [[1, 0]]; v = [[0, 0]]; }
    else { x = LAT.map(p => [...p]); const m = [0, 1].map(k => GAU.reduce((s, p) => s + p[k], 0) / N); v = GAU.map(p => [p[0] - m[0], p[1] - m[1]]); const s2 = v.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2, 0); scaleV(Math.sqrt(2.0 * NF / s2)); }
    xi = xi2 = 0; t = 0; forces(); rows = []; Tt = []; stats = { n: 0, K: 0, K2: 0, sl: 0, fast: 0 }; hist = new Float64Array(40); phase = [];
  }
  function step(record) {
    const dt = DT(), Tgt = Tg(), Nf = isOsc() ? 1 : NF, Q = isOsc() ? 1 : Nf * Tgt * st.tau ** 2, Tb = sumv2() / nf(); let info = '';
    const nhHalf = () => { if (st.th === 'nhc') { xi2 += 0.25 * dt * (Q * xi * xi - Tgt) / Q; xi *= Math.exp(-0.25 * dt * xi2); } xi += 0.5 * dt * (sumv2() - Nf * Tgt) / Q; scaleV(Math.exp(-0.5 * dt * xi)); };
    const nhHalf2 = () => { scaleV(Math.exp(-0.5 * dt * xi)); xi += 0.5 * dt * (sumv2() - Nf * Tgt) / Q; if (st.th === 'nhc') { xi *= Math.exp(-0.25 * dt * xi2); xi2 += 0.25 * dt * (Q * xi * xi - Tgt) / Q; } };
    if (st.th === 'nh' || st.th === 'nhc') nhHalf();
    v.forEach((p, i) => { p[0] += 0.5 * dt * f[i][0]; p[1] += 0.5 * dt * f[i][1]; x[i][0] += dt * p[0]; x[i][1] += dt * p[1]; if (!isOsc()) { x[i][0] -= L * Math.floor(x[i][0] / L); x[i][1] -= L * Math.floor(x[i][1] / L); } });
    forces(); v.forEach((p, i) => { p[0] += 0.5 * dt * f[i][0]; p[1] += 0.5 * dt * f[i][1]; });
    if (st.th === 'nh' || st.th === 'nhc') { nhHalf2(); info = `ξ = ${xi.toFixed(4)}`; }
    const Tm = sumv2() / nf();
    if (st.th === 'rescale') { const l = Math.sqrt(Tgt / Tm); scaleV(l); info = `λ = ${l.toFixed(5)}`; }
    if (st.th === 'ber') { const l = Math.sqrt(1 + dt / st.tau * (Tgt / Tm - 1)); scaleV(l); info = `λ = ${l.toFixed(5)}`; }
    if (st.th === 'lang') { const c = Math.exp(-st.gamma * dt), s = Math.sqrt((1 - c * c) * Tgt); v.forEach(p => { p[0] = c * p[0] + s * gauss(); if (!isOsc()) p[1] = c * p[1] + s * gauss(); }); info = `c = ${c.toFixed(5)}`; }
    if (st.th === 'bussi') { const c = Math.exp(-dt / st.tau), K = 0.5 * sumv2(), Kt = 0.5 * Nf * Tgt, R1 = gauss(); let S = 0; for (let q = 0; q < Nf - 1; q++) S += gauss() ** 2;
      const a2 = c + (1 - c) * Kt / (Nf * K) * (R1 * R1 + S) + 2 * R1 * Math.sqrt(c * (1 - c) * Kt / (Nf * K)), sg = Math.sign(R1 + Math.sqrt(c * Nf * K / ((1 - c) * Kt))) || 1; scaleV(Math.sqrt(Math.max(a2, 0)) * sg); info = `α = ${Math.sqrt(Math.max(a2, 0)).toFixed(5)}`; }
    t += dt; const Ta = sumv2() / nf();
    if (record) { rows.push([rows.length ? rows[rows.length - 1][0] + 1 : 1, Tb, info, Ta]); if (rows.length > 5) rows.shift(); }
    Tt.push(Ta); if (Tt.length > 600) Tt.shift();
    if (t > (isOsc() ? 50 : 5)) {   // statistics after a short equilibration, as running averages over ~100 time units (atoms) or ~5000 (oscillator)
      const K = 0.5 * sumv2(); stats.n++; const w = Math.max(1 / stats.n, dt / (isOsc() ? 5000 : 100)), ew = (k, val) => { stats[k] += w * (val - stats[k]); };
      ew('K', K); ew('K2', K * K); for (let k = 0; k < 40; k++) hist[k] *= 1 - w;
      if (isOsc()) { ew('fast', Math.abs(v[0][0]) > 1.5 ? 1 : 0); const b = Math.floor((v[0][0] + 4) / 8 * 40); if (b >= 0 && b < 40) hist[b]++; if (stats.n % 5 === 0) { phase.push([x[0][0], v[0][0]]); if (phase.length > 3000) phase.shift(); } }
      else { let s = 0; for (let k = 0; k < 3; k++) { let rx = 0, cy = 0; for (let q = 0; q < 3; q++) { rx += v[3 * k + q][0]; cy += v[k + 3 * q][1]; } s += (rx / 3) ** 2 + (cy / 3) ** 2; } ew('sl', 1.5 * s);
        const b = Math.floor(Ta / (3 * T0) * 40); if (b >= 0 && b < 40) hist[b]++; }
    }
  }
  const seg = (a, fn) => root.querySelectorAll(`[data-${a}]`).forEach(b => b.onclick = () => { root.querySelectorAll(`[data-${a}]`).forEach(z => z.classList.toggle('on', z === b)); fn(b.dataset[a]); reset(); draw(); });
  seg('s', val => st.sys = val); seg('t', val => { st.th = val; $('[data-tau]').hidden = val === 'lang' || val === 'nve'; $('[data-gam]').hidden = val !== 'lang'; });
  $('[data-k=tau]').oninput = e => { st.tau = +e.target.value; reset(); draw(); };
  $('[data-k=gam]').oninput = e => { st.gamma = +e.target.value; reset(); draw(); };
  $('[data-k=step]').onclick = () => { step(true); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function gammaPdf(T) { const k = nf() / 2, th = T0 / k; let lg = 0; for (let q = 1; q < k; q++) lg += Math.log(q); return Math.exp((k - 1) * Math.log(T) - T / th - lg - k * Math.log(th)); }   // distribution of T for 16 dof
  function draw() {
    if (!W || !x) return;
    cx.clearRect(0, 0, W, H);
    const Tgt = Tg(), half = H * 0.42;
    // top-left: system (atoms or phase space); top-right: T(t)
    const bw = Math.min(W * 0.42, half - 20);
    if (isOsc()) {
      const S = bw / 7, cxm = 10 + bw / 2, cym = 10 + bw / 2; cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.strokeRect(10, 10, bw, bw);
      cx.fillStyle = 'rgba(56,225,255,.35)'; phase.forEach(p => cx.fillRect(cxm + p[0] * S, cym - p[1] * S, 1.5, 1.5));
      cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(cxm + x[0][0] * S, cym - v[0][0] * S, 4, 0, 7); cx.fill();
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('phase space (x, v)', 14, 24);
    } else {
      const S = bw / L; cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.strokeRect(10, 10, bw, bw);
      x.forEach((p, i) => { const X = 10 + p[0] * S, Y = 10 + (L - p[1]) * S; cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(X, Y, 0.3 * S, 0, 7); cx.fill(); cx.strokeStyle = COL.coral; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(X, Y); cx.lineTo(X + v[i][0] * 0.2 * S, Y - v[i][1] * 0.2 * S); cx.stroke(); });
    }
    const tx0 = bw + 50, tx1 = W - 12, ty0 = 20, ty1 = half - 10, Tmax = Math.max(2.2 * Tgt, ...Tt) * 1.05, TX = k => tx0 + k / 599 * (tx1 - tx0), TY = T => ty1 - T / Tmax * (ty1 - ty0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(tx0, ty0); cx.lineTo(tx0, ty1); cx.lineTo(tx1, ty1); cx.stroke();
    cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(tx0, TY(Tgt)); cx.lineTo(tx1, TY(Tgt)); cx.stroke(); cx.setLineDash([]);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 1.5; cx.beginPath(); Tt.forEach((T, k) => k ? cx.lineTo(TX(k), TY(T)) : cx.moveTo(TX(k), TY(T))); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(isOsc() ? 'T = v² (instantaneous)' : 'instantaneous T', tx0 + 6, ty0 + 2);
    // bottom: histogram vs correct distribution
    const hx0 = 46, hx1 = W - 12, hy0 = half + 40, hy1 = H - 26, tot = hist.reduce((a, b) => a + b, 0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(hx0, hy0); cx.lineTo(hx0, hy1); cx.lineTo(hx1, hy1); cx.stroke();
    cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.fillText(W < 500 ? (isOsc() ? 'v: bars vs correct (line)' : 'T: bars vs canonical (line)') : (isOsc() ? 'distribution of v (bars) vs correct Gaussian (line)' : 'distribution of instantaneous T (bars) vs canonical (line)'), hx0 + 6, hy0 - 10);
    const lo = isOsc() ? -4 : 0, hi = isOsc() ? 4 : 3 * T0, bwid = (hi - lo) / 40, pdf = z => isOsc() ? Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI) : gammaPdf(z);
    let pmax = 0; for (let k = 0; k < 40; k++) pmax = Math.max(pmax, pdf(lo + (k + 0.5) * bwid), tot ? hist[k] / tot / bwid : 0);
    const HX = z => hx0 + (z - lo) / (hi - lo) * (hx1 - hx0), HY = p => hy1 - p / (pmax * 1.1) * (hy1 - hy0);
    if (tot) { cx.fillStyle = 'rgba(56,225,255,.45)'; for (let k = 0; k < 40; k++) { const p = hist[k] / tot / bwid; cx.fillRect(HX(lo + k * bwid) + 1, HY(p), (hx1 - hx0) / 40 - 2, hy1 - HY(p)); } }
    cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); for (let k = 0; k <= 200; k++) { const z = lo + (hi - lo) * k / 200 + 1e-6; k ? cx.lineTo(HX(z), HY(pdf(z))) : cx.moveTo(HX(z), HY(pdf(z))); } cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; (isOsc() ? [-2, 0, 2] : [0.72, 1.44]).forEach(z => cx.fillText(z, HX(z), hy1 + 14));
    // table and readouts
    o.tab.innerHTML = `<tr><th>step</th><th>T before</th><th>thermostat</th><th>T after</th></tr>` + rows.map(r => `<tr><td>${r[0]}</td><td>${r[1].toFixed(4)}</td><td>${r[2] || '—'}</td><td>${r[3].toFixed(4)}</td></tr>`).join('');
    o.tau.textContent = st.tau.toFixed(2); o.gam.textContent = st.gamma.toFixed(1);
    const n = stats.n, mK = stats.K, fl = n ? Math.sqrt(Math.max(0, stats.K2 - mK * mK)) / mK : 0, corr = Math.sqrt(2 / (isOsc() ? 1 : nf()));
    o.T.textContent = n ? `${(2 * mK / nf()).toFixed(3)} (${Tgt})` : '…';
    o.fl.textContent = n ? `${fl.toFixed(3)} (${corr.toFixed(3)})` : '…';
    o.l3.textContent = isOsc() ? 'Time with |v| > 1.5 (correct 0.134)' : `Sliding-mode KE (correct ${st.th === 'lang' ? '0.33' : '0.25'})`;
    o.sl.textContent = n ? (isOsc() ? stats.fast.toFixed(3) : (stats.sl / mK).toFixed(3)) : '…';
    o.ts.textContent = t.toFixed(1);
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let q = 0; q < (isOsc() ? 400 : 120); q++) step(false); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
