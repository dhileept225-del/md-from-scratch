// RDF lab: 256 Lennard-Jones atoms in 3D (velocity Verlet + Berendsen thermostat), with g(r) built up live.
// Top: a thin slice around one tagged atom (its first shell highlighted). Bottom: the running g(r),
// the first-shell area (coordination number) and, for the liquid, the reference curve from the code tab.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 256, DT = 0.005, RC2 = 6.25, DR = 0.025;
const STATES = { gas: { rho: 0.05, T: 2.0 }, liquid: { rho: 0.8442, T: 0.722 }, solid: { rho: 1.1, T: 0.3 } };
const REF = { r0: 0.0112, dr: 0.0224, g: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.11, 0.44, 1.06, 1.81, 2.54, 2.92, 3.06, 2.91, 2.75, 2.46, 2.14, 1.86, 1.62, 1.43, 1.23, 1.09, 0.99, 0.9, 0.81, 0.75, 0.68, 0.65, 0.62, 0.6, 0.59, 0.58, 0.57, 0.56, 0.58, 0.58, 0.6, 0.62, 0.64, 0.67, 0.72, 0.76, 0.81, 0.86, 0.92, 0.96, 1.01, 1.07, 1.11, 1.14, 1.17, 1.2, 1.25, 1.26, 1.27, 1.29, 1.29, 1.3, 1.29, 1.26, 1.24, 1.2, 1.17, 1.13, 1.1, 1.06, 1.0, 0.98, 0.93, 0.91, 0.89, 0.85, 0.84, 0.84, 0.81, 0.81, 0.81, 0.82, 0.82, 0.83, 0.86, 0.87, 0.89, 0.91, 0.93, 0.95, 0.97, 1.0, 1.02, 1.04, 1.06, 1.08, 1.09, 1.11, 1.11, 1.11, 1.11, 1.12, 1.11, 1.1, 1.1, 1.09, 1.07, 1.07, 1.06, 1.04, 1.03, 1.02, 1.0, 0.98, 0.98, 0.96, 0.95, 0.94] };   // python reference: rho 0.8442, T 0.722

export function mountRdf(root) {
  const st = { state: 'liquid', T: 0.722, running: true };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="gas">Gas</button><button data-s="liquid" class="on">Liquid</button><button data-s="solid">Solid</button></div>
    <label><span>Temperature <i class="nc">T</i>*</span><output data-o="T0"></output><input type="range" data-k="T" min="0.2" max="2.5" step="0.05" value="0.722"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset average</button></div>
    <div class="readouts">
      <div><span>Density <i class="nc">ρ</i>* / measured <i class="nc">T</i>*</span><b data-o="rt"></b></div>
      <div><span>Frames averaged</span><b data-o="nf"></b></div>
      <div><span>First peak (r, g)</span><b data-o="pk"></b></div>
      <div><span>Neighbours to 1st minimum</span><b data-o="n1"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let L = 1, x = new Float64Array(3 * N), v = new Float64Array(3 * N), f = new Float64Array(3 * N), hist = [], nfr = 0, rmax = 3, nb = 1, Tnow = 0, steps = 0;
  let rng = 12345; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());

  function build() {   // fcc lattice, 4x4x4 cells, at the chosen density
    const s = STATES[st.state], a = Math.cbrt(4 / s.rho), B = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]; L = 4 * a; let k = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let l = 0; l < 4; l++) for (const b of B) { x[3 * k] = (i + b[0]) * a; x[3 * k + 1] = (j + b[1]) * a; x[3 * k + 2] = (l + b[2]) * a; k++; }
    const m = [0, 0, 0]; for (let i = 0; i < 3 * N; i++) { v[i] = gauss() * Math.sqrt(st.T); m[i % 3] += v[i] / N; } for (let i = 0; i < 3 * N; i++) v[i] -= m[i % 3];
    rmax = Math.min(L / 2, 3.5); nb = Math.floor(rmax / DR); forces(); resetAvg(); steps = 0;
  }
  function resetAvg() { hist = new Float64Array(nb); nfr = 0; }
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
    for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
    forces(); let ke = 0;
    for (let i = 0; i < 3 * N; i++) { v[i] += 0.5 * DT * f[i]; ke += v[i] * v[i]; }
    Tnow = ke / (3 * N - 3); const lam = Math.sqrt(1 + DT / 0.1 * (st.T / Tnow - 1)); for (let i = 0; i < 3 * N; i++) v[i] *= lam;   // Berendsen, tau = 0.1
    if (++steps > 300 && steps % 5 === 0) sample();                                    // skip the first 300 steps (melting the lattice)
  }
  function sample() {
    for (let i = 0; i < N - 1; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz); if (r < rmax) hist[Math.floor(r / DR)] += 2;
    }
    nfr++;
  }
  function gofr() {   // histogram / (frames * N * rho * shell volume)
    const rho = N / L ** 3, g = new Float64Array(nb);
    for (let k = 0; k < nb; k++) { const sh = 4 / 3 * Math.PI * (((k + 1) * DR) ** 3 - (k * DR) ** 3); g[k] = nfr ? hist[k] / (nfr * N * rho * sh) : 0; }
    return g;
  }
  function shell(g) {
    let i1 = 0; for (let k = 0; k < nb && (k + 0.5) * DR < 1.6; k++) if (g[k] > g[i1]) i1 = k;
    let j = i1; for (let k = i1; k < Math.min(nb, i1 + Math.round(0.7 / DR)); k++) if (g[k] < g[j]) j = k;
    const rho = N / L ** 3; let n1 = 0; for (let k = 0; k <= j; k++) n1 += 4 * Math.PI * rho * ((k + 0.5) * DR) ** 2 * g[k] * DR;
    return { i1, j, n1 };
  }

  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { st.state = b.dataset.s; st.T = STATES[st.state].T; $('[data-k=T]').value = st.T; root.querySelectorAll('[data-s]').forEach(z => z.classList.toggle('on', z === b)); build(); });
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; resetAvg(); steps = Math.min(steps, 100); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => resetAvg();

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function draw() {
    if (!W) return;
    const g = gofr(), sh = nfr > 3 ? shell(g) : null, rMin = sh ? (sh.j + 0.5) * DR : 1.5;
    cx.clearRect(0, 0, W, H);
    // slice around tagged atom 0
    const topH = H * 0.5, S = Math.min(W, topH) / 8.6, cxv = W / 2, cyv = topH / 2 + 4; let nNow = 0;
    const pts = [];
    for (let i = 0; i < N; i++) {
      let dx = x[3 * i] - x[0], dy = x[3 * i + 1] - x[1], dz = x[3 * i + 2] - x[2];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
      const r = Math.hypot(dx, dy, dz); if (i && r < rMin) nNow++;
      if (Math.abs(dz) < 0.6 && Math.abs(dx) < 3.6 && Math.abs(dy) < 3.6) pts.push([dx, dy, dz, r, i]);
    }
    cx.save(); cx.beginPath(); cx.rect(0, 0, W, topH); cx.clip();
    cx.strokeStyle = 'rgba(56,225,255,.35)'; cx.setLineDash([5, 5]); cx.lineWidth = 1.5; cx.beginPath(); cx.arc(cxv, cyv, rMin * S, 0, 7); cx.stroke(); cx.setLineDash([]);
    pts.sort((a, b) => a[2] - b[2]).forEach(([dx, dy, dz, r, i]) => {
      cx.globalAlpha = 1 - Math.abs(dz) * 0.9; cx.fillStyle = i === 0 ? COL.amber : r < rMin ? COL.cyan : 'rgba(143,163,199,.75)';
      cx.beginPath(); cx.arc(cxv + dx * S, cyv - dy * S, 0.45 * S, 0, 7); cx.fill();
    });
    cx.restore(); cx.globalAlpha = 1;
    cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.fillStyle = COL.amber; cx.fillText('tagged atom', 10, 18); cx.fillStyle = COL.cyan; cx.fillText(`first shell (3D): ${nNow} now`, 10, 34); cx.fillStyle = COL.muted; cx.fillText('slice ±0.6 σ thick', 10, 50);
    // g(r) plot
    const x0 = 40, x1 = W - 14, yt = topH + 30, yb = H - 26, rM = 3.5, gmax = st.state === 'solid' ? 8 : st.state === 'gas' ? 2.5 : 3.6;
    const X = r => x0 + r / rM * (x1 - x0), Y = q => yb - Math.min(q, gmax * 1.05) / gmax * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(x0, Y(1)); cx.lineTo(x1, Y(1)); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 2, 3].forEach(r => cx.fillText(r, X(r), yb + 14)); cx.fillText('r (σ)', X(3.35), yb + 14);
    cx.textAlign = 'right'; [1, Math.round(gmax / 2), Math.floor(gmax)].forEach(q => cx.fillText(q, x0 - 5, Y(q) + 4));
    cx.textAlign = 'left'; cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillText('g(r)', x0 + 6, yt - 10);
    if (st.state === 'liquid' && Math.abs(st.T - 0.722) < 0.01) { cx.strokeStyle = 'rgba(242,245,250,.55)'; cx.setLineDash([6, 4]); cx.lineWidth = 1.5; cx.beginPath(); REF.g.forEach((q, k) => { const r = REF.r0 + k * REF.dr; k ? cx.lineTo(X(r), Y(q)) : cx.moveTo(X(r), Y(q)); }); cx.stroke(); cx.setLineDash([]); cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText('dashed: code tab (3000 steps)', x1, yt - 10); }
    if (nfr) {
      if (sh && st.state !== 'gas') { cx.fillStyle = 'rgba(255,194,75,.28)'; cx.beginPath(); cx.moveTo(X(0), Y(0)); for (let k = 0; k <= sh.j; k++) cx.lineTo(X((k + 0.5) * DR), Y(g[k])); cx.lineTo(X((sh.j + 0.5) * DR), Y(0)); cx.closePath(); cx.fill(); }
      cx.strokeStyle = COL.cyan; cx.lineWidth = 2.5; cx.beginPath(); for (let k = 0; k < nb; k++) { const px = X((k + 0.5) * DR), py = Y(g[k]); k ? cx.lineTo(px, py) : cx.moveTo(px, py); } cx.stroke();
    } else { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '13px "Space Grotesk", sans-serif'; cx.fillText('equilibrating…', (x0 + x1) / 2, (yt + yb) / 2); }
    // readouts
    o.T0.textContent = st.T.toFixed(2);
    o.rt.textContent = `${(N / L ** 3).toFixed(3)} / ${Tnow.toFixed(2)}`;
    o.nf.textContent = nfr; o.pk.textContent = sh ? `${((sh.i1 + 0.5) * DR).toFixed(2)} σ, ${g[sh.i1].toFixed(2)}` : '…';
    o.n1.textContent = st.state === 'gas' ? 'no clear shell in a gas' : sh ? `${sh.n1.toFixed(1)} (to r = ${rMin.toFixed(2)} σ)` : '…';
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 5; k++) step(); draw(); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
