// Lab 3.8 widget: barostats on the 9-atom 2D LJ box (Langevin thermostat on the atoms, rc = 1.5, P0 = 3.0).
// Step shows one box update; Run collects the area distribution and compares it with the correct NPT distribution
// that the lab computes from the equation of state (step 2).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 9, RC = 1.5, T0 = 0.72, P0 = 3.0, DT = 0.005, BETA = 0.03, A_OK = 11.006, S_OK = 0.564;
const GRID0 = 9.2, DG = 0.2;   // correct p(A) on A = 9.2, 9.4, ..., 13.0 (from the lab script)
const PEX = [0.0000, 0.0001, 0.0010, 0.0060, 0.0221, 0.0552, 0.0986, 0.1339, 0.1484, 0.1427, 0.1200, 0.0910, 0.0656, 0.0448, 0.0289, 0.0179, 0.0109, 0.0066, 0.0039, 0.0023];
const NAMES = { ber: 'Berendsen', cres: 'C-rescale', piston: 'Piston' };

export function mountLab3b(root) {
  const st = { b: 'ber', tau: 0.5, W: 0.1, A0: 11, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view bp-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg">${Object.entries(NAMES).map(([k, n]) => `<button data-b="${k}" class="${k === 'ber' ? 'on' : ''}">${n}</button>`).join('')}</div>
    <label data-tau><span>Coupling time τ_P</span><output data-o="tau"></output><input type="range" data-k="tau" min="0.1" max="5" step="0.1" value="0.5"></label>
    <label data-w hidden><span>Piston mass W</span><output data-o="w"></output><input type="range" data-k="w" min="-2" max="0.7" step="0.05" value="-1"></label>
    <div class="seg"><button data-a="11" class="on">Start A = 11</button><button data-a="14">Start A = 14</button></div>
    <div class="row"><button class="btn sm primary" data-k="step">Step</button><button class="btn sm ghost" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="tbl" style="font-size:12px"><table data-o="tab"></table></div>
    <div class="readouts">
      <div><span>Average area (correct 11.006)</span><b data-o="A"></b></div>
      <div><span>Area fluctuation (correct 0.564)</span><b data-o="s"></b></div>
      <div><span>Area now</span><b data-o="An"></b></div>
      <div><span>Time</span><b data-o="t"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.bp-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 7; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let x, v, f, vir, A, L, Pi, P, t, rows, At, stats, hist;
  function forces() {
    f = x.map(() => [0, 0]); vir = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      let dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r2 = dx * dx + dy * dy; if (r2 > RC * RC) continue;
      const ir6 = 1 / (r2 * r2 * r2), ff = 24 * (2 * ir6 * ir6 - ir6) / r2;
      f[i][0] -= ff * dx; f[i][1] -= ff * dy; f[j][0] += ff * dx; f[j][1] += ff * dy; vir += ff * r2;
    }
  }
  const pres = () => (v.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2, 0) + vir) / (2 * A);
  const scaleX = m => { x.forEach(p => { p[0] *= m; p[1] *= m; }); };
  function reset() {
    A = st.A0; L = Math.sqrt(A); const a = L / 3;
    x = [0, 1, 2].flatMap(j => [0, 1, 2].map(i => [(i + 0.5) * a, (j + 0.5) * a]));
    v = x.map(() => [gauss() * Math.sqrt(T0), gauss() * Math.sqrt(T0)]);
    forces(); P = pres(); Pi = 0; t = 0; rows = []; At = []; stats = { n: 0, s: 0, s2: 0 }; hist = new Float64Array(22);
  }
  function step(record) {
    const dt = DT, c = Math.exp(-dt), sn = Math.sqrt((1 - c * c) * T0), W = st.W, Ab = A, Pb = P; let info = '';
    if (st.b === 'piston') { Pi += 0.5 * dt * (P - P0); const k = Math.exp(-0.25 * dt * Pi / W / A); v.forEach(p => { p[0] *= k; p[1] *= k; }); }
    v.forEach((p, i) => { p[0] += 0.5 * dt * f[i][0]; p[1] += 0.5 * dt * f[i][1]; x[i][0] += dt * p[0]; x[i][1] += dt * p[1]; });
    if (st.b === 'piston') { const An = Math.max(A + dt * Pi / W, 8); scaleX(Math.sqrt(An / A)); A = An; L = Math.sqrt(A); }
    x.forEach(p => { p[0] -= L * Math.floor(p[0] / L); p[1] -= L * Math.floor(p[1] / L); });
    forces(); v.forEach((p, i) => { p[0] += 0.5 * dt * f[i][0]; p[1] += 0.5 * dt * f[i][1]; });
    if (st.b === 'piston') { const k = Math.exp(-0.25 * dt * Pi / W / A); v.forEach(p => { p[0] *= k; p[1] *= k; }); }
    v.forEach(p => { p[0] = c * p[0] + sn * gauss(); p[1] = c * p[1] + sn * gauss(); });
    P = pres();
    if (st.b === 'piston') { Pi += 0.5 * dt * (P - P0); info = `box speed ${(Pi / W).toFixed(3)}`; }
    if (st.b === 'ber') { const mu2 = 1 - BETA * dt / st.tau * (P0 - P); scaleX(Math.sqrt(mu2)); A *= mu2; info = `μ² = ${mu2.toFixed(6)}`; }
    if (st.b === 'cres') { const R = gauss(), de = -(BETA / st.tau) * (P0 - P - T0 / A) * dt + Math.sqrt(2 * T0 * BETA * dt / (A * st.tau)) * R;
      scaleX(Math.exp(de / 2)); v.forEach(p => { p[0] *= Math.exp(-de / 2); p[1] *= Math.exp(-de / 2); }); A *= Math.exp(de); info = `d(ln A) = ${de >= 0 ? '+' : ''}${de.toFixed(5)}`; }
    if (st.b !== 'piston') { L = Math.sqrt(A); forces(); P = pres(); }
    t += dt;
    if (record) { rows.push([rows.length ? rows[rows.length - 1][0] + 1 : 1, Pb, info, Ab, A]); if (rows.length > 5) rows.shift(); }
    if (Math.round(t / dt) % 5 === 0) { At.push(A); if (At.length > 800) At.shift(); }
    if (t > 10) { stats.n++; stats.s += A; stats.s2 += A * A; const b = Math.floor((A - 9.0) / 0.2); if (b >= 0 && b < 22) hist[b]++; }
  }
  root.querySelectorAll('[data-b]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-b]').forEach(z => z.classList.toggle('on', z === b)); st.b = b.dataset.b;
    $('[data-tau]').hidden = st.b === 'piston'; $('[data-w]').hidden = st.b !== 'piston'; reset(); draw(); });
  root.querySelectorAll('[data-a]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-a]').forEach(z => z.classList.toggle('on', z === b)); st.A0 = +b.dataset.a; reset(); draw(); });
  $('[data-k=tau]').oninput = e => { st.tau = +e.target.value; reset(); draw(); };
  $('[data-k=w]').oninput = e => { st.W = 10 ** +e.target.value; reset(); draw(); };
  $('[data-k=step]').onclick = () => { step(true); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !x) return;
    cx.clearRect(0, 0, W, H);
    const half = H * 0.42, bw = Math.min(W * 0.42, half - 20), S = bw / 3.9, ox = 10 + (bw - L * S) / 2, oy = 10 + (bw - L * S) / 2;
    cx.strokeStyle = 'rgba(143,163,199,.35)'; cx.setLineDash([3, 4]); const Lok = Math.sqrt(A_OK); cx.strokeRect(10 + (bw - Lok * S) / 2, 10 + (bw - Lok * S) / 2, Lok * S, Lok * S); cx.setLineDash([]);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.strokeRect(ox, oy, L * S, L * S);
    x.forEach(p => { cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(ox + p[0] * S, oy + (L - p[1]) * S, 0.28 * S, 0, 7); cx.fill(); });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`P now ${P.toFixed(2)}`, 12, bw + 22);
    // A(t)
    const tx0 = bw + 50, tx1 = W - 12, ty0 = 20, ty1 = half - 10, lo = 9, hi = 14.5, TX = k => tx0 + k / 799 * (tx1 - tx0), TY = a => ty1 - (Math.min(Math.max(a, lo), hi) - lo) / (hi - lo) * (ty1 - ty0);
    cx.fillStyle = 'rgba(255,194,75,.12)'; cx.fillRect(tx0, TY(A_OK + S_OK), tx1 - tx0, TY(A_OK - S_OK) - TY(A_OK + S_OK));
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(tx0, ty0); cx.lineTo(tx0, ty1); cx.lineTo(tx1, ty1); cx.stroke();
    cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(tx0, TY(A_OK)); cx.lineTo(tx1, TY(A_OK)); cx.stroke(); cx.setLineDash([]);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 1.5; cx.beginPath(); At.forEach((a, k) => k ? cx.lineTo(TX(k), TY(a)) : cx.moveTo(TX(k), TY(a))); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.fillText('box area A', tx0 + 6, ty0 + 2);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [10, 12, 14].forEach(a => cx.fillText(a, tx0 - 4, TY(a) + 4));
    // histogram
    const hx0 = 46, hx1 = W - 12, hy0 = half + 40, hy1 = H - 26, tot = hist.reduce((a, b) => a + b, 0), HX = a => hx0 + (a - 9) / 4.4 * (hx1 - hx0);
    let pmax = Math.max(...PEX) / DG; if (tot) for (let k = 0; k < 22; k++) pmax = Math.max(pmax, hist[k] / tot / 0.2);
    const HY = p => hy1 - p / (pmax * 1.1) * (hy1 - hy0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(hx0, hy0); cx.lineTo(hx0, hy1); cx.lineTo(hx1, hy1); cx.stroke();
    cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.fillText(W < 500 ? 'A: bars vs correct (line)' : 'distribution of the box area (bars) vs correct NPT (line)', hx0 + 6, hy0 - 10);
    if (tot) { cx.fillStyle = 'rgba(56,225,255,.45)'; for (let k = 0; k < 22; k++) { const p = hist[k] / tot / 0.2; cx.fillRect(HX(9 + 0.2 * k) + 1, HY(p), (hx1 - hx0) / 22 - 2, hy1 - HY(p)); } }
    cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); PEX.forEach((p, k) => { const a = GRID0 + DG * k; k ? cx.lineTo(HX(a), HY(p / DG)) : cx.moveTo(HX(a), HY(p / DG)); }); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [10, 11, 12, 13].forEach(a => cx.fillText(a, HX(a), hy1 + 14));
    // table and readouts
    o.tab.innerHTML = `<tr><th>step</th><th>P</th><th>barostat</th><th>A before → after</th></tr>` + rows.map(r => `<tr><td>${r[0]}</td><td>${r[1].toFixed(3)}</td><td>${r[2]}</td><td>${r[3].toFixed(4)} → ${r[4].toFixed(4)}</td></tr>`).join('');
    o.tau.textContent = st.tau.toFixed(1); o.w.textContent = st.W < 0.1 ? st.W.toFixed(3) : st.W.toFixed(2);
    const n = stats.n, m = n ? stats.s / n : 0;
    o.A.textContent = n ? m.toFixed(3) : 'after t = 10'; o.s.textContent = n ? Math.sqrt(Math.max(0, stats.s2 / n - m * m)).toFixed(3) : 'after t = 10';
    o.An.textContent = A.toFixed(3); o.t.textContent = t.toFixed(1);
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let q = 0; q < 80; q++) step(false); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
