// Finite-size lab widget, two modes.
// g(r) tail: ideal gas (random positions) in a 2D periodic box with N particles; g(r) -> 1 - 1/N; running Kirkwood-Buff
//            integral with and without the Ganguly-van der Vegt correction.
// Diffusion: the four 3D LJ box sizes of the lab vs Yeh-Hummer, plus a calculator for real liquids.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const LJ = [[3.420, 0.0440, 0.0014], [5.130, 0.0565, 0.0013], [6.840, 0.0597, 0.0021], [8.550, 0.0626, 0]];
const XI = 2.837297, KB = 1.380649e-23;

export function mountLab10(root) {
  const st = { mode: 'g', N: 30, gv: false, eta: 50, L: 5, D: 1.0, T: 300, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view fs-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="g" class="on">g(r) tail</button><button data-m="d">Diffusion</button></div>
    <div data-p="g">
      <label><span>Number of particles N</span><output data-o="N"></output><input type="range" data-k="N" min="5" max="200" step="1" value="30"></label>
      <div class="seg"><button data-c="0" class="on">Raw g(r)</button><button data-c="1">Ganguly–van der Vegt</button></div>
      <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    </div>
    <div data-p="d" hidden>
      <label><span>Viscosity η (mPa s)</span><output data-o="eta"></output><input type="range" data-k="eta" min="-0.5" max="3" step="0.05" value="1.7"></label>
      <label><span>Box side L (nm)</span><output data-o="L"></output><input type="range" data-k="L" min="2" max="20" step="0.5" value="5"></label>
      <label><span>D in the box (10⁻¹¹ m²/s)</span><output data-o="D"></output><input type="range" data-k="D" min="-1" max="2.5" step="0.05" value="0"></label>
    </div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.fs-view'), $ = s => root.querySelector(s);
  let seed = 77; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  const NB = 30; let hist, frames;
  const Lbox = () => Math.sqrt(st.N), RMAX = () => Lbox() / 2;
  function reset() { hist = new Float64Array(NB); frames = 0; }
  function frame() {
    const N = st.N, L = Lbox(), x = Array.from({ length: N }, () => [rand() * L, rand() * L]), dr = RMAX() / NB;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { let dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
      const r = Math.hypot(dx, dy), b = Math.floor(r / dr); if (b < NB) hist[b]++; }
    frames++;
  }
  function gofr() {
    const N = st.N, L = Lbox(), dr = RMAX() / NB, rho = N / (L * L), g = [], G = [], Gc = []; let run = 0, runc = 0;
    for (let b = 0; b < NB; b++) { const r0 = b * dr, r1 = r0 + dr, ideal = rho * Math.PI * (r1 * r1 - r0 * r0) * N / 2; g.push(frames ? hist[b] / frames / ideal : 0); }
    const gc = g.map((v, b) => { const r = (b + 0.5) * dr, x = Math.PI * r * r / (L * L), dN = rho * (run + (v - 1) * 2 * Math.PI * r * dr * 0.5);
      const f = N * (1 - x) / (N * (1 - x) - dN - 1); run += (v - 1) * 2 * Math.PI * r * dr; G.push(run); return v * f; });
    gc.forEach((v, b) => { const r = (b + 0.5) * dr; runc += (v - 1) * 2 * Math.PI * r * dr; Gc.push(runc); });
    return { g, gc, G, Gc, dr };
  }
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); draw(); });
  root.querySelectorAll('[data-c]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-c]').forEach(z => z.classList.toggle('on', z === b)); st.gv = b.dataset.c === '1'; draw(); });
  $('[data-k=N]').oninput = e => { st.N = +e.target.value; reset(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };
  $('[data-k=eta]').oninput = e => { st.eta = 10 ** +e.target.value; draw(); };
  $('[data-k=L]').oninput = e => { st.L = +e.target.value; draw(); };
  $('[data-k=D]').oninput = e => { st.D = 10 ** +e.target.value; draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x, y); };
  const hline = (y, x0, x1, c, dash = true) => { if (dash) cx.setLineDash([5, 4]); cx.strokeStyle = c; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y); cx.lineTo(x1, y); cx.stroke(); cx.setLineDash([]); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const R = $('[data-ro]'), half = H * 0.48;
    if (st.mode === 'g') {
      const { g, gc, G, Gc, dr } = gofr(), show = st.gv ? gc : g, x0 = 46, x1 = W - 12, y0 = 22, y1 = half, X = r => x0 + r / RMAX() * (x1 - x0), lo = 0.8, hi = 1.1, Y = v => y1 - (Math.min(hi, Math.max(lo, v)) - lo) / (hi - lo) * (y1 - y0);
      axes(x0, y0, x1, y1); label(`g(r) of an ideal gas, ${st.N} particles (${frames} frames)`, x0 + 6, y0 - 6);
      hline(Y(1), x0, x1, 'rgba(242,245,250,.7)'); hline(Y(1 - 1 / st.N), x0, x1, COL.coral);
      if (frames) { cx.fillStyle = st.gv ? COL.amber : COL.cyan; show.forEach((v, b) => { cx.beginPath(); cx.arc(X((b + 0.5) * dr), Y(v), 3.5, 0, 7); cx.fill(); }); }
      label('white: 1   coral: 1 − 1/N', x1 - 4, y1 - 6, COL.muted, 'right', '11px');
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [0.9, 1].forEach(v => cx.fillText(v, x0 - 4, Y(v) + 4));
      const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 28, gl = Math.PI * RMAX() ** 2 / st.N * 1.2, GY = v => (by0 + by1) / 2 - Math.max(-gl, Math.min(gl, v)) / gl * (by1 - by0) / 2;
      axes(bx0, by0, bx1, by1); label('running Kirkwood–Buff integral G(R) (true value 0)', bx0 + 6, by0 - 10); hline(GY(0), bx0, bx1, 'rgba(242,245,250,.6)');
      if (frames) { [[G, COL.cyan], [Gc, COL.amber]].forEach(([arr, c]) => { cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath(); arr.forEach((v, b) => { const xx = bx0 + (b + 0.5) / NB * (bx1 - bx0); b ? cx.lineTo(xx, GY(v)) : cx.moveTo(xx, GY(v)); }); cx.stroke(); }); }
      label('cyan: raw   amber: corrected', bx1 - 4, by1 - 6, COL.muted, 'right', '11px');
      const tail = frames ? show.slice(NB / 2).reduce((a, b) => a + b, 0) / (NB / 2) : 0;
      R.innerHTML = `<div><span>Tail of g(r) (outer half)</span><b>${frames ? tail.toFixed(4) : '…'}</b></div><div><span>1 − 1/N</span><b>${(1 - 1 / st.N).toFixed(4)}</b></div>
        <div><span>G at R = L/2, raw</span><b>${frames ? G[NB - 1].toFixed(3) : '…'}</b></div><div><span>G at R = L/2, corrected</span><b>${frames ? Gc[NB - 1].toFixed(3) : '…'}</b></div>`;
      $('[data-o=N]').textContent = st.N;
    } else {
      // LJ data vs 1/L
      const x0 = 50, x1 = W - 12, y0 = 22, y1 = half, X = q => x0 + q / 0.32 * (x1 - x0), Y = d => y1 - (d - 0.03) / 0.05 * (y1 - y0);
      axes(x0, y0, x1, y1); label('3D LJ: D vs 1/L (dots), line D = 0.0755 − 0.1054/L', x0 + 6, y0 - 6, COL.ink, 'left', W < 500 ? '600 10px' : '600 12px');
      cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(X(0), Y(0.0755)); cx.lineTo(X(0.32), Y(0.0755 - 0.1054 * 0.32)); cx.stroke();
      LJ.forEach(([L, D, e]) => { cx.strokeStyle = COL.cyan; cx.beginPath(); cx.moveTo(X(1 / L), Y(D - e)); cx.lineTo(X(1 / L), Y(D + e)); cx.stroke(); cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(X(1 / L), Y(D), 5, 0, 7); cx.fill();
        label(`${Math.round(0.8 * L ** 3)}`, X(1 / L), Y(D) - 9, COL.muted, 'center', '11px'); });
      cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(X(0), Y(0.0755), 5, 0, 7); cx.fill(); label('infinite box', X(0) + 8, Y(0.0755) + 4, COL.amber, 'left', '11px');
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0, 0.1, 0.2, 0.3].forEach(q => cx.fillText(q, X(q), y1 + 14)); cx.textAlign = 'right'; [0.04, 0.06, 0.08].forEach(d => cx.fillText(d, x0 - 4, Y(d) + 4));
      // calculator: correction (%) vs L for the chosen liquid
      const corr = L => XI * KB * st.T / (6 * Math.PI * st.eta * 1e-3 * L * 1e-9), Dsi = st.D * 1e-11;
      const bx0 = 50, bx1 = W - 12, by0 = half + 46, by1 = H - 28, BX = L => bx0 + (L - 2) / 18 * (bx1 - bx0), pmax = Math.min(200, Math.max(10, 100 * corr(2) / Dsi)), BY = p => by1 - Math.min(p, pmax) / pmax * (by1 - by0);
      axes(bx0, by0, bx1, by1); label('Yeh–Hummer correction (% of D in the box) vs box side L', bx0 + 6, by0 - 10, COL.ink, 'left', W < 500 ? '600 10px' : '600 12px');
      cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); for (let L = 2; L <= 20; L += 0.25) L > 2 ? cx.lineTo(BX(L), BY(100 * corr(L) / Dsi)) : cx.moveTo(BX(L), BY(100 * corr(L) / Dsi)); cx.stroke();
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(BX(st.L), BY(100 * corr(st.L) / Dsi), 6, 0, 7); cx.fill();
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [2, 5, 10, 15, 20].forEach(L => cx.fillText(L + ' nm', BX(L), by1 + 14)); cx.textAlign = 'right'; cx.fillText(pmax.toFixed(0) + ' %', bx0 - 4, by0 + 4);
      const c = corr(st.L);
      R.innerHTML = `<div><span>Correction (m²/s)</span><b>${c.toExponential(2)}</b></div><div><span>Relative to D in the box</span><b>+${(100 * c / Dsi).toFixed(1)} %</b></div>
        <div><span>Corrected D (m²/s)</span><b>${(Dsi + c).toExponential(3)}</b></div><div><span>Atoms in the box (~100 per nm³)</span><b>${Math.round(100 * st.L ** 3).toLocaleString()}</b></div>`;
      $('[data-o=eta]').textContent = st.eta < 10 ? st.eta.toFixed(2) : st.eta.toFixed(0); $('[data-o=L]').textContent = st.L.toFixed(1); $('[data-o=D]').textContent = st.D < 10 ? st.D.toFixed(2) : st.D.toFixed(0);
    }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running || st.mode !== 'g') return; const n = Math.max(1, Math.round(4000 / st.N)); for (let k = 0; k < n; k++) frame(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
