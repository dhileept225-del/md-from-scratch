// Dielectric lab widget: 216 unit dipoles on a 6x6x6 periodic lattice, nearest-neighbour coupling J (kT = 1),
// rotational Brownian dynamics (D_r = 1) or random jumps. Running Kirkwood factor G_K = <M^2>/N, single-dipole C1 and C2,
// and the box-dipole correlation.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const NL = 6, N = NL ** 3, DT = 0.002, EVERY = 10, FDT = DT * EVERY, KMAX = 150, BUF = 3000, NSUB = 40;

export function mountLab14(root) {
  const st = { J: 0.4, jump: false, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view di-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Neighbour coupling J (k_BT)</span><output data-o="J"></output><input type="range" data-k="J" min="-0.6" max="0.6" step="0.05" value="0.4"></label>
    <div class="seg"><button data-r="0" class="on">Small random rotations</button><button data-r="1">Large jumps</button></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span>Kirkwood factor G_K (running)</span><b data-o="gk"></b></div>
      <div><span>Water-like dielectric constant</span><b data-o="eps"></b></div>
      <div><span>One dipole: time · P1/P2 ratio</span><b data-o="t1"></b></div>
      <div><span>Box-dipole time</span><b data-o="tm"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.di-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let seed = 5; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };
  const gauss = () => Math.sqrt(-2 * Math.log(rand())) * Math.cos(2 * Math.PI * rand());
  const nbr = []; for (let i = 0; i < NL; i++) for (let j = 0; j < NL; j++) for (let k = 0; k < NL; k++) { const id = (a, b, c) => ((a + NL) % NL) * NL * NL + ((b + NL) % NL) * NL + (c + NL) % NL; nbr.push([id(i + 1, j, k), id(i - 1, j, k), id(i, j + 1, k), id(i, j - 1, k), id(i, j, k + 1), id(i, j, k - 1)]); }
  let u, t, Mbuf, Ubuf, nb, m2sum, m2n, gTrace, acf = null;
  const randUnit = () => { const v = [gauss(), gauss(), gauss()], n = Math.hypot(...v); return v.map(q => q / n); };
  function reset() { u = Array.from({ length: N }, randUnit); t = 0; Mbuf = []; Ubuf = []; nb = 0; m2sum = 0; m2n = 0; gTrace = []; acf = null; }
  function step() {
    if (st.jump) { for (let i = 0; i < N; i++) if (rand() < 2 * DT) u[i] = randUnit(); }
    else {
      const s2 = Math.sqrt(2 * DT), nu = new Array(N);
      for (let i = 0; i < N; i++) {
        const ui = u[i]; let h = [0, 0, 0]; for (const j of nbr[i]) { h[0] += u[j][0]; h[1] += u[j][1]; h[2] += u[j][2]; } h = h.map(q => q * st.J);
        const hu = h[0] * ui[0] + h[1] * ui[1] + h[2] * ui[2], xi = [gauss(), gauss(), gauss()], xu = xi[0] * ui[0] + xi[1] * ui[1] + xi[2] * ui[2];
        const v = [0, 1, 2].map(k => ui[k] + (h[k] - hu * ui[k]) * DT + s2 * (xi[k] - xu * ui[k])), n = Math.hypot(...v); nu[i] = v.map(q => q / n);
      }
      u = nu;
    }
    t++;
    if (t % EVERY === 0 && t > 3000) {
      const M = [0, 0, 0]; u.forEach(v => { M[0] += v[0]; M[1] += v[1]; M[2] += v[2]; });
      Mbuf.push(M); Ubuf.push(u.slice(0, NSUB).map(v => v.slice())); if (Mbuf.length > BUF) { Mbuf.shift(); Ubuf.shift(); }
      m2sum += M[0] ** 2 + M[1] ** 2 + M[2] ** 2; m2n++; nb++; if (nb % 20 === 0) { gTrace.push(m2sum / m2n / N); if (gTrace.length > 400) gTrace.shift(); }
    }
  }
  function computeAcf() {
    const n = Mbuf.length; if (n < KMAX * 2) return; const cM = [], c1 = [], c2 = [];
    for (let k = 0; k < KMAX; k++) { let sM = 0, s1 = 0, s2 = 0, cnt = 0, cu = 0;
      for (let a = 0; a + k < n; a += 2) { const A = Mbuf[a], B = Mbuf[a + k]; sM += A[0] * B[0] + A[1] * B[1] + A[2] * B[2]; cnt++;
        for (let i = 0; i < NSUB; i++) { const x = Ubuf[a][i], y = Ubuf[a + k][i], c = x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; s1 += c; s2 += 1.5 * c * c - 0.5; cu++; } }
      cM.push(sM / cnt); c1.push(s1 / cu); c2.push(s2 / cu); }
    const area = arr => FDT * (arr.reduce((p, q) => p + q, 0) - 0.5 * arr[0]);
    acf = { cM: cM.map(q => q / cM[0]), c1, c2, t1: area(c1), t2: area(c2), tM: area(cM.map(q => q / cM[0])) };
  }
  $('[data-k=J]').oninput = e => { st.J = +e.target.value; reset(); draw(); };
  root.querySelectorAll('[data-r]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-r]').forEach(z => z.classList.toggle('on', z === b)); st.jump = b.dataset.r === '1'; $('[data-k=J]').disabled = st.jump; reset(); draw(); });
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !u) return;
    cx.clearRect(0, 0, W, H);
    const half = H * 0.45, bw = Math.min(W * 0.42, half - 16), cell = bw / NL;
    cx.strokeStyle = 'rgba(143,163,199,.35)'; cx.strokeRect(8, 8, bw, bw);
    for (let i = 0; i < NL; i++) for (let j = 0; j < NL; j++) { const v = u[i * NL * NL + j * NL], X = 8 + (i + 0.5) * cell, Y = 8 + (j + 0.5) * cell, L = cell * 0.42;
      cx.strokeStyle = `hsl(${190 + 150 * (v[2] + 1) / 2},90%,65%)`; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(X - v[0] * L, Y + v[1] * L); cx.lineTo(X + v[0] * L, Y - v[1] * L); cx.stroke();
      cx.fillStyle = cx.strokeStyle; cx.beginPath(); cx.arc(X + v[0] * L, Y - v[1] * L, 2.5, 0, 7); cx.fill(); }
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('one layer of dipoles (colour: z part)', 10, bw + 24);
    // G_K trace
    const gx0 = bw + 50, gx1 = W - 12, gy0 = 22, gy1 = half, gmax = Math.max(2, ...gTrace) * 1.15, GX = k => gx0 + k / 399 * (gx1 - gx0), GY = v => gy1 - v / gmax * (gy1 - gy0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(gx0, gy0); cx.lineTo(gx0, gy1); cx.lineTo(gx1, gy1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.fillText('running G_K = ⟨M²⟩/N', gx0 + 6, gy0 - 6);
    cx.setLineDash([4, 4]); cx.strokeStyle = 'rgba(242,245,250,.5)'; cx.beginPath(); cx.moveTo(gx0, GY(1)); cx.lineTo(gx1, GY(1)); cx.stroke(); cx.setLineDash([]);
    cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); gTrace.forEach((v, k) => k ? cx.lineTo(GX(k), GY(v)) : cx.moveTo(GX(k), GY(v))); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [0, 1].forEach(v => cx.fillText(v, gx0 - 4, GY(v) + 4));
    // correlation functions
    const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 28, BX = k => bx0 + k / (KMAX - 1) * (bx1 - bx0), BY = c => by1 - Math.max(0, Math.min(1, c)) * (by1 - by0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(bx0, by0); cx.lineTo(bx0, by1); cx.lineTo(bx1, by1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(W < 500 ? 'C₁ (cyan), C₂ (coral), box M (amber)' : 'one dipole: C₁ (cyan), C₂ (coral);  whole box: M(0)·M(t) (amber)', bx0 + 6, by0 - 10);
    if (acf) [[acf.c1, COL.cyan], [acf.c2, COL.coral], [acf.cM, COL.amber]].forEach(([arr, c]) => { cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath(); arr.forEach((v, k) => k ? cx.lineTo(BX(k), BY(v)) : cx.moveTo(BX(k), BY(v))); cx.stroke(); });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0, 1, 2, 3].forEach(tt => cx.fillText(tt, BX(tt / FDT), by1 + 14));
    const gk = m2n ? m2sum / m2n / N : NaN;
    o.J.textContent = st.jump ? 'off (jumps)' : st.J.toFixed(2); o.gk.textContent = m2n ? gk.toFixed(3) : '…'; o.eps.textContent = m2n ? (1 + 18.65 * gk).toFixed(1) : '…';
    o.t1.textContent = acf ? `${acf.t1.toFixed(2)} · ${(acf.t1 / acf.t2).toFixed(2)}` : '…'; o.tm.textContent = acf ? acf.tM.toFixed(2) : '…';
  }
  let raf = 0, visible = true, fc = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; for (let k = 0; k < 100; k++) step(); if (++fc % 15 === 0) computeAcf(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
