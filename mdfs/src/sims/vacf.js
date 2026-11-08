// VACF lab: 256 Lennard-Jones atoms in 3D. Velocity autocorrelation C(t) = <v(t).v(0)>, averaged over atoms
// and all time origins, built up live. Bottom panel: the Green-Kubo running integral D(t) = 1/3 int C dt,
// or the vibrational density of states (cosine transform of C).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 256, DT = 0.005, RC2 = 6.25, K = 150, EQ = 800;
const STATES = { gas: { rho: 0.05, T: 2.0, every: 10 }, liquid: { rho: 0.8442, T: 0.722, every: 2 }, solid: { rho: 1.1, T: 0.3, every: 2 } };

export function mountVacf(root) {
  const st = { state: 'liquid', T: 0.722, running: true, view: 'gk' };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-s="gas">Gas</button><button data-s="liquid" class="on">Liquid</button><button data-s="solid">Solid</button></div>
    <label><span>Temperature <i class="nc">T</i>*</span><output data-o="T0"></output><input type="range" data-k="T" min="0.2" max="2.5" step="0.05" value="0.722"></label>
    <div class="seg"><button data-v="gk" class="on">Green–Kubo D(t)</button><button data-v="vdos">Spectrum</button></div>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset average</button></div>
    <div class="readouts">
      <div><span>Green–Kubo <i class="nc">D</i>*</span><b data-o="D"></b></div>
      <div><span>First zero of C(t)</span><b data-o="z"></b></div>
      <div><span>Minimum C/C(0)</span><b data-o="mn"></b></div>
      <div><span>Measured <i class="nc">T</i>* / origins</span><b data-o="T"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let L = 1, x = new Float64Array(3 * N), v = new Float64Array(3 * N), f = new Float64Array(3 * N), Tnow = 0, steps = 0, buf = [], sum, cnt;
  let rng = 777; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  const every = () => STATES[st.state].every, hh = () => every() * DT;

  function build() {
    const s = STATES[st.state], a = Math.cbrt(4 / s.rho), B = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]; L = 4 * a; let k = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let l = 0; l < 4; l++) for (const b of B) { x[3 * k] = (i + b[0]) * a; x[3 * k + 1] = (j + b[1]) * a; x[3 * k + 2] = (l + b[2]) * a; k++; }
    const m = [0, 0, 0]; for (let i = 0; i < 3 * N; i++) { v[i] = gauss() * Math.sqrt(st.T); m[i % 3] += v[i] / N; } for (let i = 0; i < 3 * N; i++) v[i] -= m[i % 3];
    forces(); steps = 0; resetAvg();
  }
  function resetAvg() { buf = []; sum = new Float64Array(K + 1); cnt = new Float64Array(K + 1); }
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
    Tnow = ke / (3 * N - 3);
    if (steps < EQ) { const lam = steps < EQ / 2 ? Math.sqrt(1 + DT / 0.1 * (st.T / Tnow - 1)) : Math.sqrt(st.T / Tnow); for (let i = 0; i < 3 * N; i++) v[i] *= lam; }   // thermostat only while equilibrating; NVE after
    steps++;
    if (steps === EQ) resetAvg();
    if (steps > EQ && steps % every() === 0) sample();
  }
  function sample() {   // new frame = end point of every lag 0..K
    buf.push(Float64Array.from(v)); if (buf.length > K + 1) buf.shift();
    const n = buf.length - 1, c = buf[n];
    for (let k = 0; k <= n; k++) { const p = buf[n - k]; let s = 0; for (let i = 0; i < 3 * N; i++) s += c[i] * p[i]; sum[k] += s / N; cnt[k]++; }
  }

  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { st.state = b.dataset.s; st.T = STATES[st.state].T; $('[data-k=T]').value = st.T; root.querySelectorAll('[data-s]').forEach(z => z.classList.toggle('on', z === b)); build(); });
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; steps = 0; };
  root.querySelectorAll('[data-v]').forEach(b => b.onclick = () => { st.view = b.dataset.v; root.querySelectorAll('[data-v]').forEach(z => z.classList.toggle('on', z === b)); draw(); });
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => resetAvg();

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function panel(yt, yb, xmax, ylo, yhi, title, xlab, yt_) {
    const x0 = 46, x1 = W - 14, X = t => x0 + t / xmax * (x1 - x0), Y = q => yb - (q - ylo) / (yhi - ylo) * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    if (ylo < 0) { cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke(); cx.setLineDash([]); }
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(title, x0 + 6, yt - 8);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText(xlab, x1, yb + 14);
    cx.textAlign = 'center'; const xs = xmax > 20 ? [10, 20, 30] : xmax > 3 ? [2, 4, 6] : [0.5, 1.0]; xs.forEach(t => { if (t < xmax * 0.9) cx.fillText(t, X(t), yb + 14); });
    cx.textAlign = 'right'; yt_.forEach(q => cx.fillText(+q.toPrecision(2), x0 - 5, Y(q) + 4));
    return { X, Y };
  }
  const line = (xs, ys, P, c, w = 2.5) => { cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); xs.forEach((t, k) => k ? cx.lineTo(P.X(t), P.Y(ys[k])) : cx.moveTo(P.X(t), P.Y(ys[k]))); cx.stroke(); };

  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const h = hh(), n = cnt ? Array.from(cnt).findIndex(c => c === 0) : 0, nk = n < 0 ? K + 1 : n, ok = nk > 5 && cnt[0] > 20;
    const C = ok ? Array.from({ length: nk }, (_, k) => sum[k] / cnt[k]) : [], t = C.map((_, k) => k * h), Cn = C.map(c => c / C[0]);
    const Dt = [0]; for (let k = 1; k < C.length; k++) Dt.push(Dt[k - 1] + 0.5 * (C[k] + C[k - 1]) * h / 3);
    const xmax = K * h, mid = H * 0.5;
    const P1 = panel(28, mid - 26, xmax, st.state === 'solid' ? -0.7 : -0.3, 1.05, 'C(t) / C(0) = ⟨v(t)·v(0)⟩ / ⟨v²⟩', 't (τ)', st.state === 'solid' ? [-0.5, 0, 0.5, 1] : [0, 0.5, 1]);
    if (!ok) { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '13px "Space Grotesk", sans-serif'; cx.fillText('equilibrating…', W / 2, mid / 2); }
    else line(t, Cn, P1, COL.cyan);
    if (st.view === 'gk') {
      const D = Dt.length ? Dt[Dt.length - 1] : 0, dmax = Math.max(0.01, ...Dt) * 1.15;
      const P2 = panel(mid + 30, H - 26, xmax, 0, dmax, 'D(t) = ⅓ ∫₀ᵗ C(s) ds  (Green–Kubo)', 'upper limit (τ)', [dmax / 2.3, dmax / 1.15].map(q => +q.toPrecision(2)));
      if (ok) { line(t, Dt, P2, COL.amber); cx.fillStyle = COL.amber; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; const rising = Dt.length > 20 && Math.abs(D - Dt[Math.floor(Dt.length * 0.8)]) > 0.03 * Math.abs(D); cx.fillText(`${rising ? 'still rising: integrate longer!' : 'plateau →'} D* = ${D.toFixed(st.state === 'gas' ? 2 : 4)}`, W - 16, P2.Y(D) - 8); }
    } else {
      const wmax = st.state === 'gas' ? 6 : 40, ws = Array.from({ length: 200 }, (_, j) => j * wmax / 199);
      const S = ok ? ws.map(w => { let s = 0; for (let k = 0; k < Cn.length; k++) { const win = 0.5 * (1 + Math.cos(Math.PI * k / Cn.length)); s += (k === 0 ? 0.5 : 1) * Cn[k] * win * Math.cos(w * k * h); } return 2 * h * s; }) : [];
      const smax = Math.max(0.05, ...S) * 1.1;
      const P2 = panel(mid + 30, H - 26, wmax, 0, smax, 'vibrational density of states S(ω)', 'ω (1/τ)', [smax / 2.2, smax / 1.1].map(q => +q.toPrecision(2)));
      if (ok) { line(ws, S, P2, COL.coral); cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(`S(0) = ${S[0].toFixed(3)}: nonzero only if atoms diffuse`, 52, mid + 46); }
    }
    o.T0.textContent = st.T.toFixed(2); o.T.textContent = `${Tnow.toFixed(2)} / ${cnt ? cnt[0] | 0 : 0}`;
    if (ok) {
      const z = Cn.findIndex(c => c < 0), mn = Math.min(...Cn), D = Dt[Dt.length - 1];
      o.D.textContent = Math.abs(D) < 1e-3 ? '≈ 0 (no diffusion)' : D.toFixed(st.state === 'gas' ? 2 : 4) + (nk <= K ? ' (filling)' : '');
      o.z.textContent = z < 0 ? 'never (in this window)' : `t = ${(z * h).toFixed(2)} τ`;
      o.mn.textContent = mn.toFixed(3);
    } else o.D.textContent = o.z.textContent = o.mn.textContent = '…';
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 6; k++) step(); draw(); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
