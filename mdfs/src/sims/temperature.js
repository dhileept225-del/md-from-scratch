// Temperature lab: 2D Lennard-Jones NVE. Instantaneous kinetic temperature, its fluctuations, equipartition,
// the drop after a lattice start, and velocity rescaling.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RHO = 0.8, RC = 2.5, DT = 0.005, SUB = 4, WIN = 300, HIST = 600;

function heat(x) {
  x = Math.min(1, Math.max(0, x));
  const S = [[56, 225, 255], [255, 194, 75], [255, 107, 107]], k = x < 0.5 ? 0 : 1, f = x < 0.5 ? x * 2 : x * 2 - 1;
  return `rgb(${S[k].map((c, i) => Math.round(c + (S[k + 1][i] - c) * f)).join(',')})`;
}
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountTemperature(root) {
  const st = { N: 64, T0: 1.0, running: true };
  let N, Lx, Ly, x, v, f, hist = [], marks = [];
  const Nf = () => 2 * N - 2;                       // 2D: 2N, minus 2 for fixed total momentum
  const temp = () => { let s = 0; for (let i = 0; i < 2 * N; i++) s += v[i] * v[i]; return s / Nf(); };
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
  function setT(T) { const s = Math.sqrt(T / temp()); for (let i = 0; i < 2 * N; i++) v[i] *= s; }
  function restart() {
    N = st.N; x = new Float64Array(2 * N); v = new Float64Array(2 * N); f = new Float64Array(2 * N);
    const n = Math.round(Math.sqrt(N)), a = Math.sqrt(2 / (Math.sqrt(3) * RHO)), h = a * Math.sqrt(3) / 2;   // triangular crystal
    Lx = n * a; Ly = n * h;
    for (let i = 0; i < N; i++) { const r = Math.floor(i / n); x[2 * i] = (i % n + 0.25 + 0.5 * (r % 2)) * a; x[2 * i + 1] = (r + 0.5) * h; }
    let px = 0, py = 0; for (let i = 0; i < N; i++) { v[2 * i] = gauss(); v[2 * i + 1] = gauss(); px += v[2 * i]; py += v[2 * i + 1]; }
    for (let i = 0; i < N; i++) { v[2 * i] -= px / N; v[2 * i + 1] -= py / N; }
    setT(st.T0); forces(); hist = [temp()]; marks = [{ k: 0, T: st.T0, txt: 'start' }];
  }
  function step() {
    for (let s = 0; s < SUB; s++) {
      for (let i = 0; i < 2 * N; i++) { const L = i % 2 ? Ly : Lx; v[i] += 0.5 * DT * f[i]; x[i] += DT * v[i]; x[i] -= L * Math.floor(x[i] / L); }
      forces(); for (let i = 0; i < 2 * N; i++) v[i] += 0.5 * DT * f[i];
    }
    hist.push(temp()); if (hist.length > HIST) { hist.shift(); marks.forEach(m => m.k--); marks = marks.filter(m => m.k > -HIST); }
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label>Atoms N <output data-o="N"></output><input type="range" data-k="N" min="0" max="5" step="1" value="1"></label>
    <label>Starting temperature <output data-o="T0"></output><input type="range" data-k="T0" min="0.2" max="2.5" step="0.1" value="1"></label>
    <div class="row"><button class="btn sm primary" data-k="restart">Restart on a lattice</button><button class="btn sm ghost" data-k="run">❚❚ Pause</button></div>
    <div class="row"><button class="btn sm ghost" data-k="hot">Heat: v × 1.1</button><button class="btn sm ghost" data-k="cold">Cool: v × 0.9</button></div>
    <div class="readouts">
      <div><span>Instant T</span><b data-o="T"></b></div>
      <div><span>Average T (recent)</span><b data-o="avg"></b></div>
      <div><span>Fluctuation std/<i class="nc">⟨T⟩</i></span><b data-o="rel"></b></div>
      <div><span>Fluctuation × √N</span><b data-o="sq"></b></div>
      <div class="wide"><span>Equipartition <i class="nc">⟨v<sub>x</sub>²⟩ : ⟨v<sub>y</sub>²⟩</i></span><b data-o="eq"></b></div>
    </div>
  </div>`;
  const Ns = [36, 64, 100, 144, 256, 400];
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let vx2 = 0, vy2 = 0, nEq = 0;
  const resetEq = () => { vx2 = vy2 = nEq = 0; };
  $('[data-k=N]').oninput = e => { st.N = Ns[+e.target.value]; restart(); resetEq(); };
  $('[data-k=T0]').oninput = e => { st.T0 = +e.target.value; o.T0.textContent = st.T0.toFixed(1); };
  $('[data-k=restart]').onclick = () => { restart(); resetEq(); };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  const kick = s => { for (let i = 0; i < 2 * N; i++) v[i] *= s; hist.push(temp()); marks.push({ k: hist.length - 1, txt: s > 1 ? 'heat' : 'cool' }); resetEq(); };
  $('[data-k=hot]').onclick = () => kick(1.1); $('[data-k=cold]').onclick = () => kick(0.9);

  let W = 600, H = 560;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  const recent = () => { const k0 = Math.max(0, hist.length - WIN, ...marks.map(m => m.k + 50)); return hist.slice(Math.min(k0, hist.length - 1)); };
  function draw() {
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(150, H * 0.34), S = Math.min((W - 24) / Lx, (H - plotH - 34) / Ly), ox = (W - Lx * S) / 2, oy = 12;
    cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.lineWidth = 1.5; cx.strokeRect(ox, oy, Lx * S, Ly * S);
    const Tn = Math.max(0.5, st.T0 * 1.5), r = Math.max(2.5, 0.4 * S);
    for (let i = 0; i < N; i++) {
      const sp2 = (v[2 * i] ** 2 + v[2 * i + 1] ** 2) / 2;
      cx.fillStyle = heat(sp2 / (2 * Tn)); cx.beginPath(); cx.arc(ox + x[2 * i] * S, oy + x[2 * i + 1] * S, r, 0, 7); cx.fill();
    }
    // temperature trace
    const px0 = 46, px1 = W - 14, py1 = H - 22, py0 = oy + Ly * S + 44;
    const hi = Math.max(1.2 * Math.max(...hist), st.T0 * 1.2, 0.5), Y = T => py1 - T / hi * (py1 - py0), X = k => px0 + (px1 - px0) * k / (HIST - 1);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(px0, py0); cx.lineTo(px0, py1); cx.lineTo(px1, py1); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right';
    [0, hi / 2, hi].forEach(t => cx.fillText(t.toFixed(1), px0 - 6, Y(t) + 4));
    cx.textAlign = 'left'; cx.fillStyle = COL.cyan; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillText('instant T', px0 + 8, py0 - 6);
    cx.fillStyle = COL.amber; cx.fillText('recent average', px0 + 90, py0 - 6);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 1.5; cx.beginPath(); hist.forEach((T, k) => k ? cx.lineTo(X(k), Y(T)) : cx.moveTo(X(k), Y(T))); cx.stroke();
    const rc = recent(), m = rc.reduce((s, q) => s + q, 0) / rc.length;
    cx.strokeStyle = COL.amber; cx.setLineDash([6, 5]); cx.lineWidth = 2; cx.beginPath(); cx.moveTo(X(hist.length - rc.length), Y(m)); cx.lineTo(X(hist.length - 1), Y(m)); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif';
    marks.forEach(mk => { if (mk.k < 0) return; cx.strokeStyle = 'rgba(255,107,107,.5)'; cx.beginPath(); cx.moveTo(X(mk.k), py0); cx.lineTo(X(mk.k), py1); cx.stroke(); cx.fillText(mk.txt, X(mk.k) + 4, py1 - 6); });
  }
  function readouts() {
    o.N.textContent = st.N; o.T0.textContent = st.T0.toFixed(1);
    const rc = recent(), m = rc.reduce((s, q) => s + q, 0) / rc.length, sd = Math.sqrt(rc.reduce((s, q) => s + (q - m) ** 2, 0) / rc.length);
    o.T.textContent = hist[hist.length - 1].toFixed(3); o.avg.textContent = m.toFixed(3);
    const enough = rc.length > 100;
    o.rel.textContent = enough ? (100 * sd / m).toFixed(1) + ' %' : 'measuring…';
    o.sq.textContent = enough ? (sd / m * Math.sqrt(N)).toFixed(2) : 'measuring…';
    o.eq.textContent = nEq > 50 ? `${(vx2 / nEq).toFixed(3)} : ${(vy2 / nEq).toFixed(3)}` : 'measuring…';
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => {
    raf = requestAnimationFrame(loop); if (!visible) return;
    if (st.running) { step(); let a = 0, b = 0; for (let i = 0; i < N; i++) { a += v[2 * i] ** 2; b += v[2 * i + 1] ** 2; } vx2 += a / N; vy2 += b / N; nEq++; }
    draw(); if (frame++ % 8 === 0) readouts();
  };
  restart(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
