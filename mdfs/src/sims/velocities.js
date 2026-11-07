// Velocity initialiser: Gaussian draw -> remove drift -> rescale, plus collision relaxation to Maxwell–Boltzmann.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 2000;

export function mountVelocities(root) {
  const st = { T: 1.0, m: 1, collide: false, coll: 0 };
  const v = new Float64Array(3 * N);
  const gauss = () => Math.sqrt(-2 * Math.log(Math.random() + 1e-12)) * Math.cos(2 * Math.PI * Math.random());

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view vel-view"><canvas></canvas><div class="sim-msg" hidden></div></div>
  <div class="sim-ctrl">
    <label>Target temperature T* <output data-o="Tt"></output><input type="range" data-k="T" min="0.2" max="3" step="0.05" value="1"></label>
    <div class="seg"><button data-m="1" class="on">Light atoms (m = 1)</button><button data-m="4">Heavy (m = 4)</button></div>
    <div class="steps-btns">
      <button class="btn sm primary" data-k="draw">1 · Draw Gaussian velocities</button>
      <button class="btn sm ghost" data-k="drift">2 · Remove the drift</button>
      <button class="btn sm ghost" data-k="scale">3 · Scale to exact T</button>
    </div>
    <div class="row"><button class="btn sm ghost" data-k="same">Bad start: all same speed</button></div>
    <label class="check"><input type="checkbox" data-k="coll"> Collisions on (energy and momentum conserved)</label>
    <div class="readouts">
      <div><span>Temperature now</span><b data-o="T"></b></div>
      <div><span>Drift speed |v_cm|</span><b data-o="cm"></b></div>
      <div><span>Mean speed</span><b data-o="mean"></b></div>
      <div><span>Theory √(8T/πm)</span><b data-o="th"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.vel-view'), msg = root.querySelector('.sim-msg');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const say = h => { msg.hidden = false; msg.innerHTML = h; clearTimeout(say.t); say.t = setTimeout(() => msg.hidden = true, 4500); };

  const temp = () => { let s = 0; for (let i = 0; i < 3 * N; i++) s += v[i] * v[i]; return st.m * s / (3 * N - 3); };
  const com = () => { const c = [0, 0, 0]; for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) c[d] += v[3 * i + d]; return c.map(q => q / N); };
  function draw() { const s = Math.sqrt(st.T / st.m); for (let i = 0; i < 3 * N; i++) v[i] = s * gauss(); for (let i = 0; i < N; i++) v[3 * i] += 0.08 * s; say('Drawn from a Gaussian with variance T/m. A real draw always has a little net drift (here exaggerated slightly). Press 2.'); }
  function drift() { const c = com(); for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) v[3 * i + d] -= c[d]; say('Centre-of-mass velocity removed: the box as a whole no longer flies off. Press 3.'); }
  function scale() { const f = Math.sqrt(st.T / temp()); for (let i = 0; i < 3 * N; i++) v[i] *= f; say(`Every velocity multiplied by √(T_target/T_now) = ${f.toFixed(4)}. The temperature is now exact.`); }
  function same() { const sp = Math.sqrt(3 * st.T / st.m); for (let i = 0; i < N; i++) { let a = gauss(), b = gauss(), c = gauss(); const n = Math.hypot(a, b, c); v[3 * i] = sp * a / n; v[3 * i + 1] = sp * b / n; v[3 * i + 2] = sp * c / n; } drift(); scale(); st.coll = 0; say('Every atom has exactly the same speed: right temperature, wrong distribution. Tick "Collisions on" and watch.'); }
  function collide() {   // random binary collisions: rotate the relative velocity, keep the pair's momentum and energy
    for (let k = 0; k < N / 4; k++) {
      const i = Math.floor(Math.random() * N), j = Math.floor(Math.random() * N); if (i === j) continue;
      const c = [0, 1, 2].map(d => 0.5 * (v[3 * i + d] + v[3 * j + d])), g = Math.hypot(v[3 * i] - v[3 * j], v[3 * i + 1] - v[3 * j + 1], v[3 * i + 2] - v[3 * j + 2]) / 2;
      let a = gauss(), b = gauss(), e = gauss(); const n = Math.hypot(a, b, e); a /= n; b /= n; e /= n;
      v[3 * i] = c[0] + g * a; v[3 * i + 1] = c[1] + g * b; v[3 * i + 2] = c[2] + g * e;
      v[3 * j] = c[0] - g * a; v[3 * j + 1] = c[1] - g * b; v[3 * j + 2] = c[2] - g * e;
    }
    st.coll += N / 4;
  }
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; };
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.m = +b.dataset.m; root.querySelectorAll('[data-m]').forEach(x => x.classList.toggle('on', x === b)); say('Mass changed but the velocities did not, so the temperature changed. Redraw (1, 2, 3) for the new mass.'); });
  $('[data-k=draw]').onclick = draw; $('[data-k=drift]').onclick = drift; $('[data-k=scale]').onclick = scale; $('[data-k=same]').onclick = same;
  $('[data-k=coll]').onchange = e => { st.collide = e.target.checked; };

  let W = 600, H = 460;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function hist(vals, lo, hi, nb) { const h = new Float64Array(nb), w = (hi - lo) / nb; for (const q of vals) { const k = Math.floor((q - lo) / w); if (k >= 0 && k < nb) h[k]++; } for (let k = 0; k < nb; k++) h[k] /= vals.length * w; return h; }
  function panel(x0, y0, w, h, data, lo, hi, theory, label, legend) {
    const nb = data.length, top = Math.max(...data, ...Array.from({ length: 60 }, (_, k) => theory(lo + (k + 0.5) / 60 * (hi - lo)))) * 1.12 || 1;
    cx.font = '500 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    let lx = x0; legend.forEach(([t, c]) => { cx.fillStyle = c; cx.fillRect(lx, y0 - 11, 14, 4); cx.fillText(t, lx + 18, y0 - 6); lx += 30 + cx.measureText(t).width; });
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y0 + h); cx.lineTo(x0 + w, y0 + h); cx.stroke();
    cx.save(); cx.beginPath(); cx.rect(x0, y0, w, h); cx.clip();
    const bw = w / nb; cx.fillStyle = 'rgba(56,225,255,.75)';
    for (let k = 0; k < nb; k++) { const bh = data[k] / top * h; cx.fillRect(x0 + k * bw + 1, y0 + h - bh, bw - 2, bh); }
    cx.strokeStyle = COL.amber; cx.lineWidth = 2.5; cx.beginPath();
    for (let k = 0; k <= 120; k++) { const q = lo + k / 120 * (hi - lo), X = x0 + k / 120 * w, Yv = y0 + h - theory(q) / top * h; k ? cx.lineTo(X, Yv) : cx.moveTo(X, Yv); } cx.stroke();
    cx.restore();
    cx.fillStyle = COL.muted; cx.textAlign = 'right'; cx.fillText(label, x0 + w, y0 + h + 16);
  }
  function render() {
    if (st.collide) collide();
    cx.clearRect(0, 0, W, H);
    const s2 = st.T / st.m, vmax = 4.6 * Math.sqrt(s2), speeds = new Float64Array(N), vx = new Float64Array(N);
    for (let i = 0; i < N; i++) { speeds[i] = Math.hypot(v[3 * i], v[3 * i + 1], v[3 * i + 2]); vx[i] = v[3 * i]; }
    const MB = q => 4 * Math.PI * Math.pow(1 / (2 * Math.PI * s2), 1.5) * q * q * Math.exp(-q * q / (2 * s2));
    const GA = q => Math.exp(-q * q / (2 * s2)) / Math.sqrt(2 * Math.PI * s2);
    const x0 = 46, w = W - 66, h1 = (H - 90) * 0.6, h2 = (H - 90) * 0.4;
    panel(x0, 24, w, h1, hist(speeds, 0, vmax, 36), 0, vmax, MB, 'speed |v|', [['your atoms', COL.cyan], ['Maxwell–Boltzmann at target T', COL.amber]]);
    panel(x0, 24 + h1 + 44, w, h2, hist(vx, -vmax * 0.85, vmax * 0.85, 36), -vmax * 0.85, vmax * 0.85, GA, 'velocity in x  (0 in the middle)', [['vx', COL.cyan], ['Gaussian', COL.amber]]);
    const T = temp(), c = com(), mean = speeds.reduce((a, q) => a + q, 0) / N;
    o.Tt.textContent = st.T.toFixed(2);
    o.T.innerHTML = `<span style="color:${Math.abs(T / st.T - 1) < 0.005 ? COL.cyan : COL.amber}">${T.toFixed(3)}</span>`;
    o.cm.textContent = Math.hypot(...c).toExponential(1); o.mean.textContent = mean.toFixed(3); o.th.textContent = Math.sqrt(8 * st.T / (Math.PI * st.m)).toFixed(3);
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (visible) render(); };
  same(); msg.hidden = true; loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); clearTimeout(say.t); };
}
