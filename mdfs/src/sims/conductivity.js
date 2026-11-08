// Conductivity lab (2D): cations and anions diffuse in an implicit solvent. A free cation that meets a free anion
// can bind into a neutral pair; the pair moves as one body and breaks after a random lifetime.
// Compares the collective charge MSD (Einstein-Helfand conductivity) with the self terms only (Nernst-Einstein).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const n = 120, L = 30, A = L * L, DT = 0.01, RB = 1.0, EVERY = 10, K = 100, H = EVERY * DT;

export function mountConductivity(root) {
  const st = { kon: 5, life: 1.0, dpair: 0.5, running: true };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Pairing rate on contact</span><output data-o="kon"></output><input type="range" data-k="kon" min="0" max="40" step="1" value="5"></label>
    <label><span>Pair lifetime</span><output data-o="life"></output><input type="range" data-k="life" min="0.2" max="5" step="0.1" value="1.0"></label>
    <label><span>Pair diffusion / free diffusion</span><output data-o="dp"></output><input type="range" data-k="dp" min="0.2" max="1" step="0.05" value="0.5"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset average</button></div>
    <div class="readouts">
      <div><span>Ions in pairs</span><b data-o="fp"></b></div>
      <div><span>Self diffusion <i class="nc">D</i></span><b data-o="D"></b></div>
      <div><span>Ionicity <i class="nc">σ</i> / <i class="nc">σ</i><sub>NE</sub></span><b data-o="ion"></b></div>
      <div><span>Haven ratio · data so far</span><b data-o="hv"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 2024; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  const x = new Float64Array(4 * n), u = new Float64Array(4 * n), partner = new Int32Array(n).fill(-1), taken = new Uint8Array(n);   // ions 0..n-1: +, n..2n-1: -
  let buf = [], sumM, sumS, cnt, steps = 0, fpAvg = 0, fpN = 0;
  function init() { for (let i = 0; i < 4 * n; i++) { x[i] = rand() * L; } u.set(x); partner.fill(-1); taken.fill(0); resetAvg(); }
  function resetAvg() { buf = []; sumM = new Float64Array(K + 1); sumS = new Float64Array(K + 1); cnt = new Float64Array(K + 1); fpAvg = 0; fpN = 0; }
  const mi = d => d - L * Math.round(d / L);
  function step() {
    const sf = Math.sqrt(2 * DT), sp = Math.sqrt(2 * st.dpair * DT);
    for (let i = 0; i < n; i++) if (partner[i] >= 0 && rand() < DT / st.life) { taken[partner[i]] = 0; partner[i] = -1; }   // pairs break
    if (st.kon > 0) for (let i = 0; i < n; i++) { if (partner[i] >= 0 || rand() > st.kon * DT) continue;                     // free cation tries to bind
      let best = -1, bd = RB; for (let j = 0; j < n; j++) { if (taken[j]) continue; const d = Math.hypot(mi(x[2 * (n + j)] - x[2 * i]), mi(x[2 * (n + j) + 1] - x[2 * i + 1])); if (d < bd) { bd = d; best = j; } }
      if (best >= 0) { partner[i] = best; taken[best] = 1; } }
    const moved = new Uint8Array(2 * n);
    for (let i = 0; i < n; i++) if (partner[i] >= 0) { const j = n + partner[i], gx = sp * gauss(), gy = sp * gauss(); for (const k of [i, j]) { x[2 * k] += gx; x[2 * k + 1] += gy; u[2 * k] += gx; u[2 * k + 1] += gy; moved[k] = 1; } }
    for (let k = 0; k < 2 * n; k++) if (!moved[k]) { const gx = sf * gauss(), gy = sf * gauss(); x[2 * k] += gx; x[2 * k + 1] += gy; u[2 * k] += gx; u[2 * k + 1] += gy; }
    for (let k = 0; k < 4 * n; k++) x[k] -= L * Math.floor(x[k] / L);
    steps++; if (steps % EVERY === 0) sample();
  }
  function sample() {
    let Mx = 0, My = 0; for (let k = 0; k < 2 * n; k++) { const q = k < n ? 1 : -1; Mx += q * u[2 * k]; My += q * u[2 * k + 1]; }
    buf.push({ u: Float64Array.from(u), M: [Mx, My] }); if (buf.length > K + 1) buf.shift();
    const m = buf.length - 1, c = buf[m];
    for (let k = 1; k <= m; k++) { const p = buf[m - k]; sumM[k] += (c.M[0] - p.M[0]) ** 2 + (c.M[1] - p.M[1]) ** 2;
      let s = 0; for (let i = 0; i < 4 * n; i++) { const d = c.u[i] - p.u[i]; s += d * d; } sumS[k] += s; cnt[k]++; }
    let np = 0; for (let i = 0; i < n; i++) if (partner[i] >= 0) np++; fpAvg += np / n; fpN++;
  }
  const fit = (s) => { let N = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; for (let k = 10; k <= 60; k++) { if (!cnt[k]) continue; const t = k * H, y = s[k] / cnt[k]; N++; sx += t; sy += y; sxx += t * t; sxy += t * y; } return N > 10 ? (N * sxy - sx * sy) / (N * sxx - sx * sx) : null; };

  const bind = (k, key) => { $(`[data-k=${k}]`).oninput = e => { st[key] = +e.target.value; resetAvg(); }; };
  bind('kon', 'kon'); bind('life', 'life'); bind('dp', 'dpair');
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => resetAvg();

  let W = 600, Hc = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; Hc = view.clientHeight; cv.width = W * d; cv.height = Hc * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, Hc);
    const topH = Hc * 0.52, S = Math.min(W - 20, topH - 10) / L, ox = (W - L * S) / 2, oy = 6;
    cx.strokeStyle = 'rgba(56,225,255,.35)'; cx.lineWidth = 1.5; cx.strokeRect(ox, oy, L * S, L * S);
    cx.strokeStyle = 'rgba(242,245,250,.8)'; cx.lineWidth = 2;
    for (let i = 0; i < n; i++) if (partner[i] >= 0) { const j = n + partner[i], dx = mi(x[2 * j] - x[2 * i]), dy = mi(x[2 * j + 1] - x[2 * i + 1]); cx.beginPath(); cx.moveTo(ox + x[2 * i] * S, oy + x[2 * i + 1] * S); cx.lineTo(ox + (x[2 * i] + dx) * S, oy + (x[2 * i + 1] + dy) * S); cx.stroke(); }
    for (let k = 0; k < 2 * n; k++) { cx.fillStyle = k < n ? COL.coral : COL.cyan; cx.beginPath(); cx.arc(ox + x[2 * k] * S, oy + x[2 * k + 1] * S, Math.max(2.5, 0.4 * S), 0, 7); cx.fill(); }
    // plot: charge MSD vs self-terms-only
    const x0 = 46, x1 = W - 14, yt = topH + 34, yb = Hc - 26, tmax = K * H;
    const MM = Array.from({ length: K + 1 }, (_, k) => cnt[k] ? sumM[k] / cnt[k] / (4 * A) : null), SS = Array.from({ length: K + 1 }, (_, k) => cnt[k] ? sumS[k] / cnt[k] / (4 * A) : null);
    const ymax = Math.max(0.1, ...MM.filter(v => v != null), ...SS.filter(v => v != null)) * 1.1, X = t => x0 + t / tmax * (x1 - x0), Y = q => yb - q / ymax * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('⟨|ΔM|²⟩ / 4A   (slope = conductivity, units e²/kT)', x0 + 6, yt - 10);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [2, 4, 6, 8].forEach(t => cx.fillText(t, X(t), yb + 14)); cx.textAlign = 'right'; cx.fillText('lag t', x1, yb + 14);
    const curve = (arr, c, dash) => { cx.setLineDash(dash || []); cx.strokeStyle = c; cx.lineWidth = 2.5; cx.beginPath(); let on = false; arr.forEach((v, k) => { if (v == null) return; on ? cx.lineTo(X(k * H), Y(v)) : cx.moveTo(X(k * H), Y(v)); on = true; }); cx.stroke(); cx.setLineDash([]); };
    if (cnt[1] > 3) { curve(SS, COL.cyan, [6, 4]); curve(MM, COL.coral); cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = COL.coral; cx.fillText('all terms (Einstein–Helfand σ)', x0 + 10, yt + 12); cx.fillStyle = COL.cyan; cx.fillText('self terms only (Nernst–Einstein σ_NE)', x0 + 10, yt + 27); }
    else { cx.fillStyle = COL.muted; cx.textAlign = 'center'; cx.font = '13px "Space Grotesk", sans-serif'; cx.fillText('collecting…', (x0 + x1) / 2, (yt + yb) / 2); }
    // readouts
    o.kon.textContent = st.kon; o.life.textContent = st.life.toFixed(1); o.dp.textContent = st.dpair.toFixed(2);
    let np = 0; for (let i = 0; i < n; i++) if (partner[i] >= 0) np++;
    const fbar = fpN ? fpAvg / fpN : np / n; o.fp.textContent = `${Math.round(100 * np / n)} % now (${Math.round(100 * fbar)} % avg)`;
    const sM = fit(sumM), sS = fit(sumS);
    if (sM != null && sS != null) {
      const D = sS / (4 * 2 * n), ionic = sM / sS, exact = (1 - fbar) / ((1 - fbar) + fbar * st.dpair);
      o.D.textContent = D.toFixed(3) + ' (free: 1)'; o.ion.textContent = `${ionic.toFixed(2)} (model: ${exact.toFixed(2)})`; o.hv.textContent = `${(1 / Math.max(0.01, ionic)).toFixed(2)} · ${(cnt[1] * H).toFixed(0)} time units`;
    } else o.D.textContent = o.ion.textContent = o.hv.textContent = 'collecting…';
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 30; k++) step(); draw(); };
  init(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
