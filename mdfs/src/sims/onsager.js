// Onsager lab (2D): a Li salt in an implicit solvent, made of free Li+, free A-, neutral LiA pairs and [LiA2]- triplets.
// Every cluster diffuses as one body. From the summed displacements R+ and R- we estimate L++, L--, L+-,
// the conductivity, and the transference number, and compare with the 'apparent' t+ from self-diffusion (as in NMR).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const NLI = 200, LBOX = 40, DT = 0.05, EVERY = 4, H = EVERY * DT, K = 30, DA = 0.8, DP = 0.5;

export function mountOnsager(root) {
  const st = { pair: 20, trip: 0, dli: 1.2, dtrip: 0.4, running: true };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Li⁺ in neutral LiA pairs</span><output data-o="pair"></output><input type="range" data-k="pair" min="0" max="80" step="5" value="20"></label>
    <label><span>Li⁺ in [LiA₂]⁻ triplets</span><output data-o="trip"></output><input type="range" data-k="trip" min="0" max="45" step="5" value="0"></label>
    <label><span><i class="nc">D</i> of free Li⁺ (free A⁻: 0.8)</span><output data-o="dli"></output><input type="range" data-k="dli" min="0.2" max="1.5" step="0.05" value="1.2"></label>
    <label><span><i class="nc">D</i> of triplets (pairs: 0.5)</span><output data-o="dtr"></output><input type="range" data-k="dtr" min="0.1" max="1.2" step="0.05" value="0.4"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">Reset average</button></div>
    <div class="readouts">
      <div><span><i class="nc">t</i>+ Onsager (exact)</span><b data-o="tp"></b></div>
      <div><span><i class="nc">t</i>+ from <i class="nc">D</i> (NMR-like)</span><b data-o="ta"></b></div>
      <div><span>L++ · L−− · L+−</span><b data-o="L"></b></div>
      <div><span><i class="nc">σ</i> / <i class="nc">σ</i><sub>NE</sub> · data</span><b data-o="sg"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 31; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let cl = [], buf = [], spp, smm, spm, cnt, warn = '';
  function spec() {   // cluster types [count, Li, A, D]
    let nT = Math.round(NLI * st.trip / 100), nP = Math.round(NLI * st.pair / 100); warn = '';
    let nF = NLI - nP - nT; if (nT > nF) { nT = Math.floor((NLI - nP) / 2); nF = NLI - nP - nT; warn = 'triplets limited: every [LiA₂]⁻ needs a free Li⁺ to keep the salt neutral'; }
    return [[nF, 1, 0, st.dli], [nF - nT, 0, 1, DA], [nP, 1, 1, DP], [nT, 1, 2, st.dtrip]];
  }
  function build() {
    cl = []; for (const [n, a, b, D] of spec()) for (let k = 0; k < n; k++) cl.push({ x: rand() * LBOX, y: rand() * LBOX, ux: 0, uy: 0, a, b, D, ph: rand() * 6.28 });
    resetAvg();
  }
  function resetAvg() { buf = []; spp = new Float64Array(K + 1); smm = new Float64Array(K + 1); spm = new Float64Array(K + 1); cnt = new Float64Array(K + 1); }
  let steps = 0;
  function step() {
    for (const c of cl) { const s = Math.sqrt(2 * c.D * DT), dx = s * gauss(), dy = s * gauss(); c.x += dx; c.y += dy; c.ux += dx; c.uy += dy; c.x -= LBOX * Math.floor(c.x / LBOX); c.y -= LBOX * Math.floor(c.y / LBOX); }
    if (++steps % EVERY === 0) sample();
  }
  function sample() {
    let Rp = [0, 0], Rm = [0, 0]; for (const c of cl) { Rp[0] += c.a * c.ux; Rp[1] += c.a * c.uy; Rm[0] += c.b * c.ux; Rm[1] += c.b * c.uy; }
    buf.push([Rp, Rm]); if (buf.length > K + 1) buf.shift();
    const m = buf.length - 1, [P, M] = buf[m];
    for (let k = 1; k <= m; k++) { const [P0, M0] = buf[m - k], dp = [P[0] - P0[0], P[1] - P0[1]], dm = [M[0] - M0[0], M[1] - M0[1]];
      spp[k] += dp[0] ** 2 + dp[1] ** 2; smm[k] += dm[0] ** 2 + dm[1] ** 2; spm[k] += dp[0] * dm[0] + dp[1] * dm[1]; cnt[k]++; }
  }
  const slope = s => { let N = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; for (let k = 5; k <= 25; k++) { if (!cnt[k]) continue; const t = k * H, y = s[k] / cnt[k]; N++; sx += t; sy += y; sxx += t * t; sxy += t * y; } return N > 8 ? (N * sxy - sx * sy) / (N * sxx - sx * sx) / 4 : null; };   // 2D: L = slope / 4
  function exact() {
    const S = spec(), L = { pp: 0, mm: 0, pm: 0 }; let Np = 0, Nm = 0, Dp = 0, Dm = 0;
    for (const [n, a, b, D] of S) { L.pp += n * a * a * D; L.mm += n * b * b * D; L.pm += n * a * b * D; Np += n * a; Nm += n * b; Dp += n * a * D; Dm += n * b * D; }
    return { L, tapp: Dp / (Dp + Dm), sNE: Dp + Dm };
  }
  const tplus = L => (L.pp - L.pm) / (L.pp + L.mm - 2 * L.pm), sigma = L => L.pp + L.mm - 2 * L.pm;

  const bindS = (k, key, out, f) => { $(`[data-k=${k}]`).oninput = e => { st[key] = +e.target.value; build(); }; };
  bindS('pair', 'pair'); bindS('trip', 'trip'); bindS('dli', 'dli'); bindS('dtr', 'dtrip');
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => resetAvg();

  let W = 600, Hc = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; Hc = view.clientHeight; cv.width = W * d; cv.height = Hc * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const ion = (x, y, r, pos) => { cx.fillStyle = pos ? COL.coral : COL.cyan; cx.beginPath(); cx.arc(x, y, r, 0, 7); cx.fill(); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, Hc);
    const topH = Hc * 0.5, S = Math.min(W - 20, topH - 12) / LBOX, ox = (W - LBOX * S) / 2, oy = 6, r = Math.max(2.5, 0.42 * S);
    cx.strokeStyle = 'rgba(56,225,255,.35)'; cx.lineWidth = 1.5; cx.strokeRect(ox, oy, LBOX * S, LBOX * S);
    for (const c of cl) { const X = ox + c.x * S, Y = oy + c.y * S;
      if (c.a && c.b) { cx.strokeStyle = c.b === 2 ? 'rgba(255,194,75,.9)' : 'rgba(242,245,250,.7)'; cx.lineWidth = 1.5; cx.beginPath(); cx.arc(X, Y, r * (c.b === 2 ? 2.6 : 2.1), 0, 7); cx.stroke(); }
      if (c.a && !c.b) ion(X, Y, r, true); else if (!c.a) ion(X, Y, r, false);
      else if (c.b === 1) { ion(X - r, Y, r, true); ion(X + r, Y, r, false); }
      else { ion(X, Y - r * 0.6, r, true); ion(X - r * 1.4, Y + r * 0.8, r, false); ion(X + r * 1.4, Y + r * 0.8, r, false); } }
    if (warn) { cx.fillStyle = COL.amber; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(warn, W / 2, topH + 2); }
    // cross-MSD curves
    const x0 = 46, x1 = W - 14, yt = topH + 34, yb = Hc - 26, tmax = K * H;
    const cur = s => Array.from({ length: K + 1 }, (_, k) => cnt[k] ? s[k] / cnt[k] / 4 : null), C1 = cur(spp), C2 = cur(smm), C3 = cur(spm);
    const vals = [...C1, ...C2, ...C3].filter(v => v != null), ymax = Math.max(1, ...vals) * 1.1, ymin = Math.min(0, ...vals) * 1.1;
    const X = t => x0 + t / tmax * (x1 - x0), Y = q => yb - (q - ymin) / (ymax - ymin) * (yb - yt);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, yt); cx.lineTo(x0, yb); cx.lineTo(x1, yb); cx.stroke();
    if (ymin < 0) { cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke(); cx.setLineDash([]); }
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('⟨ΔR_i · ΔR_j⟩ / 4   (slopes = Onsager coefficients)', x0 + 6, yt - 10);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 2, 3, 4, 5].forEach(t => cx.fillText(t, X(t), yb + 14)); cx.textAlign = 'right'; cx.fillText('lag t', x1, yb + 14);
    const curve = (arr, c) => { cx.strokeStyle = c; cx.lineWidth = 2.5; cx.beginPath(); let on = false; arr.forEach((v, k) => { if (v == null) return; on ? cx.lineTo(X(k * H), Y(v)) : cx.moveTo(X(k * H), Y(v)); on = true; }); cx.stroke(); };
    if (cnt[1] > 3) { curve(C1, COL.coral); curve(C2, COL.cyan); curve(C3, COL.amber); cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; [[COL.coral, 'Li–Li (L++)'], [COL.cyan, 'A–A (L−−)'], [COL.amber, 'Li–A (L+−)']].forEach(([c, t], j) => { cx.fillStyle = c; cx.fillText(t, x0 + 10, yt + 12 + 15 * j); }); }
    // readouts
    o.pair.textContent = st.pair + ' %'; o.trip.textContent = st.trip + ' %'; o.dli.textContent = st.dli.toFixed(2); o.dtr.textContent = st.dtrip.toFixed(2);
    const E = exact(), te = tplus(E.L);
    const L = { pp: slope(spp), mm: slope(smm), pm: slope(spm) };
    o.ta.textContent = E.tapp.toFixed(2);
    if (L.pp != null) {
      o.tp.textContent = `${tplus(L).toFixed(2)} (${te.toFixed(2)})`;
      o.L.textContent = `${L.pp.toFixed(0)} · ${L.mm.toFixed(0)} · ${L.pm.toFixed(0)}`;
      o.sg.textContent = `${(sigma(L) / E.sNE).toFixed(2)} · ${(cnt[5] * H).toFixed(0)} t`;
    } else { o.tp.textContent = `… (${te.toFixed(2)})`; o.L.textContent = o.sg.textContent = 'collecting…'; }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) for (let k = 0; k < 20; k++) step(); draw(); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
