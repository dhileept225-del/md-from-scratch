// Polarization lab: Li+ ions near polarizable sites (alpha = 1.45 A^3, water-like), in a plane.
// Two models: point induced dipoles (solved exactly) and Drude oscillators (relaxed iteratively).
// Compares the true many-body induction energy with the sum of two-body (one ion, one site) terms,
// and shows the polarization catastrophe when Thole damping is switched off and sites get too close.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const K = 332.06, ALPHA = 1.45, KD = 1000, QD = -Math.sqrt(ALPHA * KD / K), DEBYE = 4.803;
const A_DIP = 0.39, A_DRUDE = 2.6 / Math.cbrt(ALPHA);   // Thole parameters (AMOEBA style, CHARMM Drude style)

// ---------- induced point dipoles: (1/alpha - T) mu = E0, solved with Cholesky (fails = catastrophe) ----------
function solveDipoles(ions, sites, damp) {
  const n = sites.length, N = 2 * n, A = Array.from({ length: N }, () => new Float64Array(N)), b = new Float64Array(N);
  for (let i = 0; i < n; i++) {
    A[2 * i][2 * i] = A[2 * i + 1][2 * i + 1] = 1 / ALPHA;
    for (const q of ions) { const dx = sites[i][0] - q[0], dy = sites[i][1] - q[1], r = Math.hypot(dx, dy); b[2 * i] += dx / r ** 3; b[2 * i + 1] += dy / r ** 3; }
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const dx = sites[i][0] - sites[j][0], dy = sites[i][1] - sites[j][1], r = Math.hypot(dx, dy), ux = dx / r, uy = dy / r;
      let l3 = 1, l5 = 1;
      if (damp) { const u3 = r ** 3 / ALPHA, e = Math.exp(-A_DIP * u3); l3 = 1 - e; l5 = 1 - (1 + A_DIP * u3) * e; }
      const T = [[(3 * l5 * ux * ux - l3) / r ** 3, 3 * l5 * ux * uy / r ** 3], [3 * l5 * ux * uy / r ** 3, (3 * l5 * uy * uy - l3) / r ** 3]];
      for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) A[2 * i + p][2 * j + q] = -T[p][q];
    }
  }
  const L = Array.from({ length: N }, () => new Float64Array(N));   // Cholesky A = L L^T
  for (let i = 0; i < N; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
    if (i === j) { if (s <= 1e-9) return { bad: true }; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
  }
  const y = new Float64Array(N), x = new Float64Array(N);
  for (let i = 0; i < N; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; }
  for (let i = N - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < N; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; }
  let U = 0; for (let k = 0; k < N; k++) U -= 0.5 * x[k] * b[k];
  const mu = sites.map((_, i) => [x[2 * i], x[2 * i + 1]]);
  if (Math.max(...mu.map(m => Math.hypot(...m))) > 5) return { bad: true };
  return { mu, U: U * K, iters: 0 };
}

// ---------- Drude oscillators: core charge -qD at the site, Drude charge qD on a spring ----------
function screen(r, damp) { if (!damp) return [1, 0]; const ar = A_DRUDE * r, e = Math.exp(-ar); return [1 - (1 + ar / 2) * e, e * (A_DRUDE / 2 + A_DRUDE * ar / 2)]; }  // S, dS/dr
function solveDrude(ions, sites, damp) {
  const n = sites.length, d = sites.map(() => [0, 0]);
  const charges = () => { const c = []; sites.forEach((s, i) => { c.push([s[0], s[1], -QD, i]); c.push([s[0] + d[i][0], s[1] + d[i][1], QD, i]); }); return c; };
  let it = 0, conv = false;
  for (; it < 2000; it++) {
    const c = charges(); let change = 0;
    for (let i = 0; i < n; i++) {
      const px = sites[i][0] + d[i][0], py = sites[i][1] + d[i][1]; let Fx = 0, Fy = 0;
      for (const q of ions) { const dx = px - q[0], dy = py - q[1], r = Math.hypot(dx, dy); Fx += K * QD * dx / r ** 3; Fy += K * QD * dy / r ** 3; }
      for (const [x, y, q, j] of c) { if (j === i) continue; const dx = px - x, dy = py - y, r = Math.hypot(dx, dy), [S, dS] = screen(r, damp), f = K * QD * q * (S / r ** 2 - dS / r) / r; Fx += f * dx; Fy += f * dy; }
      const nx = d[i][0] + 0.7 * (Fx / KD - d[i][0]), ny = d[i][1] + 0.7 * (Fy / KD - d[i][1]);
      change = Math.max(change, Math.hypot(nx - d[i][0], ny - d[i][1])); d[i] = [nx, ny];
      if (Math.hypot(nx, ny) > 0.7) return { bad: true };
    }
    if (change < 1e-10) { conv = true; break; }
  }
  if (!conv) return { bad: true };
  const c = charges(); let U = 0;
  d.forEach(v => U += 0.5 * KD * (v[0] ** 2 + v[1] ** 2));
  for (const q of ions) for (const [x, y, qq] of c) U += K * qq / Math.hypot(x - q[0], y - q[1]);
  for (let a = 0; a < c.length; a++) for (let b = a + 1; b < c.length; b++) { if (c[a][3] === c[b][3]) continue; const r = Math.hypot(c[a][0] - c[b][0], c[a][1] - c[b][1]); U += K * c[a][2] * c[b][2] * screen(r, damp)[0] / r; }
  return { d, mu: d.map(v => [QD * v[0], QD * v[1]]), U, iters: it };
}

export function mountPolarizable(root) {
  const st = { setup: 'sites', model: 'dip', damp: true, s: 2.4, ions: [[3.2, 0.6]] };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Drag to move an ion</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-setup="sites" class="on">1 ion · 3 sites</button><button data-setup="ions">3 ions · 1 site</button></div>
    <div class="seg"><button data-model="dip" class="on">Induced dipoles</button><button data-model="drude">Drude</button></div>
    <div class="seg"><button data-damp="1" class="on">Thole on</button><button data-damp="0">Thole off</button></div>
    <label data-sp>Site spacing <output data-o="s"></output><input type="range" data-k="s" min="1.2" max="3.5" step="0.05" value="2.4"></label>
    <div class="readouts">
      <div><span>Many-body energy</span><b data-o="U"></b></div>
      <div><span>Sum of pair terms</span><b data-o="P"></b></div>
      <div><span>Largest induced dipole</span><b data-o="mu"></b></div>
      <div><span>Closest ion–site</span><b data-o="r"></b></div>
      <div class="wide"><span>Status</span><b data-o="st"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const seg = (attr, f) => root.querySelectorAll(`[data-${attr}]`).forEach(b => b.onclick = () => { root.querySelectorAll(`[data-${attr}]`).forEach(z => z.classList.toggle('on', z === b)); f(b.dataset[attr]); draw(); });
  seg('setup', v => { st.setup = v; st.ions = v === 'sites' ? [[3.2, 0.6]] : [0, 60, 120].map(t => [2.2 * Math.cos((t - 60) * Math.PI / 180), 2.2 * Math.sin((t - 60) * Math.PI / 180)]); $('[data-sp]').hidden = v !== 'sites'; push(); });
  seg('model', v => { st.model = v; push(); });
  seg('damp', v => st.damp = v === '1');
  $('[data-k=s]').oninput = e => { st.s = +e.target.value; push(); draw(); };

  const sites = () => st.setup === 'ions' ? [[0, 0]] : [0, 1, 2].map(k => { const t = Math.PI / 2 + k * 2 * Math.PI / 3, R = st.s / Math.sqrt(3); return [R * Math.cos(t) - 0.6, R * Math.sin(t)]; });
  function push() {   // keep every ion at least RMIN from every site (stands in for short-range repulsion)
    const S = sites(), RMIN = st.model === 'dip' ? 1.6 : 2.0;   // a Drude particle closer than ~1.9 A to Li+ falls onto it
    for (const q of st.ions) for (let rep = 0; rep < 4; rep++) for (const s of S) { const dx = q[0] - s[0], dy = q[1] - s[1], r = Math.hypot(dx, dy) || 1e-6; if (r < RMIN) { q[0] = s[0] + dx / r * RMIN; q[1] = s[1] + dy / r * RMIN; } }
  }
  const solve = (ions, S) => st.model === 'dip' ? solveDipoles(ions, S, st.damp) : solveDrude(ions, S, st.damp);

  let W = 600, H = 600, sc = 40, ox = 300, oy = 250, grab = -1;
  const toW = (x, y) => [(x - ox) / sc, (oy - y) / sc], toS = p => [ox + p[0] * sc, oy - p[1] * sc];
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const pos = e => { const r = cv.getBoundingClientRect(); return toW(e.clientX - r.left, e.clientY - r.top); };
  cv.addEventListener('pointerdown', e => {
    const p = pos(e); let best = 1e9; st.ions.forEach((q, k) => { const r = Math.hypot(q[0] - p[0], q[1] - p[1]); if (r < best) { best = r; grab = k; } });
    cv.setPointerCapture(e.pointerId); move(p); e.preventDefault();
  });
  const move = p => { if (grab < 0) return; const lim = [W / sc / 2 - 0.4, (H * 0.72) / sc / 2 - 0.4]; st.ions[grab] = [Math.max(-lim[0], Math.min(lim[0], p[0])), Math.max(-lim[1], Math.min(lim[1], p[1]))]; push(); draw(); };
  cv.addEventListener('pointermove', e => { if (grab >= 0) move(pos(e)); });
  const up = () => { grab = -1; }; cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);

  function draw() {
    if (!W) return;
    const S = sites(), res = solve(st.ions, S);
    let pair = 0, pairBad = false;
    for (const q of st.ions) for (const s of S) { const r1 = solve([q], [s]); if (r1.bad) pairBad = true; else pair += r1.U; }
    const plotH = H * 0.72; sc = Math.min(W / 10, plotH / 7.5); ox = W / 2; oy = plotH / 2 + 6;
    cx.clearRect(0, 0, W, H);
    // field of the ions (faint radial lines)
    cx.lineWidth = 1;
    for (const q of st.ions) { const [X, Y] = toS(q); for (let k = 0; k < 16; k++) { const t = k * Math.PI / 8; cx.strokeStyle = 'rgba(255,107,107,.13)'; cx.beginPath(); cx.moveTo(X + 22 * Math.cos(t), Y + 22 * Math.sin(t)); cx.lineTo(X + 4 * sc * Math.cos(t), Y + 4 * sc * Math.sin(t)); cx.stroke(); } }
    S.forEach((s, i) => {
      const [X, Y] = toS(s), m = res.bad ? [0, 0] : res.mu[i], R = 0.85 * sc;
      if (st.model === 'dip') {
        const shift = res.bad ? [0, 0] : [-m[0] * 0.9 * sc, m[1] * 0.9 * sc];    // electron cloud moves opposite to the dipole
        const g = cx.createRadialGradient(X + shift[0], Y + shift[1], 1, X + shift[0], Y + shift[1], R);
        g.addColorStop(0, res.bad ? 'rgba(255,107,107,.7)' : 'rgba(56,225,255,.55)'); g.addColorStop(0.7, 'rgba(56,225,255,.16)'); g.addColorStop(1, 'rgba(56,225,255,0)');
        cx.fillStyle = g; cx.beginPath(); cx.arc(X + shift[0], Y + shift[1], R, 0, 7); cx.fill();
        cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(X, Y, 5, 0, 7); cx.fill();
      } else {
        const dd = res.bad ? [0, 0] : res.d[i], D = [X + dd[0] * sc * 3, Y - dd[1] * sc * 3];
        cx.strokeStyle = COL.muted; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(X, Y); const L = Math.hypot(D[0] - X, D[1] - Y);
        if (L > 2) { const ux = (D[0] - X) / L, uy = (D[1] - Y) / L; for (let j = 1; j < 12; j++) { const f = j / 12, a = (j % 2 ? 5 : -5); cx.lineTo(X + ux * L * f - uy * a, Y + uy * L * f + ux * a); } }
        cx.lineTo(...D); cx.stroke();
        cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X, Y, 10, 0, 7); cx.fill();
        cx.fillStyle = res.bad ? COL.coral : COL.cyan; cx.beginPath(); cx.arc(D[0], D[1], 6, 0, 7); cx.fill();
      }
      if (!res.bad) { const mm = Math.hypot(...m); if (mm > 0.01) { const Lp = Math.min(2.2 * sc, mm * 2.2 * sc), ux = m[0] / mm, uy = -m[1] / mm; arrow(X, Y, ux * Lp, uy * Lp); } }
    });
    st.ions.forEach(q => { const [X, Y] = toS(q); cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X, Y, 14, 0, 7); cx.fill(); cx.fillStyle = '#060C1A'; cx.font = '700 16px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText('+', X, Y + 5.5); });
    // 1 A scale bar and legend
    cx.strokeStyle = COL.muted; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(12, plotH - 8); cx.lineTo(12 + sc, plotH - 8); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('1 Å', 12, plotH - 14);
    cx.fillStyle = COL.amber; cx.fillText(st.model === 'dip' ? 'arrows: induced dipoles' : 'arrows: dipoles · Drude shift drawn ×3', 12, 20);
    if (res.bad) { cx.fillStyle = 'rgba(6,12,26,.85)'; cx.fillRect(W / 2 - 160, 30, 320, 30); cx.fillStyle = COL.coral; cx.font = '700 15px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText('Polarization catastrophe!', W / 2, 51); }
    // energy bars
    const top = plotH + 18, bh = 20, x0 = 130, x1 = W - 70, vals = [res.bad ? null : res.U, pairBad ? null : pair], vmax = Math.max(5, ...vals.filter(v => v != null).map(Math.abs));
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.fillText('Induction energy (kcal/mol)', 12, top);
    [['many-body', COL.coral], ['sum of pairs', COL.cyan]].forEach(([nm, c], k) => {
      const y = top + 14 + k * (bh + 10), v = vals[k]; cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText(nm, x0 - 8, y + 14);
      if (v == null) { cx.fillStyle = COL.coral; cx.textAlign = 'left'; cx.fillText('diverges', x0, y + 14); return; }
      const w = Math.abs(v) / vmax * (x1 - x0); cx.fillStyle = c; cx.fillRect(x0, y, w, bh); cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.fillText(v.toFixed(2), x0 + w + 6, y + 14);
    });
    // readouts
    const rmin = Math.min(...st.ions.flatMap(q => S.map(s => Math.hypot(q[0] - s[0], q[1] - s[1]))));
    o.s.textContent = st.s.toFixed(2) + ' Å';
    o.U.textContent = res.bad ? 'diverges' : res.U.toFixed(2) + ' kcal/mol';
    o.P.textContent = pairBad ? 'diverges' : pair.toFixed(2) + ' kcal/mol';
    o.mu.textContent = res.bad ? '∞' : (Math.max(...res.mu.map(m => Math.hypot(...m))) * DEBYE).toFixed(2) + ' D';
    o.r.textContent = rmin.toFixed(2) + ' Å';
    const ratio = !res.bad && Math.abs(pair) > 0.05 ? ` · many-body / pairs = ${(res.U / pair).toFixed(2)}` : '';
    o.st.textContent = res.bad ? (st.damp ? 'no stable solution' : 'dipoles run away: switch Thole damping on, or spread the sites') : (st.model === 'dip' ? 'solved exactly' : `relaxed in ${res.iters} iterations (ions kept ≥ 2 Å away)`) + ratio;
  }
  function arrow(x, y, dx, dy) {
    const L = Math.hypot(dx, dy); if (L < 4) return; const ux = dx / L, uy = dy / L;
    cx.strokeStyle = cx.fillStyle = COL.amber; cx.lineWidth = 3; cx.beginPath(); cx.moveTo(x, y); cx.lineTo(x + dx - ux * 8, y + dy - uy * 8); cx.stroke();
    cx.beginPath(); cx.moveTo(x + dx, y + dy); cx.lineTo(x + dx - ux * 11 - uy * 6, y + dy - uy * 11 + ux * 6); cx.lineTo(x + dx - ux * 11 + uy * 6, y + dy - uy * 11 - ux * 6); cx.closePath(); cx.fill();
  }
  push();
  return () => ro.disconnect();
}
