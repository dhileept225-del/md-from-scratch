// Coulomb-sum lab: the Madelung constant of rock salt (NaCl, nearest-neighbour distance 1) computed with a plain
// spherical cutoff, neutral expanding cubes (Evjen), or the damped-shifted Wolf sum. 3D sums; the picture shows z = 0.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const M_EXACT = 1.747565, NMAX = 13;
function erfc(z) { const t = 1 / (1 + 0.3275911 * z); return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-z * z); }

// lattice around a +1 reference ion at the origin
const P = [];
for (let i = -NMAX; i <= NMAX; i++) for (let j = -NMAX; j <= NMAX; j++) for (let k = -NMAX; k <= NMAX; k++) {
  if (!i && !j && !k) continue; const r = Math.hypot(i, j, k);
  P.push({ i, j, k, r, q: (Math.abs(i + j + k) % 2) ? -1 : 1, m: Math.max(Math.abs(i), Math.abs(j), Math.abs(k)), nb: [i, j, k].filter(v => Math.abs(v) === Math.max(Math.abs(i), Math.abs(j), Math.abs(k))).length });
}
function madelung(method, R, alpha) {
  let s = 0, net = 1, n = 0;
  if (method === 'sphere') { for (const p of P) if (p.r < R) { s += p.q / p.r; net += p.q; n++; } return { M: -s, net, n }; }
  if (method === 'cube') {
    const L = Math.max(1, Math.round(R));
    for (const p of P) { if (p.m > L) continue; const w = p.m === L ? 0.5 ** p.nb : 1; s += w * p.q / p.r; net += w * p.q; n++; }
    return { M: -s, net, n };
  }
  const ec = erfc(alpha * R) / R;                     // Wolf
  for (const p of P) if (p.r < R) { s += p.q * (erfc(alpha * p.r) / p.r - ec); net += p.q; n++; }
  const e = 0.5 * s - (ec / 2 + alpha / Math.sqrt(Math.PI));
  return { M: -2 * e, net, n };
}

export function mountCoulomb(root) {
  const st = { method: 'sphere', R: 4.5, alpha: 0.3 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="sphere" class="on">Sphere cutoff</button><button data-m="cube">Neutral cubes</button><button data-m="wolf">Wolf sum</button></div>
    <label><span>Size R (in nearest-neighbour distances)</span><output data-o="R"></output><input type="range" data-k="R" min="1.5" max="12" step="0.1" value="4.5"></label>
    <label data-a hidden><span>Damping <i class="nc">α</i></span><output data-o="a"></output><input type="range" data-k="a" min="0" max="0.8" step="0.02" value="0.3"></label>
    <div class="readouts">
      <div><span>Madelung estimate</span><b data-o="M"></b></div>
      <div><span>Error (exact 1.747565)</span><b data-o="err"></b></div>
      <div><span>Net charge included</span><b data-o="net"></b></div>
      <div><span>Ions included</span><b data-o="n"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let curve = [];
  const recompute = () => { curve = []; for (let R = 1.5; R <= 12.001; R += 0.1) curve.push([R, madelung(st.method, R, st.alpha).M]); draw(); };
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.method = b.dataset.m; root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); $('[data-a]').hidden = st.method !== 'wolf'; recompute(); });
  $('[data-k=R]').oninput = e => { st.R = +e.target.value; draw(); };
  $('[data-k=a]').oninput = e => { st.alpha = +e.target.value; recompute(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function draw() {
    if (!W || !curve.length) return;
    cx.clearRect(0, 0, W, H);
    const top = H * 0.5, sp = Math.min(W, top) / 26, c0 = [W / 2, top / 2 + 4], L = Math.max(1, Math.round(st.R));
    for (let i = -12; i <= 12; i++) for (let j = -12; j <= 12; j++) {
      const x = c0[0] + i * sp, y = c0[1] + j * sp; if (Math.abs(i * sp) > W / 2 - 4 || Math.abs(j * sp) > top / 2 - 2) continue;
      const inside = st.method === 'cube' ? Math.max(Math.abs(i), Math.abs(j)) <= L : Math.hypot(i, j) < st.R;
      const edge = st.method === 'cube' && Math.max(Math.abs(i), Math.abs(j)) === L;
      cx.globalAlpha = inside ? (edge ? 0.55 : 1) : 0.15; cx.fillStyle = (Math.abs(i + j) % 2) ? COL.cyan : COL.coral;
      cx.beginPath(); cx.arc(x, y, Math.max(1.5, sp * 0.3), 0, 7); cx.fill();
    }
    cx.globalAlpha = 1; cx.fillStyle = COL.ink; cx.beginPath(); cx.arc(c0[0], c0[1], Math.max(3, sp * 0.4), 0, 7); cx.fill();
    cx.strokeStyle = COL.amber; cx.lineWidth = 2;
    if (st.method === 'cube') cx.strokeRect(c0[0] - L * sp, c0[1] - L * sp, 2 * L * sp, 2 * L * sp); else { cx.beginPath(); cx.arc(c0[0], c0[1], st.R * sp, 0, 7); cx.stroke(); }
    cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillStyle = COL.coral; cx.fillText('+ ions', 10, 18); cx.fillStyle = COL.cyan; cx.fillText('− ions', 60, 18); cx.fillStyle = COL.muted; cx.fillText('slice z = 0 of the 3D sum', 10, 34);
    // curve
    const x0 = 46, x1 = W - 12, y0 = top + 30, y1 = H - 24, lo = st.method === 'sphere' ? -10 : 1.3, hi = st.method === 'sphere' ? 10 : 2.1;
    const X = R => x0 + (R - 1.5) / 10.5 * (x1 - x0), Y = M => y1 - (Math.max(lo, Math.min(hi, M)) - lo) / (hi - lo) * (y1 - y0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [2, 4, 6, 8, 10, 12].forEach(R => cx.fillText(R, X(R), y1 + 14));
    cx.textAlign = 'right'; cx.fillText(hi, x0 - 5, y0 + 8); cx.fillText(lo, x0 - 5, y1);
    cx.setLineDash([6, 5]); cx.strokeStyle = COL.ink; cx.beginPath(); cx.moveTo(x0, Y(M_EXACT)); cx.lineTo(x1, Y(M_EXACT)); cx.stroke(); cx.setLineDash([]);
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.cyan; cx.fillText('Madelung estimate vs R', x0 + 6, y0 - 8); cx.fillStyle = COL.ink; cx.fillText('exact', x1 - 40, Y(M_EXACT) - 6);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); curve.forEach(([R, M], k) => (k ? cx.lineTo(X(R), Y(M)) : cx.moveTo(X(R), Y(M)))); cx.stroke();
    const res = madelung(st.method, st.R, st.alpha);
    cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X(st.R), Y(res.M), 6, 0, 7); cx.fill();
    o.R.textContent = st.method === 'cube' ? `${st.R.toFixed(1)} → cube half-width ${L}` : st.R.toFixed(1); o.a.textContent = st.alpha.toFixed(2);
    const err = res.M - M_EXACT;
    o.M.textContent = res.M.toFixed(5); o.err.innerHTML = `<span style="color:${Math.abs(err) < 0.01 ? COL.cyan : COL.coral}">${err >= 0 ? '+' : ''}${err.toFixed(5)}</span>`;
    o.net.textContent = (Math.abs(res.net) < 1e-9 ? 0 : +res.net.toFixed(3)) + ' e'; o.n.textContent = res.n.toLocaleString('en-US');
  }
  recompute();
  return () => ro.disconnect();
}
