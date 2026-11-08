// Lab 5.5 widget, three modes.
// Ewald: Na+/Cl- in a periodic box (L = 1 nm); alpha and k_max sliders; real/reciprocal/self bars and |error| vs alpha.
// Sums: 1D chain partial sums and 3D NaCl sphere vs cube (Evjen) sums.
// Dipoles: ion + two polarizable atoms on a line; iterate the induced dipoles; catastrophe below (2 alpha')^(1/3).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const F = 138.935, L = 1.0;
function erfc(x) { // Numerical Recipes erfc, relative error < 1.2e-7
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
function ewald(d, alpha, kmax, nreal) {
  const pos = [[0, 0, 0], [d, 0, 0]], q = [1, -1]; let Er = 0;
  for (let a = -nreal; a <= nreal; a++) for (let b = -nreal; b <= nreal; b++) for (let c = -nreal; c <= nreal; c++)
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      if (!a && !b && !c && i === j) continue;
      let dx = pos[j][0] - pos[i][0], dy = pos[j][1] - pos[i][1], dz = pos[j][2] - pos[i][2];
      if (nreal === 0) { dx -= L * Math.round(dx / L); }   // minimum image
      else { dx += a * L; dy += b * L; dz += c * L; }
      const r = Math.hypot(dx, dy, dz); Er += 0.5 * F * q[i] * q[j] * erfc(alpha * r) / r;
    }
  let Ek = 0; const tp = 2 * Math.PI / L;
  for (let a = -kmax; a <= kmax; a++) { const S2 = 2 - 2 * Math.cos(tp * a * d);   // |S|^2 depends only on kx
    if (S2 === 0 && a === 0) { continue; }
    for (let b = -kmax; b <= kmax; b++) for (let c = -kmax; c <= kmax; c++) { if (!a && !b && !c) continue;
      const k2 = tp * tp * (a * a + b * b + c * c); Ek += F * 2 * Math.PI / (L * L * L) * Math.exp(-k2 / (4 * alpha * alpha)) / k2 * S2; } }
  const Es = -F * alpha / Math.sqrt(Math.PI) * 2;
  return [Er, Ek, Es, Er + Ek + Es];
}

export function mountLab5(root) {
  const st = { mode: 'ewald', d: 0.3, alpha: 5, kmax: 4, dd: 0.25, ap: 1.5, it: 0 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view el-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="ewald" class="on">Ewald</button><button data-m="sums">Sums</button><button data-m="dip">Dipoles</button></div>
    <div data-p="ewald">
      <label><span>Ion separation (nm)</span><output data-o="d"></output><input type="range" data-k="d" min="0.15" max="0.5" step="0.01" value="0.3"></label>
      <label><span>Splitting α (1/nm)</span><output data-o="al"></output><input type="range" data-k="al" min="1" max="12" step="0.1" value="5"></label>
      <label><span>k-vectors: |nᵢ| ≤ k_max</span><output data-o="km"></output><input type="range" data-k="km" min="1" max="8" step="1" value="4"></label>
    </div>
    <div data-p="dip" hidden>
      <label><span>Distance between the two atoms (nm)</span><output data-o="dd"></output><input type="range" data-k="dd" min="0.12" max="0.40" step="0.005" value="0.25"></label>
      <label><span>Polarizability α' (Å³)</span><output data-o="ap"></output><input type="range" data-k="ap" min="0.5" max="3" step="0.1" value="1.5"></label>
      <div class="row"><button class="btn sm primary" data-k="iter">Iterate once</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
      <div class="tbl" style="font-size:12px"><table data-o="tab"></table></div>
    </div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.el-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let exact = 0, parts = [0, 0, 0, 0], curve = [], mu = [0, 0], rows = [];
  function recompute() {
    exact = ewald(st.d, 5, 10, 3)[3]; parts = ewald(st.d, st.alpha, st.kmax, 0);
    curve = []; for (let a = 1; a <= 12.001; a += 0.25) curve.push([a, Math.abs(ewald(st.d, a, st.kmax, 0)[3] - exact)]);
  }
  function resetDip() { mu = [0, 0]; rows = []; st.it = 0; }
  function iterate() {   // ion +1 at 0, A at 0.25, B at 0.25 + dd, all on x; mu in e nm
    const ap = st.ap * 1e-3, rA = 0.25, rB = 0.25 + st.dd, c = 2 / st.dd ** 3;
    const nA = ap * (1 / rA ** 2 + c * mu[1]), nB = ap * (1 / rB ** 2 + c * mu[0]); mu = [nA, nB]; st.it++;
    rows.push([st.it, nA, nB]); if (rows.length > 6) rows.shift();
  }
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); draw(); });
  $('[data-k=d]').oninput = e => { st.d = +e.target.value; recompute(); draw(); };
  $('[data-k=al]').oninput = e => { st.alpha = +e.target.value; parts = ewald(st.d, st.alpha, st.kmax, 0); draw(); };
  $('[data-k=km]').oninput = e => { st.kmax = +e.target.value; recompute(); draw(); };
  $('[data-k=dd]').oninput = e => { st.dd = +e.target.value; resetDip(); draw(); };
  $('[data-k=ap]').oninput = e => { st.ap = +e.target.value; resetDip(); draw(); };
  $('[data-k=iter]').onclick = () => { iterate(); draw(); };
  $('[data-k=reset]').onclick = () => { resetDip(); draw(); };

  // Madelung data (computed once)
  const chain = []; { let s = 0; for (let n = 1; n <= 40; n++) { s += (n % 2 ? 1 : -1) / n; chain.push(s); } }
  function mad(shape, R) {
    const n = Math.floor(R) + 1; let s = 0;
    for (let x = -n; x <= n; x++) for (let y = -n; y <= n; y++) for (let z = -n; z <= n; z++) {
      const r = Math.hypot(x, y, z); if (!r) continue; let w = 1;
      if (shape === 'sphere') { if (r > R + 1e-9) continue; } else { if (Math.abs(x) > R || Math.abs(y) > R || Math.abs(z) > R) continue; for (const A of [x, y, z]) if (Math.abs(A) === R) w *= 0.5; }
      s -= ((x + y + z) % 2 ? -1 : 1) * w / r;
    }
    return s;
  }
  let sph = null, cub = null;
  const ensureMad = () => { if (sph) return; sph = []; cub = []; for (let R = 1; R <= 10; R++) { sph.push(mad('sphere', R + 0.5)); cub.push(mad('cube', R)); } };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x, y, c = COL.ink, al = 'left', fnt = '600 12px') => { cx.fillStyle = c; cx.font = `${fnt} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x, y); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const ro_ = $('[data-ro]'), half = H * 0.42;
    if (st.mode === 'ewald') {
      // box with images (2D slice z = 0)
      const bw = Math.min(W * 0.42, half - 20), S = bw / 3;
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
        const ox = 10 + (a + 1) * S, oy = 10 + (1 - b) * S, main = !a && !b;
        cx.strokeStyle = main ? COL.cyan : 'rgba(143,163,199,.25)'; cx.lineWidth = main ? 2 : 1; cx.strokeRect(ox, oy, S, S);
        const px = (x, y) => [ox + (0.25 + x) / L * S * 0.999, oy + S * 0.5 - y * S];
        cx.globalAlpha = main ? 1 : 0.35;
        let p = px(0, 0); cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(p[0], p[1], 6, 0, 7); cx.fill();
        p = px(st.d, 0); cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(p[0], p[1], 7, 0, 7); cx.fill(); cx.globalAlpha = 1;
      }
      label('Na⁺ (amber), Cl⁻ (cyan) and images', 12, bw + 26, COL.muted, 'left', '11px');
      // error vs alpha (log)
      const tx0 = bw + 56, tx1 = W - 12, ty0 = 22, ty1 = half - 10, lo = -4, hi = 3, AX = a => tx0 + (a - 1) / 11 * (tx1 - tx0), EY = e => ty1 - (Math.min(hi, Math.max(lo, Math.log10(e + 1e-9))) - lo) / (hi - lo) * (ty1 - ty0);
      axes(tx0, ty0, tx1, ty1); label('|total − exact| vs α (log scale)', tx0 + 6, ty0 - 6);
      cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); curve.forEach(([a, e], k) => k ? cx.lineTo(AX(a), EY(e)) : cx.moveTo(AX(a), EY(e))); cx.stroke();
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(AX(st.alpha), EY(Math.abs(parts[3] - exact)), 5, 0, 7); cx.fill();
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [[-3, '0.001'], [0, '1'], [2, '100']].forEach(([e, t]) => cx.fillText(t, tx0 - 4, EY(10 ** e) + 4));
      cx.textAlign = 'center'; [2, 4, 6, 8, 10, 12].forEach(a => cx.fillText(a, AX(a), ty1 + 13));
      // bars
      const bx0 = 46, bx1 = W - 12, by0 = half + 40, by1 = H - 26, vals = [['real', parts[0]], ['reciprocal', parts[1]], ['self', parts[2]], ['total', parts[3]]];
      const vmax = Math.max(800, ...vals.map(v => Math.abs(v[1]))) * 1.3, Y0 = (by0 + by1) / 2, BY = v => Y0 - v / vmax * (by1 - by0) / 2, bwid = (bx1 - bx0) / 4;
      axes(bx0, by0, bx1, by1); cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.beginPath(); cx.moveTo(bx0, Y0); cx.lineTo(bx1, Y0); cx.stroke();
      label(`Ewald parts (kJ/mol), α = ${st.alpha.toFixed(1)}`, bx0 + 6, by0 - 10);
      vals.forEach(([n, v], k) => { const x = bx0 + k * bwid + bwid * 0.2; cx.fillStyle = k === 3 ? COL.coral : 'rgba(56,225,255,.55)'; cx.fillRect(x, Math.min(Y0, BY(v)), bwid * 0.6, Math.abs(BY(v) - Y0));
        label(n, x + bwid * 0.3, by1 + 14, COL.muted, 'center', '11px'); label(v.toFixed(1), x + bwid * 0.3, v < 0 ? BY(v) + 14 : BY(v) - 4, COL.ink, 'center', '11px'); });
      cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(bx0 + 3 * bwid, BY(exact)); cx.lineTo(bx1, BY(exact)); cx.stroke(); cx.setLineDash([]);
      o.d.textContent = st.d.toFixed(2); o.al.textContent = st.alpha.toFixed(1); o.km.textContent = `${st.kmax} (${(2 * st.kmax + 1) ** 3 - 1} vectors)`;
      ro_.innerHTML = `<div><span>Ewald total</span><b>${parts[3].toFixed(3)}</b></div><div><span>Exact periodic</span><b>${exact.toFixed(3)}</b></div>
        <div><span>Bare / minimum image</span><b>${(-F / st.d).toFixed(3)}</b></div><div><span>Ewald error</span><b>${(parts[3] - exact).toFixed(3)}</b></div>`;
    } else if (st.mode === 'sums') {
      ensureMad();
      const x0 = 46, x1 = W - 12, y0 = 24, y1 = half - 10, X = n => x0 + (n - 1) / 39 * (x1 - x0), Y = s => y1 - (s - 0.4) / 0.7 * (y1 - y0);
      axes(x0, y0, x1, y1); label('1D chain: sum after n neighbours (exact ln 2)', x0 + 6, y0 - 8);
      cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(x0, Y(Math.LN2)); cx.lineTo(x1, Y(Math.LN2)); cx.stroke(); cx.setLineDash([]);
      cx.strokeStyle = COL.cyan; cx.lineWidth = 1.5; cx.beginPath(); chain.forEach((s, k) => k ? cx.lineTo(X(k + 1), Y(s)) : cx.moveTo(X(k + 1), Y(s))); cx.stroke();
      chain.forEach((s, k) => { cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(X(k + 1), Y(s), 2.5, 0, 7); cx.fill(); });
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 10, 20, 30, 40].forEach(n => cx.fillText(n, X(n), y1 + 13));
      const bx0 = 46, bx1 = W - 12, by0 = half + 40, by1 = H - 26, RX = R => bx0 + (R - 1) / 9 * (bx1 - bx0), MY = m => by1 - (Math.max(-2.5, Math.min(3.5, m)) + 2.5) / 6 * (by1 - by0);
      axes(bx0, by0, bx1, by1); label('NaCl: sphere (coral) vs cube (cyan) around one ion', bx0 + 6, by0 - 10);
      cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(bx0, MY(1.747565)); cx.lineTo(bx1, MY(1.747565)); cx.stroke(); cx.setLineDash([]);
      [[sph, COL.coral], [cub, COL.cyan]].forEach(([arr, c]) => { cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath(); arr.forEach((m, k) => k ? cx.lineTo(RX(k + 1), MY(m)) : cx.moveTo(RX(k + 1), MY(m))); cx.stroke();
        arr.forEach((m, k) => { cx.fillStyle = c; cx.beginPath(); cx.arc(RX(k + 1), MY(m), 3.5, 0, 7); cx.fill(); }); });
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 4, 7, 10].forEach(R => cx.fillText(R, RX(R), by1 + 14));
      cx.textAlign = 'right'; [-2, 0, 1.75, 3].forEach(m => cx.fillText(m, bx0 - 4, MY(m) + 4));
      ro_.innerHTML = `<div><span>1D after 40 neighbours</span><b>${chain[39].toFixed(4)}</b></div><div><span>1D exact (ln 2)</span><b>0.6931</b></div>
        <div><span>Sphere, radius 10</span><b>${sph[9].toFixed(4)}</b></div><div><span>Cube, half-side 10</span><b>${cub[9].toFixed(4)}</b></div>`;
    } else {
      const ap = st.ap * 1e-3, factor = 2 * ap / st.dd ** 3, dc = Math.cbrt(2 * ap), y = half * 0.55, x0 = 30, sc = (W - 60) / 0.7;
      const X = x => x0 + x * sc, rA = 0.25, rB = 0.25 + st.dd;
      cx.strokeStyle = 'rgba(143,163,199,.3)'; cx.beginPath(); cx.moveTo(x0, y); cx.lineTo(W - 20, y); cx.stroke();
      cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(X(0), y, 12, 0, 7); cx.fill(); label('+', X(0), y + 4, '#0B1430', 'center');
      [[rA, mu[0], 'A'], [rB, mu[1], 'B']].forEach(([r, m, n]) => { cx.fillStyle = 'rgba(56,225,255,.3)'; cx.beginPath(); cx.arc(X(r), y, 14, 0, 7); cx.fill();
        const len = Math.min(200, m * 2500); cx.strokeStyle = COL.coral; cx.lineWidth = 3; cx.beginPath(); cx.moveTo(X(r) - len / 2, y); cx.lineTo(X(r) + len / 2, y); cx.stroke();
        cx.fillStyle = COL.coral; cx.beginPath(); cx.moveTo(X(r) + len / 2 + 7, y); cx.lineTo(X(r) + len / 2 - 2, y - 6); cx.lineTo(X(r) + len / 2 - 2, y + 6); cx.fill();
        label(n, X(r), y - 24, COL.ink, 'center'); });
      label(`ion at 0, A at 0.25 nm, B at ${rB.toFixed(3)} nm; arrows: induced dipoles`, 12, 24, COL.muted, 'left', '11px');
      label(factor < 1 ? `each round multiplies by 2α'/d³ = ${factor.toFixed(3)}: converges` : `2α'/d³ = ${factor.toFixed(3)} ≥ 1: the dipoles run away`, 12, half + 30, factor < 1 ? COL.ink : COL.coral);
      label(`catastrophe distance (2α')^(1/3) = ${dc.toFixed(3)} nm`, 12, half + 52, COL.muted, 'left', '12px');
      // convergence plot of mu_B over iterations (computed ahead)
      const bx0 = 46, bx1 = W - 12, by0 = half + 80, by1 = H - 26; let m = [0, 0]; const seq = [];
      for (let k = 0; k < 20; k++) { const c = 2 / st.dd ** 3; m = [ap * (1 / rA ** 2 + c * m[1]), ap * (1 / rB ** 2 + c * m[0])]; seq.push(m[1]); }
      const vmax = Math.max(...seq.map(Math.abs)) * 1.1 || 0.01, IX = k => bx0 + k / 19 * (bx1 - bx0), IY = v => by1 - Math.min(v, vmax) / vmax * (by1 - by0);
      axes(bx0, by0, bx1, by1); label('μ_B over 20 iterations', bx0 + 6, by0 - 8, COL.ink, 'left', '600 11px');
      cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); seq.forEach((v, k) => k ? cx.lineTo(IX(k), IY(v)) : cx.moveTo(IX(k), IY(v))); cx.stroke();
      o.dd.textContent = st.dd.toFixed(3); o.ap.textContent = st.ap.toFixed(1);
      o.tab.innerHTML = `<tr><th>iteration</th><th>μ_A (e nm)</th><th>μ_B (e nm)</th></tr>` + rows.map(r => `<tr><td>${r[0]}</td><td>${r[1].toFixed(5)}</td><td>${r[2].toFixed(5)}</td></tr>`).join('');
      const UA = -0.5 * F * ap / rA ** 4;
      ro_.innerHTML = `<div><span>Dipole of A alone (e nm)</span><b>${(ap / rA ** 2).toFixed(5)}</b></div><div><span>Dipole of B alone (e nm)</span><b>${(ap / rB ** 2).toFixed(5)}</b></div>
        <div><span>Energy of A alone (kJ/mol)</span><b>${UA.toFixed(2)}</b></div><div><span>Amplification per round</span><b>${factor.toFixed(3)}</b></div>`;
    }
  }
  recompute(); draw();
  return () => { ro.disconnect(); };
}
