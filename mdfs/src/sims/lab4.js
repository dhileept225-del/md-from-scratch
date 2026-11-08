// Lab 4.6 widget: united-atom butane (TraPPE), GROMACS units (nm, ps, g/mol, kJ/mol).
// Langevin mode: dihedral populations of 40 molecules vs the exact Boltzmann curve, with a choice of 1-4 rule.
// Time-step mode: NVE energy fluctuation with flexible bonds or rigid bonds (SHAKE/RATTLE).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const kB = 0.0083145, M = [15.035, 14.027, 14.027, 15.035], B0 = 0.154, KB = 224262.4, TH0 = 114 * Math.PI / 180, KTH = 62500 * kB;
const C1 = 355.03 * kB, C2 = -68.19 * kB, C3 = 791.32 * kB, E3 = 98 * kB, S3 = 0.375;
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], nrm = a => Math.sqrt(dot(a, a));
const sc = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
const Udih = p => C1 * (1 + Math.cos(p)) + C2 * (1 - Math.cos(2 * p)) + C3 * (1 + Math.cos(3 * p));
const LJ = r => 4 * E3 * ((S3 / r) ** 12 - (S3 / r) ** 6);
function place(phiDeg) {
  const p = phiDeg * Math.PI / 180, b = B0, th = TH0, r1 = [0, 0, 0], r2 = [b, 0, 0], r3 = [b - b * Math.cos(th), b * Math.sin(th), 0];
  const bc = sc(sub(r3, r2), 1 / b), nv0 = cross(sub(r2, r1), bc), nv = sc(nv0, 1 / nrm(nv0)), mv = cross(nv, bc);
  const d = [-b * Math.cos(th), b * Math.sin(th) * Math.cos(p), b * Math.sin(th) * Math.sin(p)];
  return [r1, r2, r3, [0, 1, 2].map(k => r3[k] + d[0] * bc[k] + d[1] * mv[k] + d[2] * nv[k])];
}
const R14 = phiDeg => { const x = place(phiDeg); return nrm(sub(x[3], x[0])); };

function forces(x, s14, bonds) {   // returns [U, f, phi] for one molecule
  const f = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]; let U = 0;
  const add = (i, g, k = 1) => { f[i][0] += k * g[0]; f[i][1] += k * g[1]; f[i][2] += k * g[2]; };
  if (bonds) for (let i = 0; i < 3; i++) { const d = sub(x[i + 1], x[i]), r = nrm(d); U += 0.5 * KB * (r - B0) ** 2; const g = sc(d, -KB * (r - B0) / r); add(i + 1, g); add(i, g, -1); }
  for (const j of [1, 2]) {
    const u = sub(x[j - 1], x[j]), v = sub(x[j + 1], x[j]), lu = nrm(u), lv = nrm(v), c = Math.max(-1, Math.min(1, dot(u, v) / lu / lv)), th = Math.acos(c), s = Math.sqrt(1 - c * c) || 1e-9, dU = KTH * (th - TH0);
    U += 0.5 * KTH * (th - TH0) ** 2;
    const fi = [0, 1, 2].map(k => dU / s * (v[k] / (lu * lv) - c / (lu * lu) * u[k])), fk = [0, 1, 2].map(k => dU / s * (u[k] / (lu * lv) - c / (lv * lv) * v[k]));
    add(j - 1, fi); add(j + 1, fk); add(j, fi, -1); add(j, fk, -1);
  }
  const b1 = sub(x[1], x[0]), b2 = sub(x[2], x[1]), b3 = sub(x[3], x[2]), m = cross(b1, b2), n = cross(b2, b3), lb2 = nrm(b2);
  const phi = Math.atan2(lb2 * dot(b1, n), dot(m, n)); U += Udih(phi);
  const dU = -C1 * Math.sin(phi) + 2 * C2 * Math.sin(2 * phi) - 3 * C3 * Math.sin(3 * phi);
  const f1 = sc(m, dU * lb2 / dot(m, m)), f4 = sc(n, -dU * lb2 / dot(n, n)), p = dot(b1, b2) / (lb2 * lb2), q = dot(b3, b2) / (lb2 * lb2);
  add(0, f1); add(3, f4); add(1, f1, -1 - p); add(1, f4, q); add(2, f4, -1 - q); add(2, f1, p);
  if (s14) { const d = sub(x[3], x[0]), r2 = dot(d, d), sr6 = (S3 * S3 / r2) ** 3; U += s14 * 4 * E3 * (sr6 * sr6 - sr6); const g = sc(d, s14 * 24 * E3 * (2 * sr6 * sr6 - sr6) / r2); add(3, g); add(0, g, -1); }
  return [U, f, phi];
}
function shake(x, xo, v, dt) {
  for (let it = 0; it < 500; it++) { let err = 0;
    for (let i = 0; i < 3; i++) { const d = sub(x[i + 1], x[i]), d0 = sub(xo[i + 1], xo[i]), diff = dot(d, d) - B0 * B0, g = diff / (2 * (1 / M[i] + 1 / M[i + 1]) * dot(d, d0));
      for (let k = 0; k < 3; k++) { x[i][k] += g / M[i] * d0[k]; x[i + 1][k] -= g / M[i + 1] * d0[k]; v[i][k] += g / M[i] / dt * d0[k]; v[i + 1][k] -= g / M[i + 1] / dt * d0[k]; }
      err = Math.max(err, Math.abs(diff) / (2 * B0 * B0)); }
    if (err < 1e-10) return; }
}
function rattleV(x, v) {
  for (let it = 0; it < 500; it++) { let err = 0;
    for (let i = 0; i < 3; i++) { const d = sub(x[i + 1], x[i]), dv = sub(v[i + 1], v[i]), k = dot(d, dv) / (B0 * B0 * (1 / M[i] + 1 / M[i + 1]));
      for (let c = 0; c < 3; c++) { v[i][c] += k / M[i] * d[c]; v[i + 1][c] -= k / M[i + 1] * d[c]; } err = Math.max(err, Math.abs(k)); }
    if (err < 1e-12) return; }
}

export function mountLab4(root) {
  const st = { mode: 'lang', T: 300, s14: 0, rigid: false, dt: 2, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view bu-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="lang" class="on">Langevin: populations</button><button data-m="nve">NVE: time step</button></div>
    <div data-lang>
      <label><span>Temperature (K)</span><output data-o="T"></output><input type="range" data-k="T" min="150" max="600" step="10" value="300"></label>
      <div class="seg"><button data-s="0" class="on">1–4 excluded</button><button data-s="0.5">1–4 × 0.5</button><button data-s="1">1–4 full</button></div>
    </div>
    <div data-nve hidden>
      <div class="seg"><button data-r="0" class="on">Flexible bonds</button><button data-r="1">Rigid (SHAKE)</button></div>
      <div class="seg" style="flex-wrap:wrap">${[1, 2, 4, 6, 8, 10, 12].map(d => `<button data-d="${d}" class="${d === 2 ? 'on' : ''}">${d} fs</button>`).join('')}</div>
    </div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts">
      <div><span data-o="l1"></span><b data-o="v1"></b></div>
      <div><span data-o="l2"></span><b data-o="v2"></b></div>
      <div><span data-o="l3"></span><b data-o="v3"></b></div>
      <div><span>Time (ps)</span><b data-o="t"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.bu-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 5; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let mol, t, hist, nh, jumps, state, E0, Es, Etr, blown, exact, exactTrans;
  const NB = () => st.mode === 'lang' ? 40 : 20;
  function computeExact() {
    const kT = kB * st.T; exact = new Float64Array(72); let tot = 0, tr = 0;
    for (let k = 0; k < 360; k++) { const ph = -179.5 + k, U = Udih(ph * Math.PI / 180) + (st.s14 ? st.s14 * LJ(R14(ph)) : 0), w = Math.exp(-U / kT); exact[Math.floor(k / 5)] += w; tot += w; if (Math.abs(ph) > 120) tr += w; }
    for (let k = 0; k < 72; k++) exact[k] /= tot * 5; exactTrans = tr / tot;
  }
  function reset() {
    const kT = kB * (st.mode === 'lang' ? st.T : 300); mol = [];
    for (let b = 0; b < NB(); b++) {
      const x = place(st.mode === 'lang' && b % 2 ? 60 : 180), v = M.map(m => [0, 1, 2].map(() => gauss() * Math.sqrt(kT / m)));
      const P = [0, 1, 2].map(k => v.reduce((s, vi, i) => s + M[i] * vi[k], 0) / 58.124); v.forEach(vi => { for (let k = 0; k < 3; k++) vi[k] -= P[k]; });
      const rig = st.mode === 'nve' && st.rigid; if (rig) rattleV(x, v);
      const [U, f, phi] = forces(x, st.mode === 'lang' ? st.s14 : 0, !rig); mol.push({ x, v, f, U, phi });
    }
    t = 0; hist = new Float64Array(72); nh = 0; jumps = 0; state = mol.map(() => 0); E0 = mol.map(energy); Es = mol.map(() => [0, 0, 0]); Etr = []; blown = false; computeExact();
  }
  const energy = m => m.U + 0.5 * m.v.reduce((s, vi, i) => s + M[i] * dot(vi, vi), 0);
  function stepAll() {
    const lang = st.mode === 'lang', dt = lang ? 0.002 : st.dt / 1000, rig = !lang && st.rigid, s14 = lang ? st.s14 : 0, kT = kB * st.T, c = Math.exp(-5 * dt);
    mol.forEach((m, b) => {
      m.v.forEach((vi, i) => { for (let k = 0; k < 3; k++) vi[k] += 0.5 * dt * m.f[i][k] / M[i]; });
      const xo = m.x.map(p => [...p]); m.x.forEach((p, i) => { for (let k = 0; k < 3; k++) p[k] += dt * m.v[i][k]; });
      if (rig) shake(m.x, xo, m.v, dt);
      [m.U, m.f, m.phi] = forces(m.x, s14, !rig);
      m.v.forEach((vi, i) => { for (let k = 0; k < 3; k++) vi[k] += 0.5 * dt * m.f[i][k] / M[i]; });
      if (rig) rattleV(m.x, m.v);
      if (lang) { m.v.forEach((vi, i) => { const sn = Math.sqrt((1 - c * c) * kT / M[i]); for (let k = 0; k < 3; k++) vi[k] = c * vi[k] + sn * gauss(); });
        const ap = Math.abs(m.phi) * 180 / Math.PI, ns = ap > 150 ? 1 : (ap > 30 && ap < 90 ? 2 : state[b]); if (ns !== state[b] && state[b]) jumps++; state[b] = ns;
        if (t > 2) { const k = Math.floor((m.phi * 180 / Math.PI + 180) / 5); if (k >= 0 && k < 72) hist[k]++; } }
      else { const dE = energy(m) - E0[b]; if (!isFinite(dE) || Math.abs(dE) > 1000) blown = true; else { Es[b][0]++; Es[b][1] += dE; Es[b][2] += dE * dE; } }
    });
    if (lang && t > 2) nh += mol.length;
    t += dt;
    if (!lang && !blown) { Etr.push(energy(mol[0]) - E0[0]); if (Etr.length > 600) Etr.shift(); }
  }
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m; $('[data-lang]').hidden = st.mode !== 'lang'; $('[data-nve]').hidden = st.mode !== 'nve'; reset(); draw(); });
  const seg = (a, fn) => root.querySelectorAll(`[data-${a}]`).forEach(b => b.onclick = () => { root.querySelectorAll(`[data-${a}]`).forEach(z => z.classList.toggle('on', z === b)); fn(b.dataset[a]); reset(); draw(); });
  seg('s', v => st.s14 = +v); seg('r', v => st.rigid = v === '1'); seg('d', v => st.dt = +v);
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; reset(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !mol) return;
    cx.clearRect(0, 0, W, H);
    const half = H * 0.42, bw = Math.min(W * 0.42, half - 20), m = mol[0];
    // Newman projection along bond 2->3 of molecule 1
    const ccx = 10 + bw / 2, ccy = 10 + bw / 2, R = bw * 0.2, b2 = sub(m.x[2], m.x[1]), e = sc(b2, 1 / nrm(b2));
    const perp = a => { const w = sub(a, sc(e, dot(a, e))); return w; }, u1 = perp(sub(m.x[0], m.x[1]));
    const ex = sc(u1, 1 / nrm(u1)), ey = cross(e, ex), u4 = perp(sub(m.x[3], m.x[2])), a4 = Math.atan2(dot(u4, ey), dot(u4, ex));
    cx.strokeStyle = 'rgba(143,163,199,.5)'; cx.lineWidth = 2; cx.beginPath(); cx.arc(ccx, ccy, R, 0, 7); cx.stroke();
    const L = bw * 0.42; cx.strokeStyle = COL.coral; cx.lineWidth = 4; cx.beginPath(); cx.moveTo(ccx + R * Math.sin(a4), ccy - R * Math.cos(a4)); cx.lineTo(ccx + L * Math.sin(a4), ccy - L * Math.cos(a4)); cx.stroke();
    cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(ccx + L * Math.sin(a4), ccy - L * Math.cos(a4), 9, 0, 7); cx.fill();
    cx.strokeStyle = COL.cyan; cx.beginPath(); cx.moveTo(ccx, ccy); cx.lineTo(ccx, ccy - L); cx.stroke(); cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(ccx, ccy - L, 9, 0, 7); cx.fill(); cx.beginPath(); cx.arc(ccx, ccy, 6, 0, 7); cx.fill();
    const pd = m.phi * 180 / Math.PI, ap = Math.abs(pd);
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(`φ = ${pd.toFixed(0)}°  ${ap > 150 ? 'trans' : ap > 30 && ap < 90 ? 'gauche' : ap <= 30 ? 'cis' : 'crossing'}`, ccx, bw + 22);
    // top right: U(phi) curve or energy trace
    const tx0 = bw + 50, tx1 = W - 12, ty0 = 22, ty1 = half - 10;
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(tx0, ty0); cx.lineTo(tx0, ty1); cx.lineTo(tx1, ty1); cx.stroke();
    cx.textAlign = 'left'; cx.fillStyle = COL.ink;
    if (st.mode === 'lang') {
      const PX = p => tx0 + (p + 180) / 360 * (tx1 - tx0), Umax = 25, PY = u => ty1 - Math.min(u, Umax) / Umax * (ty1 - ty0), Uf = p => Udih(p * Math.PI / 180) + (st.s14 ? st.s14 * LJ(R14(p)) : 0);
      cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); for (let p = -180; p <= 180; p += 3) p > -180 ? cx.lineTo(PX(p), PY(Uf(p))) : cx.moveTo(PX(p), PY(Uf(p))); cx.stroke();
      cx.setLineDash([4, 4]); cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(tx0, PY(kB * st.T)); cx.lineTo(tx1, PY(kB * st.T)); cx.stroke(); cx.setLineDash([]);
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(PX(pd), PY(Uf(pd)), 5, 0, 7); cx.fill();
      cx.fillStyle = COL.ink; cx.fillText(W < 500 ? 'U(φ), kJ/mol' : 'U(φ), kJ/mol (dashed: k_B T)', tx0 + 6, ty0 - 6);
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [-120, 0, 120].forEach(p => cx.fillText(p + '°', PX(p), ty1 + 13));
    } else {   // side view of molecule 1 (projected on its own best plane: x-y of the starting frame)
      const cxm = (tx0 + tx1) / 2, cym = (ty0 + ty1) / 2, S = Math.min(tx1 - tx0, ty1 - ty0) / 0.5, c0 = sc(m.x.reduce((a, p) => [a[0] + p[0], a[1] + p[1], a[2] + p[2]], [0, 0, 0]), 0.25);
      const P2 = p => [cxm + (p[0] - c0[0]) * S, cym - (p[1] - c0[1]) * S];
      cx.strokeStyle = COL.muted; cx.lineWidth = 3; cx.beginPath(); m.x.forEach((p, i) => { const q = P2(p); i ? cx.lineTo(q[0], q[1]) : cx.moveTo(q[0], q[1]); }); cx.stroke();
      m.x.forEach((p, i) => { const q = P2(p); cx.fillStyle = i === 0 || i === 3 ? COL.coral : COL.cyan; cx.beginPath(); cx.arc(q[0], q[1], i === 0 || i === 3 ? 9 : 8, 0, 7); cx.fill(); });
      const bl = [0, 1, 2].map(i => nrm(sub(m.x[i + 1], m.x[i])) * 10).map(v => v.toFixed(3)).join(', ');
      cx.fillStyle = COL.ink; cx.fillText(`bonds (Å): ${bl}`, tx0 + 6, ty0 - 6);
    }
    // bottom: dihedral histogram vs exact (Langevin) or energy trace (NVE)
    const hx0 = 46, hx1 = W - 12, hy0 = half + 40, hy1 = H - 26;
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(hx0, hy0); cx.lineTo(hx0, hy1); cx.lineTo(hx1, hy1); cx.stroke();
    cx.fillStyle = COL.ink; cx.textAlign = 'left'; cx.font = '600 12px "Space Grotesk", sans-serif';
    if (st.mode === 'nve') {
      const lim = Math.max(0.05, ...Etr.map(Math.abs)) * 1.1, EX = k => hx0 + k / 599 * (hx1 - hx0), EY = d => (hy0 + hy1) / 2 - d / lim * (hy1 - hy0) / 2;
      cx.strokeStyle = 'rgba(143,163,199,.4)'; cx.beginPath(); cx.moveTo(hx0, EY(0)); cx.lineTo(hx1, EY(0)); cx.stroke();
      cx.strokeStyle = COL.cyan; cx.lineWidth = 1.5; cx.beginPath(); Etr.forEach((d, k) => k ? cx.lineTo(EX(k), EY(d)) : cx.moveTo(EX(k), EY(d))); cx.stroke();
      cx.fillStyle = blown ? COL.coral : COL.ink; cx.fillText(blown ? 'energy exploded: time step too large' : `molecule 1: total energy − start (axis ±${lim.toFixed(3)} kJ/mol)`, hx0 + 6, hy0 - 10);
    } else {
      const HX = p => hx0 + (p + 180) / 360 * (hx1 - hx0); let pmax = Math.max(...exact); if (nh) for (let k = 0; k < 72; k++) pmax = Math.max(pmax, hist[k] / nh / 5);
      const HY = p => hy1 - p / (pmax * 1.1) * (hy1 - hy0);
      cx.fillText(W < 500 ? 'φ: bars vs exact (line)' : 'dihedral distribution of 40 molecules (bars) vs exact (line)', hx0 + 6, hy0 - 10);
      if (nh) { cx.fillStyle = 'rgba(56,225,255,.45)'; for (let k = 0; k < 72; k++) { const p = hist[k] / nh / 5; cx.fillRect(HX(-180 + 5 * k) + 0.5, HY(p), (hx1 - hx0) / 72 - 1, hy1 - HY(p)); } }
      cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); for (let k = 0; k < 72; k++) { const X = HX(-177.5 + 5 * k), Y = HY(exact[k]); k ? cx.lineTo(X, Y) : cx.moveTo(X, Y); } cx.stroke();
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [-180, -60, 60, 180].forEach(p => cx.fillText(p + '°', HX(p), hy1 + 14));
    }
    // readouts
    o.T.textContent = st.T; o.t.textContent = t.toFixed(1);
    if (st.mode === 'lang') {
      let tr = 0; for (let k = 0; k < 72; k++) if (Math.abs(-177.5 + 5 * k) > 120) tr += hist[k];
      o.l1.textContent = `Trans fraction (exact ${exactTrans.toFixed(3)})`; o.v1.textContent = nh ? (tr / nh).toFixed(3) : '…';
      o.l2.textContent = 'Jumps trans ↔ gauche'; o.v2.textContent = jumps;
      o.l3.textContent = 'Time per jump, one molecule (ps)'; o.v3.textContent = jumps ? (mol.length * t / jumps).toFixed(0) : '…';
    } else {
      const sd = Es.map(e => e[0] ? Math.sqrt(Math.max(0, e[2] / e[0] - (e[1] / e[0]) ** 2)) : 0), mean = sd.reduce((a, b) => a + b, 0) / sd.length;
      o.l1.textContent = 'Energy fluctuation (kJ/mol)'; o.v1.textContent = blown ? 'unstable' : t > 0 ? mean.toFixed(4) : '…';
      o.l2.textContent = 'Steps per picosecond'; o.v2.textContent = (1000 / st.dt).toFixed(0);
      o.l3.textContent = 'Bonds'; o.v3.textContent = st.rigid ? 'rigid, SHAKE' : 'flexible';
    }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return;
    if (st.mode === 'nve' && (blown || t >= 10)) { st.running = false; $('[data-k=run]').textContent = '▶ Run'; draw(); return; }
    const n = st.mode === 'lang' ? 120 : Math.max(5, Math.round(40 / st.dt)); for (let q = 0; q < n; q++) stepAll(); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
