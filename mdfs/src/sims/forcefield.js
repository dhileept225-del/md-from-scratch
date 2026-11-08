// Force-field lab: two rigid water molecules (3-site models) in the screen plane. Drag molecule B, rotate it, and
// change the model parameters. Shows LJ and Coulomb contributions, each charge pair, and an energy-vs-distance scan.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const COUL = 138.935, ROH = 0.09572, HOH = 104.52 * Math.PI / 180;
const MODELS = { tip3p: { qO: -0.834, s: 0.315061, e: 0.6364, r: 0.09572, a: 104.52 }, spce: { qO: -0.8476, s: 0.316557, e: 0.650194, r: 0.1, a: 109.47 } };

export function mountForceField(root) {
  const st = { ...MODELS.tip3p, rotB: 0, B: [0.28, 0], drag: false };
  // molecule A fixed: O at origin, one O-H pointing at +x (donor); molecule B: O at st.B, rotation rotB (deg)
  function sites(O, rot) {
    const a = st.a * Math.PI / 180, r = st.r, t = rot * Math.PI / 180;
    return [O, [O[0] + r * Math.cos(t - a / 2), O[1] + r * Math.sin(t - a / 2)], [O[0] + r * Math.cos(t + a / 2), O[1] + r * Math.sin(t + a / 2)]];
  }
  const A = () => sites([0, 0], st.a / 2);                         // first H exactly along +x
  function energy(BO, rot) {
    const a = A(), b = sites(BO, rot), q = [st.qO, -st.qO / 2, -st.qO / 2], pairs = []; let Ec = 0;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { const r = Math.hypot(a[i][0] - b[j][0], a[i][1] - b[j][1]); const e = COUL * q[i] * q[j] / r; Ec += e; pairs.push([a[i], b[j], e]); }
    const rOO = Math.hypot(BO[0], BO[1]), sr6 = (st.s / rOO) ** 6, Elj = 4 * st.e * (sr6 * sr6 - sr6);
    return { Elj, Ec, E: Elj + Ec, pairs, rOO, a, b };
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Drag the right molecule</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="tip3p" class="on">TIP3P</button><button data-m="spce">SPC/E</button></div>
    <label><span>Rotate molecule B</span><output data-o="rot"></output><input type="range" data-k="rot" min="0" max="360" step="1" value="0"></label>
    <label><span>Oxygen charge (H = −q/2)</span><output data-o="q"></output><input type="range" data-k="q" min="-1.2" max="-0.4" step="0.01" value="-0.834"></label>
    <label><span>LJ <i class="nc">σ</i> of O</span><output data-o="s"></output><input type="range" data-k="s" min="0.26" max="0.36" step="0.001" value="0.315"></label>
    <label><span>LJ <i class="nc">ε</i> of O</span><output data-o="e"></output><input type="range" data-k="e" min="0" max="2" step="0.01" value="0.636"></label>
    <div class="row"><button class="btn sm primary" data-k="opt">Find the minimum</button></div>
    <div class="readouts">
      <div><span>O–O distance</span><b data-o="r"></b></div>
      <div><span>Total energy</span><b data-o="E"></b></div>
      <div><span>Lennard-Jones part</span><b data-o="lj"></b></div>
      <div><span>Coulomb part</span><b data-o="c"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const sync = () => { $('[data-k=q]').value = st.qO; $('[data-k=s]').value = st.s; $('[data-k=e]').value = st.e; $('[data-k=rot]').value = st.rotB; draw(); };
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { Object.assign(st, MODELS[b.dataset.m]); root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); sync(); });
  const unpreset = () => root.querySelectorAll('[data-m]').forEach(z => z.classList.remove('on'));
  $('[data-k=rot]').oninput = e => { st.rotB = +e.target.value; draw(); };
  $('[data-k=q]').oninput = e => { st.qO = +e.target.value; unpreset(); draw(); };
  $('[data-k=s]').oninput = e => { st.s = +e.target.value; unpreset(); draw(); };
  $('[data-k=e]').oninput = e => { st.e = +e.target.value; unpreset(); draw(); };
  $('[data-k=opt]').onclick = () => {          // brute-force search over position and rotation (in the plane)
    let best = null;
    for (let r = 0.22; r <= 0.40; r += 0.0025) for (let ang = -60; ang <= 60; ang += 3) for (let rot = 0; rot < 360; rot += 3) {
      const BO = [r * Math.cos(ang * Math.PI / 180), r * Math.sin(ang * Math.PI / 180)], E = energy(BO, rot).E;
      if (!best || E < best.E) best = { E, BO, rot };
    }
    st.B = best.BO; st.rotB = best.rot; sync();
  };

  let W = 600, H = 600, S = 1000, ox = 200, oy = 200;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  let grab = null;
  const toWorld = (px, py) => [(px - ox) / S, -(py - oy) / S];
  cv.addEventListener('pointerdown', e => { const r = cv.getBoundingClientRect(), p = toWorld(e.clientX - r.left, e.clientY - r.top); if (Math.hypot(p[0] - st.B[0], p[1] - st.B[1]) < 0.12) { grab = [p[0] - st.B[0], p[1] - st.B[1]]; cv.setPointerCapture(e.pointerId); } });
  cv.addEventListener('pointermove', e => { if (!grab) return; const r = cv.getBoundingClientRect(), p = toWorld(e.clientX - r.left, e.clientY - r.top); let B = [p[0] - grab[0], p[1] - grab[1]]; const d = Math.hypot(B[0], B[1]); if (d < 0.2) B = B.map(v => v / d * 0.2); st.B = B; draw(); });
  const up = () => { grab = null; }; cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);

  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const plotH = Math.max(160, H * 0.36), areaH = H - plotH - 24;
    S = Math.min(W / 0.78, areaH / 0.46); ox = W * 0.3; oy = areaH * 0.5;
    const res = energy(st.B, st.rotB), P = q => [ox + q[0] * S, oy - q[1] * S];
    // charge-pair lines
    const emax = Math.max(...res.pairs.map(p => Math.abs(p[2])));
    res.pairs.forEach(([p, q, e]) => { const A_ = P(p), B_ = P(q); cx.strokeStyle = e > 0 ? `rgba(255,107,107,${0.15 + 0.6 * Math.abs(e) / emax})` : `rgba(56,225,255,${0.15 + 0.6 * Math.abs(e) / emax})`; cx.lineWidth = 1 + 4 * Math.abs(e) / emax; cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(A_[0], A_[1]); cx.lineTo(B_[0], B_[1]); cx.stroke(); });
    cx.setLineDash([]);
    [res.a, res.b].forEach(m => {
      const Pm = m.map(P); cx.strokeStyle = 'rgba(143,163,199,.9)'; cx.lineWidth = 7; cx.beginPath(); cx.moveTo(Pm[1][0], Pm[1][1]); cx.lineTo(Pm[0][0], Pm[0][1]); cx.lineTo(Pm[2][0], Pm[2][1]); cx.stroke();
      cx.fillStyle = COL.ink; [1, 2].forEach(k => { cx.beginPath(); cx.arc(Pm[k][0], Pm[k][1], 0.012 * S, 0, 7); cx.fill(); });
      cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(Pm[0][0], Pm[0][1], 0.02 * S, 0, 7); cx.fill();
      // LJ sphere (sigma / 2 radius, faint)
      cx.strokeStyle = 'rgba(255,194,75,.35)'; cx.lineWidth = 1.5; cx.beginPath(); cx.arc(Pm[0][0], Pm[0][1], st.s / 2 * S, 0, 7); cx.stroke();
    });
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.fillStyle = COL.cyan; cx.fillText('attracting charge pairs', 12, 22); cx.fillStyle = COL.coral; cx.fillText('repelling charge pairs', 12, 40); cx.fillStyle = COL.amber; cx.fillText('LJ size of O (σ/2)', 12, 58);
    cx.fillStyle = COL.muted; cx.fillText('A', P([0, 0])[0] - 30, P([0, 0])[1] - 22); cx.fillText('B', P(st.B)[0] + 22, P(st.B)[1] - 22);
    // scan along the O-O line with the current orientation
    const x0 = 46, x1 = W - 12, y0 = areaH + 30, y1 = H - 26, dir = st.B.map(v => v / res.rOO), rs = [];
    for (let r = 0.22; r <= 0.6; r += 0.005) rs.push([r, energy([dir[0] * r, dir[1] * r], st.rotB).E]);
    const lo = Math.min(-35, ...rs.map(q => q[1])), hi = 20, X = r => x0 + (r - 0.22) / 0.38 * (x1 - x0), Y = E => y1 - (Math.min(hi, E) - lo) / (hi - lo) * (y1 - y0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0.25, 0.3, 0.4, 0.5, 0.6].forEach(r => cx.fillText(r.toFixed(2), X(r), y1 + 14));
    cx.textAlign = 'right'; cx.fillText('0', x0 - 5, Y(0) + 4); cx.fillText(lo.toFixed(0), x0 - 5, y1);
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.cyan; cx.fillText('energy (kJ/mol) vs O–O distance, same orientation', x0 + 6, y0 - 8);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2.2; cx.beginPath(); rs.forEach((q, k) => (k ? cx.lineTo(X(q[0]), Y(q[1])) : cx.moveTo(X(q[0]), Y(q[1])))); cx.stroke();
    if (res.rOO <= 0.6) { cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X(res.rOO), Y(res.E), 6, 0, 7); cx.fill(); }
    o.r.textContent = res.rOO.toFixed(3) + ' nm'; o.E.textContent = res.E.toFixed(2) + ' kJ/mol'; o.lj.textContent = res.Elj.toFixed(2); o.c.textContent = res.Ec.toFixed(2);
    o.rot.textContent = st.rotB + '°'; o.q.textContent = `${st.qO.toFixed(3)} e  (H ${(-st.qO / 2).toFixed(3)})`; o.s.textContent = st.s.toFixed(4) + ' nm'; o.e.textContent = st.e.toFixed(3) + ' kJ/mol';
  }
  resize();
  return () => { ro.disconnect(); };
}
