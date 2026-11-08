// Water-model explorer: geometry, charges and the electrostatic potential of popular rigid water models,
// with a comparison of liquid properties (literature values at 298 K).
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const COUL = 138.935, DEBYE = 48.032;
const MODELS = {
  'TIP3P': { r: 0.09572, a: 104.52, qH: 0.4170, dM: 0, s: 0.315061, e: 0.6364, dimer: [-27.2, 0.275], props: { rho: 0.98, D: 5.5, eps: 94, Tm: 146, TMD: 182 } },
  'SPC/E': { r: 0.1, a: 109.47, qH: 0.4238, dM: 0, s: 0.316557, e: 0.65019, dimer: [-30.0, 0.273], props: { rho: 0.994, D: 2.5, eps: 71, Tm: 215, TMD: 241 } },
  'TIP4P/2005': { r: 0.09572, a: 104.52, qH: 0.5564, dM: 0.01546, s: 0.31589, e: 0.7749, dimer: [-28.7, 0.277], props: { rho: 0.997, D: 2.1, eps: 58, Tm: 252, TMD: 278 } },
  'OPC': { r: 0.08724, a: 103.6, qH: 0.6791, dM: 0.01594, s: 0.316655, e: 0.89036, dimer: [-30.7, 0.280], props: { rho: 0.997, D: 2.3, eps: 78, Tm: null, TMD: 272 } },
};
const EXP = { rho: 0.997, D: 2.3, eps: 78.4, Tm: 273, TMD: 277 };
const PROPS = { eps: ['Dielectric constant', ''], D: ['Diffusion (10⁻⁵ cm²/s)', ''], Tm: ['Melting point of ice', ' K'], TMD: ['Temperature of max. density', ' K'], rho: ['Density (g/cm³)', ''] };

export function mountWater(root) {
  const st = { m: 'TIP4P/2005', p: 'eps' };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg">${Object.keys(MODELS).map(k => `<button data-m="${k}" class="${k === st.m ? 'on' : ''}">${k}</button>`).join('')}</div>
    <div class="seg">${Object.entries(PROPS).map(([k, v]) => `<button data-p="${k}" class="${k === st.p ? 'on' : ''}">${{ eps: 'ε', D: 'D', Tm: 'T melt', TMD: 'T max ρ', rho: 'ρ' }[k]}</button>`).join('')}</div>
    <div class="readouts">
      <div><span>O–H / angle</span><b data-o="geo"></b></div>
      <div><span>Charges</span><b data-o="q"></b></div>
      <div><span>Dipole moment</span><b data-o="mu"></b></div>
      <div><span>Gas dimer (exp. −21 kJ/mol)</span><b data-o="dim"></b></div>
      <div class="wide"><span>LJ on oxygen</span><b data-o="lj"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.m = b.dataset.m; root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); draw(); });
  root.querySelectorAll('[data-p]').forEach(b => b.onclick = () => { st.p = b.dataset.p; root.querySelectorAll('[data-p]').forEach(z => z.classList.toggle('on', z === b)); draw(); });

  function sites(p) {   // 2D: O at origin, H atoms upwards, M on the bisector
    const h = p.a * Math.PI / 360;
    return { O: [0, 0], H1: [-p.r * Math.sin(h), p.r * Math.cos(h)], H2: [p.r * Math.sin(h), p.r * Math.cos(h)], M: [0, p.dM] };
  }
  let W = 600, H = 600;
  const off = document.createElement('canvas'), GW = 110, GH = 80; off.width = GW; off.height = GH; const octx = off.getContext('2d'), img = octx.createImageData(GW, GH);
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function draw() {
    if (!W) return;
    const p = MODELS[st.m], S_ = sites(p), qs = [[S_.H1, p.qH], [S_.H2, p.qH], [p.dM ? S_.M : S_.O, -2 * p.qH]];
    cx.clearRect(0, 0, W, H);
    // electrostatic potential map in the molecular plane
    const mapH = H * 0.52, span = 0.44, sc = Math.min(W, mapH * GW / GH) / span, cxm = W / 2, cym = mapH * 0.58;
    for (let gy = 0; gy < GH; gy++) for (let gx = 0; gx < GW; gx++) {
      const X = (gx - GW / 2) / GW * span, Y = (GH * 0.58 - gy) / GH * span * GH / GW; let V = 0;
      for (const [r, q] of qs) V += COUL * q / Math.max(0.02, Math.hypot(X - r[0], Y - r[1]));
      const t = Math.max(-1, Math.min(1, V / 600)), k = 4 * (gy * GW + gx), c = t < 0 ? [255, 107, 107] : [56, 225, 255];
      const edge = Math.min(1, gx / 12, (GW - 1 - gx) / 12, gy / 10, (GH - 1 - gy) / 10);
      img.data[k] = c[0]; img.data[k + 1] = c[1]; img.data[k + 2] = c[2]; img.data[k + 3] = Math.round(190 * Math.abs(t) ** 0.7 * edge);
    }
    octx.putImageData(img, 0, 0); cx.imageSmoothingEnabled = true; cx.drawImage(off, cxm - GW / 2 / GW * span * sc, cym - GH * 0.58 / GW * span * sc, span * sc, span * sc * GH / GW);
    const P = r => [cxm + r[0] * sc, cym - r[1] * sc];
    cx.strokeStyle = 'rgba(255,194,75,.5)'; cx.setLineDash([5, 5]); cx.lineWidth = 1.5; cx.beginPath(); cx.arc(...P(S_.O), p.s / 2 * sc, 0, 7); cx.stroke(); cx.setLineDash([]);
    cx.strokeStyle = 'rgba(242,245,250,.8)'; cx.lineWidth = 5; cx.beginPath(); cx.moveTo(...P(S_.H1)); cx.lineTo(...P(S_.O)); cx.lineTo(...P(S_.H2)); cx.stroke();
    cx.fillStyle = COL.ink; [S_.H1, S_.H2].forEach(r => { cx.beginPath(); cx.arc(...P(r), 0.018 * sc, 0, 7); cx.fill(); });
    cx.fillStyle = p.dM ? 'rgba(255,107,107,.6)' : COL.coral; cx.beginPath(); cx.arc(...P(S_.O), 0.03 * sc, 0, 7); cx.fill();
    if (p.dM) { cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(...P(S_.M), 0.009 * sc, 0, 7); cx.fill(); cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('M', P(S_.M)[0] + 8, P(S_.M)[1] + 4); }
    cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.fillStyle = COL.coral; cx.fillText('negative potential', 12, 22); cx.fillStyle = COL.cyan; cx.fillText('positive potential', 12, 40); cx.fillStyle = COL.amber; cx.fillText('LJ size (σ/2)', 12, 58);
    // property bars
    const x0 = 20, x1 = W - 20, top = mapH + 40, yb = H - 34, names = [...Object.keys(MODELS), 'experiment'];
    const vals = names.map(n => n === 'experiment' ? EXP[st.p] : MODELS[n].props[st.p]), vmax = Math.max(...vals.filter(v => v != null)) * 1.15;
    const bw = (x1 - x0) / names.length;
    cx.fillStyle = COL.ink; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(PROPS[st.p][0] + (['eps', 'D', 'rho'].includes(st.p) ? ' at 298 K' : '') + ' (literature)', x0, top - 12);
    cx.strokeStyle = 'rgba(242,245,250,.5)'; cx.setLineDash([6, 5]); cx.lineWidth = 1.5; const yE = yb - EXP[st.p] / (Math.max(...[...Object.keys(MODELS).map(n => MODELS[n].props[st.p]), EXP[st.p]].filter(v => v != null)) * 1.15) * (yb - top - 10); cx.beginPath(); cx.moveTo(x0, yE); cx.lineTo(x1, yE); cx.stroke(); cx.setLineDash([]);
    const Yv = v => yb - v / vmax * (yb - top - 10);
    names.forEach((n, k) => {
      const v = vals[k], x = x0 + k * bw + bw * 0.15, w = bw * 0.7, col = n === 'experiment' ? COL.muted : n === st.m ? COL.amber : 'rgba(56,225,255,.55)';
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; cx.fillText(n === 'experiment' ? 'exp.' : n, x + w / 2, yb + 14);
      if (v == null) { cx.fillText('n/a', x + w / 2, yb - 6); return; }
      cx.fillStyle = col; cx.fillRect(x, Yv(v), w, yb - Yv(v)); cx.fillStyle = 'rgba(6,12,26,.85)'; cx.fillRect(x + w / 2 - 26, Yv(v) - 20, 52, 17); cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.fillText(v + PROPS[st.p][1], x + w / 2, Yv(v) - 7);
    });
    // readouts
    const mu = Math.hypot(...[0, 1].map(d => qs.reduce((s, [r, q]) => s + q * r[d], 0))) * DEBYE;
    o.geo.textContent = `${p.r.toFixed(4)} nm / ${p.a}°`;
    o.q.textContent = p.dM ? `H +${p.qH.toFixed(4)}, M ${(-2 * p.qH).toFixed(4)} (${(p.dM * 10).toFixed(4)} Å from O)` : `H +${p.qH.toFixed(4)}, O ${(-2 * p.qH).toFixed(4)}`;
    o.mu.textContent = mu.toFixed(2) + ' D (gas 1.85)'; o.dim.textContent = `${p.dimer[0]} kJ/mol at ${p.dimer[1]} nm`;
    o.lj.textContent = `σ = ${p.s} nm, ε = ${p.e} kJ/mol; ${p.dM ? '4 sites (M is massless, no LJ)' : '3 sites'}`;
  }
  return () => ro.disconnect();
}
