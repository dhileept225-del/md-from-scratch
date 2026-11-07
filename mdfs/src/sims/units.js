// Reduced-units explorer: four noble gases, real vs reduced LJ curves, and reduced -> real conversion.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7', green: '#7EE787' };
const kB = 1.380649e-23, amu = 1.66053907e-27, NA = 6.02214076e23;
const GASES = { Ne: { s: 2.75, e: 36.0, m: 20.18, c: COL.cyan }, Ar: { s: 3.405, e: 119.8, m: 39.948, c: COL.amber }, Kr: { s: 3.60, e: 171.0, m: 83.80, c: COL.green }, Xe: { s: 4.10, e: 221.0, m: 131.29, c: COL.coral } };
const ENERGY = { 'kcal/mol': 4.184, 'kJ/mol': 1, 'eV': 96.48533, 'K (E/kB)': 0.0083144626 };   // in kJ/mol

export function mountUnits(root) {
  const st = { gas: 'Ar', reduced: false, k: 0, T: 0.72, rho: 0.84, dt: 0.005 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view units-view"><canvas></canvas><div class="sim-tip">Real units: r in Å, energy in K (ε/k<sub>B</sub>)</div></div>
  <div class="sim-ctrl">
    <div class="seg">${Object.keys(GASES).map(g => `<button data-g="${g}" class="${g === 'Ar' ? 'on' : ''}">${g}</button>`).join('')}</div>
    <label class="check"><input type="checkbox" data-k="red"> Show in reduced units (r/σ, V/ε)</label>
    <label>Temperature T* <output data-o="Ts"></output><input type="range" data-k="T" min="0.4" max="2.5" step="0.01" value="0.72"></label>
    <label>Density ρ* <output data-o="rs"></output><input type="range" data-k="rho" min="0.02" max="1.05" step="0.01" value="0.84"></label>
    <label>Time step Δt* <output data-o="dts"></output><input type="range" data-k="dt" min="0.001" max="0.01" step="0.0005" value="0.005"></label>
    <div class="readouts">
      <div><span>Temperature</span><b data-o="T"></b></div>
      <div><span>Density</span><b data-o="rho"></b></div>
      <div><span>Time unit τ</span><b data-o="tau"></b></div>
      <div><span>Time step</span><b data-o="dt"></b></div>
      <div><span>Pressure unit</span><b data-o="P"></b></div>
      <div><span>ε</span><b data-o="eps"></b></div>
      <div class="wide"><span>Rough state (LJ phase diagram)</span><b data-o="phase"></b></div>
    </div>
    <div class="conv"><span class="muted small">Energy converter</span>
      <div class="row"><input type="number" data-k="ev" value="1" step="any"><select data-k="eu">${Object.keys(ENERGY).map(u => `<option>${u}</option>`).join('')}</select></div>
      <div class="conv-out" data-o="conv"></div></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.units-view'), tip = root.querySelector('.sim-tip');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-g]').forEach(b => b.onclick = () => { st.gas = b.dataset.g; root.querySelectorAll('[data-g]').forEach(x => x.classList.toggle('on', x === b)); update(); });
  $('[data-k=red]').onchange = e => { st.reduced = e.target.checked; };
  [['T', 'T'], ['rho', 'rho'], ['dt', 'dt']].forEach(([k, f]) => $(`[data-k=${k}]`).oninput = e => { st[f] = +e.target.value; update(); });
  $('[data-k=ev]').oninput = conv; $('[data-k=eu]').onchange = conv;

  function units(g) {
    const s = g.s * 1e-10, e = g.e * kB, m = g.m * amu;
    return { tau: s * Math.sqrt(m / e) * 1e12, T: g.e, rho: m / s ** 3 / 1000, P: e / s ** 3 / 1e5, ekj: e * NA / 1000 };
  }
  function phase(T, r) {
    if (T > 1.31) return r < 0.3 ? 'supercritical fluid (low density)' : 'supercritical fluid';
    if (r < 0.1) return 'gas';
    if (T < 0.69 && r > 0.84) return 'solid (fcc crystal)';
    if (r > 0.65) return T < 0.69 ? 'solid or dense liquid' : 'liquid';
    return 'gas + liquid coexistence (would separate)';
  }
  function update() {
    const g = GASES[st.gas], u = units(g);
    o.Ts.textContent = st.T.toFixed(2); o.rs.textContent = st.rho.toFixed(2); o.dts.textContent = st.dt.toFixed(4);
    o.T.textContent = (st.T * u.T).toFixed(0) + ' K'; o.rho.textContent = (st.rho * u.rho).toFixed(3) + ' g/cm³';
    o.tau.textContent = u.tau.toFixed(3) + ' ps'; o.dt.textContent = (st.dt * u.tau * 1000).toFixed(1) + ' fs';
    o.P.textContent = u.P.toFixed(0) + ' bar'; o.eps.textContent = `${u.ekj.toFixed(3)} kJ/mol · ${(u.ekj / 4.184).toFixed(3)} kcal/mol`;
    o.phase.textContent = phase(st.T, st.rho);
  }
  function conv() {
    const v = +$('[data-k=ev]').value, kj = v * ENERGY[$('[data-k=eu]').value];
    o.conv.innerHTML = Object.entries(ENERGY).map(([n, f]) => `<span>${(kj / f).toPrecision(5)} <em>${n}</em></span>`).join('');
  }

  let W = 600, H = 440;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const LJ = (r, s, e) => 4 * e * ((s / r) ** 12 - (s / r) ** 6);

  function draw() {
    st.k += ((st.reduced ? 1 : 0) - st.k) * 0.08;
    const k = st.k, x0 = 56, x1 = W - 20, y0 = 30, y1 = H - 40;
    // axes ranges blend from real (Å, K) to reduced (r/σ, V/ε)
    const rlo = 2.3 + (0.85 - 2.3) * k, rhi = 8.5 + (3.0 - 8.5) * k, vlo = -240 + (-1.25 + 240) * k, vhi = 160 + (1.2 - 160) * k;
    const X = r => x0 + (r - rlo) / (rhi - rlo) * (x1 - x0), Y = v => y1 - (v - vlo) / (vhi - vlo) * (y1 - y0);
    cx.clearRect(0, 0, W, H);
    cx.strokeStyle = 'rgba(143,163,199,.55)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '500 13px "Space Grotesk", sans-serif'; cx.textAlign = 'right';
    cx.fillText(k > 0.5 ? 'r / σ' : 'r (Å)', x1, Y(0) - 8); cx.fillText(k > 0.5 ? 'V/ε' : 'V (K)', x0 - 6, y0 + 4); cx.textAlign = 'left';
    tip.innerHTML = k > 0.5 ? 'Reduced units: r in σ, energy in ε' : 'Real units: r in Å, energy in K (ε/k<sub>B</sub>)';
    cx.save(); cx.beginPath(); cx.rect(x0, y0 - 10, x1 - x0, y1 - y0 + 10); cx.clip();
    Object.entries(GASES).forEach(([name, g]) => {
      const sel = name === st.gas; cx.strokeStyle = g.c; cx.globalAlpha = sel ? 1 : 0.55; cx.lineWidth = sel ? 4 : 2; cx.beginPath(); let first = true;
      for (let r = 0.86 * g.s; r <= 9 * g.s / 3.4 + 9; r += 0.01) {
        const xr = r + (r / g.s - r) * k, yv = LJ(r, g.s, g.e) * (1 + (1 / g.e - 1) * k);
        if (xr > rhi) break; first ? cx.moveTo(X(xr), Y(yv)) : cx.lineTo(X(xr), Y(yv)); first = false;
      }
      cx.stroke();
      if (k < 0.5) { const rm = g.s * 1.122; cx.fillStyle = g.c; cx.font = '700 14px "Space Grotesk", sans-serif'; cx.fillText(name, X(rm) + 8, Y(-g.e) + 16); }
    });
    cx.restore(); cx.globalAlpha = 1;
    if (k > 0.5) { cx.fillStyle = COL.ink; cx.font = '700 15px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; cx.fillText('All four gases: one curve', x1, y0 + 24); }
  }
  let raf = 0; const loop = () => { raf = requestAnimationFrame(loop); draw(); };
  update(); conv(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); };
}
