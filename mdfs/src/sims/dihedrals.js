// Dihedral lab: build a torsion potential U(φ) = Σ k_n [1 + cos(nφ)] from three cosine terms (k_n < 0 is the same as
// δ_n = 180°), see the Boltzmann population at temperature T, and watch a Langevin rotor hop between the wells.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RB = [9.28, 12.16, -13.12, -3.06, 26.24, -31.5];            // butane, GROMACS (kJ/mol), psi = phi - 180
const Urb = ph => { const c = Math.cos(ph - Math.PI); return RB.reduce((s, C, n) => s + C * c ** n, 0); };
const fourier = f => [1, 2, 3, 4, 5].map(n => { let a = 0; for (let k = 0; k < 720; k++) { const ph = k * Math.PI / 360; a += f(ph) * Math.cos(n * ph); } return +(a / 360).toFixed(2); });
const PRESETS = { butane: fourier(Urb), ethane: [0, 0, 6, 0, 0], double: [0, -25, 0, 0, 0] };   // butane: exact RB as n = 1..5
const gauss = () => Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());

export function mountDihedrals(root) {
  const st = { k: [...PRESETS.butane], T: 300, running: true };
  let phi = Math.PI, w = 0, t = 0, hops = 0, region = 't', drag = false, trail = [];
  const U = p => st.k.reduce((s, k, n) => s + k * (1 + Math.cos((n + 1) * p)), 0);
  const dU = p => st.k.reduce((s, k, n) => s - k * (n + 1) * Math.sin((n + 1) * p), 0);
  const absDeg = p => Math.abs(((p * 180 / Math.PI + 540) % 360) - 180);
  const regionOf = p => (absDeg(p) > 120 ? 't' : 'g');
  const nameOf = p => { const d = absDeg(p); return d > 120 ? 'trans' : d > 30 ? 'gauche' : 'cis'; };
  function step() {                                   // Langevin rotor, BAOAB, I = 0.34 amu nm^2, gamma = 5/ps, dt = 1 fs
    const I = 0.34, g = 5, dt = 0.001, kT = 0.008314 * st.T, c1 = Math.exp(-g * dt), c2 = Math.sqrt((1 - c1 * c1) * kT / I);
    for (let s = 0; s < 400; s++) {
      w -= 0.5 * dt * dU(phi) / I; phi += 0.5 * dt * w; w = c1 * w + c2 * gauss(); phi += 0.5 * dt * w; w -= 0.5 * dt * dU(phi) / I;
      phi = ((phi + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI; t += dt;
      const r = regionOf(phi); if (r !== region) { hops++; region = r; }
    }
    trail.push(phi); if (trail.length > 160) trail.shift();
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view tm-view"><canvas style="touch-action:none"></canvas><div class="sim-tip" style="top:8px;right:10px;bottom:auto;left:auto">Drag the red atom to twist</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-pre="butane" class="on">Butane</button><button data-pre="ethane">3-fold</button><button data-pre="double">Double bond</button></div>
    ${[1, 2, 3].map(n => `<label><span>k<sub>${n}</sub> (n = ${n})</span><output data-o="k${n}"></output><input type="range" data-k="k${n}" min="-30" max="30" step="0.5"></label>`).join('')}
    <label>Temperature <output data-o="T"></output><input type="range" data-k="T" min="100" max="1000" step="25" value="300"></label>
    <div class="readouts">
      <div><span>Trans (<i class="nc">|φ|</i> &gt; 120°)</span><b data-o="pt"></b></div>
      <div><span>Gauche (30–120°)</span><b data-o="pg"></b></div>
      <div><span>Cis (<i class="nc">|φ|</i> &lt; 30°)</span><b data-o="pc"></b></div>
      <div><span>Barrier from trans</span><b data-o="bar"></b></div>
      <div class="wide"><span>Rotor: trans ↔ other changes per ns</span><b data-o="rate"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.tm-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const resetRate = () => { t = 0; hops = 0; region = regionOf(phi); };
  const syncSliders = () => [1, 2, 3].forEach(n => { $(`[data-k=k${n}]`).value = st.k[n - 1]; });
  root.querySelectorAll('[data-pre]').forEach(b => b.onclick = () => { st.k = [...PRESETS[b.dataset.pre]]; syncSliders(); root.querySelectorAll('[data-pre]').forEach(z => z.classList.toggle('on', z === b)); resetRate(); });
  [1, 2, 3].forEach(n => $(`[data-k=k${n}]`).oninput = e => { st.k[n - 1] = +e.target.value; st.k[3] = st.k[4] = 0; root.querySelectorAll('[data-pre]').forEach(z => z.classList.remove('on')); resetRate(); });
  $('[data-k=T]').oninput = e => { st.T = +e.target.value; resetRate(); };
  syncSliders();

  let W = 600, H = 600, nw = { x: 0, y: 0, R: 0 };
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  const setFromPointer = e => { const r = cv.getBoundingClientRect(), px = e.clientX - r.left - nw.x, py = e.clientY - r.top - nw.y; phi = Math.atan2(py, px) + Math.PI / 2; phi = ((phi + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI; w = 0; };
  cv.addEventListener('pointerdown', e => { const r = cv.getBoundingClientRect(); if (Math.hypot(e.clientX - r.left - nw.x, e.clientY - r.top - nw.y) < nw.R * 1.4) { drag = true; cv.setPointerCapture(e.pointerId); setFromPointer(e); } });
  cv.addEventListener('pointermove', e => { if (drag) setFromPointer(e); });
  const up = () => { if (drag) { drag = false; resetRate(); } }; cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);

  function draw() {
    cx.clearRect(0, 0, W, H);
    // Newman projection
    const top = H * 0.42; nw = { x: W / 2, y: top / 2 + 6, R: Math.min(W * 0.3, top * 0.38) };
    const { x: X0, y: Y0, R } = nw;
    cx.strokeStyle = COL.muted; cx.lineWidth = 3; cx.beginPath(); cx.arc(X0, Y0, R * 0.38, 0, 7); cx.stroke();
    const fa = -Math.PI / 2, ba = -Math.PI / 2 + phi;
    cx.lineCap = 'round'; cx.lineWidth = 6;
    cx.strokeStyle = COL.coral; cx.beginPath(); cx.moveTo(X0 + R * 0.38 * Math.cos(ba), Y0 + R * 0.38 * Math.sin(ba)); cx.lineTo(X0 + R * Math.cos(ba), Y0 + R * Math.sin(ba)); cx.stroke();
    cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(X0 + R * Math.cos(ba), Y0 + R * Math.sin(ba), 13, 0, 7); cx.fill();
    cx.strokeStyle = COL.amber; cx.beginPath(); cx.moveTo(X0, Y0); cx.lineTo(X0 + R * Math.cos(fa), Y0 + R * Math.sin(fa)); cx.stroke();
    cx.fillStyle = COL.amber; cx.beginPath(); cx.arc(X0 + R * Math.cos(fa), Y0 + R * Math.sin(fa), 13, 0, 7); cx.fill();
    cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(X0, Y0, 9, 0, 7); cx.fill();
    let deg = phi * 180 / Math.PI; cx.fillStyle = COL.ink; cx.font = '600 15px "Space Grotesk", sans-serif'; cx.textAlign = 'left';
    cx.fillText(`φ = ${deg.toFixed(0)}°  (${nameOf(phi)})`, 12, 22);
    cx.fillStyle = COL.muted; cx.font = '12px "Space Grotesk", sans-serif'; cx.fillText(`U = ${(U(phi) - minU).toFixed(1)} kJ/mol`, 12, 40);
    // U(phi) plot with Boltzmann distribution
    const x0 = 44, x1 = W - 12, y0 = top + 24, y1 = H - 26, Xp = d => x0 + (d + 180) / 360 * (x1 - x0);
    const N = 241, us = [], ps = []; let mx = -1e9, mn = 1e9, Z = 0; const kT = 0.008314 * st.T;
    for (let k = 0; k < N; k++) { const d = -180 + 360 * k / (N - 1), u = U(d * Math.PI / 180); us.push(u); mx = Math.max(mx, u); mn = Math.min(mn, u); }
    us.forEach(u => { const p = Math.exp(-(u - mn) / kT); ps.push(p); Z += p; }); const pmax = Math.max(...ps);
    const span = Math.max(5, mx - mn), Yp = u => y1 - (u - mn) / span * (y1 - y0);
    cx.fillStyle = 'rgba(255,194,75,.22)'; cx.beginPath(); cx.moveTo(x0, y1); ps.forEach((p, k) => cx.lineTo(Xp(-180 + 360 * k / (N - 1)), y1 - p / pmax * (y1 - y0) * 0.6)); cx.lineTo(x1, y1); cx.closePath(); cx.fill();
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [-180, -120, -60, 0, 60, 120, 180].forEach(d => cx.fillText(d + '°', Xp(d), y1 + 14));
    cx.textAlign = 'right'; cx.fillText(span.toFixed(0), x0 - 5, y0 + 8); cx.fillText('0', x0 - 5, y1);
    cx.textAlign = 'left'; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.fillStyle = COL.cyan; cx.fillText('U(φ), kJ/mol', x0 + 6, y0 - 8); cx.fillStyle = COL.amber; cx.fillText('Boltzmann population', x0 + 110, y0 - 8);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2.5; cx.beginPath(); us.forEach((u, k) => { const X = Xp(-180 + 360 * k / (N - 1)); k ? cx.lineTo(X, Yp(u)) : cx.moveTo(X, Yp(u)); }); cx.stroke();
    trail.forEach((p, k) => { cx.fillStyle = `rgba(255,107,107,${0.1 + 0.5 * k / trail.length})`; cx.beginPath(); cx.arc(Xp(p * 180 / Math.PI), Yp(U(p)), 3, 0, 7); cx.fill(); });
    cx.fillStyle = COL.coral; cx.beginPath(); cx.arc(Xp(phi * 180 / Math.PI), Yp(U(phi)), 7, 0, 7); cx.fill();
    // stash for readouts
    const frac = (a, b) => ps.reduce((s, p, k) => { const d = Math.abs(-180 + 360 * k / (N - 1)); return s + (d > a && d <= b ? p : 0); }, 0) / Z;
    stats = { pt: frac(120, 181), pg: frac(30, 120), pc: frac(-1, 30), mn };
  }
  let stats = { pt: 0, pg: 0, pc: 0, mn: 0 }, minU = 0;
  function readouts() {
    [1, 2, 3].forEach(n => o['k' + n].textContent = st.k[n - 1].toFixed(1) + ' kJ/mol' + (n === 3 && (st.k[3] || st.k[4]) ? `  (+ k₄ = ${st.k[3]}, k₅ = ${st.k[4]})` : '')); o.T.textContent = st.T + ' K';
    minU = stats.mn;
    o.pt.textContent = (100 * stats.pt).toFixed(1) + ' %'; o.pg.textContent = (100 * stats.pg).toFixed(1) + ' %'; o.pc.textContent = (100 * stats.pc).toFixed(1) + ' %';
    let barr = -1e9; for (let d = 60; d <= 180; d += 1) barr = Math.max(barr, U(d * Math.PI / 180)); o.bar.textContent = (barr - U(Math.PI)).toFixed(1) + ' kJ/mol';
    o.rate.textContent = t < 0.05 ? 'measuring…' : (hops / t * 1000).toFixed(0) + `  (${(t).toFixed(1)} ps simulated)`;
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running && !drag) step(); draw(); if (frame++ % 8 === 0) readouts(); };
  resetRate(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
