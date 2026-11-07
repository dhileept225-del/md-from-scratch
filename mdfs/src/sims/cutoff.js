// Cutoff lab: truncation schemes for LJ, jumps at rc, tail correction and cost, measured on a 500-atom configuration.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7', green: '#7EE787' };
const LJ = r => 4 * (r ** -12 - r ** -6), LJF = r => 24 * (2 * r ** -13 - r ** -7);   // F = -dV/dr
const SCHEMES = { trunc: 'Truncate', shift: 'Shift', fshift: 'Force-shift', switch: 'Switch' };

// a fixed, slightly disordered fcc configuration: 500 atoms at rho* = 0.8442
const CONF = (() => {
  let seed = 7; const R = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  const n = 5, N = 4 * n ** 3, rho = 0.8442, L = Math.cbrt(N / rho), a = L / n, P = [];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]].forEach(b => P.push([(i + b[0]) * a + (R() - .5) * .28, (j + b[1]) * a + (R() - .5) * .28, (k + b[2]) * a + (R() - .5) * .28]));
  const r = [];
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { let dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = P[i][2] - P[j][2]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L); const d = Math.hypot(dx, dy, dz); if (d < L / 2) r.push(d); }
  return { N, rho, L, r: Float64Array.from(r) };
})();
const tail = rc => 8 / 3 * Math.PI * CONF.rho * (rc ** -9 / 3 - rc ** -3);   // energy per atom beyond rc (g = 1)
const REF = (() => { const rc = CONF.L / 2 - 0.01; let e = 0; for (const d of CONF.r) if (d < rc) e += LJ(d); return e / CONF.N + tail(rc); })();

export function mountCutoff(root) {
  const st = { rc: 2.5, scheme: 'shift', tail: false, zoom: true };
  const ron = () => Math.max(1.0, st.rc - 0.4);
  function V(r) {
    const rc = st.rc; if (r >= rc) return 0;
    switch (st.scheme) {
      case 'trunc': return LJ(r);
      case 'shift': return LJ(r) - LJ(rc);
      case 'fshift': return LJ(r) - LJ(rc) + (r - rc) * LJF(rc);
      case 'switch': { const a = ron(); if (r <= a) return LJ(r); const S = ((rc * rc - r * r) ** 2 * (rc * rc + 2 * r * r - 3 * a * a)) / ((rc * rc - a * a) ** 3); return LJ(r) * S; }
    }
  }
  const F = r => { const h = 1e-5; return -(V(r + h) - V(r - h)) / (2 * h); };

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view cut-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg">${Object.entries(SCHEMES).map(([k, n]) => `<button data-s="${k}" class="${k === 'shift' ? 'on' : ''}">${n}</button>`).join('')}</div>
    <label>Cutoff rc <output data-o="rc"></output><input type="range" data-k="rc" min="1.2" max="4.1" step="0.05" value="2.5"></label>
    <label class="check"><input type="checkbox" data-k="zoom" checked> Zoom in near the cutoff</label>
    <label class="check"><input type="checkbox" data-k="tail"> Add tail correction to the energy</label>
    <div class="readouts">
      <div><span>Energy jump at rc</span><b data-o="vj"></b></div>
      <div><span>Force jump at rc</span><b data-o="fj"></b></div>
      <div><span>Neighbours per atom</span><b data-o="nb"></b></div>
      <div><span>Relative cost (rc = 2.5 is 1×)</span><b data-o="cost"></b></div>
      <div class="wide"><span>Energy per atom vs full LJ (500 atoms)</span><b data-o="E"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.cut-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { st.scheme = b.dataset.s; root.querySelectorAll('[data-s]').forEach(x => x.classList.toggle('on', x === b)); update(); });
  $('[data-k=rc]').oninput = e => { st.rc = +e.target.value; update(); };
  $('[data-k=zoom]').onchange = e => { st.zoom = e.target.checked; update(); };
  $('[data-k=tail]').onchange = e => { st.tail = e.target.checked; update(); };

  let W = 600, H = 460;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); update(); };
  const ro = new ResizeObserver(resize); ro.observe(view);

  function plot(y0, h, full, cut, label, color, yr) {
    const x0 = 50, x1 = W - 16, r0 = 1.0, r1 = Math.max(3.0, st.rc + 0.8);
    const X = r => x0 + (r - r0) / (r1 - r0) * (x1 - x0), Y = v => y0 + h / 2 - v / yr * (h / 2);
    cx.fillStyle = color; cx.font = '600 13px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillRect(x0, y0 - 12, 14, 4); cx.fillText(label, x0 + 20, y0 - 6);
    cx.fillStyle = 'rgba(242,245,250,.45)'; cx.fillRect(x0 + 150, y0 - 12, 14, 2); cx.fillText('full LJ', x0 + 170, y0 - 6);
    cx.strokeStyle = 'rgba(143,163,199,.45)'; cx.lineWidth = 1; cx.strokeRect(x0, y0, x1 - x0, h);
    cx.beginPath(); cx.moveTo(x0, Y(0)); cx.lineTo(x1, Y(0)); cx.stroke();
    cx.save(); cx.beginPath(); cx.rect(x0, y0, x1 - x0, h); cx.clip();
    [[full, 'rgba(242,245,250,.4)', 1.5], [cut, color, 3]].forEach(([f, c, w]) => {
      cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); let prev = null;
      for (let k = 0; k <= 500; k++) { const r = r0 + k / 500 * (r1 - r0), v = f(r), X1 = X(r), Y1 = Y(v); if (prev === null || (prev < st.rc) !== (r < st.rc)) cx.moveTo(X1, Y1); else cx.lineTo(X1, Y1); prev = r; }
      cx.stroke();
    });
    cx.strokeStyle = COL.amber; cx.setLineDash([5, 5]); cx.beginPath(); cx.moveTo(X(st.rc), y0); cx.lineTo(X(st.rc), y0 + h); cx.stroke(); cx.setLineDash([]);
    cx.restore();
    cx.fillStyle = COL.muted; cx.font = '500 12px "Space Grotesk", sans-serif'; cx.textAlign = 'center';
    for (let r = 1; r <= r1 + 1e-9; r += 0.5) cx.fillText(r.toFixed(1), X(r), y0 + h + 14);
    cx.fillStyle = COL.amber; cx.fillText('rc', X(st.rc), y0 + 12);
  }
  function update() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const vr = st.zoom ? Math.max(0.03, Math.abs(LJ(st.rc)) * 5) : 1.2, fr = st.zoom ? Math.max(0.06, Math.abs(LJF(st.rc)) * 5) : 2.5;
    const h = (H - 90) / 2;
    plot(24, h, LJ, V, 'V(r) with ' + SCHEMES[st.scheme].toLowerCase(), COL.cyan, vr);
    plot(24 + h + 46, h, LJF, F, 'force F(r)', COL.coral, fr);
    cx.fillStyle = COL.muted; cx.textAlign = 'right'; cx.fillText('r / σ', W - 16, H - 4);
    // numbers from the configuration
    let e = 0, n = 0; for (const d of CONF.r) if (d < st.rc) { e += V(d); n++; }
    let E = e / CONF.N;
    if (st.tail) E += tail(st.rc) + (st.scheme === 'shift' ? n * LJ(st.rc) / CONF.N : 0);
    const nb = 2 * n / CONF.N, cost = (st.rc / 2.5) ** 3, err = E - REF;
    const vj = st.scheme === 'trunc' ? LJ(st.rc) : 0, fj = (st.scheme === 'trunc' || st.scheme === 'shift') ? LJF(st.rc) : 0;
    o.rc.textContent = st.rc.toFixed(2) + ' σ' + (st.rc > CONF.L / 2 ? ' (> L/2 of this box)' : '');
    o.vj.innerHTML = Math.abs(vj) > 1e-9 ? `<span style="color:${COL.coral}">${vj.toFixed(4)} ε</span>` : `<span style="color:${COL.cyan}">0 (smooth)</span>`;
    o.fj.innerHTML = Math.abs(fj) > 1e-9 ? `<span style="color:${COL.amber}">${fj.toFixed(4)}</span>` : `<span style="color:${COL.cyan}">0 (smooth)</span>`;
    o.nb.textContent = nb.toFixed(1); o.cost.textContent = '× ' + cost.toFixed(2);
    const note = st.scheme === 'shift' && !st.tail ? ' (shifted energies are offset; tick the tail box to correct)' : st.scheme === 'fshift' || st.scheme === 'switch' ? ' (modified potential: energies differ by design)' : '';
    o.E.innerHTML = `${E.toFixed(3)} vs ${REF.toFixed(3)}: <span style="color:${Math.abs(err) < 0.02 ? COL.cyan : Math.abs(err) < 0.2 ? COL.amber : COL.coral}">error ${err >= 0 ? '+' : ''}${err.toFixed(3)} ε</span>${note}`;
  }
  resize();
  return () => { ro.disconnect(); };
}
