// Lifetime lab widget: model pair distances r(t) (overdamped Langevin in a well at r = 1.12, kT = 1) for 150 pairs,
// h(t) = r < cutoff; continuous and intermittent correlation functions with a chosen saving interval and tolerance t*.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const NP = 150, TT = 60, DTS = 0.002, D = 0.14, SAVE0 = 0.01;
const W_ = r => -1.4 * Math.exp(-(((r - 1.12) / 0.2) ** 2)), dW = r => 1.4 * 2 * (r - 1.12) / 0.04 * Math.exp(-(((r - 1.12) / 0.2) ** 2));

export function mountLab12(root) {
  const st = { rc: 1.45, save: 0.02, tstar: 0 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view lt-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Cutoff</span><output data-o="rc"></output><input type="range" data-k="rc" min="1.25" max="1.7" step="0.01" value="1.45"></label>
    <label><span>Saving interval</span><output data-o="save"></output><input type="range" data-k="save" min="0.01" max="0.5" step="0.01" value="0.02"></label>
    <label><span>Tolerance t*</span><output data-o="tstar"></output><input type="range" data-k="tstar" min="0" max="0.4" step="0.01" value="0"></label>
    <div class="row"><button class="btn sm ghost" data-k="new">New trajectories</button></div>
    <div class="readouts">
      <div><span>Continuous lifetime</span><b data-o="tc"></b></div>
      <div><span>Intermittent lifetime (baseline removed)</span><b data-o="ti"></b></div>
      <div><span>Mean bonded stretch</span><b data-o="ms"></b></div>
      <div><span>Fraction of time bonded</span><b data-o="hb"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.lt-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let seed = 99; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };
  const gauss = () => Math.sqrt(-2 * Math.log(rand())) * Math.cos(2 * Math.PI * rand());
  let R = [], res = null;
  function simulate() {   // overdamped: dr = -D dW/dr dt + sqrt(2 D dt) xi ; reflecting walls at 0.95 and 4.5
    R = []; const nsave = Math.round(TT / SAVE0), per = Math.round(SAVE0 / DTS), s2 = Math.sqrt(2 * D * DTS);
    for (let p = 0; p < NP; p++) { const a = new Float32Array(nsave); let r = 1.12 + 0.1 * gauss();
      for (let k = 0; k < nsave; k++) { for (let q = 0; q < per; q++) { r += -D * dW(r) * DTS + s2 * gauss(); if (r < 0.95) r = 1.9 - r; if (r > 4.5) r = 9 - r; } a[k] = r; }
      R.push(a); }
  }
  function analyse() {
    const sub = Math.max(1, Math.round(st.save / SAVE0)), dt = sub * SAVE0, ts = Math.round(st.tstar / dt), nf = Math.floor(R[0].length / sub), kmax = Math.min(nf - 1, Math.round(15 / dt));
    let lens = [], bonded = 0, total = 0; const num = new Float64Array(kmax), den = new Float64Array(kmax), H = [];
    R.forEach(a => {
      const h = new Uint8Array(nf); for (let k = 0; k < nf; k++) h[k] = a[k * sub] < st.rc ? 1 : 0;
      if (ts) { let i = 0; while (i < nf) { if (!h[i]) { let j = i; while (j < nf && !h[j]) j++; if (i > 0 && j < nf && j - i <= ts) for (let q = i; q < j; q++) h[q] = 1; i = j; } else i++; } }
      H.push(h); let L = 0; for (let k = 0; k <= nf; k++) { if (k < nf && h[k]) L++; else if (L) { lens.push(L); L = 0; } }
      for (let k = 0; k < nf; k++) { bonded += h[k]; total++; }
      for (let lag = 0; lag < kmax; lag++) { let s = 0, d = 0; for (let k = 0; k + lag < nf; k++) { if (h[k]) { d++; if (h[k + lag]) s++; } } num[lag] += s; den[lag] += d; }
    });
    const avg = bonded / total, sumL = lens.reduce((a, b) => a + b, 0);
    const cc = Array.from({ length: kmax }, (_, k) => lens.reduce((s, L) => s + Math.max(0, L - k), 0) / sumL);
    const ci = Array.from({ length: kmax }, (_, k) => num[k] / den[k]), cin = ci.map(v => (v - avg) / (1 - avg));
    let kz = cin.findIndex(v => v < 0.01); if (kz < 0) kz = kmax;
    res = { dt, cc, ci: cin, tc: dt * (cc.reduce((a, b) => a + b, 0) - 0.5), ti: dt * (cin.slice(0, kz).reduce((a, b) => a + b, 0) - 0.5), ms: dt * sumL / lens.length, avg, H0: H[0] };
  }
  ['rc', 'save', 'tstar'].forEach(k => $(`[data-k=${k}]`).oninput = e => { st[k] = +e.target.value; analyse(); draw(); });
  $('[data-k=new]').onclick = () => { simulate(); analyse(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  function draw() {
    if (!W || !res) return;
    cx.clearRect(0, 0, W, H);
    // r(t) of pair 1 for the first 10 time units, with the cutoff and the bonded strip
    const x0 = 46, x1 = W - 12, y0 = 24, y1 = H * 0.36, tshow = 10, X = t => x0 + t / tshow * (x1 - x0), Y = r => y1 - (Math.min(3.0, r) - 0.9) / 2.1 * (y1 - y0), a = R[0];
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(W < 500 ? 'pair 1: r(t); lit = bonded' : 'pair 1: distance r(t); lit strip = counted as bonded', x0 + 6, y0 - 8);
    const nshow = Math.round(tshow / res.dt); for (let k = 0; k < nshow && k < res.H0.length; k++) if (res.H0[k]) { cx.fillStyle = 'rgba(56,225,255,.22)'; cx.fillRect(X(k * res.dt), y0, X(res.dt) - X(0) + 0.5, y1 - y0); }
    cx.strokeStyle = COL.ink; cx.lineWidth = 1; cx.beginPath(); for (let k = 0; k <= tshow / SAVE0; k++) k ? cx.lineTo(X(k * SAVE0), Y(a[k])) : cx.moveTo(X(0), Y(a[0])); cx.stroke();
    cx.setLineDash([5, 4]); cx.strokeStyle = COL.coral; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(x0, Y(st.rc)); cx.lineTo(x1, Y(st.rc)); cx.stroke(); cx.setLineDash([]);
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'right'; [1, 1.5, 2, 3].forEach(r => cx.fillText(r, x0 - 4, Y(r) + 4)); cx.textAlign = 'center'; [0, 5, 10].forEach(t => cx.fillText(t, X(t), y1 + 13));
    // correlation functions
    const bx0 = 46, bx1 = W - 12, by0 = H * 0.36 + 50, by1 = H - 28, tmax = 10, BX = t => bx0 + t / tmax * (bx1 - bx0), BY = c => by1 - Math.max(0, Math.min(1, c)) * (by1 - by0);
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(bx0, by0); cx.lineTo(bx0, by1); cx.lineTo(bx1, by1); cx.stroke();
    cx.fillStyle = COL.ink; cx.font = '600 12px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText(W < 500 ? 'continuous (cyan), intermittent (amber)' : 'continuous (cyan) and intermittent (amber) correlation, 150 pairs', bx0 + 6, by0 - 10);
    [[res.cc, COL.cyan], [res.ci, COL.amber]].forEach(([arr, c]) => { cx.strokeStyle = c; cx.lineWidth = 2; cx.beginPath(); arr.forEach((v, k) => { const t = k * res.dt; if (t > tmax) return; k ? cx.lineTo(BX(t), BY(v)) : cx.moveTo(BX(t), BY(v)); }); cx.stroke(); });
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [0, 2, 4, 6, 8, 10].forEach(t => cx.fillText(t, BX(t), by1 + 14)); cx.textAlign = 'right'; [0, 0.5, 1].forEach(c => cx.fillText(c, bx0 - 4, BY(c) + 4));
    o.rc.textContent = st.rc.toFixed(2); o.save.textContent = res.dt.toFixed(2); o.tstar.textContent = st.tstar.toFixed(2);
    o.tc.textContent = res.tc.toFixed(2); o.ti.textContent = res.ti.toFixed(2); o.ms.textContent = res.ms.toFixed(2); o.hb.textContent = res.avg.toFixed(3);
  }
  simulate(); analyse(); draw();
  return () => ro.disconnect();
}
