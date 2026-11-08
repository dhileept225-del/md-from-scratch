// Lab 6.7 widget, two modes.
// MSD: one 2D random walker (D = 1); MSD from the first time origin vs averaged over all origins vs the exact 4Dt.
// Haven: 50 cations + 50 anions in 1D Brownian motion (D = 1); a chosen fraction is bound in neutral pairs.
//        Self MSD vs collective (charge) MSD per ion; Haven ratio sigma_NE/sigma from each run, collected over runs.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };

export function mountLab6(root) {
  const st = { mode: 'msd', n: 400, f: 0.34 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view an-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="msd" class="on">MSD: time origins</button><button data-m="hav">Haven ratio</button></div>
    <div data-p="msd"><label><span>Trajectory length (frames)</span><output data-o="n"></output><input type="range" data-k="n" min="50" max="2000" step="50" value="400"></label></div>
    <div data-p="hav" hidden><label><span>Fraction of ions in neutral pairs</span><output data-o="f"></output><input type="range" data-k="f" min="0" max="0.9" step="0.02" value="0.34"></label></div>
    <div class="row"><button class="btn sm primary" data-k="new">New run</button><button class="btn sm ghost" data-k="ten">10 runs</button><button class="btn sm ghost" data-k="clear">Clear</button></div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.an-view'), $ = s => root.querySelector(s);
  let rng = 3; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  const H = 0.01, STEPS = 2000, LAG = 100, NP = 50;
  let walk = null, msdAll = [], msd0 = [], hv = null, runs = [];
  function newWalk() {
    walk = [[0, 0]]; for (let k = 1; k < st.n; k++) { const p = walk[k - 1]; walk.push([p[0] + gauss() * Math.sqrt(2 * H), p[1] + gauss() * Math.sqrt(2 * H)]); }
    const maxlag = Math.floor(st.n * 0.9); msdAll = []; msd0 = [];
    for (let l = 1; l <= maxlag; l++) { let s = 0; for (let k = 0; k + l < st.n; k++) s += (walk[k + l][0] - walk[k][0]) ** 2 + (walk[k + l][1] - walk[k][1]) ** 2;
      msdAll.push(s / (st.n - l)); msd0.push((walk[l][0] - walk[0][0]) ** 2 + (walk[l][1] - walk[0][1]) ** 2); }
  }
  function newHaven() {
    const nb = Math.round(st.f * NP), X = Array.from({ length: 2 * NP }, () => new Float64Array(STEPS)), q = i => i < NP ? 1 : -1;
    for (let t = 1; t < STEPS; t++) for (let i = 0; i < 2 * NP; i++) {
      if (i >= NP && i < NP + nb) { X[i][t] = X[i][t - 1] + (X[i - NP][t] - X[i - NP][t - 1]); continue; }   // bound anion copies its cation
      X[i][t] = X[i][t - 1] + gauss() * Math.sqrt(2 * H);
    }
    const lags = [5, 10, 20, 40, 60, 80, 100, 150, 200], selfC = [], colC = [];
    for (const l of lags) { let s = 0, c = 0, n = 0;
      for (let t = 0; t + l < STEPS; t += 5) { let M = 0; for (let i = 0; i < 2 * NP; i++) { const d = X[i][t + l] - X[i][t]; s += d * d; M += q(i) * d; } c += M * M; n++; }
      selfC.push(s / n / (2 * NP)); colC.push(c / n / (2 * NP)); }
    const k = lags.indexOf(LAG), Ds = selfC[k] / (2 * LAG * H), Dsig = colC[k] / (2 * LAG * H);
    hv = { lags, selfC, colC, Ds, Dsig, Hr: Ds / Dsig, nb }; runs.push(hv.Hr);
  }
  const fresh = () => st.mode === 'msd' ? newWalk() : newHaven();
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); $('[data-k=ten]').hidden = st.mode === 'msd'; $('[data-k=clear]').hidden = st.mode === 'msd';
    if (st.mode === 'hav' && !hv) newHaven(); if (st.mode === 'msd' && !walk) newWalk(); draw(); });
  $('[data-k=n]').oninput = e => { st.n = +e.target.value; newWalk(); draw(); };
  $('[data-k=f]').oninput = e => { st.f = +e.target.value; runs = []; newHaven(); draw(); };
  $('[data-k=new]').onclick = () => { fresh(); draw(); };
  $('[data-k=ten]').onclick = () => { for (let k = 0; k < 10; k++) newHaven(); draw(); };
  $('[data-k=clear]').onclick = () => { runs = []; newHaven(); draw(); };
  $('[data-k=ten]').hidden = true; $('[data-k=clear]').hidden = true;

  let W = 600, Hh = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; Hh = view.clientHeight; cv.width = W * d; cv.height = Hh * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x, y); };
  const line = (pts, c, w = 2) => { cx.strokeStyle = c; cx.lineWidth = w; cx.beginPath(); pts.forEach((p, k) => k ? cx.lineTo(p[0], p[1]) : cx.moveTo(p[0], p[1])); cx.stroke(); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, Hh);
    const half = Hh * 0.42, R = $('[data-ro]');
    if (st.mode === 'msd') {
      // trajectory
      const bw = Math.min(W * 0.42, half - 20), xs = walk.map(p => p[0]), ys = walk.map(p => p[1]);
      const mn = [Math.min(...xs), Math.min(...ys)], span = Math.max(Math.max(...xs) - mn[0], Math.max(...ys) - mn[1]) || 1, S = (bw - 16) / span;
      cx.strokeStyle = 'rgba(143,163,199,.3)'; cx.strokeRect(10, 10, bw, bw);
      line(walk.map(p => [18 + (p[0] - mn[0]) * S, 10 + bw - 8 - (p[1] - mn[1]) * S]), 'rgba(56,225,255,.8)', 1.2);
      label(`2D walk, ${st.n} frames`, 12, bw + 26, COL.muted, 'left', '11px');
      // MSD plot
      const x0 = bw + 56, x1 = W - 12, y0 = 24, y1 = Hh - 30, nl = msdAll.length, ymax = Math.max(4 * nl * H, ...msdAll, ...msd0) * 1.05;
      const X = l => x0 + l / nl * (x1 - x0), Y = m => y1 - m / ymax * (y1 - y0);
      axes(x0, y0, x1, y1); label('MSD vs time lag', x0 + 6, y0 - 8);
      line([[X(0), Y(0)], [X(nl), Y(4 * nl * H)]], COL.amber, 2);
      line(msd0.map((m, k) => [X(k + 1), Y(m)]), 'rgba(255,107,107,.8)', 1.2);
      line(msdAll.map((m, k) => [X(k + 1), Y(m)]), COL.cyan, 2.2);
      label('amber: exact 4Dt   cyan: all origins   coral: origin t = 0 only', x0 + 6, y1 + 18, COL.muted, 'left', '11px');
      const fit = (arr, a, b) => { let sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0; for (let l = a; l < b; l++) { const xx = (l + 1) * H; sx += xx; sy += arr[l]; sxx += xx * xx; sxy += xx * arr[l]; n++; } return (n * sxy - sx * sy) / (n * sxx - sx * sx) / 4; };
      const a = Math.floor(nl * 0.05), b = Math.max(a + 3, Math.floor(nl * 0.3));
      R.innerHTML = `<div><span>D, all origins (first 30 % of lags)</span><b>${fit(msdAll, a, b).toFixed(3)}</b></div><div><span>D, origin t = 0 only</span><b>${fit(msd0, a, b).toFixed(3)}</b></div>
        <div><span>Exact D</span><b>1.000</b></div><div><span>Samples at the longest lag</span><b>${st.n - nl}</b></div>`;
      $('[data-o=n]').textContent = st.n;
    } else {
      const x0 = 50, x1 = W - 12, y0 = 24, y1 = half - 6, lmax = 200, ymax = 4 * lmax * H * 1.15;
      const X = l => x0 + l / lmax * (x1 - x0), Y = m => y1 - Math.min(m, ymax) / ymax * (y1 - y0);
      axes(x0, y0, x1, y1); label('per-ion MSD vs lag: self (cyan), charge / collective (coral)', x0 + 6, y0 - 8, COL.ink, 'left', W < 500 ? '600 10px' : '600 12px');
      line([[X(0), Y(0)], [X(lmax), Y(2 * lmax * H)]], 'rgba(255,194,75,.7)', 1.5);
      line([[X(0), Y(0)], ...hv.lags.map((l, k) => [X(l), Y(hv.selfC[k])])], COL.cyan, 2.2);
      line([[X(0), Y(0)], ...hv.lags.map((l, k) => [X(l), Y(hv.colC[k])])], COL.coral, 2.2);
      cx.setLineDash([3, 4]); cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.beginPath(); cx.moveTo(X(LAG), y0); cx.lineTo(X(LAG), y1); cx.stroke(); cx.setLineDash([]);
      // histogram of Haven ratios
      const hx0 = 46, hx1 = W - 12, hy0 = half + 46, hy1 = Hh - 28, lo = 0.5, hi = 4.5, nb = 32, bins = new Array(nb).fill(0);
      runs.forEach(h => { const b = Math.floor((Math.min(Math.max(h, lo), hi - 1e-9) - lo) / (hi - lo) * nb); bins[b]++; });
      const bmax = Math.max(1, ...bins), HX = v => hx0 + (v - lo) / (hi - lo) * (hx1 - hx0), HY = c => hy1 - c / bmax * (hy1 - hy0);
      axes(hx0, hy0, hx1, hy1); label(`Haven ratio from each run (${runs.length} runs); dashed: exact`, hx0 + 6, hy0 - 10);
      cx.fillStyle = 'rgba(56,225,255,.55)'; bins.forEach((c, k) => { if (c) cx.fillRect(HX(lo + k * (hi - lo) / nb) + 1, HY(c), (hx1 - hx0) / nb - 2, hy1 - HY(c)); });
      const ex = 1 / (1 - 2 * hv.nb / (2 * NP)); if (ex < hi) { cx.setLineDash([4, 4]); cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(HX(ex), hy0); cx.lineTo(HX(ex), hy1); cx.stroke(); cx.setLineDash([]); }
      cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; [1, 2, 3, 4].forEach(v => cx.fillText(v, HX(v), hy1 + 14));
      const m = runs.reduce((a, b) => a + b, 0) / runs.length, sd = Math.sqrt(runs.reduce((a, b) => a + (b - m) ** 2, 0) / runs.length);
      R.innerHTML = `<div><span>D_self (exact 1)</span><b>${hv.Ds.toFixed(3)}</b></div><div><span>D from charge, per ion (exact ${(1 - 2 * hv.nb / (2 * NP)).toFixed(2)})</span><b>${hv.Dsig.toFixed(3)}</b></div>
        <div><span>Haven, this run (exact ${isFinite(ex) ? ex.toFixed(2) : '∞'})</span><b>${hv.Hr.toFixed(2)}</b></div><div><span>Haven over runs</span><b>${runs.length > 1 ? `${m.toFixed(2)} ± ${sd.toFixed(2)}` : '…'}</b></div>`;
      $('[data-o=f]').textContent = `${(2 * Math.round(st.f * NP) / (2 * NP)).toFixed(2)} (${Math.round(st.f * NP)} pairs)`;
    }
  }
  newWalk(); draw();
  return () => { ro.disconnect(); };
}
