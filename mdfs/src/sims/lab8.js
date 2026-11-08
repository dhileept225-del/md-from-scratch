// Error-bar lab widget: correlated series x(t+1) = phi x(t) + noise (variance 1), block averaging vs the exact error,
// and the scatter of means over independent runs.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };

export function mountLab8(root) {
  const st = { phi: 0.9, n: 20000 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view er-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Correlation φ (memory of the previous value)</span><output data-o="phi"></output><input type="range" data-k="phi" min="0" max="0.99" step="0.01" value="0.9"></label>
    <label><span>Number of samples N</span><output data-o="n"></output><input type="range" data-k="n" min="1000" max="50000" step="1000" value="20000"></label>
    <div class="row"><button class="btn sm primary" data-k="new">New run</button><button class="btn sm ghost" data-k="many">20 runs</button></div>
    <div class="readouts">
      <div><span>Mean of this run (true 0)</span><b data-o="m"></b></div>
      <div><span>Naive error</span><b data-o="nv"></b></div>
      <div><span>Block plateau / exact error</span><b data-o="bp"></b></div>
      <div><span>Inefficiency · independent samples</span><b data-o="g"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.er-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  let rng = 21; const rand = () => (rng = (rng * 16807) % 2147483647) / 2147483647;
  const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
  let x, blocks, runMeans = [];
  function series() { const a = new Float64Array(st.n), e = Math.sqrt(1 - st.phi ** 2); a[0] = gauss(); for (let i = 1; i < st.n; i++) a[i] = st.phi * a[i - 1] + e * gauss(); return a; }
  const mean = a => { let s = 0; for (const v of a) s += v; return s / a.length; };
  function analyse() {
    blocks = []; for (let b = 1; st.n / b >= 8; b *= 2) { const nb = Math.floor(st.n / b), mm = []; for (let k = 0; k < nb; k++) { let s = 0; for (let i = 0; i < b; i++) s += x[k * b + i]; mm.push(s / b); }
      const m = mean(mm), v = mm.reduce((s, q) => s + (q - m) ** 2, 0) / (nb - 1); blocks.push([b, Math.sqrt(v / nb), nb]); }
  }
  function newRun() { x = series(); analyse(); }
  $('[data-k=phi]').oninput = e => { st.phi = +e.target.value; runMeans = []; newRun(); draw(); };
  $('[data-k=n]').oninput = e => { st.n = +e.target.value; runMeans = []; newRun(); draw(); };
  $('[data-k=new]').onclick = () => { newRun(); runMeans.push(mean(x)); draw(); };
  $('[data-k=many]').onclick = () => { for (let k = 0; k < 20; k++) { runMeans.push(mean(series())); } newRun(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x, y); };
  function draw() {
    if (!W || !x) return;
    cx.clearRect(0, 0, W, H);
    const half = H * 0.36, g = (1 + st.phi) / (1 - st.phi), exact = Math.sqrt(g / st.n), naive = blocks[0][1];
    // series (first 2000 points)
    const x0 = 46, x1 = W - 12, y0 = 22, y1 = half, nshow = Math.min(2000, st.n), X = i => x0 + i / (nshow - 1) * (x1 - x0), Y = v => (y0 + y1) / 2 - Math.max(-3.5, Math.min(3.5, v)) / 3.5 * (y1 - y0) / 2;
    axes(x0, y0, x1, y1); label(`first ${nshow} samples`, x0 + 6, y0 - 6);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 1; cx.beginPath(); for (let i = 0; i < nshow; i++) i ? cx.lineTo(X(i), Y(x[i])) : cx.moveTo(X(i), Y(x[i])); cx.stroke();
    // block curve (log2 x)
    const bx0 = 46, bx1 = W - 12, by0 = half + 46, by1 = H - 30, nb = blocks.length, BX = k => bx0 + k / Math.max(1, nb - 1) * (bx1 - bx0);
    const emax = Math.max(exact, ...blocks.map(b => b[1]), runMeans.length > 1 ? Math.sqrt(runMeans.reduce((s, q) => s + q * q, 0) / runMeans.length) : 0) * 1.25, BY = e => by1 - e / emax * (by1 - by0);
    axes(bx0, by0, bx1, by1); label('error of the mean vs block size', bx0 + 6, by0 - 10);
    cx.setLineDash([5, 4]); cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); cx.moveTo(bx0, BY(exact)); cx.lineTo(bx1, BY(exact)); cx.stroke();
    cx.strokeStyle = 'rgba(255,107,107,.7)'; cx.beginPath(); cx.moveTo(bx0, BY(naive)); cx.lineTo(bx1, BY(naive)); cx.stroke(); cx.setLineDash([]);
    cx.strokeStyle = COL.cyan; cx.lineWidth = 2; cx.beginPath(); blocks.forEach((b, k) => k ? cx.lineTo(BX(k), BY(b[1])) : cx.moveTo(BX(k), BY(b[1]))); cx.stroke();
    blocks.forEach((b, k) => { cx.fillStyle = b[2] < 20 ? 'rgba(56,225,255,.35)' : COL.cyan; cx.beginPath(); cx.arc(BX(k), BY(b[1]), 3.5, 0, 7); cx.fill(); });
    if (runMeans.length > 1) { const sd = Math.sqrt(runMeans.reduce((s, q) => s + q * q, 0) / runMeans.length); cx.fillStyle = 'rgba(56,225,255,.25)'; cx.fillRect(bx1 - 18, BY(sd), 14, by1 - BY(sd)); label(`runs: ${sd.toFixed(4)}`, bx1 - 22, BY(sd) - 4, COL.cyan, 'right', '11px'); }
    cx.fillStyle = COL.muted; cx.font = '11px "Space Grotesk", sans-serif'; cx.textAlign = 'center'; blocks.forEach((b, k) => { if (k % 2 === 0) cx.fillText(b[0], BX(k), by1 + 14); });
    label('amber: exact   coral: naive   faded dots: fewer than 20 blocks', bx0 + 6, by1 + 26, COL.muted, 'left', '10px');
    const plat = Math.max(...blocks.filter(b => b[2] >= 20).map(b => b[1]));
    o.phi.textContent = st.phi.toFixed(2); o.n.textContent = st.n.toLocaleString();
    o.m.textContent = mean(x).toFixed(4); o.nv.textContent = naive.toFixed(4); o.bp.textContent = `${plat.toFixed(4)} / ${exact.toFixed(4)}`;
    o.g.textContent = `${g.toFixed(1)} · ${Math.round(st.n / g).toLocaleString()}`;
  }
  newRun(); draw();
  return () => ro.disconnect();
}
