// Neighbour search: all pairs vs cell list vs Verlet list (built with cells), with live operation counts and timing.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const RC = 2.5, RHO = 0.5, STEP = 0.025;

export function mountNeighbours(root) {
  const st = { N: 800, method: 'verlet', skin: 0.4, running: true, lazy: false };
  let L, X, ref, list = [], since = 0, hist = [], rebuilds = [], ok = '–', sel = 0, lastBuild = 0;
  function reset() {
    L = Math.sqrt(st.N / RHO); X = new Float64Array(2 * st.N);
    const n = Math.ceil(Math.sqrt(st.N)), a = L / n;
    for (let i = 0; i < st.N; i++) { X[2 * i] = ((i % n) + 0.5 + (Math.random() - 0.5) * 0.6) * a; X[2 * i + 1] = (Math.floor(i / n) + 0.5 + (Math.random() - 0.5) * 0.6) * a; }
    let bd = 1e9; for (let i = 0; i < st.N; i++) { const d = Math.hypot(X[2 * i] - L / 2, X[2 * i + 1] - L / 2); if (d < bd) { bd = d; sel = i; } }
    ref = null; hist = []; rebuilds = []; since = 0; ok = '–';
  }
  const mi = d => d - L * Math.round(d / L);

  // cell list: returns pairs within r (and the number of distance tests)
  function cellPairs(r) {
    const nc = Math.max(3, Math.floor(L / r)), cs = L / nc, head = new Int32Array(nc * nc).fill(-1), next = new Int32Array(st.N), out = []; let tests = 0;
    for (let i = 0; i < st.N; i++) { const cx = Math.floor(X[2 * i] / cs) % nc, cy = Math.floor(X[2 * i + 1] / cs) % nc, c = cy * nc + cx; next[i] = head[c]; head[c] = i; }
    const r2 = r * r;
    for (let cy = 0; cy < nc; cy++) for (let cx = 0; cx < nc; cx++) {
      const c = cy * nc + cx;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const c2 = ((cy + dy + nc) % nc) * nc + ((cx + dx + nc) % nc);
        for (let i = head[c]; i >= 0; i = next[i]) for (let j = head[c2]; j >= 0; j = next[j]) {
          if (j <= i) continue; tests++;
          const ddx = mi(X[2 * i] - X[2 * j]), ddy = mi(X[2 * i + 1] - X[2 * j + 1]);
          if (ddx * ddx + ddy * ddy < r2) out.push(i, j);
        }
      }
    }
    return { pairs: out, tests, nc };
  }
  function allPairs(r) {
    const r2 = r * r; let n = 0;
    for (let i = 0; i < st.N; i++) for (let j = i + 1; j < st.N; j++) { const dx = mi(X[2 * i] - X[2 * j]), dy = mi(X[2 * i + 1] - X[2 * j + 1]); if (dx * dx + dy * dy < r2) n++; }
    return { count: n, tests: st.N * (st.N - 1) / 2 };
  }
  function stepOnce() {
    for (let i = 0; i < 2 * st.N; i++) { X[i] += STEP * (Math.random() * 2 - 1); X[i] -= L * Math.floor(X[i] / L); }
    const t0 = performance.now(); let tests = 0, found = 0;
    if (st.method === 'all') { const r = allPairs(RC); tests = r.tests; found = r.count; }
    else if (st.method === 'cell') { const r = cellPairs(RC); tests = r.tests; found = r.pairs.length / 2; }
    else {
      let maxd = 0; if (ref) for (let i = 0; i < st.N; i++) { const dx = mi(X[2 * i] - ref[2 * i]), dy = mi(X[2 * i + 1] - ref[2 * i + 1]); maxd = Math.max(maxd, dx * dx + dy * dy); }
      if (!ref || (st.lazy ? since >= 50 : Math.sqrt(maxd) > st.skin / 2)) { const r = cellPairs(RC + st.skin); list = r.pairs; tests += r.tests; ref = Float64Array.from(X); rebuilds.push(since); since = 0; lastBuild = performance.now(); if (rebuilds.length > 30) rebuilds.shift(); }
      since++;
      const r2 = RC * RC; for (let k = 0; k < list.length; k += 2) { const i = list[k], j = list[k + 1], dx = mi(X[2 * i] - X[2 * j]), dy = mi(X[2 * i + 1] - X[2 * j + 1]); if (dx * dx + dy * dy < r2) found++; }
      tests += list.length / 2;
    }
    hist.push([tests, performance.now() - t0, found]); if (hist.length > 40) hist.shift();
  }

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view nb-view"><canvas></canvas><div class="sim-tip"><b style="color:${COL.amber}">rc</b> and <b style="color:${COL.cyan}">rc + skin</b> around the yellow atom · highlighted = atoms it must test</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="all">All pairs</button><button data-m="cell">Cell list</button><button data-m="verlet" class="on">Verlet list</button></div>
    <label>Atoms N <output data-o="N"></output><input type="range" data-k="N" min="200" max="3000" step="100" value="800"></label>
    <label data-skin>Skin <output data-o="skin"></output><input type="range" data-k="skin" min="0" max="1.2" step="0.05" value="0.4"></label>
    <label class="check" data-lazy><input type="checkbox" data-k="lazy"><span style="flex:1">Ignore the skin/2 rule: rebuild every 50 steps</span></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="check">Check vs all pairs</button></div>
    <div class="readouts">
      <div><span>Distance tests per step</span><b data-o="tests"></b></div>
      <div><span>Fewer than all pairs</span><b data-o="speed"></b></div>
      <div><span>Time per step</span><b data-o="ms"></b></div>
      <div><span>Steps between rebuilds</span><b data-o="reb"></b></div>
      <div class="wide"><span>Same pairs as all-pairs search?</span><b data-o="ok"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.nb-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.method = b.dataset.m; root.querySelectorAll('[data-m]').forEach(x => x.classList.toggle('on', x === b)); $('[data-skin]').hidden = $('[data-lazy]').hidden = st.method !== 'verlet'; ref = null; hist = []; rebuilds = []; ok = '–'; });
  $('[data-k=N]').oninput = e => { st.N = +e.target.value; reset(); };
  $('[data-k=lazy]').onchange = e => { st.lazy = e.target.checked; ref = null; rebuilds = []; ok = '–'; };
  $('[data-k=skin]').oninput = e => { st.skin = +e.target.value; ref = null; rebuilds = []; };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=check]').onclick = () => {
    const a = allPairs(RC).count, f = hist.length ? hist[hist.length - 1][2] : NaN;
    ok = a === f ? `<span style="color:${COL.cyan}">yes: ${f} pairs, identical</span>` : `<span style="color:${COL.coral}">no: ${a - f} of ${a} pairs missed, the list is out of date</span>`;
  };

  let W = 600, H = 460;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function draw() {
    cx.clearRect(0, 0, W, H);
    const S = Math.min(W - 20, H - 20) / L, ox = (W - L * S) / 2, oy = (H - L * S) / 2;
    cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.lineWidth = 1.5; cx.strokeRect(ox, oy, L * S, L * S);
    const sx = X[2 * sel], sy = X[2 * sel + 1];
    const nc = Math.max(3, Math.floor(L / (st.method === 'verlet' ? RC + st.skin : RC))), cs = L / nc, scx = Math.floor(sx / cs), scy = Math.floor(sy / cs);
    if (st.method !== 'all') {
      cx.strokeStyle = 'rgba(143,163,199,.25)'; cx.lineWidth = 1;
      for (let k = 1; k < nc; k++) { cx.beginPath(); cx.moveTo(ox + k * cs * S, oy); cx.lineTo(ox + k * cs * S, oy + L * S); cx.moveTo(ox, oy + k * cs * S); cx.lineTo(ox + L * S, oy + k * cs * S); cx.stroke(); }
    }
    if (st.method === 'cell') { cx.fillStyle = 'rgba(56,225,255,.1)'; for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) cx.fillRect(ox + ((scx + dx + nc) % nc) * cs * S, oy + ((scy + dy + nc) % nc) * cs * S, cs * S, cs * S); }
    const inList = new Uint8Array(st.N);
    if (st.method === 'verlet') for (let k = 0; k < list.length; k += 2) { if (list[k] === sel) inList[list[k + 1]] = 1; if (list[k + 1] === sel) inList[list[k]] = 1; }
    const r = Math.max(1.2, 0.32 * S);
    for (let i = 0; i < st.N; i++) {
      let hl;
      if (st.method === 'all') hl = true;
      else if (st.method === 'cell') { const c = [Math.floor(X[2 * i] / cs), Math.floor(X[2 * i + 1] / cs)]; const ddx = Math.abs(c[0] - scx), ddy = Math.abs(c[1] - scy); hl = (ddx <= 1 || ddx === nc - 1) && (ddy <= 1 || ddy === nc - 1); }
      else hl = inList[i];
      cx.fillStyle = i === sel ? COL.amber : hl ? COL.cyan : '#34466B';
      cx.beginPath(); cx.arc(ox + X[2 * i] * S, oy + X[2 * i + 1] * S, i === sel ? r * 1.8 : r, 0, Math.PI * 2); cx.fill();
    }
    cx.strokeStyle = COL.amber; cx.lineWidth = 2; cx.beginPath(); cx.arc(ox + sx * S, oy + sy * S, RC * S, 0, Math.PI * 2); cx.stroke();
    if (st.method === 'verlet') { cx.strokeStyle = COL.cyan; cx.setLineDash([5, 5]); cx.beginPath(); cx.arc(ox + sx * S, oy + sy * S, (RC + st.skin) * S, 0, Math.PI * 2); cx.stroke(); cx.setLineDash([]); }
    if (st.method === 'verlet' && performance.now() - lastBuild < 250) { cx.fillStyle = COL.coral; cx.font = '700 14px "Space Grotesk", sans-serif'; cx.textAlign = 'left'; cx.fillText('list rebuilt', ox + 8, oy + 18); }
  }
  function readouts() {
    o.N.textContent = st.N; o.skin.textContent = st.skin.toFixed(2) + ' σ';
    if (!hist.length) return;
    const avg = k => hist.reduce((s, q) => s + q[k], 0) / hist.length, tests = avg(0), all = st.N * (st.N - 1) / 2;
    o.tests.textContent = Math.round(tests).toLocaleString('en-US'); o.speed.textContent = (all / tests).toFixed(1) + '×';
    o.ms.textContent = avg(1).toFixed(2) + ' ms';
    o.reb.textContent = st.method !== 'verlet' ? '— (searches every step)' : rebuilds.length > 1 ? (rebuilds.slice(1).reduce((s, q) => s + q, 0) / (rebuilds.length - 1)).toFixed(1) : 'measuring…';
    o.ok.innerHTML = ok;
  }
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) stepOnce(); draw(); if (frame++ % 8 === 0) readouts(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
