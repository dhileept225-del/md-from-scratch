// PBC playground: periodic copies, minimum-image neighbours inside a cutoff, and wrapped vs unwrapped paths.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 30;

export function mountPBC(root) {
  const st = { copies: true, mic: true, unwrap: true, rc: 0.3, running: true, sel: 0 };
  let X, V, X0, trail;
  function reset() {
    X = Array.from({ length: N }, () => [Math.random(), Math.random()]); X0 = X.map(p => [...p]);
    V = Array.from({ length: N }, () => [(Math.random() - 0.5) * 0.3, (Math.random() - 0.5) * 0.3]); trail = [X[st.sel].slice()];
  }
  const fr = q => q - Math.floor(q);

  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view pbc-view"><canvas></canvas><div class="sim-tip">Click an atom to select it</div></div>
  <div class="sim-ctrl">
    <label class="check"><input type="checkbox" data-k="copies" checked> Show periodic copies</label>
    <label class="check"><input type="checkbox" data-k="mic" checked> Minimum image (nearest copy)</label>
    <label class="check"><input type="checkbox" data-k="unwrap" checked> Path: unwrapped (off = wrapped)</label>
    <label>Cutoff rc <output data-o="rc"></output><input type="range" data-k="rc" min="0.1" max="0.75" step="0.01" value="0.3"></label>
    <div class="row"><button class="btn sm primary" data-k="run">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">↺ Reset</button></div>
    <div class="readouts">
      <div><span>Neighbours (PBC)</span><b data-o="nb"></b></div>
      <div><span>Neighbours (no PBC)</span><b data-o="nn"></b></div>
      <div><span>True distance travelled</span><b data-o="du"></b></div>
      <div><span>From wrapped coords</span><b data-o="dw"></b></div>
      <div class="wide"><span>Cutoff check</span><b data-o="chk"></b></div>
    </div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.pbc-view');
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  ['copies', 'mic', 'unwrap'].forEach(k => $(`[data-k=${k}]`).onchange = e => { st[k] = e.target.checked; });
  $('[data-k=rc]').oninput = e => { st.rc = +e.target.value; };
  const runBtn = $('[data-k=run]'); runBtn.onclick = () => { st.running = !st.running; runBtn.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = reset;

  let W = 600, H = 460, S = 200, ox = 0, oy = 0;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();
  cv.addEventListener('pointerdown', e => {
    const b = cv.getBoundingClientRect(), px = (e.clientX - b.left - ox) / S, py = (e.clientY - b.top - oy) / S;
    let best = -1, bd = 1e9;
    X.forEach((p, i) => { let dx = fr(p[0]) - px, dy = fr(p[1]) - py; dx -= Math.round(dx); dy -= Math.round(dy); const d = dx * dx + dy * dy; if (d < bd) { bd = d; best = i; } });
    if (best >= 0 && bd < 0.01) { st.sel = best; X0[best] = X[best].slice(); trail = [X[best].slice()]; }
  });

  function step(dt) {
    for (let i = 0; i < N; i++) for (let d = 0; d < 2; d++) { V[i][d] += (-1.5 * V[i][d]) * dt + 0.35 * Math.sqrt(dt) * (Math.random() * 2 - 1) * 1.7; X[i][d] += V[i][d] * dt; }
    trail.push(X[st.sel].slice()); if (trail.length > 900) trail.shift();
  }
  function draw() {
    cx.clearRect(0, 0, W, H);
    S = Math.min(W, H) / (st.copies ? 3.25 : 1.35); ox = W / 2 - S / 2; oy = H / 2 - S / 2;
    const P = (p, gx = 0, gy = 0) => [ox + (p[0] + gx) * S, oy + (p[1] + gy) * S];
    const rng = st.copies ? [-1, 0, 1] : [0];
    // boxes and atoms
    rng.forEach(gx => rng.forEach(gy => {
      const c = !gx && !gy; cx.globalAlpha = c ? 1 : 0.35;
      cx.strokeStyle = COL.cyan; cx.lineWidth = c ? 2.5 : 1; cx.strokeRect(ox + gx * S, oy + gy * S, S, S);
      X.forEach((p, i) => { const q = P([fr(p[0]), fr(p[1])], gx, gy); cx.fillStyle = i === st.sel ? COL.amber : i % 3 ? COL.cyan : COL.coral; cx.beginPath(); cx.arc(q[0], q[1], i === st.sel ? 7 : 5, 0, Math.PI * 2); cx.fill(); });
    }));
    cx.globalAlpha = 1;
    // path of the selected atom
    cx.strokeStyle = st.unwrap ? COL.amber : COL.coral; cx.lineWidth = 2; cx.beginPath();
    trail.forEach((p, k) => {
      const q = st.unwrap ? P(p) : P([fr(p[0]), fr(p[1])]);
      const jump = k > 0 && (Math.floor(p[0]) !== Math.floor(trail[k - 1][0]) || Math.floor(p[1]) !== Math.floor(trail[k - 1][1]));
      if (k === 0 || (!st.unwrap && jump)) cx.moveTo(q[0], q[1]); else cx.lineTo(q[0], q[1]);
    });
    cx.stroke();
    // neighbours of the selected atom
    const s = [fr(X[st.sel][0]), fr(X[st.sel][1])], sp = P(s); let nb = 0, nn = 0;
    X.forEach((p, i) => {
      if (i === st.sel) return; const q = [fr(p[0]), fr(p[1])];
      let dx = q[0] - s[0], dy = q[1] - s[1];
      if (Math.hypot(dx, dy) < st.rc) nn++;
      if (st.mic) { dx -= Math.round(dx); dy -= Math.round(dy); }
      if (Math.hypot(dx, dy) < st.rc) { nb++; cx.strokeStyle = 'rgba(255,194,75,.75)'; cx.lineWidth = 1.5; cx.beginPath(); cx.moveTo(sp[0], sp[1]); cx.lineTo(sp[0] + dx * S, sp[1] + dy * S); cx.stroke(); }
    });
    const bad = st.rc > 0.5;
    cx.strokeStyle = bad ? COL.coral : COL.amber; cx.lineWidth = 2; cx.setLineDash([6, 5]); cx.beginPath(); cx.arc(sp[0], sp[1], st.rc * S, 0, Math.PI * 2); cx.stroke(); cx.setLineDash([]);
    // readouts
    const u = X[st.sel], u0 = X0[st.sel], dU = Math.hypot(u[0] - u0[0], u[1] - u0[1]), dW = Math.hypot(fr(u[0]) - fr(u0[0]), fr(u[1]) - fr(u0[1]));
    o.rc.textContent = st.rc.toFixed(2) + ' L'; o.nb.textContent = st.mic ? nb : '(turn on minimum image)'; o.nn.textContent = nn;
    o.du.textContent = dU.toFixed(2) + ' L'; o.dw.textContent = dW.toFixed(2) + ' L';
    o.chk.innerHTML = bad ? `<span style="color:${COL.coral}">rc > L/2: an atom can meet its own copy, and pairs get counted twice</span>` : `<span style="color:${COL.cyan}">rc ≤ L/2: safe</span>`;
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; if (st.running) step(1 / 60); draw(); };
  reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
