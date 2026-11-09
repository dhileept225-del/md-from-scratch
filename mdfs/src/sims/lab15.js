// MC vs MD lab widget.
// Step size: 9-atom 2D LJ box (kT = 0.72) sampled by Metropolis with step delta; acceptance, energy trace, inefficiency g.
// Barrier: double well U = 5 (s^2 - 1)^2 (kT = 1): Langevin MD vs MC small steps vs MC with 10 % flip moves.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const L = 3.6, N = 9, RC = 1.8, T0 = 0.72, USH = 4 * (RC ** -12 - RC ** -6);
const U2 = s => 5 * (s * s - 1) ** 2, dU2 = s => 20 * s * (s * s - 1);

export function mountLab15(root) {
  const st = { mode: 'step', delta: 0.15, running: false };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view mc-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-m="step" class="on">Step size</button><button data-m="bar">Barrier</button></div>
    <div data-p="step"><label><span>Maximum step δ</span><output data-o="d"></output><input type="range" data-k="d" min="-2" max="0" step="0.05" value="-0.82"></label></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.mc-view'), $ = s => root.querySelector(s);
  let seed = 17; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };
  const gauss = () => Math.sqrt(-2 * Math.log(rand())) * Math.cos(2 * Math.PI * rand());
  // ---- step-size mode ----
  let x, Us, acc, tries, sweeps;
  function uAtom(i, p) { let u = 0; for (let j = 0; j < N; j++) { if (j === i) continue; let dx = x[j][0] - p[0], dy = x[j][1] - p[1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); const r2 = dx * dx + dy * dy; if (r2 < RC * RC) { const s6 = 1 / (r2 * r2 * r2); u += 4 * (s6 * s6 - s6) - USH; } } return u; }
  function uTot() { // plain LJ over pairs inside rc (the lab's observable; the shift only enters the moves)
    let u = 0; for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { let dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); const r2 = dx * dx + dy * dy; if (r2 < RC * RC) { const s6 = 1 / (r2 * r2 * r2); u += 4 * (s6 * s6 - s6); } } return u; }
  function resetStep() { x = [0, 1, 2].flatMap(j => [0, 1, 2].map(i => [(i + 0.5) * 1.2, (j + 0.5) * 1.2])); Us = []; acc = 0; tries = 0; sweeps = 0; }
  function sweep() {
    for (let m = 0; m < N; m++) { const i = Math.floor(rand() * N), old = x[i], nw = [((old[0] + st.delta * (2 * rand() - 1)) % L + L) % L, ((old[1] + st.delta * (2 * rand() - 1)) % L + L) % L];
      const dU = uAtom(i, nw) - uAtom(i, old); tries++; if (rand() < Math.exp(-dU / T0)) { x[i] = nw; acc++; } }
    sweeps++; if (sweeps > 200) { Us.push(uTot() / N); if (Us.length > 6000) Us.shift(); }
  }
  function ineff() {
    const n = Us.length; if (n < 400) return NaN; const m = Us.reduce((a, b) => a + b, 0) / n, v = Us.reduce((a, b) => a + (b - m) ** 2, 0) / n; let g = 1;
    for (let k = 1; k < Math.min(500, n / 4); k++) { let c = 0; for (let i = 0; i + k < n; i++) c += (Us[i] - m) * (Us[i + k] - m); const r = c / (n - k) / v; if (r < 0) break; g += 2 * r; }
    return g;
  }
  // ---- barrier mode ----
  let W3, H3, nb3;
  function resetBar() { W3 = { md: { s: -1, v: 0, right: 0, cross: 0, side: -1, n: 0 }, mc: { s: -1, right: 0, cross: 0, side: -1, n: 0 }, fl: { s: -1, right: 0, cross: 0, side: -1, n: 0 } }; nb3 = 40; H3 = { md: new Float64Array(nb3), mc: new Float64Array(nb3), fl: new Float64Array(nb3) }; }
  const tally = (w, h) => { if (w.s * w.side < 0 && Math.abs(w.s) > 0.8) { w.cross++; w.side = -w.side; } w.right += w.s > 0; w.n++; const b = Math.floor((w.s + 2) / 4 * nb3); if (b >= 0 && b < nb3) h[b]++; };
  function barStep() {
    const md = W3.md, c = Math.exp(-0.005), sn = Math.sqrt(1 - c * c);
    md.v += -0.0005 * dU2(md.s); md.s += 0.001 * md.v; md.v += -0.0005 * dU2(md.s); md.v = c * md.v + sn * gauss(); tally(md, H3.md);
    for (const [key, pf] of [['mc', 0], ['fl', 0.1]]) { const w = W3[key], sn2 = rand() < pf ? -w.s : w.s + 0.3 * (2 * rand() - 1); if (rand() < Math.exp(-(U2(sn2) - U2(w.s)))) w.s = sn2; tally(w, H3[key]); }
  }
  root.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { root.querySelectorAll('[data-m]').forEach(z => z.classList.toggle('on', z === b)); st.mode = b.dataset.m;
    root.querySelectorAll('[data-p]').forEach(p => p.hidden = p.dataset.p !== st.mode); draw(); });
  $('[data-k=d]').oninput = e => { st.delta = 10 ** +e.target.value; resetStep(); draw(); };
  $('[data-k=run]').onclick = e => { st.running = !st.running; e.target.textContent = st.running ? '❚❚ Pause' : '▶ Run'; };
  $('[data-k=reset]').onclick = () => { st.mode === 'step' ? resetStep() : resetBar(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const axes = (x0, y0, x1, y1) => { cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke(); };
  const label = (t, x_, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x_, y); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const R = $('[data-ro]'), half = H * 0.45;
    if (st.mode === 'step') {
      const narrow = W < 520, bw = narrow ? Math.min(W * 0.5, H * 0.4) : Math.min(W * 0.42, half - 16), S = bw / L;
      cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.strokeRect(8, 8, bw, bw);
      x.forEach(p => { cx.fillStyle = COL.cyan; cx.beginPath(); cx.arc(8 + p[0] * S, 8 + (L - p[1]) * S, 0.3 * S, 0, 7); cx.fill(); });
      cx.strokeStyle = 'rgba(255,194,75,.5)'; cx.setLineDash([3, 3]); const p0 = x[4]; cx.strokeRect(8 + (p0[0] - st.delta) * S, 8 + (L - p0[1] - st.delta) * S, 2 * st.delta * S, 2 * st.delta * S); cx.setLineDash([]);
      label('dashed square: trial range of one atom', 10, bw + 24, COL.muted, 'left', '11px');
      const x0 = narrow ? 40 : bw + 50, x1 = W - 12, y0 = narrow ? bw + 52 : 22, y1 = H - 20, show = Us.slice(-400), lo = Math.min(-2.1, ...show) - 0.05, hi = Math.max(-1.5, ...show) + 0.05, X = k => x0 + k / 399 * (x1 - x0), Y = u => y1 - (u - lo) / (hi - lo) * (y1 - y0);
      axes(x0, y0, x1, y1); label('U/N per sweep (last 400 sweeps)', x0 + 6, y0 - 6);
      cx.strokeStyle = COL.amber; cx.lineWidth = 1.5; cx.beginPath(); show.forEach((u, k) => k ? cx.lineTo(X(k), Y(u)) : cx.moveTo(X(k), Y(u))); cx.stroke();
      const g = ineff(), a = tries ? acc / tries : NaN;
      R.innerHTML = `<div><span>Acceptance</span><b>${tries ? a.toFixed(3) : '…'}</b></div><div><span>Inefficiency g (sweeps)</span><b>${isFinite(g) ? g.toFixed(1) : 'collecting…'}</b></div>
        <div><span>Independent samples / 1000 sweeps</span><b>${isFinite(g) ? (1000 / g).toFixed(0) : '…'}</b></div><div><span>Average U/N (MD: −1.809)</span><b>${Us.length ? (Us.reduce((p, q) => p + q, 0) / Us.length).toFixed(3) : '…'}</b></div>`;
      $('[data-o=d]').textContent = st.delta.toFixed(3);
    } else {
      const x0 = 46, x1 = W - 12, y0 = 22, y1 = H * 0.6, X = s => x0 + (s + 2) / 4 * (x1 - x0), Z = Array.from({ length: nb3 }, (_, b) => Math.exp(-U2(-2 + (b + 0.5) * 4 / nb3))), zs = Z.reduce((a, b) => a + b, 0), ex = Z.map(z => z / zs);
      let pmax = Math.max(...ex); for (const k of ['md', 'mc', 'fl']) { const n = W3[k].n || 1; H3[k].forEach(c => pmax = Math.max(pmax, c / n)); }
      const Y = p => y1 - p / (pmax * 1.1) * (y1 - y0); axes(x0, y0, x1, y1); label('where each walker spends its time; white: exact', x0 + 6, y0 - 6);
      cx.strokeStyle = 'rgba(242,245,250,.8)'; cx.lineWidth = 2; cx.beginPath(); ex.forEach((p, b) => { const xx = X(-2 + (b + 0.5) * 4 / nb3); b ? cx.lineTo(xx, Y(p)) : cx.moveTo(xx, Y(p)); }); cx.stroke();
      cx.setLineDash([5, 4]);
      [['md', COL.coral], ['mc', COL.amber], ['fl', COL.cyan]].forEach(([k, c]) => { const n = W3[k].n || 1; cx.strokeStyle = c; cx.lineWidth = 1.8; cx.beginPath(); H3[k].forEach((h, b) => { const xx = X(-2 + (b + 0.5) * 4 / nb3); b ? cx.lineTo(xx, Y(h / n)) : cx.moveTo(xx, Y(h / n)); }); cx.stroke(); });
      cx.setLineDash([]); cx.strokeStyle = 'rgba(242,245,250,.9)'; cx.lineWidth = 1.5; cx.beginPath(); ex.forEach((p, b) => { const xx = X(-2 + (b + 0.5) * 4 / nb3); b ? cx.lineTo(xx, Y(p)) : cx.moveTo(xx, Y(p)); }); cx.stroke();
      label('coral: MD   amber: MC small steps   cyan: MC + flips', x0 + 6, y1 + 18, COL.muted, 'left', '11px');
      const row_ = (k, name) => `<div><span>${name}: right fraction · crossings</span><b>${W3[k].n ? (W3[k].right / W3[k].n).toFixed(3) : '…'} · ${W3[k].cross}</b></div>`;
      R.innerHTML = row_('md', 'MD') + row_('mc', 'MC small steps') + row_('fl', 'MC + flips') + `<div><span>Steps / moves each</span><b>${W3.md.n.toLocaleString()}</b></div>`;
    }
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return; if (st.mode === 'step') for (let k = 0; k < 40; k++) sweep(); else for (let k = 0; k < 4000; k++) barStep(); draw(); };
  resetStep(); resetBar(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
