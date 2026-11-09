// Packing and equilibrating lab widget: 64-atom 2D LJ box at rho = 0.70, kT = 0.72.
// 1 Insert (random with tolerance d, or a lattice) -> 2 Minimize (steepest descent, GROMACS step rule) -> 3 NVT (Bussi) with automatic t0 detection.
const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7' };
const N = 64, RHO = 0.70, L = Math.sqrt(N / RHO), RC = 2.5, T0 = 0.72, DT = 0.005, EMTOL = 3.42, USH = 4 * (RC ** -12 - RC ** -6);

export function mountLab16(root) {
  const st = { d: 0.9, stage: 'empty', running: false, msg: '' };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view pk-view"><canvas></canvas></div>
  <div class="sim-ctrl">
    <label><span>Insertion tolerance d (0 = fully random)</span><output data-o="d"></output><input type="range" data-k="d" min="0" max="1.1" step="0.05" value="0.9"></label>
    <div class="row"><button class="btn sm primary" data-k="ins">1 Insert</button><button class="btn sm ghost" data-k="lat">or lattice</button></div>
    <div class="row"><button class="btn sm primary" data-k="min">2 Minimize</button><button class="btn sm primary" data-k="md">3 NVT run</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
    <div class="readouts" data-ro></div>
  </div>`;
  const cv = root.querySelector('canvas'), cx = cv.getContext('2d'), view = root.querySelector('.pk-view'), $ = s => root.querySelector(s);
  let seed = 7; const rand = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296 + 1e-12; };
  const gauss = () => Math.sqrt(-2 * Math.log(rand())) * Math.cos(2 * Math.PI * rand());
  let x, v, F, U, tries, fails, trace, h, iter, tm;
  const mi = a => a - L * Math.round(a / L);
  function force() {
    F = x.map(() => [0, 0]); U = 0;
    for (let i = 0; i < x.length; i++) for (let j = i + 1; j < x.length; j++) {
      const dx = mi(x[i][0] - x[j][0]), dy = mi(x[i][1] - x[j][1]), r2 = Math.max(dx * dx + dy * dy, 1e-6);
      if (r2 < RC * RC) { const s6 = 1 / (r2 * r2 * r2), f = 24 * (2 * s6 * s6 - s6) / r2; U += 4 * (s6 * s6 - s6) - USH; F[i][0] += f * dx; F[i][1] += f * dy; F[j][0] -= f * dx; F[j][1] -= f * dy; }
    }
  }
  const fmax = () => Math.max(...F.map(f => Math.hypot(f[0], f[1])));
  const close = () => { let c = 0; for (let i = 0; i < x.length; i++) for (let j = i + 1; j < x.length; j++) if (mi(x[i][0] - x[j][0]) ** 2 + mi(x[i][1] - x[j][1]) ** 2 < 0.64) c++; return c; };
  const temp = () => v.reduce((a, p) => a + p[0] * p[0] + p[1] * p[1], 0) / (2 * N - 2);
  function reset() { x = []; v = null; F = null; U = 0; tries = 0; fails = 0; trace = []; h = 0.01; iter = 0; tm = 0; st.stage = 'empty'; st.running = false; st.msg = ''; }
  function insertStep() {
    for (let k = 0; k < 60 && x.length < N && fails < 4000; k++) {
      const c = [rand() * L, rand() * L]; tries++;
      if (x.some(p => mi(p[0] - c[0]) ** 2 + mi(p[1] - c[1]) ** 2 < st.d * st.d)) { fails++; continue; }
      x.push(c); fails = 0; trace.push([tries, x.length]);
    }
    if (x.length >= N || fails >= 4000) { st.running = false; force(); st.stage = 'inserted'; st.msg = x.length < N ? `jammed: only ${x.length} of ${N} fit at d = ${st.d.toFixed(2)}` : ''; }
  }
  function lattice() { reset(); const nx = 8, ny = 8, a = L / nx; for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) x.push([(i + 0.5 * (j % 2) + 0.25) * a, (j + 0.5) * (L / ny)]); force(); st.stage = 'inserted'; trace = [[0, N]]; }
  function minStep() {
    for (let k = 0; k < 3; k++) {
      const fm = fmax(); if (fm < EMTOL || iter > 3000) { st.running = false; st.stage = 'minimized'; return; }
      const x0 = x.map(p => [...p]), F0 = F, U0 = U;
      x = x.map((p, i) => [((p[0] + h * F0[i][0] / fm) % L + L) % L, ((p[1] + h * F0[i][1] / fm) % L + L) % L]); force();
      if (U < U0) h *= 1.2; else { x = x0; F = F0; U = U0; h *= 0.2; }
      iter++; trace.push([iter, Math.log10(fmax())]);
    }
  }
  function bussi() {
    const nf = 2 * N - 2, K = 0.5 * nf * temp(), Kt = 0.5 * nf * T0, c = Math.exp(-DT / 0.5), R = gauss(); let S = 0; for (let k = 0; k < nf - 1; k++) S += gauss() ** 2;
    const a2 = c + Kt / (nf * K) * (1 - c) * (R * R + S) + 2 * R * Math.sqrt(c * Kt / (nf * K) * (1 - c)); const a = Math.sqrt(Math.max(a2, 0)); v.forEach(p => { p[0] *= a; p[1] *= a; });
  }
  function mdStep() {
    for (let k = 0; k < 10; k++) {
      v.forEach((p, i) => { p[0] += 0.5 * DT * F[i][0]; p[1] += 0.5 * DT * F[i][1]; });
      x = x.map((p, i) => [((p[0] + DT * v[i][0]) % L + L) % L, ((p[1] + DT * v[i][1]) % L + L) % L]); force();
      v.forEach((p, i) => { p[0] += 0.5 * DT * F[i][0]; p[1] += 0.5 * DT * F[i][1]; });
      const T = temp(); if (!(T < 1e4)) { st.running = false; st.msg = `blew up: T = ${T.toExponential(1)} after ${Math.round(tm / DT) + k + 1} steps (target 0.72)`; return; }
      bussi(); tm += DT;
    }
    trace.push([tm, U / N]);
  }
  function ineff(a) { const n = a.length, m = a.reduce((p, q) => p + q, 0) / n, va = a.reduce((p, q) => p + (q - m) ** 2, 0) / n; if (!va) return 1; let g = 1;
    for (let k = 1; k < n / 3; k++) { let c = 0; for (let i = 0; i + k < n; i++) c += (a[i] - m) * (a[i + k] - m); const r = c / (n - k) / va; if (r < 0) break; g += 2 * r; } return g; }
  function pickT0() { const a = trace.map(p => p[1]), n = a.length; if (n < 60) return null; let best = null;
    for (let t = 0; t < n * 0.75; t += Math.max(1, Math.floor(n / 30))) { const ne = (n - t) / ineff(a.slice(t)); if (!best || ne > best[1]) best = [t, ne]; } return best; }

  $('[data-k=d]').oninput = e => { st.d = +e.target.value; $('[data-o=d]').textContent = st.d.toFixed(2); };
  $('[data-k=ins]').onclick = () => { reset(); st.stage = 'inserting'; st.running = true; };
  $('[data-k=lat]').onclick = () => { lattice(); draw(); };
  $('[data-k=min]').onclick = () => { if (!x.length || st.stage === 'inserting') return; if (st.stage === 'md') force(); trace = []; iter = 0; h = 0.01; st.msg = ''; st.stage = 'minimizing'; st.running = true; };
  $('[data-k=md]').onclick = () => { if (x.length < N) { st.msg = x.length ? 'not all atoms are in the box yet' : 'insert atoms first'; draw(); return; }
    v = x.map(() => [gauss(), gauss()]); const m0 = v.reduce((a, p) => [a[0] + p[0] / N, a[1] + p[1] / N], [0, 0]); v.forEach(p => { p[0] -= m0[0]; p[1] -= m0[1]; }); const s = Math.sqrt(T0 / temp()); v.forEach(p => { p[0] *= s; p[1] *= s; });
    force(); trace = []; tm = 0; st.msg = ''; st.stage = 'md'; st.running = true; };
  $('[data-k=reset]').onclick = () => { reset(); draw(); };

  let W = 600, H = 600;
  const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
  const ro = new ResizeObserver(resize); ro.observe(view);
  const label = (t, x_, y, c = COL.ink, al = 'left', f = '600 12px') => { cx.fillStyle = c; cx.font = `${f} "Space Grotesk", sans-serif`; cx.textAlign = al; cx.fillText(t, x_, y); };
  function draw() {
    if (!W) return;
    cx.clearRect(0, 0, W, H);
    const narrow = W < 520, bw = narrow ? Math.min(W * 0.55, H * 0.42) : Math.min(W * 0.45, H - 40), S = bw / L, bx = narrow ? (W - bw) / 2 : 10;
    cx.strokeStyle = 'rgba(56,225,255,.5)'; cx.lineWidth = 1; cx.strokeRect(bx, 10, bw, bw);
    const bad = new Set(); for (let i = 0; i < x.length; i++) for (let j = i + 1; j < x.length; j++) if (mi(x[i][0] - x[j][0]) ** 2 + mi(x[i][1] - x[j][1]) ** 2 < 0.64) { bad.add(i); bad.add(j); }
    x.forEach((p, i) => { cx.fillStyle = bad.has(i) ? COL.coral : COL.cyan; cx.globalAlpha = 0.85; cx.beginPath(); cx.arc(bx + p[0] * S, 10 + (L - p[1]) * S, 0.5 * S, 0, 7); cx.fill(); });
    cx.globalAlpha = 1; label('coral: pairs closer than 0.8', bx, bw + 26, COL.muted, 'left', '11px');
    const x0 = narrow ? 44 : bw + 56, x1 = W - 12, y0 = narrow ? bw + 58 : 26, y1 = H - 22;
    cx.strokeStyle = 'rgba(143,163,199,.6)'; cx.lineWidth = 1.2; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x0, y1); cx.lineTo(x1, y1); cx.stroke();
    const md = st.stage === 'md', mn = st.stage === 'minimizing' || st.stage === 'minimized';
    label(md ? 'U/N vs time; dashed: detected t0' : mn ? 'log10 |F|max vs iteration; dashed: emtol' : 'atoms placed vs trials', x0 + 6, y0 - 8);
    if (trace.length > 1) {
      const xs = trace.map(p => p[0]), ys = trace.map(p => p[1]), xa = Math.min(...xs), xb = Math.max(...xs) || 1;
      let ya = Math.min(...ys), yb = Math.max(...ys); if (mn) { ya = Math.min(ya, 0); } if (!md && !mn) { ya = 0; yb = N; } if (yb - ya < 0.2) { ya -= 0.1; yb += 0.1; }
      const X = q => x0 + (q - xa) / (xb - xa || 1) * (x1 - x0), Y = q => y1 - (q - ya) / (yb - ya) * (y1 - y0);
      cx.strokeStyle = md ? COL.amber : mn ? COL.coral : COL.cyan; cx.lineWidth = 1.6; cx.beginPath(); trace.forEach((p, k) => k ? cx.lineTo(X(p[0]), Y(p[1])) : cx.moveTo(X(p[0]), Y(p[1]))); cx.stroke();
      cx.setLineDash([4, 4]); cx.strokeStyle = 'rgba(242,245,250,.7)';
      if (mn) { const yl = Y(Math.log10(EMTOL)); cx.beginPath(); cx.moveTo(x0, yl); cx.lineTo(x1, yl); cx.stroke(); }
      if (md) { const b = pickT0(); if (b) { const xt = X(trace[b[0]][0]); cx.beginPath(); cx.moveTo(xt, y0); cx.lineTo(xt, y1); cx.stroke(); } }
      cx.setLineDash([]);
      label(`${ya.toFixed(md ? 2 : mn ? 0 : 0)}`, x0 - 4, y1, COL.muted, 'right', '10px'); label(`${yb.toFixed(md ? 2 : mn ? 0 : 0)}`, x0 - 4, y0 + 8, COL.muted, 'right', '10px');
    }
    const R = $('[data-ro]'); let html = '';
    const it = (k, val) => { html += `<div><span>${k}</span><b>${val}</b></div>`; };
    if (st.stage === 'empty') it('Next', 'press 1 Insert');
    else if (st.stage === 'inserting' || st.stage === 'inserted') { it('Placed / trials', `${x.length} / ${tries.toLocaleString()}`); it('Density reached (target 0.70)', (x.length / L / L).toFixed(3)); if (F) { it('Pairs closer than 0.8', close()); it('U/N  ·  max force', `${(U / x.length).toExponential(2)} · ${fmax().toExponential(1)}`); } }
    else if (mn) { it('Iteration  ·  step h', `${iter} · ${h.toExponential(1)}`); it('Max force (stop below 3.4)', fmax().toExponential(2)); it('U/N', (U / N).toFixed(3)); it('Pairs closer than 0.8', close()); }
    else if (md) { const b = pickT0(), a = trace.map(p => p[1]); const fmt = q => Math.abs(q) < 1e4 ? q.toFixed(Math.abs(q) < 10 ? 2 : 0) : q.toExponential(1); it('Time  ·  T (target 0.72)', `${tm.toFixed(1)} · ${v ? fmt(temp()) : '…'}`); it('U/N now', fmt(U / N));
      if (b) { const k = a.slice(b[0]), m = k.reduce((p, q) => p + q, 0) / k.length; it('Detected t0', trace[b[0]][0].toFixed(1)); it('U/N after t0  ·  independent samples', `${m.toFixed(3)} · ${b[1].toFixed(0)}`); } else it('Detected t0', 'collecting…'); }
    if (st.msg) it('Note', st.msg);
    R.innerHTML = html;
  }
  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible || !st.running) return;
    if (st.stage === 'inserting') insertStep(); else if (st.stage === 'minimizing') minStep(); else if (st.stage === 'md') mdStep(); draw(); };
  $('[data-o=d]').textContent = st.d.toFixed(2); reset(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); };
}
