// Live Lennard-Jones MD in the browser: velocity Verlet + optional Berendsen thermostat.
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

const CYAN = new THREE.Color('#38E1FF'), AMBER = new THREE.Color('#FFC24B'), CORAL = new THREE.Color('#FF6B6B');

class LJ {
  constructor({ n = 5, L = 8.5, rc = 2.5 } = {}) {
    this.N = n ** 3; this.L = L; this.rc = rc; this.n = n;
    this.vc = rc ** -12 - rc ** -6;
    this.x = new Float64Array(this.N * 3); this.v = new Float64Array(this.N * 3); this.f = new Float64Array(this.N * 3);
  }
  reset(T) {
    const { n, L, N, x, v } = this, a = L / n;
    let k = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let l = 0; l < n; l++, k++) x.set([(i + 0.5) * a, (j + 0.5) * a, (l + 0.5) * a], k * 3);
    const g = () => Math.sqrt(-2 * Math.log(Math.random() + 1e-12)) * Math.cos(2 * Math.PI * Math.random());
    for (let i = 0; i < N * 3; i++) v[i] = g() * Math.sqrt(T);
    for (let d = 0; d < 3; d++) { let m = 0; for (let i = 0; i < N; i++) m += v[i * 3 + d]; m /= N; for (let i = 0; i < N; i++) v[i * 3 + d] -= m; }
    this.scaleTo(T); this.step_n = 0; this.forces();
  }
  ke() { let s = 0; for (const u of this.v) s += u * u; return 0.5 * s; }
  temp() { return 2 * this.ke() / (3 * this.N - 3); }
  scaleTo(T) { const s = Math.sqrt(T / Math.max(this.temp(), 1e-9)); for (let i = 0; i < this.v.length; i++) this.v[i] *= s; }
  forces() {
    const { N, L, x, f, rc, vc } = this, rc2 = rc * rc; f.fill(0); let pe = 0;
    for (let i = 0; i < N; i++) {
      const xi = x[i * 3], yi = x[i * 3 + 1], zi = x[i * 3 + 2];
      for (let j = i + 1; j < N; j++) {
        let dx = xi - x[j * 3], dy = yi - x[j * 3 + 1], dz = zi - x[j * 3 + 2];
        dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > rc2) continue;
        const i2 = 1 / r2, i6 = i2 * i2 * i2, ff = 24 * i2 * i6 * (2 * i6 - 1);
        f[i * 3] += ff * dx; f[i * 3 + 1] += ff * dy; f[i * 3 + 2] += ff * dz;
        f[j * 3] -= ff * dx; f[j * 3 + 1] -= ff * dy; f[j * 3 + 2] -= ff * dz;
        pe += 4 * (i6 * i6 - i6 - vc);
      }
    }
    this.pe = pe;
  }
  step(dt, T0, thermo) {
    const { x, v, f, L } = this, h = 0.5 * dt;
    for (let i = 0; i < x.length; i++) { v[i] += h * f[i]; x[i] += dt * v[i]; x[i] -= L * Math.floor(x[i] / L); }
    this.forces();
    for (let i = 0; i < v.length; i++) v[i] += h * f[i];
    if (thermo) { const T = this.temp(), lam = Math.sqrt(1 + (dt / 0.4) * (T0 / Math.max(T, 1e-9) - 1)); for (let i = 0; i < v.length; i++) v[i] *= lam; }
    this.step_n++;
  }
}

function controlsHTML() {
  return `
  <div class="sim-ctrl">
    <div class="row"><button class="btn sm primary" data-k="play">❚❚ Pause</button><button class="btn sm ghost" data-k="reset">↺ Reset</button></div>
    <label>Target temperature <output data-o="T">1.00</output><input type="range" data-k="T" min="0.1" max="2.5" step="0.05" value="1.0"></label>
    <label>Time step Δt <output data-o="dt">0.005</output><input type="range" data-k="dt" min="0.001" max="0.03" step="0.001" value="0.005"></label>
    <label class="check"><input type="checkbox" data-k="thermo" checked> Thermostat on <span class="muted">(off = energy should stay constant)</span></label>
    <div class="readouts">
      <div><span>Temperature</span><b data-o="Tnow">–</b></div>
      <div><span>Total energy</span><b data-o="E">–</b></div>
      <div><span>Steps</span><b data-o="steps">0</b></div>
      <div><span>Time (argon)</span><b data-o="time">0 ps</b></div>
    </div>
    <canvas class="echart" width="560" height="200"></canvas>
    <div class="legend"><i style="--c:#FFC24B"></i>Total <i style="--c:#38E1FF"></i>Kinetic <i style="--c:#FF6B6B"></i>Potential</div>
  </div>`;
}

export function mountLJ(root, { mini = false } = {}) {
  const md = new LJ();
  const state = { T: mini ? 0.7 : 1.0, dt: 0.005, thermo: true, playing: true, blown: false };
  md.reset(state.T);

  root.classList.add('sim');
  root.innerHTML = `<div class="sim-view"><div class="sim-msg" hidden></div>${mini ? '' : '<div class="sim-tip">Colour = speed: <b style="color:#38E1FF">slow</b> → <b style="color:#FF6B6B">fast</b> · drag to rotate</div>'}</div>${mini ? '' : controlsHTML()}`;
  const view = root.querySelector('.sim-view'), msg = root.querySelector('.sim-msg');

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  view.prepend(renderer.domElement);
  const scene = new THREE.Scene();
  const cam = new THREE.PerspectiveCamera(40, 1, 0.1, 200); cam.position.set(14, 9, 16);
  scene.add(new THREE.HemisphereLight(0xcdefff, 0x1b2550, 1.7));
  const dl = new THREE.DirectionalLight(0xffffff, 2.2); dl.position.set(6, 10, 8); scene.add(dl);
  const atoms = new THREE.InstancedMesh(new THREE.SphereGeometry(0.5, 24, 16), new THREE.MeshStandardMaterial({ roughness: 0.35, metalness: 0.05 }), md.N);
  const box = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(md.L, md.L, md.L)), new THREE.LineBasicMaterial({ color: 0x38e1ff, transparent: true, opacity: 0.5 }));
  scene.add(atoms, box);
  const orbit = new OrbitControls(cam, renderer.domElement);
  orbit.enableDamping = true; orbit.autoRotate = true; orbit.autoRotateSpeed = mini ? 1.2 : 0.6; orbit.enablePan = false;
  orbit.minDistance = 8; orbit.maxDistance = 40;
  if (mini) { orbit.enableZoom = false; cam.position.set(12, 8, 14); }

  const resize = () => { const w = view.clientWidth, h = view.clientHeight; renderer.setSize(w, h, false); cam.aspect = w / h; cam.updateProjectionMatrix(); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  // ui
  const $ = s => root.querySelector(s), outs = {};
  root.querySelectorAll('[data-o]').forEach(o => outs[o.dataset.o] = o);
  const hist = [];
  if (!mini) {
    $('[data-k=T]').oninput = e => { state.T = +e.target.value; outs.T.textContent = state.T.toFixed(2); };
    $('[data-k=dt]').oninput = e => { state.dt = +e.target.value; outs.dt.textContent = state.dt.toFixed(3); };
    $('[data-k=thermo]').onchange = e => { state.thermo = e.target.checked; hist.length = 0; };
    $('[data-k=play]').onclick = e => { state.playing = !state.playing; e.target.textContent = state.playing ? '❚❚ Pause' : '▶ Play'; };
    $('[data-k=reset]').onclick = () => { md.reset(state.T); hist.length = 0; state.blown = false; msg.hidden = true; state.playing = true; $('[data-k=play]').textContent = '❚❚ Pause'; };
  }
  const chart = $('.echart'), cx = chart?.getContext('2d');
  function drawChart() {
    const w = chart.width, h = chart.height; cx.clearRect(0, 0, w, h);
    if (hist.length < 2) return;
    let lo = Infinity, hi = -Infinity; for (const p of hist) for (const k of p) { lo = Math.min(lo, k); hi = Math.max(hi, k); }
    const pad = (hi - lo) * 0.1 + 1e-6; lo -= pad; hi += pad;
    cx.strokeStyle = 'rgba(143,163,199,.25)'; cx.lineWidth = 1; cx.beginPath(); cx.moveTo(0, h - 1); cx.lineTo(w, h - 1); cx.stroke();
    ['#FFC24B', '#38E1FF', '#FF6B6B'].forEach((c, k) => {
      cx.strokeStyle = c; cx.lineWidth = 3; cx.beginPath();
      hist.forEach((p, i) => { const X = i / 299 * w, Yv = h - (p[k] - lo) / (hi - lo) * h; i ? cx.lineTo(X, Yv) : cx.moveTo(X, Yv); });
      cx.stroke();
    });
  }

  const dummy = new THREE.Object3D(), c = new THREE.Color();
  let raf = 0, visible = true, frame = 0;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  function loop() {
    raf = requestAnimationFrame(loop);
    if (!visible) return;
    if (state.playing && !state.blown) {
      for (let s = 0; s < 6; s++) md.step(state.dt, state.T, state.thermo);
      const E = md.ke() + md.pe;
      if (!isFinite(E) || md.temp() > 50) {
        state.blown = true; msg.hidden = false;
        msg.innerHTML = '<b>The simulation blew up!</b><br>The time step was too big: atoms jumped into each other and the forces exploded. Lower Δt and press Reset.';
      } else if (!mini) {
        hist.push([E / md.N, md.ke() / md.N, md.pe / md.N]); if (hist.length > 300) hist.shift();
        if (frame++ % 4 === 0) {
          outs.Tnow.textContent = md.temp().toFixed(2); outs.E.textContent = (E / md.N).toFixed(3) + ' /atom';
          outs.steps.textContent = md.step_n.toLocaleString('en-US'); outs.time.textContent = (md.step_n * state.dt * 2.16).toFixed(1) + ' ps';
          drawChart();
        }
      }
    }
    const vref = 4.2, h = md.L / 2;   // fixed scale: colour shows how hot it is
    for (let i = 0; i < md.N; i++) {
      dummy.position.set(md.x[i * 3] - h, md.x[i * 3 + 1] - h, md.x[i * 3 + 2] - h); dummy.updateMatrix(); atoms.setMatrixAt(i, dummy.matrix);
      const sp = Math.hypot(md.v[i * 3], md.v[i * 3 + 1], md.v[i * 3 + 2]) / vref;
      sp < 0.5 ? c.copy(CYAN).lerp(AMBER, sp * 2) : c.copy(AMBER).lerp(CORAL, Math.min(1, (sp - 0.5) * 2));
      atoms.setColorAt(i, c);
    }
    atoms.instanceMatrix.needsUpdate = true; atoms.instanceColor.needsUpdate = true;
    orbit.update(); renderer.render(scene, cam);
  }
  loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); orbit.dispose(); renderer.dispose(); };
}
