// Lattice builder: SC / BCC / FCC / random starting structures in a periodic box, with overlap and energy checks.
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA' };
const BASIS = { sc: [[0, 0, 0]], bcc: [[0, 0, 0], [0.5, 0.5, 0.5]], fcc: [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]] };
const MAXN = 4 * 6 ** 3;

export function mountLattice(root) {
  const st = { type: 'fcc', n: 3, rho: 0.84, reject: false, seed: 1 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view"><div class="sim-tip">Drag to rotate · <b style="color:${COL.coral}">red</b> = atoms closer than 0.8σ</div></div>
  <div class="sim-ctrl">
    <div class="seg"><button data-t="sc">Simple cubic</button><button data-t="bcc">BCC</button><button data-t="fcc" class="on">FCC</button><button data-t="rand">Random</button></div>
    <label>Unit cells per side <output data-o="n"></output><input type="range" data-k="n" min="2" max="6" step="1" value="3"></label>
    <label>Density ρ* <output data-o="rho"></output><input type="range" data-k="rho" min="0.2" max="1.2" step="0.01" value="0.84"></label>
    <label class="check" data-rand hidden><input type="checkbox" data-k="rej"> Reject overlaps (keep atoms ≥ 0.8σ apart)</label>
    <button class="btn sm ghost" data-k="again" hidden>↻ New random placement</button>
    <div class="readouts">
      <div><span>Atoms N</span><b data-o="N"></b></div>
      <div><span>Box length L</span><b data-o="L"></b></div>
      <div><span>Closest pair</span><b data-o="dmin"></b></div>
      <div><span>Overlapping atoms</span><b data-o="ov"></b></div>
      <div class="wide"><span>LJ energy per atom (<i class="nc">ε</i>)</span><b data-o="E"></b></div>
    </div>
  </div>`;
  const view = root.querySelector('.sim-view'), $ = s => root.querySelector(s), o = {};
  root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2)); renderer.outputColorSpace = THREE.SRGBColorSpace; view.prepend(renderer.domElement);
  const scene = new THREE.Scene(), cam = new THREE.PerspectiveCamera(40, 1, 0.1, 500);
  scene.add(new THREE.HemisphereLight(0xcdefff, 0x1b2550, 1.7)); const dl = new THREE.DirectionalLight(0xffffff, 2.2); dl.position.set(6, 10, 8); scene.add(dl);
  const atoms = new THREE.InstancedMesh(new THREE.SphereGeometry(0.5, 20, 14), new THREE.MeshStandardMaterial({ roughness: 0.35 }), MAXN);
  const boxMat = new THREE.LineBasicMaterial({ color: 0x38e1ff, transparent: true, opacity: 0.6 });
  let box = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(1, 1, 1)), boxMat);
  scene.add(atoms, box);
  const orbit = new OrbitControls(cam, renderer.domElement); orbit.enableDamping = true; orbit.autoRotate = true; orbit.autoRotateSpeed = 0.8; orbit.enablePan = false;
  const resize = () => { const w = view.clientWidth, h = view.clientHeight; renderer.setSize(w, h, false); cam.aspect = w / h; cam.updateProjectionMatrix(); };
  const ro = new ResizeObserver(resize); ro.observe(view); resize();

  function rand(seed) { return () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function build() {
    const n = st.n, perCell = st.type === 'rand' ? 4 : BASIS[st.type].length, N = perCell * n ** 3, L = Math.cbrt(N / st.rho), a = L / n;
    const P = [], kind = []; let failed = false;
    if (st.type !== 'rand') {
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) BASIS[st.type].forEach((b, q) => { P.push([(i + b[0] + 0.25) * a, (j + b[1] + 0.25) * a, (k + b[2] + 0.25) * a]); kind.push(q); });
    } else {
      const R = rand(st.seed); let tries = 0;
      while (P.length < N && tries < 200000) {
        tries++; const p = [R() * L, R() * L, R() * L];
        if (st.reject && P.some(q => { let dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L); return dx * dx + dy * dy + dz * dz < 0.64; })) continue;
        P.push(p); kind.push(0);
      }
      failed = P.length < N;
    }
    // checks: closest pair, overlaps, LJ energy (rc = 2.5, minimum image)
    const M = P.length, over = new Uint8Array(M); let dmin = Infinity, E = 0;
    for (let i = 0; i < M; i++) for (let j = i + 1; j < M; j++) {
      let dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = P[i][2] - P[j][2];
      dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
      const r2 = dx * dx + dy * dy + dz * dz; dmin = Math.min(dmin, r2);
      if (r2 < 0.64) { over[i] = over[j] = 1; }
      if (r2 < 6.25) { const i6 = 1 / (r2 * r2 * r2); E += 4 * (i6 * i6 - i6); }
    }
    // draw
    const d = new THREE.Object3D(), c = new THREE.Color();
    for (let i = 0; i < MAXN; i++) {
      if (i < M) { d.position.set(P[i][0] - L / 2, P[i][1] - L / 2, P[i][2] - L / 2); d.scale.setScalar(over[i] ? 1.05 : 0.85); c.set(over[i] ? COL.coral : kind[i] ? COL.amber : COL.cyan); }
      else d.scale.setScalar(0);
      d.updateMatrix(); atoms.setMatrixAt(i, d.matrix); if (i < M) atoms.setColorAt(i, c);
    }
    atoms.instanceMatrix.needsUpdate = true; if (atoms.instanceColor) atoms.instanceColor.needsUpdate = true;
    scene.remove(box); box.geometry.dispose(); box = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(L, L, L)), boxMat); scene.add(box);
    if (!build.placed || build.lastL === undefined || Math.abs(build.lastL - L) > 0.5) { cam.position.set(L * 1.5, L * 1.0, L * 1.8); orbit.target.set(0, 0, 0); build.lastL = L; build.placed = true; }
    const nOver = over.reduce((s, q) => s + q, 0), ev = E / Math.max(1, M);
    o.n.textContent = st.n; o.rho.textContent = st.rho.toFixed(2); o.N.textContent = failed ? `${M} of ${N}` : N;
    o.L.textContent = L.toFixed(2) + ' σ'; o.dmin.textContent = Math.sqrt(dmin).toFixed(3) + ' σ';
    o.ov.innerHTML = nOver ? `<span style="color:${COL.coral}">${nOver}</span>` : '0';
    o.E.innerHTML = failed ? `<span style="color:${COL.coral}">Too dense to insert atoms at random</span>`
      : Math.abs(ev) > 1e4 ? `<span style="color:${COL.coral}">${ev.toExponential(1)}: the run would explode</span>`
      : ev > 0 ? `<span style="color:${COL.amber}">${ev.toFixed(2)}: positive, minimise first</span>` : `<span style="color:${COL.cyan}">${ev.toFixed(2)}: a good start</span>`;
  }
  let pending = false; const rebuild = () => { if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; build(); }); } };
  root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => {
    st.type = b.dataset.t; root.querySelectorAll('[data-t]').forEach(x => x.classList.toggle('on', x === b));
    $('[data-rand]').hidden = $('[data-k=again]').hidden = st.type !== 'rand'; rebuild();
  });
  $('[data-k=n]').oninput = e => { st.n = +e.target.value; rebuild(); };
  $('[data-k=rho]').oninput = e => { st.rho = +e.target.value; rebuild(); };
  $('[data-k=rej]').onchange = e => { st.reject = e.target.checked; rebuild(); };
  $('[data-k=again]').onclick = () => { st.seed++; rebuild(); };

  let raf = 0, visible = true;
  const io = new IntersectionObserver(([e]) => { visible = e.isIntersecting; }); io.observe(root);
  const loop = () => { raf = requestAnimationFrame(loop); if (!visible) return; orbit.update(); renderer.render(scene, cam); };
  build(); loop();
  return () => { cancelAnimationFrame(raf); ro.disconnect(); io.disconnect(); orbit.dispose(); renderer.dispose(); };
}
