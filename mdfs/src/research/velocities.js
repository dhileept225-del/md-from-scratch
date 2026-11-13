// Research · A.6: velocities. The four steps of velocity create, observed live with velocity histograms and checked by hand.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/src/velocity.cpp';

const script = o => IN_MELT
  .replace('velocity        all create 3.0 87287 loop geom', `velocity        all create ${o.T} ${o.seed} loop geom dist ${o.dist} mom ${o.mom} rot ${o.rot}`)
  .replace('thermo          50\nrun             250', `variable        vx equal vcm(all,x)
variable        vy equal vcm(all,y)
thermo_style    custom step temp ke v_vx v_vy
thermo_modify   norm no
thermo          100
dump            1 all custom 100 v.dump vx
run             100
undump          1`);

function readDump(text) {
  const snaps = []; const L = text.split('\n'); let i = 0;
  while (i < L.length) {
    if (L[i].startsWith('ITEM: TIMESTEP')) { const step = +L[i + 1]; const n = +L[i + 3]; const k = L.findIndex((l, j) => j > i && l.startsWith('ITEM: ATOMS')); snaps.push({ step, v: L.slice(k + 1, k + 1 + n).map(Number) }); i = k + 1 + n; } else i++;
  }
  return snaps;
}

function histo(canvas, snaps, T, T2) {
  const d = Math.min(devicePixelRatio || 1, 2), W = canvas.clientWidth, H = canvas.clientHeight; canvas.width = W * d; canvas.height = H * d;
  const c = canvas.getContext('2d'); c.setTransform(d, 0, 0, d, 0, 0); c.clearRect(0, 0, W, H);
  const vmax = 6, nb = 48, x0 = 40, x1 = W - 10, y0 = 14, y1 = H - 28, X = v => x0 + (v + vmax) / (2 * vmax) * (x1 - x0);
  const sets = snaps.map(s => { const h = new Array(nb).fill(0); s.v.forEach(v => { const b = Math.floor((v + vmax) / (2 * vmax) * nb); if (b >= 0 && b < nb) h[b]++; }); return h.map(x => x / (s.v.length * 2 * vmax / nb)); });
  const gT = t => v => Math.exp(-v * v / (2 * t)) / Math.sqrt(2 * Math.PI * t), g = gT(T);
  const ymax = Math.max(...sets.flat(), g(0), T2 ? gT(T2)(0) : 0) * 1.4, Y = p => y1 - p / ymax * (y1 - y0);
  c.strokeStyle = 'rgba(143,163,199,.6)'; c.beginPath(); c.moveTo(x0, y0); c.lineTo(x0, y1); c.lineTo(x1, y1); c.stroke();
  c.fillStyle = '#8FA3C7'; c.font = '11px "Space Grotesk", sans-serif'; c.textAlign = 'center';
  for (let v = -6; v <= 6; v += 2) c.fillText(String(v), X(v), y1 + 14);
  c.fillText('vx', (x0 + x1) / 2, y1 + 26);
  const cols = ['#FFC24B', '#38E1FF'];
  sets.forEach((h, k) => { c.strokeStyle = cols[k]; c.lineWidth = 2; c.beginPath(); h.forEach((p, b) => { const xa = X(-vmax + b * 2 * vmax / nb), xb = X(-vmax + (b + 1) * 2 * vmax / nb); b ? c.lineTo(xa, Y(p)) : c.moveTo(xa, Y(p)); c.lineTo(xb, Y(p)); }); c.stroke(); });
  c.strokeStyle = 'rgba(242,245,250,.85)'; c.setLineDash([4, 4]); c.lineWidth = 1.5; c.beginPath();
  for (let k = 0; k <= 200; k++) { const v = -vmax + k * 2 * vmax / 200; k ? c.lineTo(X(v), Y(g(v))) : c.moveTo(X(v), Y(g(v))); } c.stroke();
  if (T2) { const g2 = gT(T2); c.strokeStyle = 'rgba(56,225,255,.9)'; c.beginPath(); for (let k = 0; k <= 200; k++) { const v = -vmax + k * 2 * vmax / 200; k ? c.lineTo(X(v), Y(g2(v))) : c.moveTo(X(v), Y(g2(v))); } c.stroke(); }
  c.setLineDash([]);
  c.textAlign = 'left'; c.fillStyle = cols[0]; c.fillText(`step ${snaps[0].step}`, x0 + 8, y0 + 4); if (snaps[1]) { c.fillStyle = cols[1]; c.fillText(`step ${snaps[1].step}`, x0 + 70, y0 + 4); }
  c.fillStyle = '#F2F5FA'; const lab = `dashed: Maxwell, T = ${T}`; c.fillText(lab, x0 + 8, y0 + 18); if (T2) { c.fillStyle = cols[1]; c.fillText(`and T = ${T2.toFixed(2)}`, x0 + 14 + c.measureText(lab).width, y0 + 18); }
}

export default {
  id: 'velocities',
  title: 'Velocities',
  minutes: 20,
  summary: 'One line, velocity all create 3.0 87287, gives 4000 atoms their starting velocities in four steps: random numbers, no net drift, optional no rotation, and an exact rescale to the target temperature. You read those steps in the source, check the kinetic energy by hand, and watch the velocity distribution change shape in the first 100 steps.',
  learn: [
    'The four steps of velocity create, in velocity.cpp',
    'Why the kinetic energy at step 0 is exactly (3N − 3)/2 × T, and what mom no does',
    'Watch uniform random velocities become a Maxwell distribution through collisions',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: four steps in one line</h2>
    <div class="tbl wide"><table>${head(['Step', 'What happens', 'Keyword (default)', 'velocity.cpp'])}
${row(['1. random', 'each atom gets three random numbers, from a flat range −½…½ or a bell curve, divided by √m', '<code>dist uniform</code> or <code>gaussian</code> (uniform)', '<a href="' + GH + '#L328-L353" target="_blank" rel="noopener">lines 328–353</a>'])}
${row(['2. no drift', 'subtract the average velocity, so the whole crystal does not fly away', '<code>mom yes</code> (yes)', 'line 357'])}
${row(['3. no spin', 'remove any overall rotation', '<code>rot yes</code> (no)', 'line 358'])}
${row(['4. rescale', 'measure the temperature, multiply every velocity by √(T<sub>target</sub>/T<sub>measured</sub>)', '', 'lines 364–368'])}</table></div>
    <p>With <code>loop geom</code> (as in melt) the random numbers of each atom are seeded from the seed and the atom's own position (<code>random->reset(seed,x[i])</code>, line 334). The same atom therefore gets the same velocity whatever the number of processors, but a change of units changes its velocity (page A.4).</p>
    <div class="calc">  if (momentum_flag) zero_momentum();       <b>step 2</b>
  if (rotation_flag) zero_rotation();       <b>step 3</b>
  t = temperature->compute_scalar();        <b>step 4: measure T</b>
  rescale(t,t_desired);                     <b>multiply v by √(t_desired/t)</b></div>
  </section>

  <section class="block prose">
    <h2>2 · Predict by hand: the kinetic energy at step 0</h2>
    <p>LAMMPS defines the temperature from the kinetic energy and the number of degrees of freedom. For N atoms in 3D with the total momentum removed, there are 3N − 3 of them:</p>
    <div class="calc">T = 2 KE / (dof × k_B)        dof = 3N − 3 = 3 × 4000 − 3 = 11997
after step 4, T = 3 exactly, so
KE = dof × T / 2 = 11997 × 3 / 2 = <b>17995.5</b>          (total, in ε)
the velocity of the centre of mass after step 2: <b>0</b> (to rounding, about 10⁻¹⁷)</div>
    <p>With <code>mom no</code>, step 2 is skipped: the centre of mass keeps a small random velocity. No internal force can change it (forces between atoms cancel in pairs), so the prediction is that it stays exactly the same at every later step.</p>
  </section>

  <section class="block">
    <h2>3 · Modify and observe</h2>
    <p class="muted">Choose the options, run 100 steps of melt, and compare the vx distribution at step 0 with the Maxwell bell curve at T = 3, and at step 100 with the Maxwell curve at the temperature LAMMPS prints for step 100.</p>
    <div class="lmp" id="vel">
      <div class="row">
        <label class="rlab">dist <select class="r-in o-dist"><option>uniform</option><option>gaussian</option></select></label>
        <label class="rlab">mom <select class="r-in o-mom"><option>yes</option><option>no</option></select></label>
        <label class="rlab">rot <select class="r-in o-rot"><option>no</option><option>yes</option></select></label>
        <label class="rlab">seed <input class="r-in o-seed" value="87287"></label>
        <button class="btn sm primary" data-k="run">▶ Run 100 steps</button><span class="lmp-status muted v-st"></span>
      </div>
      <canvas class="vhist"></canvas>
      <div class="lmp-compare v-cmp"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>What the runs show</h2>
    <ul>
      <li><b>Every option gives KE = 17995.5 at step 0:</b> step 4 always rescales to the exact target.</li>
      <li><b>uniform</b> starts flat, clearly not the Maxwell curve. After 100 steps (0.5 τ) collisions have already turned it into a bell curve. It is narrower than the T = 3 curve, because melting has turned kinetic energy into potential energy: the temperature has dropped to about 1.67, and the step-100 histogram follows the Maxwell curve for that temperature. The starting distribution is forgotten quickly; the temperature it sets is not.</li>
      <li><b>mom no</b>: the centre of mass moves at a small constant velocity, identical at steps 0 and 100 to every digit. The whole system drifts; with a thermostat this drift can grow into the "flying ice cube" of Lab 3.7. Keep the default.</li>
      <li><b>A different seed</b>: the same KE at step 0, different numbers at step 100. Runs with several seeds are the honest way to estimate error bars (Lab 6.8).</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> What KE does LAMMPS print at step 0 for 500 atoms created at T = 1.5?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>dof = 3 × 500 − 3 = 1497, KE = 1497 × 1.5 / 2 = <b>1122.75</b> (total, with <code>thermo_modify norm no</code>).</div></div>
    <div class="ex"><p><b>2.</b> In a 2D run with 9 atoms (page A.3) at T = 0.72, what is the total KE at step 0?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>dof = 2N − 2 = 16, KE = 16 × 0.72 / 2 = <b>5.76</b>, or 0.64 per atom.</div></div>
    <div class="ex"><p><b>3.</b> Why are the random numbers divided by √m before rescaling?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>So that light and heavy atoms start with the same average kinetic energy ½mv² (equipartition), which is what a temperature means. Without it, light atoms (H) would start too cold and heavy atoms too hot.</div></div>
  </section>`,

  mount(root) {
    const V = root.querySelector('#vel'), cv = V.querySelector('.vhist'), st = V.querySelector('.v-st');
    let last = null;
    const ro = new ResizeObserver(() => { if (last) histo(cv, last.snaps, 3, last.T2); }); ro.observe(cv);
    V.querySelector('[data-k=run]').onclick = async e => {
      const o = { T: 3.0, seed: parseInt(V.querySelector('.o-seed').value) || 87287, dist: V.querySelector('.o-dist').value, mom: V.querySelector('.o-mom').value, rot: V.querySelector('.o-rot').value };
      const b = e.currentTarget; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(script(o), { read: ['v.dump'] }); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = `dist ${o.dist} · mom ${o.mom} · rot ${o.rot} · seed ${o.seed}`;
      const snaps = res.files['v.dump'] ? readDump(res.files['v.dump']) : [];
      last = { snaps }; 
      const { cols, rows } = res.thermo, r0 = rows.find(w => w[0] === '0'), r1 = rows.find(w => w[0] === '100'), C = k => cols.indexOf(k);
      last.T2 = r1 ? +r1[C('Temp')] : 0; if (snaps.length) histo(cv, snaps, 3, last.T2);
      const m = matchMarkup(g8(17995.5), r0[C('KinEng')]);
      const drift = o.mom === 'no';
      const mv = drift ? matchMarkup(r0[C('v_vx')], r1[C('v_vx')]) : null;
      V.querySelector('.v-cmp').innerHTML = `<div class="tbl"><table><tr><th></th><th>step 0</th><th>step 100</th><th>by hand</th></tr>
        <tr><td>Temp</td><td class="mono">${r0[C('Temp')]}</td><td class="mono">${r1[C('Temp')]}</td><td>3 at step 0</td></tr>
        <tr><td>KinEng (total)</td><td class="mono">${r0[C('KinEng')]}</td><td class="mono">${r1[C('KinEng')]}</td><td class="mono">${m.html} ${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>
        <tr><td>centre-of-mass vx</td><td class="mono">${r0[C('v_vx')]}</td><td class="mono">${r1[C('v_vx')]}</td><td>${drift ? `unchanged: <span class="mono">${mv.html}</span> ${mv.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}` : '0 (rounding only)'}</td></tr>
        <tr><td>centre-of-mass vy</td><td class="mono">${r0[C('v_vy')]}</td><td class="mono">${r1[C('v_vy')]}</td><td>${drift ? 'unchanged' : '0 (rounding only)'}</td></tr></table></div>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => ro.disconnect();
  },
};
