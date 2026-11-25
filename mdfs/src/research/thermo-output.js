// Research · D.12: thermo output, temperature and pressure. Four temperature computes and their degrees of freedom, the pressure scalar and all six tensor components, and thermo_modify temp, matched by hand on ten drifting LJ atoms.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const N = 10, L = 20.0;
const X = Array.from({ length: N }, (_, i) => [8.6 + 1.12 * (i % 3) + 0.05 * Math.sin(i), 8.8 + 1.12 * Math.floor(i / 3) % 3 + 0.04 * Math.cos(2 * i), 9.5 + 1.1 * Math.floor(i / 9) + 0.03 * i].map(v => +v.toFixed(4)));
const TYPE = Array.from({ length: N }, (_, i) => (i % 2) + 1), MASS = [0, 1.0, 2.0];
const V = Array.from({ length: N }, (_, i) => [0, 1, 2].map(k => +(0.9 * Math.sin(2.1 * i + 1.3 * k) + 0.2 * Math.cos(i * k)).toFixed(4)));
const DATA = drift => `ten atoms

${N} atoms
2 atom types

0 ${L} xlo xhi
0 ${L} ylo yhi
0 ${L} zlo zhi

Masses

1 1.0
2 2.0

Atoms # atomic

${X.map((x, i) => `${i + 1} ${TYPE[i]} ${x.join(' ')}`).join('\n')}

Velocities

${V.map((v, i) => `${i + 1} ${+(v[0] + drift).toFixed(4)} ${v[1]} ${v[2]}`).join('\n')}
`;
const script = tm => `units           lj
atom_style      atomic
read_data       ten.data
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0
compute         tcom all temp/com
compute         tpart all temp/partial 0 1 1
compute         tfull all temp
compute_modify  tfull extra/dof 0
compute         pcom all pressure tcom
${tm === 2 ? 'thermo_modify   temp tcom\n' : ''}thermo_style    custom step temp c_tcom c_tpart c_tfull press c_pcom pxx pyy pzz pxy pxz pyz
thermo_modify   ${tm === 1 ? 'temp tcom ' : ''}norm no format float %.17g
run             0`;

function hand(drift, tm) {
  const v = V.map(q => [+(q[0] + drift).toFixed(4), q[1], q[2]]), m = TYPE.map(t => MASS[t]), x = X;
  // pair forces and fdotr virial
  const f = x.map(() => [0, 0, 0]);
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { const d = [0, 1, 2].map(k => x[i][k] - x[j][k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; if (rsq < 6.25) { const r2 = 1 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48 * r6 - 24) * r2; for (let k = 0; k < 3; k++) { f[i][k] += d[k] * fp; f[j][k] -= d[k] * fp; } } }
  const W = [0, 0, 0, 0, 0, 0]; for (let i = 0; i < N; i++) { W[0] += f[i][0] * x[i][0]; W[1] += f[i][1] * x[i][1]; W[2] += f[i][2] * x[i][2]; W[3] += f[i][1] * x[i][0]; W[4] += f[i][2] * x[i][0]; W[5] += f[i][2] * x[i][1]; }
  let mt = 0; for (const q of m) mt += q;
  const p = [0, 0, 0]; v.forEach((q, i) => { for (let k = 0; k < 3; k++) p[k] += q[k] * m[i]; }); const vcm = p.map(q => q / mt);
  const sum = g => { let t = 0; v.forEach((q, i) => { t += g(q, i); }); return t; };
  const dof = { temp: 3 * N - 3, tcom: 3 * N - 3, tpart: 2 * N - (1.0 * 2 / 3) * 3, tfull: 3 * N };
  const T = {
    temp: sum((q, i) => (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) * m[i]) * (1 / (dof.temp * 1)),
    tcom: sum((q, i) => { const a = q[0] - vcm[0], b = q[1] - vcm[1], c = q[2] - vcm[2]; return (a * a + b * b + c * c) * m[i]; }) * (1 / (dof.tcom * 1)),
    tpart: sum((q, i) => (0 * q[0] * q[0] + 1 * q[1] * q[1] + 1 * q[2] * q[2]) * m[i]) * (1 / (dof.tpart * 1)),
    tfull: sum((q, i) => (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) * m[i]) * (1 / (dof.tfull * 1)),
  };
  const ke = [0, 0, 0, 0, 0, 0], kec = [0, 0, 0, 0, 0, 0];
  v.forEach((q, i) => { const mo = m[i]; ke[0] += mo * q[0] * q[0]; ke[1] += mo * q[1] * q[1]; ke[2] += mo * q[2] * q[2]; ke[3] += mo * q[0] * q[1]; ke[4] += mo * q[0] * q[2]; ke[5] += mo * q[1] * q[2];
    const w = [q[0] - vcm[0], q[1] - vcm[1], q[2] - vcm[2]]; kec[0] += mo * w[0] * w[0]; kec[1] += mo * w[1] * w[1]; kec[2] += mo * w[2] * w[2]; kec[3] += mo * w[0] * w[1]; kec[4] += mo * w[0] * w[2]; kec[5] += mo * w[1] * w[2]; });
  const inv = 1.0 / (L * L * L), Tth = tm ? T.tcom : T.temp, kt = tm ? kec : ke;
  const press = (dof.temp * 1 * Tth + W[0] + W[1] + W[2]) / 3.0 * inv * 1;
  const pcom = (dof.tcom * 1 * T.tcom + W[0] + W[1] + W[2]) / 3.0 * inv * 1;
  const tens = kt.map((q, i) => (q + W[i]) * inv * 1);
  return { T, dof, press, pcom, tens, vcm, W };
}
const close = (h, l) => h === l ? 'exact' : Math.abs(h - l) <= 1e-13 * Math.max(Math.abs(h), Math.abs(l), 1e-6) ? 'near' : 'no';

export default {
  id: 'thermo-output',
  title: 'Thermo output: temperature and pressure',
  minutes: 20,
  summary: 'The temp and press columns are computes, chosen by thermo_style and changeable with thermo_modify. On ten LJ atoms of two masses, flowing along x, you compute four temperatures (temp, temp/com, temp/partial, and temp without the 3 removed degrees of freedom), the pressure scalar, a second pressure built on temp/com, and all six tensor components by hand, then watch a drift of the whole cluster heat the temp column but not temp/com, and move the pressure with it.',
  learn: [
    'How each temperature compute counts degrees of freedom, and extra/dof',
    'temp/com removes the centre-of-mass velocity; temp/partial drops directions',
    'Pressure tensor = kinetic tensor + virial tensor; the scalar is its trace / 3',
    'thermo_modify temp: the thermo pressure uses that temperature too',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the computes behind the columns</h2>
    <div class="prose"><p>thermo creates <code>thermo_temp</code> (compute temp) and <code>thermo_press</code> (compute pressure using thermo_temp). Degrees of freedom, ${gh('src/compute_temp.cpp', 58)}:</p>
<pre class="lmp-mini">temp, temp/com:   dof = 3N − extra_dof − fix_dof            // extra_dof = 3 by default
temp/partial 0 1 1: dof = 2N − (2/3)(extra_dof + fix_dof)      // compute_temp_partial.cpp
T = mvv2e Σ m v² / (dof k_B)       temp/com: v → v − v_cm (mass-weighted)
P = (dof k_B T + Wxx + Wyy + Wzz) / 3V · nktv2p                // compute_pressure.cpp:256
P_ab = (Σ m v_a v_b · mvv2e + W_ab) / V · nktv2p               // compute_vector: ab = xx yy zz xy xz yz</pre>
<p>W is the virial Σ x·f (page C.4), <code>W_xy = Σ f_y x_x</code> (${gh('src/pair.cpp', 1765)}). Fixes such as shake or rigid add to fix_dof; <code>compute_modify ID extra/dof 0</code> stops the 3 from being removed.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: ten atoms, eleven numbers</h2>
    <div class="prose"><p>Ten LJ atoms, masses 1 and 2 alternating, close enough to interact, <code>run 0</code>. Add a drift to every vx, and choose whether thermo uses temp or temp/com.</p></div>
    <div class="row">
      <label class="lc-in">drift added to vx<select class="lc-u to-d"><option>0</option><option selected>1.5</option><option>5</option></select></label>
      <label class="lc-in">thermo temperature<select class="lc-u to-t"><option value="0">temp (default)</option><option value="1">thermo_modify temp tcom, after thermo_style</option><option value="2">thermo_modify temp tcom, before thermo_style</option></select></label>
      <button class="btn sm primary" data-k="a">▶ run 0</button><span class="lmp-status muted to-sta"></span></div>
    <div class="calc to-pred"></div>
    <div class="tbl to-out"></div>
    <div class="prose to-sum"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe</h2>
    <ul>
      <li><b>A drift is not heat.</b> Adding 1.5 to every vx raises <code>temp</code> a lot and leaves <code>temp/com</code> unchanged: a flowing liquid or a moving droplet reads hot unless you remove the bias (temp/com, temp/profile, temp/partial).</li>
      <li><b>The pressure follows the chosen temperature.</b> With <code>thermo_modify temp tcom</code> the press column and pxx … pyz use the drift-free kinetic part. The virial part does not change.</li>
      <li><b>Order matters.</b> Put <code>thermo_modify temp</code> <i>before</i> <code>thermo_style</code> and you get a mixed state: the temp column back on compute temp, but thermo_press still on tcom. LAMMPS prints "New thermo_style command, previous thermo_modify settings will be lost"; the pressure part was not lost.</li>
      <li><b>(pxx + pyy + pzz)/3 = press</b> (to round-off): the scalar is the trace. The off-diagonal pxy, pxz, pyz are what shear viscosity (Green–Kubo, lesson 6) is computed from.</li>
      <li><b>Small systems feel extra/dof.</b> With N = 10, removing 3 of 30 degrees of freedom changes T by 11 %; at N = 4000 by 0.025 %.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Ten atoms, temp/partial 0 1 1: how many degrees of freedom?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>2 × 10 − (2/3) × 3 = <b>18</b>.</div></div>
    <div class="ex"><p><b>2.</b> A water box with fix shake on all O–H bonds and H–O–H angles: what does the temperature compute subtract?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>3 per molecule (two bonds and one angle constraint) as fix_dof, plus extra_dof = 3. For 1000 molecules: 9000 − 3000 − 3 = <b>5997</b>.</div></div>
    <div class="ex"><p><b>3.</b> Why can the pressure of a small cluster in a big empty box be negative?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The virial of attracting atoms (r &gt; 2^(1/6)σ) is negative and can outweigh the kinetic term; divided by the whole box volume it is a small negative number. Pressure of a cluster in vacuum is not a property of the cluster.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pr = () => { const d = +$('.to-d').value, tm = +$('.to-t').value, H = hand(d, tm); $('.to-pred').textContent = `v_cm = (${H.vcm.map(q => q.toFixed(6)).join(', ')})   total mass 15\ndof: temp ${H.dof.temp}, temp/com ${H.dof.tcom}, temp/partial ${H.dof.tpart}, temp with extra/dof 0 ${H.dof.tfull}\nvirial W = (${H.W.slice(0, 3).map(q => q.toFixed(6)).join(', ')} | ${H.W.slice(3).map(q => q.toFixed(6)).join(', ')})`; $('.to-out').innerHTML = ''; $('.to-sum').innerHTML = ''; return H; };
    root.querySelectorAll('.to-d,.to-t').forEach(e => e.onchange = pr); pr();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.to-sta'), d = +$('.to-d').value, tm = +$('.to-t').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(script(tm), { files: { 'ten.data': DATA(d) } }); b.disabled = false;
      if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = pr(), c = res.thermo.cols, r0 = res.thermo.rows[0], g = k => parseFloat(r0[c.indexOf(k)]);
      const items = [['temp', H.T[tm === 1 ? 'tcom' : 'temp'], g('Temp')], ['c_tcom', H.T.tcom, g('c_tcom')], ['c_tpart', H.T.tpart, g('c_tpart')], ['c_tfull', H.T.tfull, g('c_tfull')], ['press', H.press, g('Press')], ['c_pcom', H.pcom, g('c_pcom')], ...['Pxx', 'Pyy', 'Pzz', 'Pxy', 'Pxz', 'Pyz'].map((k, i) => [k.toLowerCase(), H.tens[i], g(k)])];
      let ex = 0, nr = 0;
      $('.to-out').innerHTML = `<table>${head(['column', 'by hand', 'LAMMPS', ''])}${items.map(([k, h, l]) => { const s = close(h, l); if (s === 'exact') ex++; else if (s === 'near') nr++; return row([k, `<span class="mono">${s17(h)}</span>`, `<span class="mono">${s17(l)}</span>`, s === 'exact' ? `${OK} bit for bit` : s === 'near' ? `${OK} to 13 digits` : BAD]); }).join('')}</table>`;
      const tr = (H.tens[0] + H.tens[1] + H.tens[2]) / 3;
      $('.to-sum').innerHTML = `<p>${ex + nr} of ${items.length} match (${ex} bit for bit${nr ? `; ${nr} differ in the last digits, where LAMMPS adds the pair forces in neighbour-list order and the hand code in i &lt; j order` : ''}). (pxx + pyy + pzz)/3 = ${s17(tr)} vs press ${s17(H.press)}.${tm === 2 ? ' <b>Look at temp and press:</b> the temp column is the plain temperature, but press uses temp/com. thermo_style made a new Thermo and dropped the earlier thermo_modify (LAMMPS warns), while thermo_press, which already existed, kept tcom.' : ''}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
