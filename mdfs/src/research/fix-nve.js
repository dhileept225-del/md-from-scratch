// Research · D.9: the integrator as a Fix. FixNVE's constructor, init, reset_dt and group test; fix nve/limit and fix dt/reset reproduced bit for bit through a hard collision of two atoms.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const MODES = {
  nve: { fix: 'fix             1 all nve', desc: 'plain velocity Verlet' },
  limit: { fix: 'fix             1 all nve/limit 0.01', desc: 'no atom moves more than 0.01 per step' },
  group: { fix: 'group           one id 1\nfix             1 one nve', desc: 'only atom 1 is integrated' },
  noforce: { fix: 'fix             1 all nve/noforce', desc: 'x += dt·v, forces ignored' },
  dtreset: { fix: 'fix             1 all nve\nfix             2 all dt/reset 1 NULL NULL 0.01', desc: 'dt chosen each step so no atom moves more than 0.01' },
};
const script = (m, x2, v) => `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10 units box
create_box      1 box
create_atoms    1 single 5.0 5.0 5.0 units box
create_atoms    1 single ${x2} 5.0 5.0 units box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
set             atom 1 vx ${v}
set             atom 2 vx ${-v}
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
${MODES[m].fix}
thermo_style    custom step dt time pe ke etotal${m === 'limit' ? ' f_1' : ''}
thermo_modify   norm no format float %.17g
thermo          1
dump            d all custom 1 c.dump id x vx fx
dump_modify     d format float %.17g sort id
run             12`;

function force(x) {
  const f = [0, 0], d = x[0] - x[1], rsq = d * d + 0 + 0; let e = 0;
  if (rsq < 6.25) { const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48.0 * r6 - 24.0) * r2; f[0] += d * fp; f[1] -= d * fp; e = r6 * (4.0 * r6 - 4.0); }
  return { f, e };
}
function hand(m, x2, v0, n = 12) {
  const x = [5.0, x2], v = [v0, -v0]; let F = force(x), dt = 0.005, atime = 0, atimestep = 0, ncount = 0;
  const xmax = 0.01, out = [];
  const integ = i => m !== 'group' || i === 0;
  const dtreset = step => {   // FixDtReset::end_of_step, also called in setup
    let dtmin = 1e20;
    for (let i = 0; i < 2; i++) {
      const massinv = 1.0 / 1.0, vsq = v[i] * v[i] + 0 + 0, fsq = F.f[i] * F.f[i] + 0 + 0; let dtv = 1e20, dtf = 1e20;
      if (vsq > 0.0) dtv = xmax / Math.sqrt(vsq);
      if (fsq > 0.0) dtf = Math.sqrt(2.0 * xmax / (1.0 * Math.sqrt(fsq) * massinv));
      let d = Math.min(dtv, dtf); const dtsq = d * d;
      const delx = d * v[i] + 0.5 * dtsq * massinv * F.f[i] * 1.0, dely = d * 0 + 0.5 * dtsq * massinv * 0 * 1.0, delz = dely;
      const delr = Math.sqrt(delx * delx + dely * dely + delz * delz); if (delr > xmax) d *= xmax / delr;
      dtmin = Math.min(dtmin, d);
    }
    if (dtmin !== dt) { atime += (step - atimestep) * dt; atimestep = step; dt = dtmin; }
  };
  const rec = s => { let t = 0; for (let i = 0; i < 2; i++) t += (v[i] * v[i] + 0 + 0) * 1.0; const T = t * (1.0 / (3.0 * 1.0)); let ke = T; ke *= 0.5 * 3.0 * 1.0; out.push({ s, x: [...x], v: [...v], f: [...F.f], pe: 0.0 + (F.e + 0.0), ke, et: (0.0 + (F.e + 0.0)) + ke, dt, time: atime + (s - atimestep) * dt, nc: ncount }); };
  if (m === 'dtreset') dtreset(0);
  const vl2 = (xmax / dt) * (xmax / dt);
  rec(0);
  for (let s = 1; s <= n; s++) {
    const dtf = 0.5 * dt * 1.0, dtfm = dtf / 1.0;
    for (let i = 0; i < 2; i++) if (integ(i)) {
      if (m !== 'noforce') v[i] += dtfm * F.f[i];
      if (m === 'limit') { const vsq = v[i] * v[i] + 0 + 0; if (vsq > vl2) { ncount++; const sc = Math.sqrt(vl2 / vsq); v[i] *= sc; } }
      x[i] += dt * v[i];
    }
    F = force(x);
    for (let i = 0; i < 2; i++) if (integ(i) && m !== 'noforce') {
      v[i] += dtfm * F.f[i];
      if (m === 'limit') { const vsq = v[i] * v[i] + 0 + 0; if (vsq > vl2) { ncount++; const sc = Math.sqrt(vl2 / vsq); v[i] *= sc; } }
    }
    if (m === 'dtreset') dtreset(s);
    rec(s);
  }
  return out;
}

export default {
  id: 'fix-nve',
  title: 'The integrator: fix nve',
  minutes: 20,
  summary: 'fix nve is a Fix like any other: a constructor, a mask, init() that turns dt into dtv and dtf, a group test inside each loop, and reset_dt() for when somebody changes dt. You drive two atoms into a hard collision and reproduce, bit for bit, five integrators: plain nve, nve/limit (and its count of capped atoms), nve on one atom only, nve/noforce, and fix dt/reset choosing a new time step every step.',
  learn: [
    'FixNVE: setmask, init (dtv, dtf), mass or rmass, the group bit',
    'fix nve/limit: velocities capped at xmax/dt, before the drift and after the second kick',
    'fix dt/reset: dt from velocity, force and xmax, and reset_dt() for every fix',
    'Which integrators conserve energy, and which cannot',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a Fix that moves atoms</h2>
    <div class="prose"><p>${gh('src/fix_nve.cpp', 29)}: the constructor only sets <code>time_integrate = 1</code> (so LAMMPS can warn about two integrators on one atom). ${gh('src/fix_nve.cpp', 56)}:</p>
<pre class="lmp-mini">void FixNVE::init()
{
  dtv = update-&gt;dt;
  dtf = 0.5 * update-&gt;dt * force-&gt;ftm2v;
}
void FixNVE::initial_integrate(int)
{
  ...
  if (igroup == atom-&gt;firstgroup) nlocal = atom-&gt;nfirst;
  if (rmass) { ... dtfm = dtf / rmass[i]; ... }         // per-atom mass
  else for (int i = 0; i &lt; nlocal; i++)
    if (mask[i] &amp; groupbit) {                            // only atoms in the fix's group
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];  ...
      x[i][0] += dtv * v[i][0];   ...
    }
}
void FixNVE::reset_dt() { dtv = update-&gt;dt; dtf = 0.5 * update-&gt;dt * force-&gt;ftm2v; }   // 172</pre>
<p>${gh('src/fix_nve_limit.cpp', 83)} adds one test after each kick, with vlimitsq = (xmax/dt)² (line 66):</p>
<pre class="lmp-mini">vsq = v·v;
if (vsq &gt; vlimitsq) { ncount++; scale = sqrt(vlimitsq/vsq); v *= scale; }</pre>
<p>${gh('src/fix_dt_reset.cpp', 136)} runs at end_of_step (and in setup): for each atom</p>
<pre class="lmp-mini">dtv = xmax / |v|;   dtf = sqrt(2 xmax / (ftm2v |f| / m));   dt = min(dtv, dtf);
del = dt v + ½ dt² f/m · ftm2v;   if (|del| &gt; xmax) dt *= xmax / |del|;
dt = min over atoms;  update-&gt;dt = dt;  every fix-&gt;reset_dt();</pre></div>
  </section>

  <section class="block">
    <h2>2 · Match: a hard collision, five integrators</h2>
    <div class="prose"><p>Two LJ atoms on the x axis, coming together at ±vx. Start them close and fast: the repulsion is steep, and dt = 0.005 is too long for it.</p></div>
    <div class="row">
      <label class="lc-in">integrator<select class="lc-u fn-m">${Object.entries(MODES).map(([k, v]) => `<option value="${k}"${k === 'limit' ? ' selected' : ''}>${esc(v.fix.split('\n').pop().replace(/^fix\s+\d+\s+\S+\s+/, ''))}</option>`).join('')}</select></label>
      <label class="lc-in">x₂<select class="lc-u fn-x"><option>6.0</option><option selected>5.95</option><option>5.9</option></select></label>
      <label class="lc-in">vx<select class="lc-u fn-v"><option>1.0</option><option selected>3.0</option><option>6.0</option></select></label>
      <button class="btn sm primary" data-k="a">▶ Run 12 steps</button><span class="lmp-status muted fn-sta"></span></div>
    <pre class="lmp-mini fn-scr"></pre>
    <div class="tbl fn-out"></div>
    <div class="prose fn-sum"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe</h2>
    <ul>
      <li><b>nve/limit is not energy-conserving.</b> Each cap throws kinetic energy away, so TotEng drops; it is meant for relaxing overlaps at the start, not for production. Its <code>f_1</code> counts caps since the run began, both kicks included.</li>
      <li><b>One atom in the group.</b> Atom 2 keeps its x and v to the last bit while its force grows: forces are computed for everyone, only the fix's group moves. thermo temp still counts atom 2.</li>
      <li><b>dt/reset keeps the physics, changes the clock.</b> The trajectory is velocity Verlet with a different dt each step; "time" is summed step by step (<code>update_time</code>), so it is not step × dt any more. In the hard collisions (vx = 3, 6) it picks dt well below 0.005 and TotEng drifts several times less than with plain nve; in the slow one (vx = 1) it picks a <i>larger</i> dt (0.01) and does worse. xmax sets the accuracy, not dt.</li>
      <li><b>nve/noforce</b> moves atoms in straight lines: here they run into each other and PotEng explodes. It is for atoms whose velocities another fix or command sets.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> With nve/limit 0.01 and dt = 0.005, what is the speed cap?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>vlimit = xmax/dt = 0.01 / 0.005 = <b>2</b> (LJ units).</div></div>
    <div class="ex"><p><b>2.</b> At step 0 the atoms are at rest-ish with |f| = 2000. What dt does dt/reset 0.01 choose from the force alone?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>dtf = √(2 × 0.01 / 2000) = <b>0.00316</b>; the velocity term and the combined displacement test can make it smaller.</div></div>
    <div class="ex"><p><b>3.</b> You run <code>fix 1 all nve</code> and <code>fix 2 water nve</code>. What happens to water atoms, and does LAMMPS stop you?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>They are moved twice per step (both fixes' kicks and drifts). LAMMPS only warns ("One or more atoms are time integrated more than once"), using the time_integrate flag set in the constructor.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const show = () => { $('.fn-scr').textContent = MODES[$('.fn-m').value].fix; $('.fn-out').innerHTML = ''; $('.fn-sum').innerHTML = ''; };
    $('.fn-m').onchange = show; show();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.fn-sta'), m = $('.fn-m').value, x2 = +$('.fn-x').value, v = +$('.fn-v').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(script(m, x2, v), { read: ['c.dump'] }); b.disabled = false;
      const t = res.files?.['c.dump']; if (res.error || !t || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = hand(m, x2, v), c = res.thermo.cols, g = (r, k) => parseFloat(r[c.indexOf(k)]);
      const D = {}; for (const f of t.split('ITEM: TIMESTEP').slice(1)) { const ls = f.trim().split('\n'); D[+ls[0]] = ls.slice(-2).map(l => l.trim().split(/\s+/).map(Number)); }
      let n = 0, ok = 0;
      const rows = res.thermo.rows.map(r => {
        const s = +r[0], h = H[s], d = D[s], vals = [[h.x[0], d[0][1]], [h.x[1], d[1][1]], [h.v[0], d[0][2]], [h.v[1], d[1][2]], [h.f[0], d[0][3]], [h.dt, g(r, 'Dt')], [h.time, g(r, 'Time')], [h.pe, g(r, 'PotEng')], [h.ke, g(r, 'KinEng')]];
        if (m === 'limit') vals.push([h.nc, g(r, 'f_1')]);
        const same = vals.map(([a, bb]) => a === bb); same.forEach(x => { n++; if (x) ok++; });
        return row([s, s17(h.x[0]), s17(h.v[0]), s17(h.f[0]), m === 'dtreset' ? s17(h.dt) : h.dt, h.et.toPrecision(10), m === 'limit' ? h.nc : '', same.every(Boolean) ? `${OK} ${vals.length} of ${vals.length}` : `${BAD} ${same.filter(x => !x).length} differ`]);
      });
      $('.fn-out').innerHTML = `<table>${head(['step', 'x₁', 'vx₁', 'fx₁', 'dt', 'TotEng', m === 'limit' ? 'caps' : '', 'all printed numbers'])}${rows.join('')}</table>`;
      const e0 = H[0].et, e12 = H[12].et;
      $('.fn-sum').innerHTML = `<p>${ok} of ${n} printed numbers identical to your hand values (x and v of both atoms, fx, dt, time, PotEng, KinEng${m === 'limit' ? ', the cap count' : ''}). TotEng: ${e0.toPrecision(8)} → ${e12.toPrecision(8)} (${MODES[m].desc}).</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
