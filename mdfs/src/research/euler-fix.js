// Research · E.5: a broken integrator on purpose. A new fix nve/euler (explicit and symplectic Euler, 60 lines,
// derived from FixNVE); three atoms for ten steps by hand, bit for bit; an 864-atom liquid that explicit Euler
// heats until atoms are lost, while symplectic Euler and velocity Verlet stay bounded with errors O(dt) and O(dt²).
import { NATIVE } from './euler-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const DTS = ['0.005', '0.0025', '0.00125'];
const NAMES = { nve: 'velocity Verlet (fix nve)', nve_euler: 'explicit Euler', nve_euler_symplectic: 'symplectic Euler' };

// three LJ atoms (page E.1 positions), pair order 1-2, 1-3, 2-3 as in the neighbour list
const X0 = [[4.0, 5.0, 5.0], [5.2, 5.1, 5.0], [4.5, 6.0, 5.4]];
function forces(x) {
  const f = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; let e = 0.0;
  for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) {
    const d = [0, 1, 2].map(k => x[i][k] - x[j][k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (rsq < 6.25) { const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = 1.0 * (r6 * (48.0 * r6 - 24.0)) * r2; for (let k = 0; k < 3; k++) { f[i][k] += d[k] * fp; f[j][k] -= d[k] * fp; } let ev = r6 * (4.0 * r6 - 4.0) - 0.0; ev *= 1.0; e += ev; }
  }
  return { f, e: 0.0 + e };
}
function integrate(kind) {
  const x = X0.map(a => [...a]), v = [[0.3, 0, 0], [0, 0, 0], [0, 0, 0]], dt = 0.005, dtf = 0.5 * dt * 1.0, dtv = dt, out = [];
  let { f } = forces(x);
  for (let s = 1; s <= 10; s++) {
    if (kind === 'nve') {                                            // FixNVE: half kick, drift / forces / half kick
      for (let i = 0; i < 3; i++) { const m = dtf / 1.0; for (let k = 0; k < 3; k++) v[i][k] += m * f[i][k]; for (let k = 0; k < 3; k++) x[i][k] += dtv * v[i][k]; }
      ({ f } = forces(x));
      for (let i = 0; i < 3; i++) { const m = dtf / 1.0; for (let k = 0; k < 3; k++) v[i][k] += m * f[i][k]; }
    } else {                                                          // FixNVEEuler::initial_integrate, then forces
      const full = 2.0 * dtf;
      for (let i = 0; i < 3; i++) { const m = full / 1.0; for (let k = 0; k < 3; k++) { if (kind === 'nve_euler_symplectic') { v[i][k] += m * f[i][k]; x[i][k] += dtv * v[i][k]; } else { x[i][k] += dtv * v[i][k]; v[i][k] += m * f[i][k]; } } }
      ({ f } = forces(x));
    }
    out.push({ x: x.map(a => [...a]), v: v.map(a => [...a]), pe: forces(x).e });
  }
  return out;
}
const frames = t => { const L = t.split('\n'), F = []; for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const n = +L[i + 3]; F.push(L.slice(i + 9, i + 9 + n).map(l => l.trim().split(/\s+/).map(Number))); i += 8 + n; } return F; };

// small line chart, colours from the site palette
function chart(series, { w = 640, h = 220, ylab = '', xmax = 10 } = {}) {
  const ys = series.flatMap(s => s.y), ymin = Math.min(...ys), ymax = Math.max(...ys), pad = 40;
  const X = t => pad + (w - pad - 10) * t / xmax, Y = y => h - 25 - (h - 40) * (y - ymin) / (ymax - ymin || 1);
  const lines = series.map(s => `<polyline fill="none" stroke="${s.c}" stroke-width="1.8" points="${s.t.map((t, k) => `${X(t).toFixed(1)},${Y(s.y[k]).toFixed(1)}`).join(' ')}"/>`).join('');
  const leg = series.map((s, k) => `<text x="${pad + 8}" y="${16 + 14 * k}" fill="${s.c}" font-size="12">${esc(s.name)}</text>`).join('');
  return `<svg viewBox="0 0 ${w} ${h}" style="width:100%;max-width:${w}px" role="img"><line x1="${pad}" y1="${h - 25}" x2="${w - 10}" y2="${h - 25}" stroke="var(--muted)"/><line x1="${pad}" y1="10" x2="${pad}" y2="${h - 25}" stroke="var(--muted)"/>
<text x="${w - 10}" y="${h - 8}" fill="var(--muted)" font-size="11" text-anchor="end">t (τ)</text><text x="${pad - 4}" y="${Y(ymax) + 4}" fill="var(--muted)" font-size="11" text-anchor="end">${+ymax.toPrecision(4)}</text><text x="${pad - 4}" y="${Y(ymin) + 4}" fill="var(--muted)" font-size="11" text-anchor="end">${+ymin.toPrecision(4)}</text><text x="${pad + 4}" y="${h - 30}" fill="var(--muted)" font-size="11">${esc(ylab)}</text>${lines}${leg}</svg>`;
}

export default {
  id: 'euler-fix',
  title: 'A broken integrator on purpose',
  minutes: 20,
  summary: 'Write a new fix that integrates the equations of motion the textbook way, Euler, and watch it fail. fix nve/euler is 60 lines derived from fix nve, with two variants: explicit Euler (x and v both from the old values) and symplectic Euler (new velocity first). Ten steps of three atoms match a hand calculation bit for bit for all three integrators; then an 864-atom LJ liquid shows explicit Euler heating until LAMMPS loses atoms, at every time step tried, while symplectic Euler and velocity Verlet keep the energy bounded, with fluctuations that shrink as dt and dt².',
  learn: [
    'A fix that integrates: setmask, initial_integrate, final_integrate',
    'Deriving from an existing fix instead of copying it',
    'Explicit vs symplectic Euler vs velocity Verlet: one line of order',
    'Energy error O(dt) or O(dt²), and drift vs bounded fluctuation',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three ways to take a step</h2>
    <div class="prose"><p>${gh('src/fix_nve.cpp')} is velocity Verlet in two halves around the force computation (page C.2): half kick and drift in <code>initial_integrate</code>, half kick in <code>final_integrate</code>. The two Euler methods take one whole step before the forces:</p>
<pre class="lmp-mini">velocity Verlet   v += ½dt f/m;  x += dt v;   [forces];  v += ½dt f/m
explicit Euler    x += dt v;     v += dt f/m;  [forces]           (x and v from old values)
symplectic Euler  v += dt f/m;   x += dt v;    [forces]           (x uses the new v)</pre>
<p>The only difference between the two Euler methods is the order of two lines, and it decides everything. Explicit Euler multiplies the energy of a harmonic oscillator by (1 + ω²dt²) every step: it always heats. Symplectic Euler preserves phase-space volume, like Verlet, so its energy error stays bounded, but it is only first order: the error is O(dt), Verlet's O(dt²).</p></div>
  </section>

  <section class="block">
    <h2>2 · Change: a fix derived from FixNVE</h2>
    <div class="prose"><p>Nothing to copy this time: the new class inherits the constructor, init() (which sets dtv = dt and dtf = ½dt·ftm2v) and reset_dt() from FixNVE, and overrides three functions. <code>setmask</code> returns only INITIAL_INTEGRATE, so Verlet::run never calls final_integrate for it, and rRESPA is not offered.</p></div>
    <pre class="lmp-mini">${esc(NATIVE.hdr.split('\n').slice(4).join('\n').replace(/\n{3,}/g, '\n\n'))}</pre>
    <pre class="lmp-mini">${esc(NATIVE.src.split('\n').slice(5).join('\n'))}</pre>
    <div class="prose"><p>Usage: <code>fix 1 all nve/euler</code> or <code>fix 1 all nve/euler symplectic</code>. Native rebuild 13 s (the new file, modify.cpp for the regenerated fix list, lammps.cpp).</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: three atoms, ten steps, three integrators</h2>
    <div class="prose"><p>The three atoms of page E.1 (atom 1 moving at 0.3σ/τ), dt = 0.005, positions and velocities dumped every step with 17 digits. The hand code repeats each integrator line by line, with the LJ forces in the neighbour list's pair order.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Integrate by hand</button></div>
    <div class="tbl eu-a"></div>
    <div class="prose eu-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: an LJ liquid, 10 τ</h2>
    <div class="prose"><p>864 atoms at ρ = 0.8442, equilibrated with fix nve (T ≈ 0.7), written to a restart file. From it, each integrator runs 10τ at dt = 0.005, 0.0025 and 0.00125, with the total energy printed every 0.1τ.</p>
<pre class="lmp-mini">read_restart    liq.bin
timestep        0.005
fix             1 all nve/euler
thermo_style    custom step time etotal temp
run             2000</pre></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Show the runs</button></div>
    <div class="eu-ch1"></div>
    <div class="tbl eu-b"></div>
    <div class="prose eu-bsum"></div>
    <div class="eu-ch2"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> A harmonic oscillator, explicit Euler. Show that each step multiplies the energy by 1 + ω²dt².</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>x' = x + dt v, v' = v − dt ω² x. Then ω²x'² + v'² = ω²x² + v² + dt²ω²(v² + ω²x²) (the cross terms cancel) = (1 + ω²dt²)(ω²x² + v²). Over a time t that is (1 + ω²dt²)^(t/dt) ≈ e^(ω²dt·t): the heating rate is proportional to dt, so halving dt doubles the survival time, roughly what the liquid showed.</div></div>
    <div class="ex"><p><b>2.</b> Symplectic Euler's positions are those of leapfrog Verlet. Why then is its energy error O(dt), not O(dt²)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Its velocity after the step is v(t + ½dt) in leapfrog terms, but it is reported together with x(t + dt) as if both were at the same time. The kinetic energy is off by about m v·a dt/2, first order in dt. Start it with v(0) − ½dt f/m and average neighbouring velocities, and you have Verlet again. The trajectory is fine; the bookkeeping is half a step out.</div></div>
    <div class="ex"><p><b>3.</b> Why does LAMMPS report "Lost atoms" instead of an energy error?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Nothing checks energy conservation. Once atoms are hot enough to move farther than the neighbour skin and the ghost-atom cutoff between rebuilds, or across a processor boundary between exchanges, some are no longer owned by any processor; the atom count drops and LAMMPS stops. That is the usual symptom of a time step or integrator that is too aggressive.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = () => {
      let all = 0; const rows = [];
      for (const tag of ['nve', 'nve_euler', 'nve_euler_symplectic']) {
        const H = integrate(tag), F = frames(NATIVE.two[tag].dump), T = NATIVE.two[tag].thermo;
        let ok = 0, n = 0, pk = 0;
        for (let s = 1; s <= 10; s++) { for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { n += 2; ok += (F[s][i][1 + k] === H[s - 1].x[i][k]) + (F[s][i][4 + k] === H[s - 1].v[i][k]); } if (Number(T[s][1]) === H[s - 1].pe) pk++; }
        all += ok + pk;
        rows.push(row([NAMES[tag], `<span class="mono">${esc(String(F[10][0][1]))}</span>`, `<span class="mono">${esc(String(F[10][0][4]))}</span>`, Number(T[10][3]).toFixed(9), `${ok === n ? OK : BAD} ${ok}/${n}`, `${pk === 10 ? OK : BAD} ${pk}/10`]));
      }
      $('.eu-a').innerHTML = `<table>${head(['integrator', 'x of atom 1, step 10', 'vx of atom 1, step 10', 'TotEng, step 10', 'x, v (all atoms, 10 steps)', 'PotEng'])}${rows.join('')}</table>`;
      const E0 = Number(NATIVE.two.nve.thermo[0][3]);
      $('.eu-asum').innerHTML = `<p>${all} of 570 numbers identical to the last bit. Already after ten steps the total energy (started at ${E0.toFixed(9)}) shows the pattern: Verlet changed it by ${(Number(NATIVE.two.nve.thermo[10][3]) - E0).toExponential(1)}, symplectic Euler by ${(Number(NATIVE.two.nve_euler_symplectic.thermo[10][3]) - E0).toExponential(1)}, explicit Euler by +${(Number(NATIVE.two.nve_euler.thermo[10][3]) - E0).toExponential(1)}. The three atoms drift apart at different rates even though all three methods have the same force law: the integrator is part of the model.</p>`;
    };
    $('[data-k=b]').onclick = () => {
      const D = NATIVE.drift, cols = { '0.005': 'var(--coral)', '0.0025': 'var(--amber)', '0.00125': 'var(--cyan)' };
      $('.eu-ch1').innerHTML = chart(DTS.map(dt => ({ name: `explicit Euler, dt = ${dt}`, c: cols[dt], t: D[`nve_euler_${dt}`].t, y: D[`nve_euler_${dt}`].e })), { ylab: 'TotEng per atom', xmax: 4 });
      const stats = (tag) => { const d = D[tag], n = d.e.length, mt = d.t.reduce((s, x) => s + x, 0) / n, me = d.e.reduce((s, x) => s + x, 0) / n; let sxy = 0, sxx = 0; d.t.forEach((t, k) => { sxy += (t - mt) * (d.e[k] - me); sxx += (t - mt) ** 2; }); const sl = sxy / sxx; const rms = Math.sqrt(d.e.reduce((s, e, k) => s + (e - (me + sl * (d.t[k] - mt))) ** 2, 0) / n); return { sl, rms, end: d.t[n - 1], err: d.err }; };
      const rows = [];
      for (const k of ['nve', 'nve_euler_symplectic', 'nve_euler']) for (const dt of DTS) { const s = stats(`${k}_${dt}`); rows.push(row([dt === DTS[0] ? NAMES[k] : '', dt, k === 'nve_euler' ? `lost atoms just after ${s.end.toFixed(1)}τ` : `${s.end.toFixed(0)}τ, no error`, k === 'nve_euler' ? '–' : s.rms.toExponential(2), k === 'nve_euler' ? '–' : s.sl.toExponential(1)])); }
      $('.eu-b').innerHTML = `<table>${head(['integrator', 'dt', 'run', 'σ(TotEng) per atom', 'drift per τ'])}${rows.join('')}</table>`;
      const r = DTS.map(dt => stats(`nve_${dt}`).rms), q = DTS.map(dt => stats(`nve_euler_symplectic_${dt}`).rms), e = DTS.map(dt => stats(`nve_euler_${dt}`).end);
      $('.eu-bsum').innerHTML = `<p>Explicit Euler heats from the first step and LAMMPS stops with "Lost atoms" after ${e.map(x => x.toFixed(1)).join(', ')}τ: a smaller step only buys time, about ${(e[1] / e[0]).toFixed(1)}× and ${(e[2] / e[1]).toFixed(1)}× per halving, as exercise 1 predicts. The other two never drift measurably (the fitted drift is noise). Halving dt divides symplectic Euler's fluctuation by ${(q[0] / q[1]).toFixed(1)} and ${(q[1] / q[2]).toFixed(1)}: first order. Verlet's falls by ${(r[0] / r[1]).toFixed(1)} and ${(r[1] / r[2]).toFixed(1)}: second order, and 30 to 100 times smaller to begin with. Same forces, same cost per step.</p>`;
      $('.eu-ch2').innerHTML = chart([{ name: 'symplectic Euler, dt = 0.005', c: 'var(--amber)', t: D['nve_euler_symplectic_0.005'].t, y: D['nve_euler_symplectic_0.005'].e }, { name: 'velocity Verlet, dt = 0.005', c: 'var(--cyan)', t: D['nve_0.005'].t, y: D['nve_0.005'].e }], { ylab: 'TotEng per atom' });
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
