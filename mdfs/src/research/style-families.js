// Research · B.4: the families of styles. One butane molecule where pair, bond, angle, dihedral, fix and compute each report one number, all matched by hand.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
// C's %.<p>g
const gP = (x, p) => {
  if (x === 0) return '0';
  const ex = Math.floor(Math.log10(Math.abs(Number(x.toPrecision(p)))));
  if (ex < -4 || ex >= p) { let [m, e] = x.toExponential(p - 1).split('e'); m = m.replace(/\.?0+$/, ''); const n = Math.abs(+e); return `${m}e${+e < 0 ? '-' : '+'}${n < 10 ? '0' + n : n}`; }
  return String(Number(x.toPrecision(p)));
};
const g10 = x => gP(x, 10);
const D2R = Math.PI / 180;

// united-atom butane: CH3–CH2–CH2–CH3, TraPPE-like parameters (LAMMPS form, no ½)
const FF = { kb: 260.0, b0: 1.54, ka: 62.1, a0: 114.0, K: [1.411, -0.271, 3.145, 0.0], eps: 0.195, sig: 3.75 };
const GEO = { b: [1.55, 1.53, 1.56], t: [110, 116] };

const sub = (a, b) => a.map((x, i) => x - b[i]), dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(...a), unit = a => { const n = norm(a); return a.map(x => x / n); };

function build(phi) {
  const [b1, b2, b3] = GEO.b, t1 = GEO.t[0] * D2R, t2 = GEO.t[1] * D2R, p = phi * D2R;
  const A = [0, 0, 0], B = [b1, 0, 0], C = [b1 - b2 * Math.cos(t1), b2 * Math.sin(t1), 0];
  const bc = unit(sub(C, B)), n = unit(cross(sub(B, A), bc)), m = cross(n, bc);
  const d = [-b3 * Math.cos(t2), b3 * Math.sin(t2) * Math.cos(p), b3 * Math.sin(t2) * Math.sin(p)];
  const Dp = C.map((c, i) => c + d[0] * bc[i] + d[1] * m[i] + d[2] * n[i]);
  return [A, B, C, Dp].map(v => v.map(x => (Math.abs(x) < 5e-7 ? 0 : x).toFixed(6)));
}
const dataFile = X => `butane, united atoms

4 atoms
3 bonds
2 angles
1 dihedrals

2 atom types
1 bond types
1 angle types
1 dihedral types

-20 20 xlo xhi
-20 20 ylo yhi
-20 20 zlo zhi

Masses

1 15.035
2 14.027

Atoms # full

${X.map((x, i) => `${i + 1} 1 ${i === 0 || i === 3 ? 1 : 2} 0.0 ${x.join(' ')}`).join('\n')}

Bonds

1 1 1 2
2 1 2 3
3 1 3 4

Angles

1 1 1 2 3
2 1 2 3 4

Dihedrals

1 1 1 2 3 4
`;
const script = s => `units           real
atom_style      full
special_bonds   lj 0.0 0.0 ${s.w}
read_data       butane.data
pair_style      lj/cut 14.0                       # family: pair
pair_coeff      1 1 ${FF.eps} ${FF.sig}
pair_coeff      2 2 0.0914 3.95
bond_style      harmonic                          # family: bond
bond_coeff      1 ${FF.kb} ${FF.b0}
angle_style     harmonic                          # family: angle
angle_coeff     1 ${FF.ka} ${FF.a0}
dihedral_style  opls                              # family: dihedral
dihedral_coeff  1 ${FF.K.join(' ')}
${s.fix ? `fix             r all restrain dihedral 1 2 3 4 ${s.K} ${s.K} ${s.p0}   # family: fix
fix_modify      r energy yes
` : ''}compute         b all bond/local dist             # family: compute
compute         a all angle/local theta
compute         d all dihedral/local phi
dump            1 all local 1 b.dump c_b
dump            2 all local 1 a.dump c_a
dump            3 all local 1 d.dump c_d
dump_modify     1 format float %.10g
dump_modify     2 format float %.10g
dump_modify     3 format float %.10g
thermo_style    custom step ebond eangle edihed evdwl${s.fix ? ' f_r' : ''} pe
thermo_modify   format float %.10g
run             0
`;

function hand(X, s) {
  const P = X.map(v => v.map(Number));
  const r = [0, 1, 2].map(i => norm(sub(P[i + 1], P[i])));
  const ang = [0, 1].map(i => { const u = sub(P[i], P[i + 1]), v = sub(P[i + 2], P[i + 1]); return Math.acos(dot(u, v) / (norm(u) * norm(v))); });
  const b1 = sub(P[1], P[0]), b2 = sub(P[2], P[1]), b3 = sub(P[3], P[2]);
  const n1 = cross(b1, b2), n2 = cross(b2, b3);
  const phi = Math.atan2(norm(b2) * dot(b1, n2), dot(n1, n2));
  const Eb = r.reduce((t, x) => t + FF.kb * (x - FF.b0) ** 2, 0);
  const Ea = ang.reduce((t, x) => t + FF.ka * (x - FF.a0 * D2R) ** 2, 0);
  const [K1, K2, K3, K4] = FF.K;
  const Ed = 0.5 * K1 * (1 + Math.cos(phi)) + 0.5 * K2 * (1 - Math.cos(2 * phi)) + 0.5 * K3 * (1 + Math.cos(3 * phi)) + 0.5 * K4 * (1 - Math.cos(4 * phi));
  const r14 = norm(sub(P[3], P[0])), q6 = (FF.sig / r14) ** 6;
  const Ev = s.w * 4 * FF.eps * (q6 * q6 - q6);
  const Er = s.fix ? s.K * (1 + Math.cos(phi - s.p0 * D2R)) : 0;
  const ErDoc = s.fix ? s.K * (1 + Math.cos(phi - (s.p0 + 180) * D2R)) : 0;
  return { r, ang: ang.map(x => x / D2R), phi: phi / D2R, Eb, Ea, Ed, r14, Ev, Er, ErDoc, pe: Eb + Ea + Ed + Ev + Er };
}

export default {
  id: 'style-families',
  title: 'The families of styles',
  minutes: 25,
  summary: 'Pair, bond, angle, dihedral, improper, kspace, fix and compute styles each have their own moment in a timestep and their own number in the output. One united-atom butane molecule, built from bond lengths, angles and a dihedral you choose, lets every family speak once: you predict each energy and each geometric value by hand and LAMMPS reports the same. On the way, a hand calculation shows that fix restrain dihedral holds the molecule 180° away from the angle the manual promises.',
  learn: [
    'The families, the base class of each, and where in a timestep each one is called',
    'Forces come from pair, bonded and kspace styles; fixes change the motion; computes only measure',
    'Build a molecule from internal coordinates and predict E_bond, E_angle, E_dihed, E_vdwl and the fix energy by hand',
    'Special bonds: why the 1–4 LJ weight matters',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: who speaks when in one timestep</h2>
    <p>Page B.3 showed one family (pair) in detail. Every family works the same way: a base class with sockets, many plugs. What differs is <b>when</b> the main loop calls them. The order inside one step of <code>Verlet::run()</code> (${gh('src/verlet.cpp', 258)}):</p>
<pre class="lmp-mini">258  modify->initial_integrate(vflag);        // FIX: first half of velocity Verlet (nve, nvt …)
306  force_clear();                           // all forces f = 0
316  force->pair->compute(eflag,vflag);       // PAIR
321  if (force->bond) force->bond->compute(eflag,vflag);          // BOND
322  if (force->angle) force->angle->compute(eflag,vflag);        // ANGLE
323  if (force->dihedral) force->dihedral->compute(eflag,vflag);  // DIHEDRAL
324  if (force->improper) force->improper->compute(eflag,vflag);  // IMPROPER
329  force->kspace->compute(eflag,vflag);     // KSPACE (long-range Coulomb)
347  if (n_post_force_any) modify->post_force(vflag);   // FIX: add or change forces (restrain, addforce, langevin …)
348  modify->final_integrate();               // FIX: second half of velocity Verlet
349  if (n_end_of_step) modify->end_of_step();          // FIX: averages, printing, thermostats …
356  output->write(ntimestep);                // asks COMPUTEs for the numbers it prints</pre>
    <div class="tbl wide"><table>${head(['Family', 'Base class', 'Its job', 'Thermo keyword', 'In this browser'])}
${row(['pair', 'Pair', 'forces between nearby atoms', '<code>evdwl</code>, <code>ecoul</code>', 'lj/cut, morse, buck …'])}
${row(['bond', 'Bond', 'force along each listed bond', '<code>ebond</code>', 'harmonic, fene, morse …'])}
${row(['angle', 'Angle', 'force on each listed angle', '<code>eangle</code>', 'harmonic, cosine, charmm …'])}
${row(['dihedral', 'Dihedral', 'twist around the middle bond of 4 atoms', '<code>edihed</code>', 'opls, harmonic, charmm …'])}
${row(['improper', 'Improper', 'keeps 4 atoms flat or a centre chiral', '<code>eimp</code>', 'harmonic, cvff, umbrella'])}
${row(['kspace', 'KSpace', 'the long-range part of Coulomb', '<code>elong</code>', 'none here (page D.6)'])}
${row(['fix', 'Fix', 'anything that changes the motion: integrate, thermostat, add forces', '<code>f_ID</code>', '79 styles'])}
${row(['compute', 'Compute', 'measures, never changes anything', '<code>c_ID</code>', '68 styles'])}</table></div>
    <p>A useful rule: <b>forces</b> come from pair, bonded and kspace styles; <b>fixes</b> move atoms and may add forces; <b>computes</b> only read.</p>
  </section>

  <section class="block">
    <h2>2 · Modify: build a butane molecule</h2>
    <div class="prose"><p>Four united atoms, CH₃–CH₂–CH₂–CH₃, with TraPPE-like parameters. The page places the atoms from three bond lengths (1.55, 1.53, 1.56 Å), two angles (110°, 116°) and the dihedral φ you choose, writes the data file (page A.8), and gives it to LAMMPS.</p></div>
    <div class="lc-in">
      <label>dihedral φ (°) <input type="number" class="lmp-in sm" data-s="phi" value="65" step="5" min="-180" max="180"></label>
      <label>1–4 LJ weight <select class="lc-u" data-s="w"><option>0.0</option><option>0.5</option><option>1.0</option></select></label>
      <label class="rd-sp-lab" style="justify-content:flex-end"><span><input type="checkbox" data-s="fix" checked> add fix restrain</span></label>
      <label>restrain K <input type="number" class="lmp-in sm" data-s="K" value="2.0" step="0.5" min="0"></label>
      <label>restrain φ₀ (°) <input type="number" class="lmp-in sm" data-s="p0" value="120" step="10" min="-180" max="180"></label>
    </div>
    <details><summary class="muted">butane.data and the script</summary><pre class="lmp-mini sf-data"></pre><pre class="lmp-mini sf-script"></pre></details>
  </section>

  <section class="block">
    <h2>3 · Predict by hand, then run</h2>
    <div class="calc sf-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run with LAMMPS</button><span class="lmp-status muted sf-st"></span></div>
    <div class="tbl sf-cmp"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe: what each family said</h2>
    <ul>
      <li><b>bond, angle, dihedral</b>: each one sums over its own list from the data file (3 bonds, 2 angles, 1 dihedral). LAMMPS's harmonic forms have no ½ (page A.8). The OPLS dihedral is ½K₁(1 + cos φ) + ½K₂(1 − cos 2φ) + ½K₃(1 + cos 3φ) + ½K₄(1 − cos 4φ) (${gh('src/MOLECULE/dihedral_opls.cpp', 177)}).</li>
      <li><b>pair</b>: LJ is computed for every pair within the cutoff, then <code>special_bonds</code> multiplies 1–2, 1–3 and 1–4 pairs by its three weights. In butane only atoms 1 and 4 are a 1–4 pair, and with weight 0 (TraPPE) E_vdwl is exactly 0. At φ = 65° with weight 1 they add over 4 kcal/mol, more than all bonded terms together: the dihedral parameters were fitted assuming that weight, so the weight is part of the force field.</li>
      <li><b>fix</b>: <code>restrain</code> adds a force in <code>post_force</code> and, with <code>fix_modify energy yes</code>, its energy to PotEng.</li>
      <li><b>compute</b>: bond/local, angle/local and dihedral/local measured the geometry and changed nothing. Their numbers equal the hand geometry.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>5 · A surprise: where does fix restrain hold the dihedral?</h2>
    <p>The manual for fix restrain says E = K[1 + cos(nφ − d)] with d = φ₀ + 180°, so that the energy is <b>lowest at φ = φ₀</b>. The code (${gh('src/fix_restrain.cpp', 135)} and ${gh('src/fix_restrain.cpp', 642)}) uses d = φ₀:</p>
<pre class="lmp-mini">135  target[nrestrain] = DEG2RAD * utils::numeric(FLERR,arg[iarg+7],false,lmp);   // φ₀, no + 180
642  p = p*cos_target[m] + df1*sin_target[m];      // cos(φ − φ₀)
645  p += 1.0;
647  edihed += k * p;                              // E = K[1 + cos(φ − φ₀)]: HIGHEST at φ = φ₀</pre>
    <p>Section 3 already showed it: the hand value with the code's formula matches f_r; the manual's formula does not. To see what it does to a molecule, minimize with the restraint as the only dihedral term and look where φ ends:</p>
    <div class="row"><button class="btn sm primary" data-k="min">▶ Minimize with only the restraint</button><span class="lmp-status muted sf-mst"></span></div>
    <div class="sf-min"></div>
    <p>The same lines are in the current development version and in versions back to at least 2018. If you use fix restrain dihedral (for example to scan a torsion profile), check φ after the run, or give φ₀ − 180°.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Which φ gives E_dihed = 0? Predict from the OPLS formula, then set it.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>φ = <b>180°</b> (trans): cos φ = −1, cos 2φ = 1, cos 3φ = −1, so every bracket is 0. LAMMPS, like IUPAC, calls trans 180°.</div></div>
    <div class="ex"><p><b>2.</b> Set φ = 0 (cis) with 1–4 weight 1. Why does E_vdwl explode?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Atoms 1 and 4 come to about 2.6 Å, well inside σ = 3.75 Å, where (σ/r)¹² is huge. This is one reason force fields exclude or scale the 1–4 pair and put the cis barrier in the dihedral term instead.</div></div>
    <div class="ex"><p><b>3.</b> Which family would you use to keep atom 1 from moving, and why not a compute?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A fix: <code>fix 2 one setforce 0 0 0</code> (with v = 0) or <code>fix spring/self</code>. A compute can only measure; it is never called in a way that changes x, v or f.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const S = () => {
      const v = k => root.querySelector(`[data-s=${k}]`);
      const num = (k, d) => { const x = parseFloat(v(k).value); return isFinite(x) ? x : d; };
      return { phi: Math.max(-180, Math.min(180, num('phi', 65))), w: v('w').value, fix: v('fix').checked, K: num('K', 2), p0: num('p0', 120) };
    };
    let cur = null;
    const draw = () => {
      const s = S(), X = build(s.phi), data = dataFile(X), H = hand(X, { ...s, w: +s.w });
      cur = { s, data, H };
      $('.sf-data').textContent = data; $('.sf-script').textContent = script(s);
      const deg = x => g10(x) + '°';
      $('.sf-hand').innerHTML = esc(`geometry (from the coordinates in butane.data)
  bonds   r = ${H.r.map(g10).join(', ')} Å
  angles  θ = ${H.ang.map(deg).join(', ')}
  φ = ${deg(H.phi)},   r₁₄ = ${g10(H.r14)} Å
bond      E = 260 Σ(r − 1.54)²                   = `) + `<b>${g10(H.Eb)}</b>` + esc(`
angle     E = 62.1 Σ(θ − 114°)²   (θ in radians)  = `) + `<b>${g10(H.Ea)}</b>` + esc(`
dihedral  E = ½K₁(1+cos φ) + ½K₂(1−cos 2φ) + ½K₃(1+cos 3φ) = `) + `<b>${g10(H.Ed)}</b>` + esc(`
pair      E = ${s.w} × 4·0.195[(3.75/r₁₄)¹² − (3.75/r₁₄)⁶]   = `) + `<b>${g10(H.Ev)}</b>` + (s.fix ? esc(`
fix       E = K[1 + cos(φ − φ₀)]  (the code)       = `) + `<b>${g10(H.Er)}</b>` + esc(`
          E = K[1 + cos(φ − φ₀ − 180°)]  (the manual) = ${g10(H.ErDoc)}`) : '') + esc(`
PotEng    = sum                                  = `) + `<b>${g10(H.pe)}</b>`;
      $('.sf-cmp').innerHTML = ''; $('.sf-st').textContent = '';
    };
    root.querySelectorAll('[data-s]').forEach(i => { i.oninput = draw; i.onchange = draw; }); draw();
    const val = txt => (txt || '').split('\n').filter(l => /^-?[\d.]+(e[-+]\d+)?\s*$/i.test(l.trim()) && !/^ITEM/.test(l)).map(l => l.trim());
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.sf-st'), { s, data, H } = cur; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(script(s), { files: { 'butane.data': data }, read: ['b.dump', 'a.dump', 'd.dump'] }); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const tl = res.lines.findIndex(l => /^\s+Step\s+E_bond/.test(l)), cols = res.lines[tl].trim().split(/\s+/), vals = res.lines[tl + 1].trim().split(/\s+/);
      const T = k => vals[cols.indexOf(k)];
      // dump local: the numbers after "ITEM: ENTRIES"
      const ent = name => { const t = res.files[name] || ''; return t.split('ITEM: ENTRIES')[1]?.split('\n').slice(1).map(x => x.trim()).filter(Boolean) || []; };
      const B = ent('b.dump'), A = ent('a.dump'), Dd = ent('d.dump');
      const R = [
        ['compute bond/local', 'r₁₂', g10(H.r[0]), B[0]], ['', 'r₂₃', g10(H.r[1]), B[1]], ['', 'r₃₄', g10(H.r[2]), B[2]],
        ['compute angle/local', 'θ₁₂₃', g10(H.ang[0]), A[0]], ['', 'θ₂₃₄', g10(H.ang[1]), A[1]],
        ['compute dihedral/local', 'φ', g10(H.phi), Dd[0]],
        ['bond', 'E_bond', g10(H.Eb), T('E_bond')], ['angle', 'E_angle', g10(H.Ea), T('E_angle')], ['dihedral', 'E_dihed', g10(H.Ed), T('E_dihed')],
        ['pair', 'E_vdwl', g10(H.Ev), T('E_vdwl')], ...(s.fix ? [['fix', 'f_r', g10(H.Er), T('f_r')]] : []), ['all', 'PotEng', g10(H.pe), T('PotEng')],
      ];
      const close = (h, l) => h === l || (l !== undefined && Math.abs(+h - +l) <= 2e-9 * Math.max(1, Math.abs(+l)));
      $('.sf-cmp').innerHTML = `<table>${head(['family', '', 'by hand', 'LAMMPS', ''])}${R.map(([f, k, h, l]) => row([f, k, `<span class="mono">${h}</span>`, `<span class="mono">${l ?? '—'}</span>`, close(h, l) ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'])).join('')}</table>`;
    };
    $('[data-k=min]').onclick = async e => {
      const b = e.currentTarget, st = $('.sf-mst'), { s, data } = cur; b.disabled = true; st.textContent = 'Minimizing…';
      const K = s.K > 0 ? s.K : 2, p0 = s.p0;
      const sc = `units           real
atom_style      full
read_data       butane.data
bond_style      harmonic
bond_coeff      1 ${FF.kb} ${FF.b0}
angle_style     harmonic
angle_coeff     1 ${FF.ka} ${FF.a0}
dihedral_style  zero
dihedral_coeff  1
pair_style      zero 10.0
pair_coeff      * *
fix             r all restrain dihedral 1 2 3 4 ${K * 10} ${K * 10} ${p0}
fix_modify      r energy yes
compute         d all dihedral/local phi
compute         dm all reduce max c_d inputs local
thermo_style    custom step f_r c_dm
minimize        0 1e-10 2000 20000
`;
      const res = await runLammps(sc, { files: { 'butane.data': data } }); b.disabled = false; st.textContent = '';
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      const tl = res.lines.findIndex(l => /^\s+Step\s+f_r/.test(l)), rows = [];
      for (let i = tl + 1; i < res.lines.length && /^\s+\d+\s/.test(res.lines[i]); i++) rows.push(res.lines[i].trim().split(/\s+/));
      const first = rows[0], last = rows[rows.length - 1], wrap = x => ((x + 540) % 360) - 180;
      $('.sf-min').innerHTML = `<div class="tbl"><table>${head(['', 'step', 'f_r (kcal/mol)', 'φ'])}${row(['start', first[0], first[1], first[2] + '°'])}${row(['end', last[0], last[1], `<b>${(+last[2]).toFixed(2)}°</b>`])}</table></div>
        <p>You asked for φ₀ = <b>${p0}°</b>. The manual predicts the end at ${p0}°; the code's formula predicts φ₀ − 180° = <b>${wrap(p0 - 180)}°</b>. ${Math.abs(wrap(+last[2] - (p0 - 180))) < 1 ? '<span class="m-ok">LAMMPS followed the code.</span>' : ''}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
