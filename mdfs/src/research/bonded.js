// Research · D.8: bonds, angles and dihedrals. bond_harmonic, angle_harmonic and dihedral_opls reproduced operation by operation for a united-atom butane: energies and the force on every atom, a full torsion scan, and the K conventions that differ between codes.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const D2R = Math.PI / 180.0;

// force field: united-atom butane (OPLS-UA style numbers), LAMMPS conventions
const FF = { kb: 260.0, r0: 1.54, ka: 62.1, th0: 114.0, K: [1.411, -0.271, 3.145, 0.0] };

// one molecule: bond length r, angle θ (deg), dihedral φ (deg; 180 = trans), centred near c
function butane(r, th, ph, c = [10, 10, 10]) {
  const t = th * D2R, p = ph * D2R;
  const a2 = [0, 0, 0], a3 = [r, 0, 0], a1 = [r * Math.cos(t), r * Math.sin(t), 0];
  const a4 = [r - r * Math.cos(t), r * Math.sin(t) * Math.cos(p), r * Math.sin(t) * Math.sin(p)];
  return [a1, a2, a3, a4].map(v => v.map((q, k) => q + c[k]));
}
const DATA = mols => `united-atom butane

${4 * mols.length} atoms
${3 * mols.length} bonds
${2 * mols.length} angles
${mols.length} dihedrals
1 atom types
1 bond types
1 angle types
1 dihedral types

0 ${20 * mols.length} xlo xhi
0 20 ylo yhi
0 20 zlo zhi

Masses

1 14.027

Atoms # molecular

${mols.flatMap((m, k) => m.map((x, i) => `${4 * k + i + 1} ${k + 1} 1 ${x.map(String).join(' ')}`)).join('\n')}

Bonds

${mols.flatMap((m, k) => [1, 2, 3].map(i => `${3 * k + i} 1 ${4 * k + i} ${4 * k + i + 1}`)).join('\n')}

Angles

${mols.flatMap((m, k) => [1, 2].map(i => `${2 * k + i} 1 ${4 * k + i} ${4 * k + i + 1} ${4 * k + i + 2}`)).join('\n')}

Dihedrals

${mols.map((m, k) => `${k + 1} 1 ${4 * k + 1} ${4 * k + 2} ${4 * k + 3} ${4 * k + 4}`).join('\n')}
`;
const script = (extra) => `units           real
atom_style      molecular
read_data       b.data
pair_style      zero 5.0
pair_coeff      * *
bond_style      harmonic
bond_coeff      1 ${FF.kb} ${FF.r0}
angle_style     harmonic
angle_coeff     1 ${FF.ka} ${FF.th0}
dihedral_style  opls
dihedral_coeff  1 ${FF.K.join(' ')}
${extra}
thermo_style    custom step ebond eangle edihed
thermo_modify   format float %.17g
run             0`;

// --- the three compute() functions, in the C++ order ---
function bond(x, f, i1, i2) {
  const dx = x[i1][0] - x[i2][0], dy = x[i1][1] - x[i2][1], dz = x[i1][2] - x[i2][2], rsq = dx * dx + dy * dy + dz * dz, r = Math.sqrt(rsq), dr = r - FF.r0, rk = FF.kb * dr;
  const fb = r > 0 ? -2.0 * rk / r : 0.0;
  f[i1][0] += dx * fb; f[i1][1] += dy * fb; f[i1][2] += dz * fb; f[i2][0] -= dx * fb; f[i2][1] -= dy * fb; f[i2][2] -= dz * fb;
  return rk * dr;
}
function angle(x, f, i1, i2, i3) {
  const d1 = [0, 1, 2].map(k => x[i1][k] - x[i2][k]), d2 = [0, 1, 2].map(k => x[i3][k] - x[i2][k]);
  const rsq1 = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2], r1 = Math.sqrt(rsq1), rsq2 = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2], r2 = Math.sqrt(rsq2);
  let c = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]; c /= r1 * r2; if (c > 1.0) c = 1.0; if (c < -1.0) c = -1.0;
  let s = Math.sqrt(1.0 - c * c); if (s < 0.001) s = 0.001; s = 1.0 / s;
  const dth = Math.acos(c) - FF.th0 * D2R, tk = FF.ka * dth, e = tk * dth, a = -2.0 * tk * s, a11 = a * c / rsq1, a12 = -a / (r1 * r2), a22 = a * c / rsq2;
  const f1 = [0, 1, 2].map(k => a11 * d1[k] + a12 * d2[k]), f3 = [0, 1, 2].map(k => a22 * d2[k] + a12 * d1[k]);
  for (let k = 0; k < 3; k++) { f[i1][k] += f1[k]; f[i2][k] -= f1[k] + f3[k]; f[i3][k] += f3[k]; }
  return { e, theta: Math.acos(c) / D2R };
}
function opls(x, f, i1, i2, i3, i4) {
  const [k1, k2, k3, k4] = FF.K.map(v => 0.5 * v);
  const vb1 = [0, 1, 2].map(k => x[i1][k] - x[i2][k]), vb2 = [0, 1, 2].map(k => x[i3][k] - x[i2][k]), vb2m = vb2.map(v => -v), vb3 = [0, 1, 2].map(k => x[i4][k] - x[i3][k]);
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const sb1 = 1.0 / dot(vb1, vb1), sb2 = 1.0 / dot(vb2, vb2), sb3 = 1.0 / dot(vb3, vb3), rb1 = Math.sqrt(sb1), rb3 = Math.sqrt(sb3);
  const c0 = dot(vb1, vb3) * rb1 * rb3, b1mag = Math.sqrt(dot(vb1, vb1)), b2mag = Math.sqrt(dot(vb2, vb2)), b3mag = Math.sqrt(dot(vb3, vb3));
  const r12c1 = 1.0 / (b1mag * b2mag), c1mag = dot(vb1, vb2) * r12c1, r12c2 = 1.0 / (b2mag * b3mag), c2mag = dot(vb2m, vb3) * r12c2;
  let sc1 = Math.sqrt(Math.max(1.0 - c1mag * c1mag, 0.0)); if (sc1 < 0.001) sc1 = 0.001; sc1 = 1.0 / sc1;
  let sc2 = Math.sqrt(Math.max(1.0 - c2mag * c2mag, 0.0)); if (sc2 < 0.001) sc2 = 0.001; sc2 = 1.0 / sc2;
  const s1 = sc1 * sc1, s2 = sc2 * sc2; let s12 = sc1 * sc2, c = (c0 + c1mag * c2mag) * s12;
  const cx = vb1[1] * vb2[2] - vb1[2] * vb2[1], cy = vb1[2] * vb2[0] - vb1[0] * vb2[2], cz = vb1[0] * vb2[1] - vb1[1] * vb2[0], cmag = Math.sqrt(cx * cx + cy * cy + cz * cz);
  const dxx = (cx * vb3[0] + cy * vb3[1] + cz * vb3[2]) / cmag / b3mag;
  if (c > 1.0) c = 1.0; if (c < -1.0) c = -1.0;
  let phi = Math.acos(c); if (dxx < 0.0) phi *= -1.0;
  let si = Math.sin(phi); if (Math.abs(si) < 0.00001) si = 0.00001; const siinv = 1.0 / si;
  const p = k1 * (1.0 + c) + k2 * (1.0 - Math.cos(2.0 * phi)) + k3 * (1.0 + Math.cos(3.0 * phi)) + k4 * (1.0 - Math.cos(4.0 * phi));
  const pd = k1 - 2.0 * k2 * Math.sin(2.0 * phi) * siinv + 3.0 * k3 * Math.sin(3.0 * phi) * siinv - 4.0 * k4 * Math.sin(4.0 * phi) * siinv;
  const a = pd; c = c * a; s12 = s12 * a;
  const a11 = c * sb1 * s1, a22 = -sb2 * (2.0 * c0 * s12 - c * (s1 + s2)), a33 = c * sb3 * s2, a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12), a13 = -rb1 * rb3 * s12, a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);
  const sx2 = [0, 1, 2].map(k => a12 * vb1[k] + a22 * vb2[k] + a23 * vb3[k]);
  const f1 = [0, 1, 2].map(k => a11 * vb1[k] + a12 * vb2[k] + a13 * vb3[k]), f4 = [0, 1, 2].map(k => a13 * vb1[k] + a23 * vb2[k] + a33 * vb3[k]);
  const f2 = f1.map((v, k) => -sx2[k] - v), f3 = f4.map((v, k) => sx2[k] - v);
  for (let k = 0; k < 3; k++) { f[i1][k] += f1[k]; f[i2][k] += f2[k]; f[i3][k] += f3[k]; f[i4][k] += f4[k]; }
  return { e: p, phi: phi / D2R };
}
function hand(x) {
  const f = x.map(() => [0, 0, 0]);
  const eb = bond(x, f, 0, 1) + bond(x, f, 1, 2) + bond(x, f, 2, 3);
  const A1 = angle(x, f, 0, 1, 2), A2 = angle(x, f, 1, 2, 3), D = opls(x, f, 0, 1, 2, 3);
  return { eb, ea: A1.e + A2.e, ed: D.e, phi: D.phi, th: [A1.theta, A2.theta], f };
}
const close = (h, l, tol = 2e-14) => h === l || Math.abs(h - l) <= tol * Math.max(Math.abs(h), Math.abs(l), 1e-3);

export default {
  id: 'bonded',
  title: 'Bonds, angles and dihedrals',
  minutes: 25,
  summary: 'Bonded terms are three short loops over lists built from the data file. You rebuild a united-atom butane at any bond length, angle and torsion, reproduce bond_harmonic, angle_harmonic and dihedral_opls line by line, and match the three energies and all 12 force components LAMMPS prints. Then a full torsion scan in one run, and the factor of 2 that catches everyone moving angle and bond constants between codes.',
  learn: [
    'Topology lists: bondlist, anglelist, dihedrallist and newton_bond',
    'angle_harmonic: E = K(θ − θ₀)², no ½, and its two forces',
    'dihedral_opls: φ from three bond vectors, trans = 180°, forces on four atoms',
    'Converting K between LAMMPS, GROMACS and AMBER conventions',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: lists and three loops</h2>
    <div class="prose"><p>read_data stores each bond, angle and dihedral with its atoms. On every neighbour rebuild, ${gh('src/ntopo.cpp')} and its subclasses turn them into flat lists of local indices (<code>neighbor->bondlist[n] = {i1, i2, type}</code>, and so on); with <code>newton_bond on</code> each term is computed once, on one processor, and forces on ghosts go home in reverse_comm. The force loop of Verlet (page C.2) then calls bond, angle, dihedral, improper in that order. ${gh('src/MOLECULE/angle_harmonic.cpp', 56)}:</p>
<pre class="lmp-mini">c = (d1 · d2) / (r1 r2);                 // cos θ, clamped to [-1, 1]
s = 1 / max(sqrt(1 − c²), 0.001);
dtheta = acos(c) − theta0;   tk = k * dtheta;
eangle = tk * dtheta;                    // K (θ − θ0)²: K already holds the ½
a = −2 tk s;   a11 = a c / r1²;   a12 = −a / (r1 r2);   a22 = a c / r2²
f1 = a11 d1 + a12 d2;   f3 = a22 d2 + a12 d1;   f2 = −(f1 + f3)</pre>
<p>${gh('src/MOLECULE/dihedral_opls.cpp', 59)} gets cos φ from the three bond vectors, the sign of φ from (b1 × b2) · b3, and</p>
<pre class="lmp-mini">E = ½K1(1 + cos φ) + ½K2(1 − cos 2φ) + ½K3(1 + cos 3φ) + ½K4(1 − cos 4φ)     // coeff stores ½K (line 283)</pre>
<p>so trans (φ = 180°) has E = 0.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: one butane, every number</h2>
    <div class="prose"><p>Bond 260 (r − 1.54)², angle 62.1 (θ − 114°)², OPLS torsion K = 1.411, −0.271, 3.145, 0 kcal/mol. <code>pair_style zero</code>, so only bonded forces act.</p></div>
    <div class="row">
      <label class="lc-in">r (Å)<input type="number" class="lmp-in sm bd-r" value="1.58" step="0.01"></label>
      <label class="lc-in">θ (°)<input type="number" class="lmp-in sm bd-t" value="109.5" step="0.5"></label>
      <label class="lc-in">φ (°)<input type="number" class="lmp-in sm bd-p" value="65" step="5"></label>
      <button class="btn sm primary" data-k="a">▶ run 0</button><span class="lmp-status muted bd-sta"></span></div>
    <div class="calc bd-pred"></div>
    <div class="tbl bd-outa"></div>
    <div class="prose bd-suma"></div>
  </section>

  <section class="block">
    <h2>3 · Modify: the torsion scan</h2>
    <div class="prose"><p>Thirteen butanes, φ = −180° … 180° in steps of 30°, in one box, <code>run 0</code>; the torsion energy of each comes from <code>compute pe/atom dihedral</code> summed over its four atoms.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Scan</button><span class="lmp-status muted bd-stb"></span></div>
    <div class="tbl bd-outb"></div>
    <div class="prose bd-sumb"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe: the factor of 2</h2>
    <div class="tbl bd-conv"></div>
    <ul>
      <li><b>LAMMPS bond and angle harmonic have no ½.</b> GROMACS writes ½k_b(r − r₀)² in kJ/mol/nm²; AMBER files list K with no ½, like LAMMPS. From GROMACS: K_LAMMPS = k/2 ÷ 4.184 (and ÷ 100 for nm² → Å²).</li>
      <li><b>OPLS dihedral K's are the published ones</b> (the ½ is inside the formula). GROMACS uses Ryckaert–Bellemans for the same thing; converting needs the RB ↔ Fourier relation, not a factor.</li>
      <li><b>φ convention.</b> LAMMPS, like IUPAC, puts trans at 180° and cis at 0°; some papers measure from trans. A torsion that looks wrong by a sign or 180° usually is just this.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> E_angle for θ = 109.5° with K = 62.1, θ₀ = 114°?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Δθ = −4.5° = −0.078540 rad; E = 62.1 × 0.078540² = <b>0.38307 kcal/mol</b> per angle.</div></div>
    <div class="ex"><p><b>2.</b> E_OPLS at φ = 60° (gauche)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>½ [1.411 (1 + ½) − 0.271 (1 + ½) + 3.145 (1 − 1)] = ½ (2.1165 − 0.4065) = <b>0.855 kcal/mol</b> above trans.</div></div>
    <div class="ex"><p><b>3.</b> A GROMACS topology has angle k = 519.7 kJ/mol/rad². What do you put in angle_coeff (real units)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>519.7 / 2 / 4.184 = <b>62.1 kcal/mol/rad²</b>: the value used on this page.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pA = () => { const g = (s, d, lo, hi) => { const v = +$(s).value; return isFinite(v) && v >= lo && v <= hi ? v : d; }; return { r: g('.bd-r', 1.58, 1.2, 2.0), t: g('.bd-t', 109.5, 80, 170), p: g('.bd-p', 65, -180, 180) }; };
    const predA = () => { const p = pA(), x = butane(p.r, p.t, p.p), H = hand(x); $('.bd-pred').textContent = `E_bond = 3 × 260 × (${p.r} − 1.54)² = ${s17(H.eb)}\nE_angle = 2 × 62.1 × (${H.th[0].toFixed(6)}° − 114°)² in rad = ${s17(H.ea)}\nE_dihedral at φ = ${H.phi.toFixed(6)}° = ${s17(H.ed)}`; $('.bd-outa').innerHTML = ''; $('.bd-suma').innerHTML = ''; return { x, H }; };
    root.querySelectorAll('.bd-r,.bd-t,.bd-p').forEach(e => e.oninput = predA); predA();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.bd-sta'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const { x, H } = predA();
      const res = await runLammps(script('dump            d all custom 1 f.dump id fx fy fz\ndump_modify     d format float %.17g sort id'), { files: { 'b.data': DATA([x]) }, read: ['f.dump'] }); b.disabled = false;
      const t = res.files?.['f.dump']; if (res.error || !t || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const c = res.thermo.cols, r0 = res.thermo.rows[0], L = ['E_bond', 'E_angle', 'E_dihed'].map(k => parseFloat(r0[c.indexOf(k)]));
      const F = t.split('ITEM: ATOMS')[1].trim().split('\n').slice(1).map(l => l.trim().split(/\s+/).map(Number));
      let ok = 0, n = 0, exact = 0; const rows = [];
      [['E_bond', H.eb, L[0]], ['E_angle', H.ea, L[1]], ['E_dihed', H.ed, L[2]]].forEach(([k, h, l]) => { n++; const s = close(h, l); if (s) ok++; if (h === l) exact++; rows.push(row([k, `<span class="mono">${s17(h)}</span>`, `<span class="mono">${s17(l)}</span>`, h === l ? `${OK} bit for bit` : s ? `${OK} to ${(Math.abs(h - l) / Math.abs(l)).toExponential(0)}` : BAD])); });
      F.forEach(w => { for (let k = 0; k < 3; k++) { const h = H.f[w[0] - 1][k], l = w[1 + k]; n++; const s = close(h, l); if (s) ok++; if (h === l) exact++; rows.push(row([`f${'xyz'[k]} atom ${w[0]}`, `<span class="mono">${s17(h)}</span>`, `<span class="mono">${s17(l)}</span>`, h === l ? `${OK} bit for bit` : s ? `${OK} to ${(Math.abs(h - l) / Math.max(Math.abs(l), 1e-3)).toExponential(0)}` : BAD])); } });
      $('.bd-outa').innerHTML = `<table>${head(['', 'by hand', 'LAMMPS', ''])}${rows.join('')}</table>`;
      const fs = [0, 1, 2].map(k => H.f.reduce((s, q) => s + q[k], 0));
      $('.bd-suma').innerHTML = `<p>${ok} of ${n} match${exact === n ? ', every one bit for bit' : ` (${exact} bit for bit; the others differ only in the last digits, where the browser\'s acos, sin and cos round differently from JavaScript\'s)`}. Σ forces = (${fs.map(v => v.toExponential(1)).join(', ')}): bonded forces never push the molecule as a whole.</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.bd-stb'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const PH = Array.from({ length: 13 }, (_, i) => -180 + 30 * i), mols = PH.map((p, i) => butane(1.54, 114, p, [10 + 20 * i, 10, 10]));
      const res = await runLammps(script('compute         pd all pe/atom dihedral\ndump            d all custom 1 p.dump id mol c_pd\ndump_modify     d format float %.17g sort id'), { files: { 'b.data': DATA(mols) }, read: ['p.dump'] }); b.disabled = false;
      const t = res.files?.['p.dump']; if (res.error || !t) return err(st, res); st.textContent = '';
      const sum = Array(13).fill(0); for (const l of t.split('ITEM: ATOMS')[1].trim().split('\n').slice(1)) { const w = l.trim().split(/\s+/); sum[+w[1] - 1] += +w[2]; }
      let ok = 0;
      $('.bd-outb').innerHTML = `<table>${head(['φ (°)', 'E by hand', 'LAMMPS', '', ''])}${PH.map((p, i) => { const H = hand(mols[i]), s = Math.abs(H.ed - sum[i]) <= 1e-12 * Math.max(1, Math.abs(H.ed)); if (s) ok++; const bar = '█'.repeat(Math.round(H.ed / 0.25)); return row([p, H.ed.toFixed(10), sum[i].toFixed(10), s ? OK : BAD, `<span class="mono">${bar}</span>`]); }).join('')}</table>`;
      $('.bd-sumb').innerHTML = `<p>${ok} of 13 match. Trans (±180°) is 0, gauche (±60°) 0.855, the cis barrier (0°) ${hand(mols[6]).ed.toFixed(3)} kcal/mol. The quarter split of compute pe/atom adds the last-digit noise.</p>`;
    };

    $('.bd-conv').innerHTML = `<table>${head(['term', 'LAMMPS', 'GROMACS', 'from GROMACS to LAMMPS (real units)'])}
      ${row(['bond harmonic', 'K (r − r₀)², kcal/mol/Å²', '½ k (r − r₀)², kJ/mol/nm²', 'K = k / 2 / 4.184 / 100'])}
      ${row(['angle harmonic', 'K (θ − θ₀)², kcal/mol/rad²', '½ k (θ − θ₀)², kJ/mol/rad²', 'K = k / 2 / 4.184'])}
      ${row(['dihedral opls', '½ΣKₙ(1 ± cos nφ)', 'Ryckaert–Bellemans Cₙ cosⁿψ, ψ = φ − 180°', 'K₁ = −2C₁ − 3C₃/2, K₂ = −C₂ − C₄, K₃ = −C₃/2, K₄ = −C₄/4 (÷ 4.184)'])}</table>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
