// Research · F.9: machine-learning potentials. A SNAP model fitted to Stillinger–Weber silicon with LAMMPS's own
// compute snap: the linear fit redone in the browser (QR least squares, energy/force weight live), the per-atom
// energy rebuilt as β·B to the last bit, and the deployed model tested natively in MD and outside its training range.
const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const nat = () => import('./ml-native.js').then(m => m.NATIVE);
const NAT = 64, K = 30;

const STATES = [['diamond, a = 5.431 Å', '300 K'], ['diamond, a = 5.431 Å', '1000 K'], ['diamond, a = 5.431 Å', '1600 K'], ['a = 5.431 Å, melted', '3000 K'], ['compressed, a = 5.30 Å', '800 K'], ['expanded, a = 5.55 Å', '800 K']];
const CSNAP = `# for each stored configuration (64 Si atoms), with the teacher as the active pair style:
pair_style sw
pair_coeff * * Si.sw Si
read_dump cfg.s2.dump 4400 x y z
compute sn all snap 3.77118 0.99363 6 0.5 1.0     # rcutfac rfac0 2J_max radius weight
run 0
# c_sn: 1 + 3N + 6 rows x (30 + 1) columns
#   row 0       : sum over atoms of B_k            | last column: teacher energy
#   rows 1..3N  : d(sum B_k)/d(r_j,a)              | last column: teacher force F_j,a
#   last 6 rows : virial of B_k                    | last column: teacher stress`;
const DEPLOY = `pair_style snap
pair_coeff * * lin.snapcoeff lin.snapparam Si`;

// Householder QR least squares: minimize |X b - y|, X m x n (row-major Float64Array)
function lstsq(X, y, m, n) {
  const A = Float64Array.from(X), b = Float64Array.from(y);
  for (let k = 0; k < n; k++) {
    let s = 0; for (let i = k; i < m; i++) s += A[i * n + k] ** 2;
    const nrm = Math.sqrt(s), alpha = A[k * n + k] > 0 ? -nrm : nrm;
    const v = new Float64Array(m - k); for (let i = k; i < m; i++) v[i - k] = A[i * n + k]; v[0] -= alpha;
    let vv = 0; for (let i = 0; i < v.length; i++) vv += v[i] * v[i]; if (vv === 0) continue;
    for (let j = k; j < n; j++) { let d = 0; for (let i = k; i < m; i++) d += v[i - k] * A[i * n + j]; d = 2 * d / vv; for (let i = k; i < m; i++) A[i * n + j] -= d * v[i - k]; }
    let d = 0; for (let i = k; i < m; i++) d += v[i - k] * b[i]; d = 2 * d / vv; for (let i = k; i < m; i++) b[i] -= d * v[i - k];
  }
  const x = new Float64Array(n); for (let k = n - 1; k >= 0; k--) { let s = b[k]; for (let j = k + 1; j < n; j++) s -= A[k * n + j] * x[j]; x[k] = s / A[k * n + k]; }
  return x;
}

export default {
  id: 'ml-potentials',
  title: 'Machine-learning potentials',
  minutes: 16,
  summary: 'A machine-learning potential is a fit: descriptors of each atom\'s neighbourhood, a model from descriptors to energy, and training data from a more expensive "teacher". This page builds the simplest complete example with LAMMPS alone. SNAP descriptors come from compute snap, the teacher is Stillinger–Weber silicon, and a least-squares fit is deployed as pair_style snap. You redo the fit in your browser and move the energy/force weight. You check that the deployed model computes exactly β·B, then see what MD with it gets right, what it gets wrong, and how it fails outside its data.',
  learn: [
    'Descriptors, model, data: SNAP as linear regression on bispectrum components',
    'compute snap: the design matrix LAMMPS builds for you',
    'Energy vs force weighting, linear vs quadratic, test-set errors',
    'Validation beyond errors: pressure, structure, extrapolation, cost',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: descriptors, model, data</h2>
    <div class="prose"><p>SNAP writes the energy of atom i as E<sub>i</sub> = β₀ + Σ<sub>k</sub> β<sub>k</sub> B<sub>ik</sub>, where B<sub>ik</sub> are the bispectrum components of i's neighbour density within r<sub>c</sub>: rotation-invariant numbers built from its expansion in 4-d hyperspherical harmonics. With 2J<sub>max</sub> = 6 there are 30 of them (${gh('src/ML-SNAP/pair_snap.cpp', 206)}). The model is linear in β, so the forces are too: F<sub>j</sub> = −Σ<sub>k</sub> β<sub>k</sub> ∂(Σ<sub>i</sub>B<sub>ik</sub>)/∂r<sub>j</sub>. Fitting is ordinary least squares. ACE, MACE, NequIP and Allegro change the descriptors (more body orders, learned message passing) and the model (nonlinear), but not the logic: descriptors, a model, and data from a teacher.</p></div>
    <pre class="lmp-mini">${esc(CSNAP)}</pre>
    <div class="prose"><p>The teacher here is Stillinger–Weber silicon (a three-body potential, cutoff 3.77 Å). We know it exactly, so every error below is the model's, not the data's. Training data: 66 configurations of 64 atoms, 11 frames from each of six states; 9 frames per state for training, 2 held out for testing.</p></div>
    <div class="tbl ml-states"></div>
  </section>

  <section class="block">
    <h2>2 · Fit it yourself</h2>
    <div class="prose"><p>Rows: one energy row per configuration, [1, ΣB/N] → E/N, multiplied by a weight w<sub>E</sub>; and force rows [0, ∂ΣB/∂r] → F. To keep the download small, your browser gets the force rows of 8 atoms per configuration (24 rows); native numpy solved the same problem. Solve by QR (Householder), predict the held-out configurations:</p></div>
    <div class="row"><label class="lc-in">w<sub>E</sub> <select class="lc-u ml-w"><option>1</option><option>10</option><option selected>100</option><option>1000</option><option>10000</option></select></label><button class="btn sm primary" data-k="a">▶ Fit 31 coefficients</button><span class="ml-sta"></span></div>
    <div class="prose ml-asum"></div>
    <div class="tbl ml-atab"></div>
  </section>

  <section class="block">
    <h2>3 · Deploy and check: the model is β·B, exactly</h2>
    <div class="prose"><p>The coefficients go into a .snapcoeff file (β₀ first, then β₁…β₃₀), the descriptor settings into a .snapparam file, and the model runs like any pair style:</p></div>
    <pre class="lmp-mini">${esc(DEPLOY)}</pre>
    <div class="row"><button class="btn sm ghost" data-k="f">Show the two files</button></div>
    <pre class="lmp-mini ml-files" hidden></pre>
    <div class="prose"><p>For a held-out configuration, native LAMMPS dumped each atom's 30 B<sub>ik</sub> (compute sna/atom) and its energy from pair snap (compute pe/atom). Rebuild every atom's energy from the coefficients, in the order pair_snap.cpp:209–217 adds them:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Rebuild 64 atomic energies</button></div>
    <div class="prose ml-bsum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: MD with the learned potential</h2>
    <div class="prose"><p>216 atoms of crystalline Si at 1000 K: 1 ps NVT to equilibrate, 4 ps NVT to average, then 2 ps NVE, natively, with the teacher and both models (the quadratic one adds all 465 products B<sub>ik</sub>B<sub>il</sub>: 496 coefficients).</p></div>
    <div class="tbl ml-md"></div>
    <div class="prose ml-mdsum"></div>
    <div class="prose"><p>Outside the data: energy per atom of the perfect diamond lattice versus lattice constant. The training states span 5.30–5.55 Å.</p></div>
    <div class="tbl ml-ev"></div>
    <div class="prose ml-evsum"></div>
  </section>

  <section class="block prose">
    <h2>5 · What carries over to MACE, NequIP and friends</h2>
    <ul><li><b>Test errors are necessary, not sufficient.</b> The linear model's force error is 14 % of the force spread, and its pressure is 36 % off. Validate the properties you will publish (pressure, density, diffusion, structure) in MD against the teacher or experiment.</li>
    <li><b>Fit what you will use.</b> The virial rows of compute snap were left out here on purpose. Include stresses in training if you will run NPT.</li>
    <li><b>Data defines the domain.</b> A model is an interpolator. Sample the states MD will visit (high T, compression, defects, surfaces), and watch for configurations outside the data. Active learning automates that.</li>
    <li><b>Weights move errors around.</b> No weight makes both errors small, and too much weight on few energies overfits them. Capacity and data do the real work (quadratic, more descriptors, a neural network; more configurations).</li>
    <li><b>Cost.</b> Even this small SNAP is 20–35× slower per step than SW. Equivariant networks are slower again and usually run on GPUs (KOKKOS, or their own LAMMPS pair styles: pair_nequip/pair_allegro from the NequIP developers, MACE through the ML-IAP interface).</li></ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Raising w<sub>E</sub> from 1 to 100 lowers the test energy error; raising it to 10 000 makes both errors worse. Why?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The 31 coefficients cannot represent SW exactly, so there is a residual, and the weights decide where it goes. With w<sub>E</sub> = 1 the 54 energy rows hardly count against 1296 force rows. With 10 000 they dominate, and 31 coefficients fitted mainly to 54 energies overfit them: the training energies improve, the held-out ones and the forces get worse. More capacity (the quadratic model, 496 coefficients, fitted to all forces) shrinks the residual itself: 0.87 meV/atom and 0.052 eV/Å at once.</div></div>
    <div class="ex"><p><b>2.</b> Both models agree with SW within a few meV near 5.43 Å but are off by 1–2 eV/atom at 4.6 Å. Would more training data at 5.43 Å help?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. At 4.6 Å neighbours sit at distances the training set never contained, and the model is extrapolating a polynomial in B. Only data in that region helps, or a physical baseline that is right there by construction (SNAP is often combined with ZBL repulsion for short distances).</div></div>
    <div class="ex"><p><b>3.</b> Your fitted coefficients agree with numpy's only to about 10 digits. Does that matter for MD?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. A least-squares solve amplifies round-off by the problem's condition number, so LAPACK and the QR here differ in the last digits of β; the predicted energies differ by far less than the fit error (meV). What matters is that the coefficients written to the file are the ones you validated, written with enough digits (%.17g here) that LAMMPS reads back exactly those doubles.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    $('.ml-states').innerHTML = `<table>${head(['state', 'structure', 'T', 'frames'])}${STATES.map((s, i) => row([`s${i + 1}`, s[0], s[1], '11 (9 train, 2 test)'])).join('')}</table>`;

    const decode = N => { const u = Uint8Array.from(atob(N.R), c => c.charCodeAt(0)); return new Float64Array(u.buffer); };
    const build = (R, N, ids, wE) => {
      const per = 1 + 3 * N.nsub, cols = N.ncol, X = [], y = [], isE = [];
      for (const c of ids) for (let r = 0; r < per; r++) {
        const off = (c * per + r) * cols;
        if (r === 0) { X.push(wE, ...Array.from({ length: K }, (_, k) => wE * R[off + k] / NAT)); y.push(wE * R[off + K] / NAT); isE.push(true); }
        else { X.push(0, ...Array.from({ length: K }, (_, k) => R[off + k])); y.push(R[off + K]); isE.push(false); }
      }
      return { X: Float64Array.from(X), y: Float64Array.from(y), isE, m: y.length };
    };
    const ids = Array.from({ length: 66 }, (_, i) => i), train = ids.filter(i => i % 11 < 9), test = ids.filter(i => i % 11 >= 9);
    const rmse = (R, N, b) => { const T = build(R, N, test, 1), e = [], f = []; for (let i = 0; i < T.m; i++) { let p = 0; for (let k = 0; k <= K; k++) p += T.X[i * (K + 1) + k] * b[k]; (T.isE[i] ? e : f).push(p - T.y[i]); } const r = a => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length); return [1000 * r(e), r(f)]; };
    const seen = new Map();

    $('[data-k=a]').onclick = async e => {
      const btn = e.currentTarget, st = $('.ml-sta'); btn.disabled = true; st.textContent = ' loading data…';
      const N = await nat(); if (!alive) return; const R = decode(N), wE = Number($('.ml-w').value);
      st.textContent = ' solving…'; await new Promise(r => setTimeout(r, 20));
      const T = build(R, N, train, wE), t0 = performance.now(), b = lstsq(T.X, T.y, T.m, K + 1), ms = performance.now() - t0;
      const [eE, eF] = rmse(R, N, b); seen.set(wE, [eE, eF]); btn.disabled = false; st.textContent = '';
      let cmp = '';
      if (wE === 100) { let mx = 0; for (let k = 0; k <= K; k++) mx = Math.max(mx, Math.abs(b[k] / N.bs[k] - 1)); cmp = ` Native numpy on the same rows: largest relative difference in a coefficient ${mx.toExponential(1)} ${mx < 1e-6 ? OK : BAD}.`; }
      $('.ml-asum').innerHTML = `<p>${T.m} rows × 31 columns solved in ${ms.toFixed(0)} ms. β₀ = ${b[0].toFixed(6)} eV, β₁ = ${b[1].toExponential(4)}.${cmp} Held-out configurations: energy RMSE ${eE.toFixed(2)} meV/atom, force RMSE ${eF.toFixed(3)} eV/Å (forces spread ${N.err.fstd.toFixed(2)} eV/Å). Natively, with all 192 force rows per configuration: linear ${N.err.lin[0].toFixed(2)} meV/atom and ${N.err.lin[1].toFixed(3)} eV/Å; quadratic ${N.err.quad[0].toFixed(2)} meV/atom and ${N.err.quad[1].toFixed(3)} eV/Å. Try other weights.</p>`;
      $('.ml-atab').innerHTML = `<table>${head(['w<sub>E</sub>', 'energy RMSE (meV/atom)', 'force RMSE (eV/Å)'])}${[...seen.entries()].sort((p, q) => p[0] - q[0]).map(([w, v]) => row([w, v[0].toFixed(2), v[1].toFixed(3)])).join('')}</table>`;
    };

    $('[data-k=f]').onclick = async e => { const btn = e.currentTarget, p = $('.ml-files'); const N = await nat(); p.textContent = `# lin.snapparam\n${N.param}\n# lin.snapcoeff (first lines; 31 coefficients in all)\n${N.coeffhead}…`; p.hidden = !p.hidden; btn.textContent = p.hidden ? 'Show the two files' : 'Hide the two files'; };

    $('[data-k=b]').onclick = async () => {
      const N = await nat(); if (!alive) return; const b = N.bl;
      let same = 0, mx = 0, tot = 0;
      for (const w of N.pa) { const B = w.slice(2).map(Number); let E = b[0]; for (let k = 0; k < K; k++) E += b[k + 1] * B[k]; const pe = Number(w[1]); if (E === pe) same++; mx = Math.max(mx, Math.abs(E - pe)); tot += E; }
      $('.ml-bsum').innerHTML = `<p>All ${N.pa.length} atomic energies agree with LAMMPS's pe/atom to ${mx.toExponential(1)} eV ${mx < 1e-13 ? OK : BAD}, ${same} of them to the last bit. The rest is round-off in the B's themselves: pair snap computes its bispectrum in its own pass (${gh('src/ML-SNAP/pair_snap.cpp', 110)}), separately from compute sna/atom, and the two need not add the neighbour contributions in the same order. Their sum, ${tot.toFixed(10)} eV, is the configuration's energy. The same number came out natively as β₀·N + β·ΣB from the fitting matrix: the model you fitted and the model LAMMPS runs are the same function.</p>`;
    };

    nat().then(N => {
      if (!alive) return;
      const M = N.md, nm = { sw: 'SW (teacher)', lin: 'SNAP linear', quad: 'SNAP quadratic' };
      $('.ml-md').innerHTML = `<table>${head(['', 'PE (eV/atom)', 'pressure (bar)', 'NVE drift, 2 ps (meV)', 'ms per step'])}${['sw', 'lin', 'quad'].map(p => row([nm[p], M[p].pe.toFixed(5), M[p].press.toFixed(0), (1000 * (Number(M[p].e1) - Number(M[p].e0))).toFixed(1), M[p].ms.toFixed(2)])).join('')}</table>`;
      const pk = p => { const g = N.rdf[p]; let i = 0; for (let k = 0; k < g.length; k++) if (g[k] > g[i]) i = k; return [N.rr[i], g[i]]; };
      $('.ml-mdsum').innerHTML = `<p>Energy: both within ${(1000 * Math.abs(M.lin.pe - M.sw.pe)).toFixed(1)} and ${(1000 * Math.abs(M.quad.pe - M.sw.pe)).toFixed(2)} meV/atom of SW. Pressure: the linear model is ${(100 * (M.lin.press / M.sw.press - 1)).toFixed(0)} % off, the quadratic ${(100 * (M.quad.press / M.sw.press - 1)).toFixed(1)} %. The stresses were not in the fit, and a 0.18 eV/Å force error is enough to spoil them. Structure: the first RDF peak is at ${pk('sw')[0].toFixed(3)} Å with height ${pk('sw')[1].toFixed(2)} (SW), ${pk('lin')[1].toFixed(2)} (linear), ${pk('quad')[1].toFixed(2)} (quadratic). Energy is conserved in NVE by all three, because SNAP is a smooth function of positions. The price per step: ${(M.lin.ms / M.sw.ms).toFixed(0)}× and ${(M.quad.ms / M.sw.ms).toFixed(0)}× SW.</p>`;
      $('.ml-ev').innerHTML = `<table>${head(['a (Å)', 'SW', 'linear', 'quadratic', 'in training range?'])}${N.ev.map(r => row([r[0], r[1].toFixed(4), r[2].toFixed(4), r[3].toFixed(4), r[0] >= 5.3 && r[0] <= 5.55 ? 'yes' : ''])).join('')}</table>`;
      const at = a => N.ev.find(r => r[0] === a);
      $('.ml-evsum').innerHTML = `<p>Inside the range the models are within a few meV of SW; the quadratic one stays close out to 6.4 Å. At 4.6 Å (15 % compression, neighbours at 2.0 Å instead of 2.35 Å) both are wrong by electron-volts and in opposite directions: linear ${at(4.6)[2].toFixed(2)}, quadratic ${at(4.6)[3].toFixed(2)}, SW ${at(4.6)[1].toFixed(2)} eV/atom. The quadratic model would bind atoms far too strongly at short distance, a hole in the potential that an energetic collision in MD could fall into. Nothing in its test error warned of it.</p>`;
    });

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
