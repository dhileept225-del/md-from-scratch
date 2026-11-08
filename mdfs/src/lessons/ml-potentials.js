import { tex } from '../tex.js';

export default {
  id: 'ml-potentials',
  summary: 'Machine-learning interatomic potentials (MLIPs) learn the energy and forces of a quantum-mechanical method such as DFT and reproduce them at a tiny fraction of the cost. They write the energy as a sum of atomic contributions that depend on each atom\'s neighbourhood, are fitted to energies and forces, and must be trained on data that cover every configuration the simulation will visit.',
  minutes: 18,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/ml-potentials.mp4', poster: 'videos/ml-potentials.jpg' },
  sim: 'mlpot',
  learn: [
    'Why MLIPs exist: near-DFT accuracy at a cost close to classical force fields',
    'The atomic-energy decomposition, descriptors and symmetry, and fitting to energies and forces',
    'Extrapolation, uncertainty and active learning, and the current toolbox (MACE, NequIP, ACE, foundation models)',
  ],
  simHint: 'The map is a 2D energy surface and the dots are the reference ("DFT") calculations, all in one basin. Watch the walker running on the model: it leaves the basin and finds regions where the model is wrong (⚠ in the readout; check the Error view). Press "Active-learning step" a few times: each one computes the reference where the model is least sure along the path, and the RMS error drops. Compare with tapping points yourself, and try a shorter or longer kernel length.',

  idea: `
<p>DFT gives accurate energies and forces, but its cost grows roughly as N³ and a single MD step for a few hundred atoms can take minutes. Classical force fields are fast but have a fixed functional form that often misses many-body effects (polarization, bond breaking, unusual coordination). <b>Machine-learning interatomic potentials</b> (MLIPs) aim for both: train a flexible model on DFT data, then run MD with the model.</p>
<p><b>Atomic energies.</b> Almost all MLIPs write the total energy as a sum over atoms, E = Σ<sub>i</sub> ε<sub>i</sub>, where ε<sub>i</sub> depends only on the positions of neighbours within a cutoff (typically 5–6 Å). This makes the cost linear in N, like a classical force field, and lets a model trained on small cells be used for large ones.</p>
<p><b>Descriptors and symmetry.</b> The model must give the same energy if the system is translated, rotated or if identical atoms are swapped. Early models (Behler–Parrinello neural networks, GAP/SOAP, SNAP, ACE) feed each atom's neighbourhood through hand-built invariant <b>descriptors</b> (radial and angular functions) into a neural network or a linear/kernel model. Modern <b>equivariant message-passing networks</b> (NequIP, Allegro, MACE) learn the descriptors themselves, passing geometric information between neighbouring atoms with built-in rotational symmetry; they are among the most accurate and data-efficient.</p>
<p><b>Fitting to forces.</b> Forces are the negative gradient of the model energy, so they come for free and are energy-conserving. Each configuration gives one energy but 3N force components, so training on forces is far more informative (code tab, Part 1: a linear model with angular descriptors fits a three-body reference 4× better than a pair potential).</p>
<p><b>Extrapolation is the main danger.</b> An MLIP has no physics outside its training data. In unseen configurations (higher temperature, compression, a new chemical environment) errors can grow by an order of magnitude (Part 2), and a simulation can wander into unphysical "holes" where the model predicts a low energy. Two defences: training data that cover what MD will visit, and an <b>uncertainty estimate</b> (committee spread, Gaussian-process variance, or MACE/NequIP ensembles). <b>Active learning</b> runs MD with the current model, picks the configurations where it is least certain, computes DFT there, and retrains (Part 3 and the widget).</p>
<p><b>Foundation models</b> (for example MACE-MP-0 and similar models trained on large materials databases) work out of the box for many systems and are a good starting point for fine-tuning on your own DFT data.</p>
<div class="callout"><b>Key idea.</b> E = Σ atomic energies of local environments, learned from DFT energies and forces. Validate on unseen data, watch the uncertainty, and grow the training set where the model is unsure.</div>`,

  algorithm: `
<pre class="pseudo">1. initial data : DFT on a few hundred diverse configurations (AIMD snippets, distorted cells, different T/ρ)
2. train        : minimize  Σ (E_model - E_DFT)²/N² · w_E  +  Σ |F_model - F_DFT|² · w_F   (+ stress)
3. validate     : held-out test set; parity plots for E, F; also physical tests (g(r), density, D, phonons)
4. explore      : run MD with the model at the conditions you care about
5. select       : configurations with high uncertainty (committee spread / model std)
6. label        : DFT on the selected configurations; add to the training set; go to 2
7. production   : long MD with the final model; monitor the uncertainty during the run</pre>
<ol class="steps">
  <li>DFT settings (functional, dispersion correction, basis, k-points) define what the model learns; keep them identical for all data.</li>
  <li>Training codes: <code>mace</code>, <code>nequip</code>/<code>allegro</code>, <code>pacemaker</code> (ACE), <code>FitSNAP</code>, DeePMD-kit. Active-learning helpers: FLARE and committee scripts.</li>
  <li>Running: LAMMPS (<code>pair_style mace</code>, <code>pair_style allegro</code>, ML-PACE, ML-SNAP), ASE calculators, OpenMM-ML, and GROMACS through its neural-network potential interface.</li>
  <li>Watch for long-range physics: plain local models miss electrostatics beyond the cutoff; for ionic systems consider models with explicit long-range terms or charges.</li>
</ol>`,

  code: {
    note: 'A complete small MLIP pipeline in NumPy. The reference ("DFT") is a pair plus three-body potential; 64-atom liquid snapshots are labelled with energies and forces; linear models on Behler–Parrinello-style radial and angular descriptors are fitted to both. Part 1 compares a pair-only model with the many-body model, Part 2 tests extrapolation, and Part 3 shows that a committee flags the unseen configurations. Run <code>pip install numpy</code>, then <code>python ml_potential.py</code> (about 30 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(8)
RC, A3 = 2.0, 1.8                  # cutoff of everything; range of the three-body term

# ---------- the 'expensive' reference (stands in for DFT): pair term + Stillinger-Weber-like angular term ----------
def fc(r, rc=RC):                  # smooth cutoff: goes to zero with zero slope at rc
    return np.where(r < rc, 0.5 * (np.cos(np.pi * r / rc) + 1), 0.0)

def neighbours(x, L):
    d = x[None, :, :] - x[:, None, :]; d -= L * np.round(d / L)            # d[i, j] = r_j - r_i (minimum image)
    r = np.sqrt((d**2).sum(-1)); np.fill_diagonal(r, np.inf)
    return d, r

def ref_energy(x, L, lam=10.0):
    d, r = neighbours(x, L)
    iu = np.triu_indices(len(x), 1); rr = r[iu]
    E = np.sum(4 * (rr**-12 - rr**-6) * fc(rr))                            # pair part
    for i in range(len(x)):                                                # three-body part: prefers ~tetrahedral angles
        nb = np.where(r[i] < A3)[0]
        if len(nb) < 2: continue
        v = d[i, nb]; rn = r[i, nb]; f = np.exp(1.2 / (rn - A3))
        cos = (v @ v.T) / np.outer(rn, rn); ju, ku = np.triu_indices(len(nb), 1)
        E += lam * np.sum(f[ju] * f[ku] * (cos[ju, ku] + 1 / 3) ** 2)
    return E

# ---------- descriptors: what an atom 'sees' within the cutoff (Behler-Parrinello style) ----------
MU = np.linspace(0.9, 1.9, 10); GAM = [0.8, 1.6, 3.2]
def descriptors(x, L):
    """Sum over atoms of [r^-12, r^-6, radial Gaussians (10) | angular sums g(r_ij) g(r_ik) cos^n (3 widths x n = 0, 1, 2)].
    The two power laws act as a physical 'baseline' for the steep repulsion, as many MLIPs add one."""
    d, r = neighbours(x, L)
    iu = np.triu_indices(len(x), 1); rr = r[iu]
    rad = 2 * np.r_[np.sum(rr**-12 * fc(rr)), np.sum(rr**-6 * fc(rr)),
                    (np.exp(-10 * (rr[:, None] - MU[None, :]) ** 2) * fc(rr)[:, None]).sum(0)]   # each pair seen by both atoms
    ang = np.zeros(9)
    for i in range(len(x)):
        nb = np.where(r[i] < A3)[0]
        if len(nb) < 2: continue
        v = d[i, nb]; rn = r[i, nb]; cos = (v @ v.T) / np.outer(rn, rn); ju, ku = np.triu_indices(len(nb), 1)
        for m, g in enumerate(GAM):
            f = np.exp(g / (rn - A3)); w = f[ju] * f[ku]
            ang[3 * m:3 * m + 3] += [w.sum(), (w * cos[ju, ku]).sum(), (w * cos[ju, ku] ** 2).sum()]
    return np.r_[rad, ang]

def fd_grad(fun, x, L, atoms, h=1e-4):
    """Central finite differences of fun (scalar or vector) w.r.t. the coordinates of the chosen atoms.
    (Real MLIP codes use analytic gradients or automatic differentiation; finite differences keep this short.)"""
    out = []
    for i in atoms:
        for a in range(3):
            xp, xm = x.copy(), x.copy(); xp[i, a] += h; xm[i, a] -= h
            out.append((np.asarray(fun(xp, L)) - np.asarray(fun(xm, L))) / (2 * h))
    return np.array(out)                                                    # (3 * len(atoms), ...)

# ---------- training and test data: snapshots of an LJ liquid, labelled with the reference ----------
def lj_snapshots(n_snap, rho, T, N=64, every=40):
    n = round((N / 4) ** (1 / 3)); a = (4 / rho) ** (1 / 3); L = n * a
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]); g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    x = ((g[:, None, :] + b[None]).reshape(-1, 3) + 0.1) * a; v = rng.normal(0, np.sqrt(T), x.shape); out = []
    def frc(x):
        dd, r = neighbours(x, L); m = r < 2.5; ir6 = np.where(m, r**-6.0, 0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r**2, 0)
        return -(f[:, :, None] * dd).sum(1)
    f = frc(x); dt = 0.004
    for s in range(400 + n_snap * every):
        v += 0.5 * dt * f; x = (x + dt * v) % L; f = frc(x); v += 0.5 * dt * f
        v *= np.sqrt(T * (3 * N - 3) / (v**2).sum())
        if s >= 400 and (s - 400) % every == 0: out.append(x.copy())
    return out, L

def label(snaps, L, n_force_atoms=6):
    data = []
    for x in snaps:
        atoms = rng.choice(len(x), n_force_atoms, replace=False)           # forces on a few atoms per snapshot (enough here)
        data.append(dict(E=ref_energy(x, L), F=-fd_grad(ref_energy, x, L, atoms), D=descriptors(x, L),
                         dD=fd_grad(descriptors, x, L, atoms), N=len(x)))
    return data

t0 = time.time()
train = label(*lj_snapshots(30, 0.90, 1.0)); test = label(*lj_snapshots(10, 0.90, 1.0))
hot = label(*lj_snapshots(8, 0.90, 2.5)); dense = label(*lj_snapshots(8, 1.10, 1.0))
print(f"Data: 64-atom snapshots labelled with the reference: 30 train, 10 test, 8 hot (T* 2.5), 8 compressed (rho* 1.1)  ({time.time() - t0:.0f} s)")

# ---------- Part 1: fit a linear model E = w . D (+ N b) to energies AND forces ----------
def fit(data, cols, wF=1.0, ridge=1e-6):
    rows, y = [], []
    for s in data:
        rows.append(np.r_[s["D"][cols], s["N"]] / s["N"]); y.append(s["E"] / s["N"])          # energy per atom
        for k in range(len(s["F"])): rows.append(wF * np.r_[-s["dD"][k][cols], 0]); y.append(wF * s["F"][k])
    A, y = np.array(rows), np.array(y)
    return np.linalg.solve(A.T @ A + ridge * np.eye(A.shape[1]), A.T @ y)

def errors(w, data, cols):
    eE = [(np.r_[s["D"][cols], s["N"]] @ w - s["E"]) / s["N"] for s in data]
    eF = np.concatenate([-s["dD"][:, cols] @ w[:-1] - s["F"] for s in data])
    return np.sqrt(np.mean(np.square(eE))), np.sqrt(np.mean(eF**2)), np.sqrt(np.mean(np.concatenate([s["F"] for s in data])**2))

pair_cols, all_cols = np.arange(12), np.arange(21)
print("\\nPart 1: RMS errors (energy per atom and force components, in epsilon and epsilon/sigma)")
print("   model                         test E     test F     (typical |F| on test set)")
for name, cols in [("pair only (classical FF)", pair_cols), ("pair + angular descriptors", all_cols)]:
    w = fit(train, cols); eE, eF, Fs = errors(w, test, cols)
    print(f"   {name:28s}  {eE:8.4f}   {eF:8.4f}     ({Fs:.1f})")
w = fit(train, all_cols)

# ---------- Part 2: extrapolation: configurations unlike the training data ----------
print("\\nPart 2: the many-body model on data it was not trained on")
for name, data in [("test (same T, rho)", test), ("hot liquid, T* = 2.5", hot), ("compressed, rho* = 1.1", dense)]:
    eE, eF, _ = errors(w, data, all_cols); print(f"   {name:24s} E error {eE:7.4f}   F error {eF:7.4f}")
w2 = fit(train + dense[:4], all_cols); eE, eF, _ = errors(w2, dense[4:], all_cols)
print(f"   after adding 4 compressed snapshots to training: compressed E error {eE:7.4f}, F error {eF:7.4f}")

# ---------- Part 3: a committee flags where the model is unsure (the basis of active learning) ----------
print("\\nPart 3: committee of 5 models, each trained on a random half of the training set")
committee = [fit([train[k] for k in rng.choice(len(train), 15, replace=False)], all_cols) for _ in range(5)]
for name, data in [("test", test), ("hot", hot), ("compressed", dense)]:
    spread = np.mean([np.std([(np.r_[s["D"], s["N"]] @ c) / s["N"] for c in committee]) for s in data])
    print(f"   {name:10s}: committee disagreement in E/N = {spread:.4f}")
print("   large disagreement -> run the reference calculation there and add it to the training set")
`,
  },

  maths: () => `
<h4>1. Atomic decomposition</h4>
${tex(String.raw`E = \sum_i \varepsilon\big(\{\mathbf r_{ij}\}_{r_{ij} < r_c}\big), \qquad \mathbf F_k = -\nabla_{\mathbf r_k} E`, true)}
<h4>2. Behler–Parrinello symmetry functions</h4>
${tex(String.raw`G^{\text{rad}}_i = \sum_j e^{-\eta(r_{ij} - \mu)^2} f_c(r_{ij}), \qquad G^{\text{ang}}_i = \sum_{j<k} g(r_{ij})\,g(r_{ik})\,\cos^n\theta_{jik}`, true)}
<p>with a smooth cutoff f<sub>c</sub>(r) = ½[cos(πr/r<sub>c</sub>) + 1]. They are invariant to translation, rotation and permutation of neighbours.</p>
<h4>3. Loss function (energies and forces)</h4>
${tex(String.raw`\mathcal L = w_E\sum_{\text{configs}}\Big(\frac{E^{\text{model}} - E^{\text{ref}}}{N}\Big)^2 + w_F\sum_{\text{configs}}\sum_k\big|\mathbf F_k^{\text{model}} - \mathbf F_k^{\text{ref}}\big|^2`, true)}
<h4>4. Message passing (MACE/NequIP idea)</h4>
${tex(String.raw`\mathbf h_i^{(t+1)} = U\Big(\mathbf h_i^{(t)},\; \sum_{j\in\mathcal N(i)} M\big(\mathbf h_i^{(t)}, \mathbf h_j^{(t)}, \mathbf r_{ij}\big)\Big), \qquad \varepsilon_i = R\big(\mathbf h_i^{(T)}\big)`, true)}
<p>Features h carry scalars, vectors and higher-order tensors that rotate with the system (equivariance), so geometry is used without hand-made descriptors. T layers enlarge the receptive field to T × r<sub>c</sub>.</p>
<h4>5. Gaussian-process uncertainty (used in the widget)</h4>
${tex(String.raw`\mu(\mathbf x) = m + \mathbf k_*^{\top}K^{-1}(\mathbf y - m), \qquad \sigma^2(\mathbf x) = k(\mathbf x, \mathbf x) - \mathbf k_*^{\top}K^{-1}\mathbf k_*`, true)}`,

  quiz: [
    {
      q: 'Why do MLIPs write the total energy as a sum of atomic energies that depend on local neighbourhoods?',
      options: ['Because DFT does the same', 'It makes the cost linear in the number of atoms and lets models trained on small cells run on large ones', 'To include long-range electrostatics exactly', 'It is required for energy conservation'],
      answer: 1,
      why: 'Locality (a cutoff) gives O(N) cost and transferability across system sizes. The price is that physics beyond the cutoff, such as long-range electrostatics, needs extra care.',
    },
    {
      q: 'Your MLIP has a tiny test error, but an MD run at a higher temperature suddenly produces atoms that collapse together. The most likely cause?',
      options: ['A bug in the integrator', 'The run reached configurations outside the training data, where the model extrapolates badly', 'The test set was too large', 'The time step was too small'],
      answer: 1,
      why: 'Low test error only holds for configurations like the training data. At higher temperature atoms visit new regions (e.g. short distances) where the model may predict unphysically low energies. Uncertainty monitoring and active learning address this.',
    },
    {
      q: 'In active learning, which configurations should get new DFT calculations?',
      options: ['Random ones', 'Those where the model (or committee) is most uncertain', 'Those with the lowest energy', 'Only the starting structure'],
      answer: 1,
      why: 'High uncertainty marks where the model lacks information. Labelling those configurations improves the model fastest and keeps the number of expensive DFT calculations small.',
    },
  ],
};
