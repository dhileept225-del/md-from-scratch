import { tex } from '../tex.js';

export default {
  id: 'pme',
  summary: 'Particle-mesh Ewald makes Ewald summation fast enough for millions of atoms. Charges are spread onto a grid, the Fourier sum is done with fast Fourier transforms, and forces are interpolated back. The cost drops from N^1.5 (or worse) to N log N.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/pme.mp4', poster: 'videos/pme.jpg' },
  sim: 'pme',
  learn: [
    'The four PME steps: spread, FFT, multiply by the influence function, transform back',
    'B-spline interpolation, the grid spacing and the interpolation order',
    'Why PME scales as N log N, and the standard parameter choices (and PPPM)',
  ],
  simHint: 'Start with a 16³ grid and order 4, and read the relative error. Then try 8³ and 32³, and orders 2 to 8. Which change helps most? Compare the PME time with the direct Ewald time. The left picture shows how each ion\'s charge is smeared over nearby grid points; the right one is the smooth long-range potential that PME computes.',

  idea: `
<p>In plain Ewald (Lesson 5.2) the reciprocal sum needs the structure factor S(<b>k</b>) = Σ<sub>j</sub> q<sub>j</sub> e<sup>i<b>k</b>·<b>r</b><sub>j</sub></sup> for every wave vector: N atoms × N<sub>k</sub> vectors. For large systems that is the bottleneck.</p>
<p><b>The PME idea</b> (Darden, York and Pedersen 1993; smooth PME by Essmann et al. 1995): if the charges sat on a regular grid, S(<b>k</b>) for <i>all</i> k at once would be a single fast Fourier transform, costing K³ log K³. So:</p>
<ol>
<li><b>Spread</b>: assign each charge to the nearby grid points with smooth B-spline weights (order 4 means 4×4×4 points per atom).</li>
<li><b>FFT</b> the charge grid.</li>
<li><b>Multiply</b> by the Ewald influence function e<sup>−k²/4α²</sup>/k², corrected by a factor |b(<b>k</b>)|² that undoes the smoothing of the B-splines. The energy is the sum of this product with |FFT|².</li>
<li><b>Inverse FFT</b> to get the smooth potential on the grid, and <b>interpolate</b> the forces back to the atoms with the derivatives of the same B-splines.</li>
</ol>
<p>The real-space part stays as in Ewald (cutoff, neighbour lists). Because the FFT is cheap, α can be large and the real-space cutoff short (about 1 nm), and the total cost scales as N log N.</p>
<p><b>Accuracy</b> is controlled by the grid spacing (finer is better) and the interpolation order (higher is better). The defaults in most codes, 0.10–0.12 nm and order 4 with a relative tolerance of 10⁻⁵, give forces accurate to about 10⁻³ relative, which is plenty for MD. PPPM (Hockney and Eastwood, used by LAMMPS) is a closely related method with an optimised influence function.</p>
<div class="callout"><b>Key idea.</b> PME = Ewald with the Fourier sum done on a grid by FFTs. Spread → FFT → multiply → inverse FFT → interpolate. Finer grids and higher orders give higher accuracy; the cost grows only as N log N.</div>`,

  algorithm: `
<p>Smooth PME, reciprocal energy (forces follow from the same grids):</p>
<pre class="pseudo">u_i = K * r_i / L                         # scaled fractional coordinates
Q = zeros(K, K, K)
for each atom i:                           # spread with B-splines of order n
    for a, b, c in 0..n-1:
        Q[floor(u_i) - (a, b, c)] += q_i M_n(w_x + a) M_n(w_y + b) M_n(w_z + c)
F = FFT(Q)
for each grid frequency m ≠ 0:
    G(m) = exp(-π² |m/L|² / α²) / |m/L|² * |b_x|² |b_y|² |b_z|²
E_recip = Σ_m G(m) |F(m)|² / (2π V)
φ = inverse FFT(G * F)                     # smooth potential on the grid
F_i = -q_i Σ (gradient of the B-spline weights) * φ   # interpolate back</pre>
<ol class="steps">
  <li>Choose K as a product of small primes (2, 3, 5) for fast FFTs; codes do this automatically from the spacing.</li>
  <li>Spreading and interpolation cost N × n³; the FFT costs K³ log K³ with K³ ∝ N. Both are linear (or N log N).</li>
  <li>In parallel runs the FFT needs all-to-all communication; GROMACS can put PME on separate ranks or the GPU.</li>
  <li>GROMACS: <code>coulombtype = PME</code>, <code>fourierspacing = 0.12</code>, <code>pme-order = 4</code>, <code>ewald-rtol = 1e-5</code>. LAMMPS: <code>kspace_style pppm 1e-5</code> or <code>pppm/tip4p</code>.</li>
</ol>`,

  code: {
    note: 'A working smooth PME (reciprocal energy) with B-splines and FFTs. Part 1 compares it with the exact Ewald Fourier sum for 500 ions at several grid sizes and orders. Part 2 times both as N grows. Run <code>pip install numpy</code>, then <code>python pme.py</code> (about 15 s, mostly the slow direct sum).',
    src: `import numpy as np, time
from math import pi

def ewald_recip_direct(pos, q, L, alpha, kmax):
    """Exact reciprocal-space Ewald energy (Lesson 5.2), summed over all k with |m_i| <= kmax."""
    m = np.arange(-kmax, kmax + 1); M = np.stack(np.meshgrid(m, m, m, indexing="ij"), -1).reshape(-1, 3)
    M = M[(M != 0).any(1)]; k = 2 * pi * M / L; k2 = (k**2).sum(1)
    E = 0.0
    for s in range(0, len(k), 2000):                                 # chunks keep memory small
        S = (q[None, :] * np.exp(1j * k[s:s + 2000] @ pos.T)).sum(1)
        E += np.sum(np.exp(-k2[s:s + 2000] / (4 * alpha**2)) / k2[s:s + 2000] * np.abs(S)**2)
    return 2 * pi / L**3 * E

def bspline(w, n):
    """Cardinal B-spline weights M_n(w + j), j = 0..n-1, for fractional offsets w (array)."""
    M = np.zeros((len(w), n)); M[:, 0] = w; M[:, 1] = 1 - w          # order 2: M_2(w) = w, M_2(w + 1) = 1 - w
    for p in range(3, n + 1):                                        # recursion up to order n
        Mn = np.zeros_like(M)
        for j in range(p):
            x = w + j
            left = M[:, j] if j < p - 1 else 0.0
            right = M[:, j - 1] if j >= 1 else 0.0
            Mn[:, j] = (x * left + (p - x) * right) / (p - 1)
        M = Mn
    return M

def pme_recip(pos, q, L, alpha, K, n):
    """Smooth particle-mesh Ewald (Essmann et al. 1995), reciprocal energy only."""
    u = pos / L * K                                                  # scaled fractional coordinates
    Q = np.zeros((K, K, K))
    fl = np.floor(u).astype(int); w = u - fl
    W = [bspline(w[:, d], n) for d in range(3)]                     # weights per dimension
    for a in range(n):                                               # spread charges onto the grid
        for b in range(n):
            for c in range(n):
                idx = ((fl[:, 0] - a) % K, (fl[:, 1] - b) % K, (fl[:, 2] - c) % K)
                np.add.at(Q, idx, q * W[0][:, a] * W[1][:, b] * W[2][:, c])
    FQ = np.fft.fftn(Q)
    m = np.fft.fftfreq(K) * K                                        # integer frequencies
    # |b(m)|^2 factors that correct for the B-spline interpolation
    Mint = bspline(np.zeros(1), n)[0][1:]                            # M_n(1), ..., M_n(n-1)
    kk = np.arange(n - 1)
    denom = np.array([np.sum(Mint * np.exp(2j * pi * mi * kk / K)) for mi in m])
    b2 = 1.0 / np.abs(denom)**2
    B = b2[:, None, None] * b2[None, :, None] * b2[None, None, :]
    mx, my, mz = np.meshgrid(m, m, m, indexing="ij"); m2 = (mx**2 + my**2 + mz**2) / L**2
    m2[0, 0, 0] = 1.0
    G = np.exp(-pi**2 * m2 / alpha**2) / m2; G[0, 0, 0] = 0.0
    return np.sum(G * B * np.abs(FQ)**2) / (2 * pi * L**3)

rng = np.random.default_rng(1)
def system(N, rho=0.5):
    L = (N / rho) ** (1 / 3)
    return rng.uniform(0, L, (N, 3)), np.array([1.0, -1.0] * (N // 2)), L

# ---------- Part 1: accuracy vs grid size and interpolation order ----------
pos, q, L = system(500); alpha = 0.7
Eref = ewald_recip_direct(pos, q, L, alpha, kmax=12)
print(f"Part 1: 500 ions, box L = {L:.2f}, alpha = {alpha}. Exact reciprocal energy {Eref:.6f}")
print("   grid K    spacing    order 4 error    order 6 error    order 8 error")
for K in [8, 12, 16, 24, 32]:
    errs = [abs(pme_recip(pos, q, L, alpha, K, n) - Eref) / abs(Eref) for n in (4, 6, 8)]
    print(f"   {K:4d}      {L / K:.3f}      " + "     ".join(f"{e:10.1e}" for e in errs))

# ---------- Part 2: cost ----------
print("\\nPart 2: time for the reciprocal part (same accuracy target, density fixed)")
print("        N     direct Ewald (s)    PME (s)")
for N in [250, 1000, 4000, 16000]:
    pos, q, L = system(N); K = int(2 ** np.ceil(np.log2(L / 0.4))); kmax = int(np.ceil(L / 1.0))
    t0 = time.time(); pme_recip(pos, q, L, alpha, K, 4); tp = time.time() - t0
    if N <= 4000:
        t0 = time.time(); ewald_recip_direct(pos, q, L, alpha, kmax); td = time.time() - t0; tds = f"{td:10.3f}"
    else: tds = "  (too slow)"
    print(f"   {N:6d}       {tds}         {tp:7.3f}")
print("Direct Ewald grows roughly as N^2 here (k vectors x atoms); PME as N log N.")
`,
  },

  maths: () => `
<h4>1. Ewald reciprocal energy</h4>
${tex(String.raw`E_{\text{rec}} = \frac{1}{2\pi V}\sum_{\mathbf m\ne 0}\frac{e^{-\pi^2 m^2/\alpha^2}}{m^2}\,|S(\mathbf m)|^2, \qquad S(\mathbf m) = \sum_j q_j\, e^{2\pi i\,\mathbf m\cdot\mathbf r_j}`, true)}
<p>Here <b>m</b> are reciprocal lattice vectors (k = 2πm).</p>
<h4>2. B-spline interpolation of the exponentials</h4>
${tex(String.raw`e^{2\pi i m u/K} \approx b(m)\sum_{k} M_n(u - k)\, e^{2\pi i m k/K}, \qquad b(m) = \frac{e^{2\pi i (n-1) m/K}}{\sum_{k=0}^{n-2} M_n(k+1)\, e^{2\pi i m k/K}}`, true)}
<p>M<sub>n</sub> is the cardinal B-spline of order n (piecewise polynomial of degree n − 1). With this, S(<b>m</b>) ≈ b<sub>x</sub>b<sub>y</sub>b<sub>z</sub> · FFT(Q)(<b>m</b>).</p>
<h4>3. PME energy</h4>
${tex(String.raw`E_{\text{rec}} \approx \frac{1}{2\pi V}\sum_{\mathbf m\ne 0}\frac{e^{-\pi^2 m^2/\alpha^2}}{m^2}\,B(\mathbf m)\,\big|\hat Q(\mathbf m)\big|^2, \qquad B = |b_x|^2|b_y|^2|b_z|^2`, true)}
<h4>4. Cost</h4>
${tex(String.raw`\underbrace{N n^3}_{\text{spread + interpolate}} + \underbrace{K^3\log K^3}_{\text{FFT}} \;\propto\; N\log N \quad (K^3 \propto V \propto N)`, true)}`,

  quiz: [
    {
      q: 'What is the main reason PME is faster than plain Ewald summation?',
      options: ['It ignores long-range interactions', 'It computes the structure factor for all wave vectors at once with an FFT on a grid', 'It uses a smaller cutoff only', 'It uses single precision'],
      answer: 1,
      why: 'Plain Ewald needs N × N_k exponentials. PME spreads the charges to a grid, then one FFT gives all Fourier components in K³ log K³ operations.',
    },
    {
      q: 'Your PME forces are not accurate enough. Which two changes help most directly?',
      options: ['Increase the time step and the cutoff', 'Use a finer grid spacing and a higher interpolation order', 'Remove the self term', 'Use a smaller box'],
      answer: 1,
      why: 'The PME error comes from the grid interpolation. Smaller spacing and higher B-spline order both reduce it (see the code tab and the widget).',
    },
    {
      q: 'Which part of the calculation does PME change compared with Ewald?',
      options: ['Only the real-space part', 'Only the reciprocal-space (Fourier) part', 'The self-energy only', 'The Lennard-Jones part'],
      answer: 1,
      why: 'The real-space erfc sum and the self term are the same as in Ewald. PME replaces the reciprocal sum with a grid-based FFT calculation.',
    },
  ],
};
