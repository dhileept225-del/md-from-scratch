// D.14 by hand: MinCG::iterate + MinLineSearch::linemin_quadratic for a few LJ atoms in a 10σ box,
// with LAMMPS's own neighbour-list pair order (bins, half stencil, setup sort), so every number matches bit for bit.
export const BOX = 10.0, CUT = 2.5, SKIN = 0.3, CN = CUT + SKIN;
const CUTSQ = CUT * CUT, CNSQ = CN * CN, TRIG = (0.5 * SKIN) * (0.5 * SKIN);
// neighbour bins (NBinStandard::setup_bins) and the half stencil (NStencilBin<HALF,3D,ortho>::create)
const BS = 0.5 * CN, NB = Math.trunc(BOX * (1.0 / BS)), BSX = BOX / NB, BINV = 1.0 / BSX;
let SX = Math.trunc(CN * BINV); if (SX * BSX < CN) SX++;
const bd = i => (i > 0 ? (i - 1) * BSX : i === 0 ? 0.0 : (i + 1) * BSX);
export const STENCIL = [[0, 0, 0]];
for (let k = 0; k <= SX; k++) for (let j = -SX; j <= SX; j++) for (let i = -SX; i <= SX; i++) {
  if (k <= 0 && j <= 0 && (j !== 0 || i <= 0)) continue;
  if (bd(i) * bd(i) + bd(j) * bd(j) + bd(k) * bd(k) < CNSQ) STENCIL.push([i, j, k]);
}
const bin = x => [0, 1, 2].map(k => Math.min(Math.trunc(x[k] * BINV), NB - 1));
// Atom::sort at setup: its own bins, bininv = nbin/len
const SINV = Math.trunc(BOX * (1.0 / BS)) / BOX;
const sbin = x => { const b = [0, 1, 2].map(k => Math.max(0, Math.min(Math.trunc(x[k] * SINV), NB - 1))); return (b[2] * NB + b[1]) * NB + b[0]; };

export function makeSystem(P) {
  const n = P.length, order = P.map((p, i) => [sbin(p), i]).sort((a, b) => a[0] - b[0] || a[1] - b[1]).map(a => a[1]);
  const S = { n, N: 3 * n, ids: order.map(i => i + 1), x: order.flatMap(i => P[i].slice()), pairs: [], xhold: null, builds: 0 };
  build(S); return S;
}
function build(S) {
  const B = []; for (let i = 0; i < S.n; i++) B.push(bin(S.x.slice(3 * i, 3 * i + 3)));
  S.pairs = [];
  for (let i = 0; i < S.n; i++) for (let s = 0; s < STENCIL.length; s++) {
    const t = [0, 1, 2].map(k => B[i][k] + STENCIL[s][k]);
    for (let j = s === 0 ? i + 1 : 0; j < S.n; j++) {
      if (B[j][0] !== t[0] || B[j][1] !== t[1] || B[j][2] !== t[2]) continue;
      const dx = S.x[3 * i] - S.x[3 * j], dy = S.x[3 * i + 1] - S.x[3 * j + 1], dz = S.x[3 * i + 2] - S.x[3 * j + 2];
      if (dx * dx + dy * dy + dz * dz < CNSQ) S.pairs.push([i, j]);
    }
  }
  S.xhold = S.x.slice(); S.builds++;
}
// Min::energy_force: neighbor->decide() (every 1, delay 0, check yes), then PairLJCut::compute
export function energyForce(S) {
  for (let i = 0; i < S.n; i++) {
    const dx = S.x[3 * i] - S.xhold[3 * i], dy = S.x[3 * i + 1] - S.xhold[3 * i + 1], dz = S.x[3 * i + 2] - S.xhold[3 * i + 2];
    if (dx * dx + dy * dy + dz * dz > TRIG) { build(S); break; }
  }
  const f = new Array(S.N).fill(0.0); let e = 0.0;
  for (const [i, j] of S.pairs) {
    const dx = S.x[3 * i] - S.x[3 * j], dy = S.x[3 * i + 1] - S.x[3 * j + 1], dz = S.x[3 * i + 2] - S.x[3 * j + 2], rsq = dx * dx + dy * dy + dz * dz;
    if (rsq < CUTSQ) {
      const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48.0 * r6 - 24.0) * r2;
      f[3 * i] += dx * fp; f[3 * i + 1] += dy * fp; f[3 * i + 2] += dz * fp;
      f[3 * j] -= dx * fp; f[3 * j + 1] -= dy * fp; f[3 * j + 2] -= dz * fp;
      e += r6 * (4.0 * r6 - 4.0);
    }
  }
  S.f = f; S.e = 0.0 + e; return S.e;
}
const dot = (a, b) => { let s = 0.0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
// thermo keyword fnorm (Thermo::compute_fnorm): per-atom sums
export const thermoFnorm = f => { let s = 0.0; for (let i = 0; i < f.length; i += 3) s += f[i] * f[i] + f[i + 1] * f[i + 1] + f[i + 2] * f[i + 2]; return Math.sqrt(s); };

const ALPHA_MAX = 1.0, ALPHA_REDUCE = 0.5, BACKTRACK_SLOPE = 0.4, QUADRATIC_TOL = 0.1, EMACH = 1.0e-8, EPS_QUAD = 1.0e-28, DMAX = 0.1, EPS_ENERGY = 1.0e-8;

function lineminQuadratic(S, h, eoriginal, log, strict) {
  const fdothall = dot(S.f, h);
  if (fdothall <= 0.0) return 'DOWNHILL';
  let hmax = 0.0; for (const v of h) hmax = Math.max(hmax, Math.abs(v));
  const alphamax = Math.min(ALPHA_MAX, DMAX / hmax);
  if (hmax === 0.0) return 'ZEROFORCE';
  const x0 = S.x.slice();
  const step = (a, neval) => { if (a > 0.0) for (let i = 0; i < S.N; i++) S.x[i] = x0[i] + a * h[i]; else for (let i = 0; i < S.N; i++) S.x[i] = x0[i]; S.neval++; return energyForce(S); };
  let alpha = alphamax, fhprev = fdothall, engprev = eoriginal, alphaprev = 0.0, ecurrent;
  const trials = [];
  for (;;) {
    ecurrent = step(alpha);
    const fh = dot(S.f, h), delfh = fh - fhprev;
    if (Math.abs(fh) < EPS_QUAD || Math.abs(delfh) < EPS_QUAD) { step(0.0); return 'ZEROQUAD'; }
    const relerr = Math.abs(1.0 - (0.5 * (alpha - alphaprev) * (fh + fhprev) + ecurrent) / engprev), alpha0 = alpha - (alpha - alphaprev) * fh / delfh;
    const t = { alpha, de: ecurrent - eoriginal, fh, fhprev, relerr, alpha0 }; trials.push(t);
    if (relerr <= QUADRATIC_TOL && alpha0 > 0.0 && alpha0 < alphamax && !(strict === 'bracket' && fh < 0.0 && fhprev < 0.0)) {
      ecurrent = step(alpha0); t.de0 = ecurrent - eoriginal;
      if (strict === true ? ecurrent - eoriginal < 0.0 : ecurrent - eoriginal < EMACH) { t.took = 'alpha0'; if (log) log.push({ fdothall, alphamax, trials, x0, h: h.slice() }); return 0; }
    }
    const de_ideal = -BACKTRACK_SLOPE * alpha * fdothall, de = ecurrent - eoriginal;
    if (de <= de_ideal) { t.took = 'backtrack'; if (log) log.push({ fdothall, alphamax, trials, x0, h: h.slice() }); return 0; }
    fhprev = fh; engprev = ecurrent; alphaprev = alpha; alpha *= ALPHA_REDUCE;
    if (alpha <= 0.0 || de_ideal >= -EMACH) { step(0.0); if (log) log.push({ fdothall, alphamax, trials, x0, h: h.slice() }); return 'ZEROALPHA'; }
  }
}
// MinCG::iterate (Polak-Ribiere, reset every ndoftotal iterations). Returns rows [iter, pe, thermo fnorm] and the stop reason.
export function minimizeCG(P, { etol = 0.0, ftol = 1e-10, maxiter = 1000, maxeval = 10000, keepLog = false, strict = false } = {}) {
  const S = makeSystem(P); S.neval = 0; energyForce(S);
  const N = S.N, nlimit = N, log = keepLog ? [] : null;
  let h = S.f.slice(), g = S.f.slice(), gg = dot(S.f, S.f), niter = 0, stop = 'MAXITER';
  const rows = [[0, S.e, thermoFnorm(S.f)]], e0 = S.e, f0 = Math.sqrt(gg);
  for (let iter = 0; iter < maxiter; iter++) {
    niter++; const eprevious = S.e;
    const fail = lineminQuadratic(S, h, S.e, log, strict);
    if (fail) { stop = fail; break; }
    rows.push([niter, S.e, thermoFnorm(S.f)]);
    if (S.neval >= maxeval) { stop = 'MAXEVAL'; break; }
    if (Math.abs(S.e - eprevious) < etol * 0.5 * (Math.abs(S.e) + Math.abs(eprevious) + EPS_ENERGY)) { stop = 'ETOL'; break; }
    const d0 = dot(S.f, S.f), d1 = dot(S.f, g);
    if (ftol > 0.0 && d0 < ftol * ftol) { stop = 'FTOL'; break; }
    let beta = Math.max(0.0, (d0 - d1) / gg);
    if ((niter + 1) % nlimit === 0) beta = 0.0;
    gg = d0;
    for (let i = 0; i < N; i++) { g[i] = S.f[i]; h[i] = g[i] + beta * h[i]; }
    if (dot(g, h) <= 0.0) h = g.slice();
  }
  if (rows[rows.length - 1][0] !== niter) rows.push([niter, S.e, thermoFnorm(S.f)]);
  return { rows, stop, niter, neval: S.neval, e0, f0, e: S.e, fnorm2: Math.sqrt(dot(S.f, S.f)), builds: S.builds, log, S };
}
export const STOP_TEXT = { ETOL: 'energy tolerance', FTOL: 'force tolerance', MAXITER: 'max iterations', MAXEVAL: 'max force evaluations', ZEROALPHA: 'linesearch alpha is zero', ZEROQUAD: 'quadratic factors are zero', DOWNHILL: 'search direction is not downhill', ZEROFORCE: 'forces are zero' };
// seeded starts: n atoms within ±0.8σ of the box centre, no pair closer than 0.9σ, coordinates on a 0.001σ grid
export function randomStarts(count, seed = 1, n = 4, half = 0.8) {
  let a = seed >>> 0; const rnd = () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  const out = [];
  while (out.length < count) {
    const P = []; for (let i = 0; i < n; i++) P.push([0, 1, 2].map(() => Math.round((5.0 + (rnd() * 2 * half - half)) * 1000) / 1000));
    let ok = true; for (let i = 0; i < n && ok; i++) for (let j = i + 1; j < n; j++) { const d = [0, 1, 2].map(k => P[i][k] - P[j][k]); if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < 0.81) { ok = false; break; } }
    if (ok) out.push(P);
  }
  return out;
}
export const lmpScript = (P, { style = 'cg', line = '', etol = '0.0', ftol = '1.0e-10', maxiter = 1000, maxeval = 10000, thermo = 1 } = {}) => `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10 units box
create_box      1 box
${P.map(p => `create_atoms    1 single ${p.map(String).join(' ')} units box`).join('\n')}
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
thermo_style    custom step pe fnorm
thermo_modify   norm no format float %.17g
thermo          ${thermo}
${style === 'cg' ? '' : `min_style       ${style}\n`}${line ? `min_modify      line ${line}\n` : ''}minimize        ${etol} ${ftol} ${maxiter} ${maxeval}`;
