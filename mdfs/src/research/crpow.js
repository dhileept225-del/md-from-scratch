// Correctly rounded x^(p/q) for x > 0 and small integers p, q (q > 0), with BigInt arithmetic.
// glibc's pow() is correctly rounded in practice; browsers' Math.pow is not always, so hand checks use this.
const decompose = x => {   // x = M · 2^E, M an integer below 2^53
  const b = new DataView(new ArrayBuffer(8)); b.setFloat64(0, x);
  const hi = b.getUint32(0), lo = b.getUint32(4), ex = (hi >>> 20) & 0x7ff;
  let M = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo), E;
  if (ex === 0) E = -1074; else { M |= 1n << 52n; E = ex - 1075; }
  return [M, E];
};
const bitlen = n => n.toString(2).length;
// round the positive rational num/den to the nearest double (ties to even)
export function ratioToDouble(num, den) {
  let sh = 0;
  const target = 60;                               // keep 60 significant bits, the rest as a sticky bit
  const lb = bitlen(num) - bitlen(den);
  sh = target - lb;
  const N = sh >= 0 ? num << BigInt(sh) : num, D = sh >= 0 ? den : den << BigInt(-sh);
  let q = N / D; const sticky = N % D !== 0n;
  const ql = bitlen(q), drop = ql - 53;
  let mant = q >> BigInt(drop); const rem = q & ((1n << BigInt(drop)) - 1n), half = 1n << BigInt(drop - 1);
  if (rem > half || (rem === half && (sticky || (mant & 1n)))) mant += 1n;
  else if (rem === half && !sticky && !(mant & 1n)) { /* tie, even: keep */ }
  return Number(mant) * Math.pow(2, drop - sh);   // exact: mant < 2^54 and the power of two is exact
}
function iroot(n, q) {                             // floor(n^(1/q))
  if (n < 2n) return n;
  let x = 1n << BigInt(Math.ceil(bitlen(n) / q) + 1);
  for (;;) { const y = ((BigInt(q) - 1n) * x + n / x ** BigInt(q - 1)) / BigInt(q); if (y >= x) return x; x = y; }
}
export function crpow(x, p, q = 1) {
  if (x === 1) return 1;
  if (p < 0) { const [M, E] = decompose(x); return crpowParts(M, E, -p, q, true); }
  const [M, E] = decompose(x); return crpowParts(M, E, p, q, false);
}
function crpowParts(M, E, p, q, inv) {
  // y = (M^p 2^(E p))^(1/q); scale by 2^(q S) to get S extra bits from the integer root
  const S = 200 + Math.max(0, Math.ceil(-E * p / q)), num = M ** BigInt(p), e = E * p + q * S;   // e >= 0
  const v = num << BigInt(e), r = iroot(v, q), sticky = r ** BigInt(q) !== v;
  // y ≈ (r + sticky·ε) / 2^S ; one extra low bit makes the sticky information survive rounding
  const R = (r << 1n) | (sticky ? 1n : 0n), den = 1n << BigInt(S + 1);
  return inv ? ratioToDouble(den, R) : ratioToDouble(R, den);
}

// Correctly rounded natural logarithm (x > 0): ln(M/2^52) by the atanh series in 200-bit fixed point, plus e·ln 2.
const P = 200n, ONE = 1n << P;
function lnFixed(num, den) {                       // ln(num/den) · 2^200, for num/den near 1..2
  const z = ((num - den) << P) / (num + den), z2 = (z * z) >> P;
  let term = z, sum = 0n, k = 1n;
  while (term !== 0n) { sum += term / k; term = (term * z2) >> P; k += 2n; }
  return 2n * sum;
}
const LN2 = lnFixed(2n, 1n);
export function crlog(x) {
  if (x === 1) return 0;
  let [M, E] = decompose(x);
  while (M < (1n << 52n)) { M <<= 1n; E -= 1; }     // subnormals
  const L = lnFixed(M, 1n << 52n) + BigInt(E + 52) * LN2;
  if (L === 0n) return 0;
  return L < 0n ? -ratioToDouble(-L, ONE) : ratioToDouble(L, ONE);
}
