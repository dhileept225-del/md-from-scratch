// RanMars (src/random_mars.cpp), Marsaglia's universal generator, as LAMMPS seeds and uses it.
// gaussian() is the polar Box-Muller method with the saved second value; it calls log(), so a browser whose
// Math.log is not correctly rounded can differ from native LAMMPS in the last bit of a draw.
export class RanMars {
  constructor(seed) {
    let ij = Math.trunc((seed - 1) / 30082), kl = (seed - 1) - 30082 * ij, i = Math.trunc(ij / 177) % 177 + 2, j = ij % 177 + 2, k = Math.trunc(kl / 169) % 178 + 1, l = kl % 169;
    this.u = new Float64Array(98);
    for (let ii = 1; ii <= 97; ii++) { let s = 0.0, t = 0.5; for (let jj = 1; jj <= 24; jj++) { const m = ((i * j) % 179) * k % 179; i = j; j = k; k = m; l = (53 * l + 1) % 169; if ((l * m) % 64 >= 32) s = s + t; t = 0.5 * t; } this.u[ii] = s; }
    this.c = 362436.0 / 16777216.0; this.cd = 7654321.0 / 16777216.0; this.cm = 16777213.0 / 16777216.0; this.i97 = 97; this.j97 = 33; this.save = 0; this.second = 0; this.uniform();
  }
  uniform() { let uni = this.u[this.i97] - this.u[this.j97]; if (uni < 0.0) uni += 1.0; this.u[this.i97] = uni; if (--this.i97 === 0) this.i97 = 97; if (--this.j97 === 0) this.j97 = 97; this.c -= this.cd; if (this.c < 0.0) this.c += this.cm; uni -= this.c; if (uni < 0.0) uni += 1.0; return uni; }
  gaussian(log = Math.log) {
    if (!this.save) {
      let v1, v2, rsq;
      do { v1 = 2.0 * this.uniform() - 1.0; v2 = 2.0 * this.uniform() - 1.0; rsq = v1 * v1 + v2 * v2; } while (rsq >= 1.0 || rsq === 0.0);
      const fac = Math.sqrt(-2.0 * log(rsq) / rsq);
      this.second = v1 * fac; this.save = 1; return v2 * fac;
    }
    this.save = 0; return this.second;
  }
}
