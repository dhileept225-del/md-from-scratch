// Research · A.5: building the system with lattice, region, create_box and create_atoms; a lattice builder checked against LAMMPS.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

const LAT = {
  sc: { dim: 3, n: 1, basis: [[0, 0, 0]], nn: a => a, nnTxt: 'a', cell: [1, 1, 1] },
  bcc: { dim: 3, n: 2, basis: [[0, 0, 0], [.5, .5, .5]], nn: a => a * Math.sqrt(3) / 2, nnTxt: 'a√3/2', cell: [1, 1, 1] },
  fcc: { dim: 3, n: 4, basis: [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]], nn: a => a / Math.SQRT2, nnTxt: 'a/√2', cell: [1, 1, 1] },
  diamond: { dim: 3, n: 8, basis: [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0], [.25, .25, .25], [.25, .75, .75], [.75, .25, .75], [.75, .75, .25]], nn: a => a * Math.sqrt(3) / 4, nnTxt: 'a√3/4', cell: [1, 1, 1] },
  sq: { dim: 2, n: 1, basis: [[0, 0, 0]], nn: a => a, nnTxt: 'a', cell: [1, 1, 1] },
  hex: { dim: 2, n: 2, basis: [[0, 0, 0], [.5, .5, 0]], nn: a => a, nnTxt: 'a', cell: [1, Math.sqrt(3), 1] },
};

function predict(style, rho, cells) {
  const L = LAT[style];
  const a = L.dim === 3 ? Math.cbrt(L.n / rho) : Math.sqrt(L.n / (rho * L.cell[1]));
  const sp = [a * L.cell[0], a * L.cell[1], L.dim === 3 ? a : a];
  const box = [cells[0] * sp[0], cells[1] * sp[1], L.dim === 3 ? cells[2] * sp[2] : null];
  const cnt = c => { const out = []; for (let i = 0; i <= Math.ceil(c) + 1; i++) out.push(i); return out; };
  let N = 0; const eps = 1e-10;
  const nz = L.dim === 3 ? cells[2] : 1;
  for (const i of cnt(cells[0])) for (const j of cnt(cells[1])) for (const k of (L.dim === 3 ? cnt(nz) : [0])) for (const b of L.basis) {
    const x = i + b[0], y = j + b[1], z = k + b[2];
    if (x >= -eps && x < cells[0] - eps && y >= -eps && y < cells[1] - eps && (L.dim === 2 ? z === 0 : z >= -eps && z < nz - eps)) N++;
  }
  return { a, sp, box, N, nn: L.nn(a) };
}

const script = (style, rho, cells, dim) => `units           lj
${dim === 2 ? 'dimension       2\n' : ''}atom_style      atomic
lattice         ${style} ${rho}
region          box block 0 ${cells[0]} 0 ${cells[1]} ${dim === 2 ? '-0.5 0.5' : `0 ${cells[2]}`}
create_box      1 box
create_atoms    1 box
`;

export default {
  id: 'building-the-system',
  title: 'Building the system',
  minutes: 20,
  summary: 'Four lines turn a number into a crystal: lattice, region, create_box and create_atoms. You derive the lattice spacing, box size and atom count by hand for any lattice, then let LAMMPS build it and compare its log line by line, including a box size that does not fit the crystal.',
  learn: [
    'What lattice, region, create_box and create_atoms each do, and in which units',
    'Predict the lattice spacing, box size and number of atoms for sc, bcc, fcc, diamond and 2D lattices',
    'See what happens when the box is not a whole number of lattice cells',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: four lines, one crystal</h2>
    <div class="tbl wide"><table>${head(['Line', 'What it does', 'In melt'])}
${row(['<code>lattice fcc 0.8442</code>', 'describes a pattern of points: the FCC cube with 4 atoms; in LJ units the number is the density, so the cube edge is a = (4/0.8442)<sup>1/3</sup>', 'a = 1.6795962'])}
${row(['<code>region box block 0 10 0 10 0 10</code>', 'a named shape; numbers are counted in lattice spacings', 'a cube 16.795962 on each side'])}
${row(['<code>create_box 1 box</code>', 'makes the simulation box with that shape and 1 atom type', 'box (0,0,0) to (16.795962, …)'])}
${row(['<code>create_atoms 1 box</code>', 'puts an atom of type 1 on every lattice point inside the region', '4 × 10³ = 4000 atoms'])}</table></div>
    <p>In LJ units the lattice number is a <b>density</b>; in every other unit style it is the <b>lattice constant</b> itself (for copper in metal units, <code>lattice fcc 3.615</code> means a = 3.615 Å). LAMMPS confirms what it built in three log lines:</p>
    <div class="calc">Lattice spacing in x,y,z = 1.6795962 1.6795962 1.6795962
Created orthogonal box = (0 0 0) to (16.795962 16.795962 16.795962)
Created 4000 atoms</div>
  </section>

  <section class="block prose">
    <h2>2 · The rules by hand</h2>
    <div class="tbl wide"><table>${head(['Lattice', 'atoms per cell n', 'spacing a (LJ units, density ρ)', 'nearest neighbour'])}
${row(['sc (simple cubic)', '1', '(1/ρ)<sup>1/3</sup>', 'a'])}
${row(['bcc', '2', '(2/ρ)<sup>1/3</sup>', 'a√3/2'])}
${row(['fcc', '4', '(4/ρ)<sup>1/3</sup>', 'a/√2'])}
${row(['diamond', '8', '(8/ρ)<sup>1/3</sup>', 'a√3/4'])}
${row(['sq (2D square)', '1', '(1/ρ)<sup>1/2</sup>', 'a'])}
${row(['hex (2D)', '2 in an a × √3a cell', '(2/(ρ√3))<sup>1/2</sup>', 'a'])}</table></div>
    <div class="calc">box length  = number of cells × spacing
atoms       = n × cells_x × cells_y × cells_z      (whole numbers of cells)
fcc, ρ = 0.8442, 3 cells:  a = (4/0.8442)^(1/3) = 1.6795962,  L = 5.0387886,  N = 4 × 27 = 108</div>
  </section>

  <section class="block">
    <h2>3 · Modify, predict, observe: the lattice builder</h2>
    <p class="muted">Pick a lattice, a density and the number of cells. The page predicts the three log lines; LAMMPS then builds the system and the numbers are compared.</p>
    <div class="lmp" id="builder">
      <div class="row">
        <label class="rlab">lattice <select class="r-in b-style">${Object.keys(LAT).map(k => `<option${k === 'fcc' ? ' selected' : ''}>${k}</option>`).join('')}</select></label>
        <label class="rlab">density <input class="r-in b-rho" value="0.8442"></label>
        <label class="rlab">cells <input class="r-in b-cells" value="3"></label>
      </div>
      <div class="row quick"><span class="muted small">Try:</span>${[['fcc', '0.8442', '10'], ['bcc', '0.8442', '3'], ['diamond', '0.8442', '3'], ['hex', '0.8442', '3'], ['fcc', '0.8442', '2.5']].map(([s, r, c]) => `<button class="btn sm ghost" data-try="${s},${r},${c}">${s} ${c}</button>`).join('')}</div>
      <p><b>Predict by hand</b></p>
      <div class="calc b-pred"></div>
      <pre class="lmp-mini b-script"></pre>
      <div class="row"><button class="btn sm primary" data-k="build">▶ Build it in LAMMPS</button><span class="lmp-status muted b-st"></span></div>
      <div class="lmp-compare b-cmp"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>What the builder shows</h2>
    <ul>
      <li>For whole numbers of cells, three lines of arithmetic predict exactly what LAMMPS builds.</li>
      <li>With 2.5 cells LAMMPS still makes a box, and fills it with only the 63 lattice points that fit inside. Nothing warns you: across the periodic boundary the next copy starts half a cell too early, so atoms meet at the wrong distance and the crystal is broken. Choose whole numbers of cells for crystals.</li>
      <li>In 2D the hex lattice has a rectangular cell, √3 times longer in y: its spacing line shows two different numbers.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> How many atoms does <code>lattice bcc 0.8442</code> with <code>region box block 0 5 0 5 0 5</code> give, and what is the box length?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>N = 2 × 5³ = <b>250</b>. a = (2/0.8442)<sup>1/3</sup> = 1.3330964, so L = 5a = <b>6.6654819</b>. Check it in the builder.</div></div>
    <div class="ex"><p><b>2.</b> Count the 63 atoms of the 2.5-cell fcc box by hand.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>FCC points have coordinates that are multiples of ½ (in cells) with an even sum of the doubled coordinates. Inside [0, 2.5) each coordinate can be 0, ½, 1, 1½ or 2: 5 × 5 × 5 = 125 combinations, of which those with an even sum of (2x + 2y + 2z) count. That is (125 + 1)/2 = <b>63</b>.</div></div>
    <div class="ex"><p><b>3.</b> In metal units you write <code>lattice fcc 0.8442</code>. What goes wrong?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Outside LJ units the number is the lattice constant: you get a cube edge of 0.8442 Å, so neighbours 0.6 Å apart. Atoms that close have enormous repulsion and the run blows up. Use the real lattice constant (for example 5.26 Å for solid argon).</div></div>
  </section>`,

  mount(root) {
    const B = root.querySelector('#builder'), sel = B.querySelector('.b-style'), rin = B.querySelector('.b-rho'), cin = B.querySelector('.b-cells');
    const read = () => { const style = sel.value, rho = parseFloat(rin.value), c = parseFloat(cin.value), dim = LAT[style].dim; return { style, rho, c, dim, cells: [c, c, c] }; };
    function draw() {
      const { style, rho, c, dim, cells } = read();
      if (!(rho > 0) || !(c > 0) || c > 12) { B.querySelector('.b-pred').textContent = 'Use a density above 0 and between 0 and 12 cells.'; return null; }
      const p = predict(style, rho, cells), L = LAT[style];
      const whole = Number.isInteger(c);
      B.querySelector('.b-pred').innerHTML = esc(`n = ${L.n} atom${L.n > 1 ? 's' : ''} per cell${style === 'hex' ? ' (cell a × √3a)' : ''}
a = ${dim === 3 ? `(${L.n}/${rho})^(1/3)` : style === 'hex' ? `(2/(${rho}·√3))^(1/2)` : `(1/${rho})^(1/2)`} = ${g8(p.a)}
spacing x,y${dim === 3 ? ',z' : ''}: ${p.sp.slice(0, dim).map(g8).join(' ')}
box = (${dim === 3 ? '0 0 0' : '0 0'}) to (${p.box.slice(0, dim).map(g8).join(' ')})
atoms = ${whole ? `${L.n} × ${c}^${dim} = ` : 'lattice points inside the box: '}${p.N}
nearest neighbour = ${L.nnTxt} = ${g8(p.nn)}`);
      B.querySelector('.b-script').textContent = script(style, rho, cells, dim);
      B.querySelector('.b-cmp').innerHTML = '';
      return p;
    }
    [sel, rin, cin].forEach(x => x.oninput = draw);
    B.querySelectorAll('[data-try]').forEach(b => b.onclick = () => { const [s, r, c] = b.dataset.try.split(','); sel.value = s; rin.value = r; cin.value = c; draw(); B.querySelector('[data-k=build]').click(); });
    draw();
    B.querySelector('[data-k=build]').onclick = async e => {
      const p = draw(); if (!p) return; const { style, rho, cells, dim } = read();
      const btn = e.currentTarget, st = B.querySelector('.b-st'); btn.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(script(style, rho, cells, dim)); btn.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const sp = (res.lines.find(l => l.includes('Lattice spacing')) || '').split('=')[1]?.trim().split(/\s+/) || [];
      const bx = (res.lines.find(l => l.includes('Created orthogonal box')) || '').match(/to \(([^)]+)\)/)?.[1].split(/\s+/) || [];
      const n = (res.lines.find(l => /Created \d+ atoms/.test(l)) || '').match(/Created (\d+) atoms/)?.[1] || '?';
      const lines = [];
      for (let k = 0; k < dim; k++) lines.push([`spacing ${'xyz'[k]}`, p.sp[k], sp[k]]);
      for (let k = 0; k < dim; k++) lines.push([`box ${'xyz'[k]}`, p.box[k], bx[k]]);
      const trs = lines.map(([lab, h, l]) => { const m = matchMarkup(g8(h), l || ''); return `<tr><td>${lab}</td><td class="mono">${m.html}</td><td class="mono">${l}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')
        + `<tr><td>atoms</td><td class="mono">${p.N}</td><td class="mono">${n}</td><td>${String(p.N) === n ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`;
      B.querySelector('.b-cmp').innerHTML = `<div class="tbl"><table><tr><th></th><th>by hand</th><th>LAMMPS</th><th></th></tr>${trs}</table></div>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
