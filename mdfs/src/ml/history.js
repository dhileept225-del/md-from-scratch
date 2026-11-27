// ML potentials from scratch · M0.2 Who started it: 1995 to today.
// Every date was checked (Nov 2026) against the arXiv listing (v1 submission date) and the journal DOI record.
const FAM = {
  nn: { name: 'Neural network on a descriptor', col: '#38E1FF' },
  kernel: { name: 'Kernel / Gaussian process', col: '#B79CFF' },
  linear: { name: 'Linear basis', col: '#7EE787' },
  mp: { name: 'Message passing', col: '#FFC24B' },
  uni: { name: 'Universal (foundation) model', col: '#FF6B6B' },
  data: { name: 'Dataset or workflow', col: '#8FA3C7' },
};
// first: first public appearance (arXiv v1, or the journal if there was no preprint). pub: journal or conference year.
const E = [
  { n: 'UMA', who: 'Meta FAIR', fam: 'uni', first: '2025-06-30', ax: '2506.23971', pub: null, add: 'One model for molecules, materials, catalysts and more, trained on OMol25, OMat24 and other sets together.' },
  { n: 'OMol25', who: 'Meta FAIR et al.', fam: 'data', first: '2025-05-13', ax: '2505.08762', pub: null, add: 'A very large molecular DFT dataset (ωB97M-V), with released models trained on it, including eSEN.' },
  { n: 'PET-MAD', who: 'Mazitov, Ceriotti et al.', fam: 'uni', first: '2025-03-18', ax: '2503.14118', pub: 2025, venue: 'Nat. Commun. 16, 10653', doi: '10.1038/s41467-025-65662-7', add: 'A light universal model whose rotational symmetry is learned, not built in.' },
  { n: 'eSEN', who: 'Fu, Zitnick et al. (Meta FAIR)', fam: 'uni', first: '2025-02-17', ax: '2502.12147', pub: 2025, venue: 'ICML 2025', add: 'Test a model by running MD and checking energy conservation; a smooth model built to pass it.' },
  { n: 'MACE-Osaka24', who: 'Shiota, Mizukami et al.', fam: 'uni', first: '2024-12-17', ax: '2412.13088', pub: null, add: 'One dataset and one MACE model for both molecules and crystals.' },
  { n: 'Orb', who: 'Orbital Materials', fam: 'uni', first: '2024-10-29', ax: '2410.22570', pub: null, add: 'A fast universal materials model; forces predicted directly, not as the slope of the energy.' },
  { n: 'OMat24', who: 'Barroso-Luque, Ulissi et al. (Meta FAIR)', fam: 'data', first: '2024-10-16', ax: '2410.12771', pub: null, add: 'Over 110 million DFT calculations of inorganic materials, far from equilibrium, released openly.' },
  { n: 'MatterSim', who: 'Microsoft Research', fam: 'uni', first: '2024-05-08', ax: '2405.04967', pub: null, add: 'A universal materials model trained across elements, temperatures and pressures.' },
  { n: 'SevenNet', who: 'Park, Han et al.', fam: 'uni', first: '2024-02-06', ax: '2402.03789', pub: 2024, venue: 'J. Chem. Theory Comput. 20, 4857', doi: '10.1021/acs.jctc.4c00190', add: 'A way to run NequIP-type message passing in parallel in LAMMPS, plus a pretrained universal model.' },
  { n: 'MACE-MP-0', who: 'Batatia, Csányi et al.', fam: 'uni', first: '2023-12-29', ax: '2401.00096', pub: null, add: 'The first widely used foundation model for materials: one MACE trained on the Materials Project.' },
  { n: 'MACE-OFF23', who: 'Kovács, Csányi et al.', fam: 'uni', first: '2023-12-23', ax: '2312.15211', pub: 2025, venue: 'J. Am. Chem. Soc.', doi: '10.1021/jacs.4c07099', add: 'Transferable MACE force fields for organic molecules, short-ranged.' },
  { n: 'GRACE', who: 'Bochkarev, Lysogorskiy, Drautz', fam: 'mp', first: '2023-11', ax: '2311.16326', pub: 2024, venue: 'Phys. Rev. X 14, 021036', doi: '10.1103/PhysRevX.14.021036', add: 'Graph atomic cluster expansion: ACE extended to graphs, which contains message passing as a special case.' },
  { n: 'CHGNet', who: 'Deng, Ceder et al.', fam: 'uni', first: '2023-02-28', ax: '2302.14231', pub: 2023, venue: 'Nat. Mach. Intell. 5, 1031', doi: '10.1038/s42256-023-00716-3', add: 'A universal materials model that also predicts magnetic moments, as a proxy for charge state.' },
  { n: 'MACE', who: 'Batatia, Kovács, Simm, Ortner, Csányi', fam: 'mp', first: '2022-06-15', ax: '2206.07697', pub: 2022, venue: 'NeurIPS 2022', add: 'Many-body equivariant messages: two message-passing layers are enough.' },
  { n: 'Allegro', who: 'Musaelian, Batzner, Kozinsky et al.', fam: 'mp', first: '2022-04-11', ax: '2204.05249', pub: 2023, venue: 'Nat. Commun. 14, 579', doi: '10.1038/s41467-023-36329-y', add: 'Equivariant but strictly local, no message passing: scales to 100 million atoms.' },
  { n: 'M3GNet', who: 'Chen, Ong', fam: 'uni', first: '2022-02-05', ax: '2202.02450', pub: 2022, venue: 'Nat. Comput. Sci. 2, 718', doi: '10.1038/s43588-022-00349-3', add: 'A universal potential for the periodic table, with three-body interactions, trained on Materials Project relaxations.' },
  { n: 'NequIP', who: 'Batzner, Kozinsky et al.', fam: 'mp', first: '2021-01-08', ax: '2101.03164', pub: 2022, venue: 'Nat. Commun. 13, 2453', doi: '10.1038/s41467-022-29939-5', add: 'E(3)-equivariant messages: accurate from a few hundred training structures.' },
  { n: '4G-HDNNP', who: 'Ko, Finkler, Goedecker, Behler', fam: 'nn', first: '2020-09', ax: '2009.06484', pub: 2021, venue: 'Nat. Commun. 12, 398', doi: '10.1038/s41467-020-20427-2', add: 'Charges that can move across the whole system: long-range electrostatics and charge transfer.' },
  { n: 'DP-GEN', who: 'Zhang, Wang, E et al.', fam: 'data', first: '2019-10', ax: '1910.12690', pub: 2020, venue: 'Comput. Phys. Commun. 253, 107206', doi: '10.1016/j.cpc.2020.107206', add: 'Active learning as software: explore with MD, pick uncertain structures, label with DFT, retrain.' },
  { n: 'ACE', who: 'Drautz', fam: 'linear', first: '2019-01', pub: 2019, venue: 'Phys. Rev. B 99, 014104', doi: '10.1103/PhysRevB.99.014104', add: 'Atomic cluster expansion: a complete basis, systematic in body order, linear cost in neighbours.' },
  { n: 'DeePMD', who: 'Zhang, Han, Wang, Car, E', fam: 'nn', first: '2017-07-30', ax: '1707.09571', pub: 2018, venue: 'Phys. Rev. Lett. 120, 143001', doi: '10.1103/PhysRevLett.120.143001', add: 'Deep Potential: descriptors and network trained together; the DeePMD-kit software.' },
  { n: 'SchNet', who: 'Schütt, Müller et al.', fam: 'mp', first: '2017-06-26', ax: '1706.08566', pub: 2017, venue: 'NeurIPS 2017', add: 'Continuous-filter convolutions: atoms pass messages, no hand-made descriptor.' },
  { n: 'ANI-1', who: 'Smith, Isayev, Roitberg', fam: 'nn', first: '2016-10-27', ax: '1610.08935', pub: 2017, venue: 'Chem. Sci. 8, 3192', doi: '10.1039/C6SC05720A', add: 'A transferable network potential for organic molecules (H, C, N, O), trained on millions of DFT points.' },
  { n: 'MTP', who: 'Shapeev', fam: 'linear', first: '2015-12-18', ax: '1512.06054', pub: 2016, venue: 'Multiscale Model. Simul. 14, 1153', doi: '10.1137/15M1054183', add: 'Moment tensor potentials: a linear model on moment tensors, improvable step by step.' },
  { n: 'SNAP', who: 'Thompson et al. (Sandia)', fam: 'linear', first: '2014-09', ax: '1409.3880', pub: 2015, venue: 'J. Comput. Phys. 285, 316', doi: '10.1016/j.jcp.2014.12.018', add: 'A linear least-squares fit on bispectrum descriptors, built into LAMMPS.' },
  { n: 'SOAP', who: 'Bartók, Kondor, Csányi', fam: 'kernel', first: '2012-09-14', ax: '1209.3140', pub: 2013, venue: 'Phys. Rev. B 87, 184115', doi: '10.1103/PhysRevB.87.184115', add: 'Smooth overlap of atomic positions: a standard measure of how alike two environments are.' },
  { n: 'GAP', who: 'Bartók, Payne, Kondor, Csányi', fam: 'kernel', first: '2009-10-06', ax: '0910.1019', pub: 2010, venue: 'Phys. Rev. Lett. 104, 136403', doi: '10.1103/PhysRevLett.104.136403', add: 'Gaussian approximation potentials: kernel regression of atomic energies, with error bars.' },
  { n: 'Behler–Parrinello', who: 'Behler, Parrinello', fam: 'nn', first: '2007-04-02', pub: 2007, venue: 'Phys. Rev. Lett. 98, 146401', doi: '10.1103/PhysRevLett.98.146401', add: 'Energy as a sum of atomic energies, one network per element, symmetry-function inputs: any number of atoms. Bulk silicon.' },
  { n: 'Blank et al.', who: 'Blank, Brown, Calhoun, Doren', fam: 'nn', first: '1995-09-08', pub: 1995, venue: 'J. Chem. Phys. 103, 4129', doi: '10.1063/1.469597', add: 'The first neural-network fit of a potential energy surface from electronic-structure data: H₂ on a silicon cluster, a few degrees of freedom.' },
];
const MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const fmtD = s => { const [y, m, d] = s.split('-'); return `${d ? +d + ' ' : ''}${MON[+m - 1]} ${y}`; };
const yr = (e, mode) => mode === 'pub' && e.pub ? e.pub : +e.first.slice(0, 4);

export default {
  id: 'history',
  minutes: 10,
  summary: 'Machine-learned potentials are older than they look. The first neural-network fit of a potential energy surface appeared in 1995; the idea that made them work for any number of atoms came in 2007; the foundation models of 2023–2025 are the latest step. This page puts the main papers on one timeline, newest first, with the date each first appeared and where it was published, every date checked against its arXiv record and DOI.',
  learn: [
    'The four steps: small fits (1995), atomic energies (2007), message passing (2017), foundation models (2023)',
    'Which family each method belongs to: networks, kernels, linear bases, message passing',
    'What each step added that the one before lacked',
    'Where the methods taught later in this course come from',
  ],
  simHint: 'Tap a family to filter. Switch between the date a paper first appeared (usually arXiv) and its journal year: watch which entries move. Tap an entry for what it added and its links.',

  idea: `
  <p><b>1995: a fit, not a potential.</b> Blank, Brown, Calhoun and Doren fitted a small neural network to the energy of H₂ near a silicon cluster. It worked, but the network took a fixed list of coordinates as input, so it was tied to that one system and its few degrees of freedom.</p>
  <p><b>2007: atomic energies.</b> Behler and Parrinello wrote the total energy as a sum of atomic energies, each from a network that sees only the atom's neighbours within a cutoff, described by <i>symmetry functions</i> that do not change when the system is rotated or atoms are swapped. One network per element, any number of atoms. Almost every model since keeps this sum.</p>
  <p><b>2009–2019: better descriptors, other fits.</b> GAP replaced the network by a Gaussian process (a kernel fit, with error bars); SOAP gave a smooth way to compare environments; SNAP, MTP and ACE showed that a <i>linear</i> fit on a rich enough basis is as accurate and very fast. ANI-1 and DeePMD scaled the network route to large datasets.</p>
  <p><b>2017–2022: message passing.</b> SchNet let atoms exchange learned messages instead of using a fixed descriptor. NequIP made the messages equivariant (they rotate with the molecule), which cut the data needed by orders of magnitude; MACE made them many-body, so two layers suffice; Allegro kept equivariance but dropped messages to stay strictly local and parallel.</p>
  <p><b>2022–2025: one model for everything.</b> Trained on large public DFT datasets (Materials Project, then OMat24 and OMol25), M3GNet, CHGNet, MACE-MP-0 and their successors run on almost any composition without new training. They are the starting point you fine-tune in part M5, and part M6 compares them.</p>`,

  footnote: `<p><b>About the dates.</b> The year shown first is the paper's first public appearance: the arXiv v1 date, or the journal date when there was no preprint (Blank 1995, Behler–Parrinello 2007, ACE 2019). Conference papers (SchNet, MACE, eSEN) carry the conference year. Papers marked "preprint" had no journal version when this page was written.</p>
  <p>Journal years are often later: ANI-1 by 104 days, DeePMD by 248, NequIP by 481 (8 Jan 2021 → 4 May 2022). The plan for this course listed journal years; the toggle above shows both. Every date was checked against its arXiv record and DOI in November 2026.</p>`,

  quiz: [
    { q: 'What did Behler and Parrinello (2007) add that let a neural network handle any number of atoms?', options: ['A larger network', 'The energy as a sum of atomic energies, each from the atom\'s local neighbourhood', 'Training on forces', 'Periodic boundary conditions'], answer: 1, why: 'Each atom\'s energy depends on neighbours within a cutoff, through symmetry functions; summing them works for 10 atoms or 10 million.' },
    { q: 'Which appeared first as a preprint?', options: ['MACE', 'NequIP', 'Allegro', 'They appeared the same month'], answer: 1, why: 'NequIP: 8 Jan 2021. Allegro: 11 Apr 2022. MACE: 15 Jun 2022.' },
    { q: 'SNAP and GAP both use bispectrum-type descriptors. What is the main difference?', options: ['SNAP is a linear least-squares fit; GAP is a Gaussian process (kernel) fit', 'SNAP is a neural network', 'GAP only works for molecules', 'SNAP has error bars, GAP does not'], answer: 0, why: 'SNAP assumes the atomic energy is linear in the descriptors, which makes it fast and simple to fit; GAP compares environments through a kernel and gives an uncertainty.' },
    { q: 'What did NequIP (2021) change in message passing?', options: ['It removed the cutoff', 'Its messages are equivariant: they rotate with the structure', 'It used a Gaussian process instead of a network', 'It trained on energies only'], answer: 1, why: 'Equivariant messages carry directions, not just distances, which made it accurate from a few hundred training structures.' },
  ],

  mount(root) {
    const el = root.querySelector('#sim');
    el.innerHTML = `
    <div class="hist-view"><canvas></canvas></div>
    <div class="hist-ctrl">
      <div class="seg mini"><button data-m="first" class="on">First appeared</button><button data-m="pub">Journal year</button></div>
      <div class="chips">${Object.entries(FAM).map(([k, f]) => `<button class="chip on" data-f="${k}"><span style="color:${f.col}">●</span> ${f.name}</button>`).join('')}</div>
      <p class="muted small" data-o="count"></p>
    </div>
    <div class="tbl wide hist-tbl"><table><thead><tr><th>Year</th><th>Method</th><th>What it added</th></tr></thead><tbody></tbody></table></div>`;
    const cv = el.querySelector('canvas'), cx = cv.getContext('2d'), view = el.querySelector('.hist-view'), body = el.querySelector('tbody'), count = el.querySelector('[data-o=count]');
    const st = { mode: 'first', on: new Set(Object.keys(FAM)), sel: null };
    let W = 600, H = 220, hits = [];
    const shown = () => E.filter(e => st.on.has(e.fam)).sort((a, b) => yr(b, st.mode) - yr(a, st.mode) || b.first.localeCompare(a.first));
    const link = e => [e.ax && `<a href="https://arxiv.org/abs/${e.ax}" target="_blank" rel="noopener">arXiv ${e.ax}</a>`, e.doi && `<a href="https://doi.org/${e.doi}" target="_blank" rel="noopener">${e.venue}</a>`, !e.doi && e.venue && e.venue].filter(Boolean).join(' · ');
    function list() {
      const L = shown(), moved = L.filter(e => e.pub && e.pub !== +e.first.slice(0, 4)).length;
      count.textContent = st.mode === 'pub'
        ? `${L.length} shown. ${moved} of them carry a later year when dated by the journal; preprints with no journal version keep their arXiv year.`
        : `${L.length} shown, newest first, dated by first public appearance.`;
      body.innerHTML = L.map(e => {
        const y = yr(e, st.mode), f = FAM[e.fam], shift = st.mode === 'pub' && e.pub && e.pub !== +e.first.slice(0, 4);
        const open = st.sel === e.n;
        return `<tr data-e="${e.n}" style="cursor:pointer${open ? ';background:rgba(56,225,255,.06)' : ''}"><td><b${shift ? ' style="color:#FFC24B"' : ''}>${y}</b>${shift ? `<br><small class="muted">was ${e.first.slice(0, 4)}</small>` : ''}</td>
          <td><span style="color:${f.col}">●</span> <b>${e.n}</b><br><small class="muted">${e.who}</small></td>
          <td>${e.add}${open ? `<br><small class="muted">First appeared ${fmtD(e.first)}${e.pub ? '' : ' · preprint'}${e.ax || e.venue ? ' · ' + link(e) : ''}</small>` : ''}</td></tr>`;
      }).join('');
      body.querySelectorAll('[data-e]').forEach(r => r.onclick = () => { st.sel = st.sel === r.dataset.e ? null : r.dataset.e; list(); draw(); });
    }
    function draw() {
      cx.clearRect(0, 0, W, H); hits = [];
      const y0 = 1995, y1 = 2025, pad = 22, x = y => pad + (y - y0) / (y1 - y0) * (W - 2 * pad), base = H - 34;
      cx.strokeStyle = 'rgba(143,163,199,.35)'; cx.fillStyle = '#8FA3C7'; cx.font = '12px system-ui'; cx.lineWidth = 1;
      cx.beginPath(); cx.moveTo(pad, base); cx.lineTo(W - pad, base); cx.stroke();
      for (let y = y0; y <= y1; y += W < 480 ? 10 : 5) { const X = x(y); cx.beginPath(); cx.moveTo(X, base); cx.lineTo(X, base + 5); cx.stroke(); cx.fillText(y, X - 14, base + 19); }
      const stack = {}, r = Math.max(4, Math.min(7, (W - 2 * pad) / 31 / 2.4));
      [...shown()].reverse().forEach(e => {
        const y = yr(e, st.mode), k = stack[y] = (stack[y] || 0) + 1, X = x(y), Y = base - 6 - (k - 0.5) * (2 * r + 3);
        cx.fillStyle = FAM[e.fam].col; cx.globalAlpha = st.sel && st.sel !== e.n ? 0.35 : 1;
        cx.beginPath(); cx.arc(X, Y, r, 0, 7); cx.fill(); cx.globalAlpha = 1;
        if (st.sel === e.n) { cx.strokeStyle = '#F2F5FA'; cx.lineWidth = 2; cx.stroke(); cx.fillStyle = '#F2F5FA'; cx.font = '600 12px system-ui'; const t = `${e.n} ${y}`, w = cx.measureText(t).width; cx.fillText(t, Math.min(Math.max(X - w / 2, 2), W - w - 2), Y - r - 6); cx.font = '12px system-ui'; }
        hits.push([X, Y, r + 4, e.n]);
      });
      cx.fillStyle = '#8FA3C7'; cx.fillText(st.mode === 'pub' ? 'dated by journal year' : 'dated by first appearance', pad, 14);
    }
    const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
    const ro = new ResizeObserver(resize); ro.observe(view);
    cv.onclick = ev => { const b = cv.getBoundingClientRect(), px = ev.clientX - b.left, py = ev.clientY - b.top; const h = hits.find(([X, Y, R]) => (px - X) ** 2 + (py - Y) ** 2 < R * R); st.sel = h ? (st.sel === h[3] ? null : h[3]) : null; list(); draw(); };
    el.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.mode = b.dataset.m; el.querySelectorAll('[data-m]').forEach(x => x.classList.toggle('on', x === b)); list(); draw(); });
    el.querySelectorAll('[data-f]').forEach(b => b.onclick = () => { const k = b.dataset.f; st.on.has(k) ? st.on.delete(k) : st.on.add(k); if (!st.on.size) Object.keys(FAM).forEach(f => st.on.add(f)); el.querySelectorAll('[data-f]').forEach(x => x.classList.toggle('on', st.on.has(x.dataset.f))); list(); draw(); });
    list(); resize();

    return () => ro.disconnect();
  },
};
