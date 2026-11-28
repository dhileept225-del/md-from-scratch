// ML potentials from scratch: 36 lessons in seven parts plus six labs (plan of 16 Nov 2026).
// Set `ready: true` when a page exists in src/ml/<id>.js and is listed in src/ml/pages.js.
export const mlParts = [
  {
    key: 'M0', title: 'Why learn a potential', blurb: 'What machine-learned potentials are for, where they came from, and the rules any potential must obey.',
    pages: [
      { id: 'cost-ladder', title: 'The cost ladder: force field, DFT, ML', ready: true },
      { id: 'history', title: 'Who started it: 1995 to today', ready: true },
      { id: 'invariances', title: 'What a potential must obey', ready: true },
    ],
  },
  {
    key: 'M1', title: 'Learning, by hand', blurb: 'Least squares, a tiny neural network, backpropagation and forces, all on the argon dimer.',
    pages: [
      { id: 'least-squares', title: 'Fitting is learning: least squares' },
      { id: 'overfitting', title: 'Overfitting, training and test sets' },
      { id: 'neuron', title: 'A neuron and a tiny network' },
      { id: 'backprop', title: 'Gradient descent and backpropagation' },
      { id: 'forces-slope', title: 'Forces are the slope of the energy' },
      { id: 'energy-force-training', title: 'Training on energies and forces' },
      { id: 'lab-train-by-hand', title: 'Lab: train a network by hand', lab: 'ML-L1' },
    ],
  },
  {
    key: 'M2', title: "Describing an atom's surroundings", blurb: 'Cutoffs, symmetry, and the descriptors that turn neighbours into numbers.',
    pages: [
      { id: 'cutoff-function', title: 'Locality and the cutoff function' },
      { id: 'raw-coordinates', title: 'Why raw coordinates fail' },
      { id: 'symmetry-functions', title: 'Behler–Parrinello symmetry functions G2 and G4' },
      { id: 'soap-snap', title: 'SOAP and the bispectrum (SNAP)' },
      { id: 'body-order-ace', title: 'Body order and ACE' },
      { id: 'lab-fingerprint-water', title: 'Lab: fingerprint a water molecule', lab: 'ML-L2' },
    ],
  },
  {
    key: 'M3', title: 'Classical ML potentials', blurb: 'Behler–Parrinello networks, Gaussian processes, linear models and Deep Potential.',
    pages: [
      { id: 'behler-parrinello', title: 'The Behler–Parrinello network: E = Σ Eᵢ' },
      { id: 'gaussian-process', title: 'Kernels and Gaussian processes (GAP)' },
      { id: 'linear-models', title: 'Linear models: SNAP, MTP, ACE' },
      { id: 'deep-potential', title: 'Deep Potential (DeePMD)' },
      { id: 'lab-bp-argon', title: 'Lab: a Behler–Parrinello potential for argon', lab: 'ML-L3' },
      { id: 'lab-gp-error-bars', title: 'Lab: a Gaussian-process potential with error bars', lab: 'ML-L4' },
    ],
  },
  {
    key: 'M4', title: 'Message passing and equivariant models', blurb: 'Graphs, vectors that rotate, spherical harmonics, NequIP, Allegro, MACE and long range.',
    pages: [
      { id: 'message-passing', title: 'Atoms as a graph; message passing (SchNet)' },
      { id: 'equivariant-features', title: 'Invariant vs equivariant features' },
      { id: 'spherical-harmonics', title: 'Spherical harmonics and tensor products' },
      { id: 'nequip-allegro', title: 'NequIP and Allegro' },
      { id: 'mace', title: 'MACE: many-body messages' },
      { id: 'long-range', title: 'Long range: charges and dispersion' },
      { id: 'lab-equivariance', title: 'Lab: equivariance, checked', lab: 'ML-L5' },
    ],
  },
  {
    key: 'M5', title: 'Data, training and trust', blurb: 'Reference data, validation beyond parity plots, active learning, fine-tuning and MD tests.',
    pages: [
      { id: 'dft-data', title: 'Where data comes from: DFT and its settings' },
      { id: 'validation', title: 'Validation: parity plots are not enough' },
      { id: 'active-learning', title: 'Active learning and committees (DP-GEN)' },
      { id: 'fine-tuning', title: 'Fine-tuning a foundation model' },
      { id: 'md-testing', title: 'Testing in MD: stability, RDF, diffusion' },
      { id: 'lab-when-it-breaks', title: 'Lab: when an ML potential breaks', lab: 'ML-L6' },
    ],
  },
  {
    key: 'M6', title: 'Foundation models and using them', blurb: 'The pretrained models, their data, how to choose one and how to run it in MD.',
    pages: [
      { id: 'foundation-datasets', title: 'Foundation models and their datasets' },
      { id: 'materials-models', title: 'Materials models' },
      { id: 'molecular-models', title: 'Molecular models' },
      { id: 'universal-models', title: 'Both: UMA and MACE-Osaka24' },
      { id: 'conservative-forces', title: 'Conservative vs direct forces' },
      { id: 'choosing-models', title: 'Choosing a model, benchmarks and licences' },
      { id: 'speed-md-codes', title: 'Speed and running in MD codes' },
    ],
  },
];

// lessons are numbered M1.1, M1.2 … within a part; labs carry their own label (ML-L1 …)
export const allMl = mlParts.flatMap(p => { let n = 0; return p.pages.map(pg => ({ ...pg, part: p.key, num: pg.lab || `${p.key}.${++n}` })); });
