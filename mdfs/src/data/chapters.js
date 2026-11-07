// The whole course map. Set `ready: true` when a lesson file exists in src/lessons/.
// level: 'core' = undergraduate path, 'adv' = graduate / research path.
export const chapters = [
  {
    n: 1, title: 'Basics', blurb: 'What MD is, and the smallest loop that makes atoms move.',
    lessons: [
      { id: 'what-is-md', title: 'What is molecular dynamics?', level: 'core', ready: true },
      { id: 'lennard-jones', title: 'Atoms and forces: the Lennard-Jones potential', level: 'core', ready: true },
      { id: 'newton-in-code', title: "Newton's law in code", level: 'core', ready: true },
      { id: 'time-step', title: 'The time step', level: 'core', ready: true },
      { id: 'verlet', title: 'Euler vs velocity Verlet', level: 'core' },
      { id: 'energy-check', title: 'Energy conservation: checking your simulation', level: 'core' },
    ],
  },
  {
    n: 2, title: 'Building a system', blurb: 'Units, starting structures and the periodic box.',
    lessons: [
      { id: 'units', title: 'Units in MD', level: 'core' },
      { id: 'lattice-start', title: 'Starting positions: lattices', level: 'core' },
      { id: 'maxwell-boltzmann', title: 'Starting velocities: Maxwell–Boltzmann', level: 'core' },
      { id: 'pbc', title: 'Periodic boundary conditions', level: 'core' },
      { id: 'cutoffs', title: 'Cutoffs and the minimum image', level: 'core' },
      { id: 'neighbour-lists', title: 'Neighbour lists', level: 'adv' },
    ],
  },
  {
    n: 3, title: 'Ensembles and temperature', blurb: 'Controlling temperature and pressure.',
    lessons: [
      { id: 'temperature', title: 'What temperature means in MD', level: 'core' },
      { id: 'ensembles', title: 'NVE, NVT and NPT', level: 'core' },
      { id: 'berendsen', title: 'The Berendsen thermostat', level: 'core' },
      { id: 'nose-hoover', title: 'Nosé–Hoover thermostats', level: 'adv' },
      { id: 'langevin', title: 'Langevin dynamics', level: 'adv' },
      { id: 'barostats', title: 'Pressure and barostats', level: 'adv' },
    ],
  },
  {
    n: 4, title: 'Molecules', blurb: 'Bonds, angles and force fields.',
    lessons: [
      { id: 'bonds-angles', title: 'Bonds and angles', level: 'core' },
      { id: 'dihedrals', title: 'Dihedrals', level: 'core' },
      { id: 'force-fields', title: 'What is a force field?', level: 'core' },
      { id: 'shake', title: 'Constraints and SHAKE', level: 'adv' },
      { id: 'water-models', title: 'Water models', level: 'core' },
    ],
  },
  {
    n: 5, title: 'Electrostatics', blurb: 'Long-range charges in a periodic world.',
    lessons: [
      { id: 'coulomb-pbc', title: "Coulomb's law in a periodic box", level: 'core' },
      { id: 'ewald', title: 'Ewald summation', level: 'adv' },
      { id: 'pme', title: 'Particle-mesh Ewald', level: 'adv' },
      { id: 'polarizable', title: 'Polarizable models', level: 'adv' },
    ],
  },
  {
    n: 6, title: 'Analysis', blurb: 'Turning trajectories into numbers you can compare with experiment.',
    lessons: [
      { id: 'rdf', title: 'Radial distribution function', level: 'core' },
      { id: 'msd', title: 'Mean squared displacement and diffusion', level: 'core' },
      { id: 'vacf', title: 'Velocity autocorrelation', level: 'adv' },
      { id: 'van-hove', title: 'Van Hove functions', level: 'adv' },
      { id: 'conductivity', title: 'Ionic conductivity', level: 'adv' },
      { id: 'onsager', title: 'Onsager coefficients and transference numbers', level: 'adv' },
    ],
  },
  {
    n: 7, title: 'Beyond', blurb: 'From your own code to research tools.',
    lessons: [
      { id: 'real-codes', title: 'From toy code to LAMMPS and GROMACS', level: 'core' },
      { id: 'ml-potentials', title: 'Machine-learning potentials', level: 'adv' },
      { id: 'metadynamics', title: 'Enhanced sampling: metadynamics', level: 'adv' },
    ],
  },
];

export const allLessons = chapters.flatMap(c => c.lessons.map((l, i) => ({ ...l, chapter: c.n, num: `${c.n}.${i + 1}` })));
