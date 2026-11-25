// Links between the theory lessons and the LAMMPS pages that show the same idea inside LAMMPS.
// Used by the app and by the prerender script, so search engines see the links too.
const MAP = {
  'lennard-jones': ['pair-lj-cut', 'styles', 'forces-in-script'],
  'time-step': ['log-and-errors', 'other-loops'],
  'energy-check': ['log-and-errors', 'energy-bookkeeping'],
  units: ['units', 'units-constants'],
  'lattice-start': ['building-the-system'],
  'maxwell-boltzmann': ['velocities'],
  temperature: ['velocities', 'finding-your-way'],
  verlet: ['verlet-loop', 'two-atom-step'],
  cutoffs: ['pair-lj-cut', 'forces-in-script', 'energy-bookkeeping'],
  'nose-hoover': ['energy-bookkeeping'],
  berendsen: ['hooks', 'energy-bookkeeping'],
  langevin: ['hooks', 'energy-bookkeeping'],
  'neighbour-lists': ['neighbor-lists', 'forces-in-script', 'multiple-runs', 'timing-table'],
  'bonds-angles': ['read-data', 'style-families'],
  dihedrals: ['style-families'],
  'force-fields': ['style-families', 'styles'],
  'water-models': ['read-data'],
  'real-codes': ['lammps-start', 'source-tree'],
  msd: ['atom-storage', 'domain-pbc'],
  pbc: ['atom-storage', 'domain-pbc'],
  'lab-molecules': ['read-data', 'style-families'],
  'lab-building': ['building-the-system'],
  'lab-basics': ['styles'],
};
export const lammpsFor = lessonId => MAP[lessonId] || [];
export const lessonsFor = lammpsId => Object.keys(MAP).filter(k => MAP[k].includes(lammpsId));

// html for the box; `find` returns {id, num, title, ready} for an id
export function lessonBox(lessonId, findLammps) {
  const L = lammpsFor(lessonId).map(findLammps).filter(p => p && p.ready);
  if (!L.length) return '';
  return `<section class="block xlink"><h2>See it inside LAMMPS</h2><p class="muted">The same idea in the real code, run live in your browser and checked by hand.</p><ul>${L.map(p => `<li><a href="/lammps/${p.id}">${p.num} · ${p.title}</a></li>`).join('')}</ul></section>`;
}
export function lammpsBox(lammpsId, findLesson) {
  const L = lessonsFor(lammpsId).map(findLesson).filter(l => l && l.ready);
  if (!L.length) return '';
  return `<section class="block xlink"><h2>The theory behind this page</h2><ul>${L.map(l => `<li><a href="/lessons/${l.id}">${l.lab ? 'Lab' : 'Lesson'} ${l.num} · ${l.title}</a></li>`).join('')}</ul></section>`;
}
