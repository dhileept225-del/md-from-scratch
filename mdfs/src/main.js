import '@fontsource/space-grotesk/400.css';
import { SITE, NAME, HOME_TITLE, HOME_DESC } from './site.js';
import { lessonBox, lammpsBox } from './crosslinks.js';
import '@fontsource/space-grotesk/500.css';
import '@fontsource/space-grotesk/700.css';
import './style.css';
import { inject } from '@vercel/analytics';
inject(); // Vercel Web Analytics: counts page views, including /lessons/<id> navigation
import { chapters, allLessons } from './data/chapters.js';
import { researchParts, allResearch } from './research/index.js';
import lammpsStart from './research/lammps-start.js';
import sourceCompileRun from './research/source-compile-run.js';
import terminalBasics from './research/terminal-basics.js';
import getLammps from './research/get-lammps.js';
import buildLammps from './research/build-lammps.js';
import readCpp from './research/read-cpp.js';
import gitDiff from './research/git-diff.js';
import inputScript from './research/input-script.js';
import readingALine from './research/reading-a-line.js';
import fourBlocks from './research/four-blocks.js';
import unitsPage from './research/units.js';
import buildingTheSystem from './research/building-the-system.js';
import velocitiesPage from './research/velocities.js';
import forcesInScript from './research/forces-in-script.js';
import readData from './research/read-data.js';
import logAndErrors from './research/log-and-errors.js';
import sourceTree from './research/source-tree.js';
import lammpsClass from './research/lammps-class.js';
import stylesPage from './research/styles.js';
import styleFamilies from './research/style-families.js';
import packagesPage from './research/packages.js';
import atomStorage from './research/atom-storage.js';
import { mountLJ } from './sims/lj.js';
import { mountPair } from './sims/pair.js';
import { mountNewton } from './sims/newton.js';
import { mountTimestep } from './sims/timestep.js';
import { mountIntegrators } from './sims/integrators.js';
import { mountEnergyLab } from './sims/energylab.js';
import { mountUnits } from './sims/units.js';
import { mountLattice } from './sims/lattice.js';
import { mountVelocities } from './sims/velocities.js';
import { mountPBC } from './sims/pbc.js';
import { mountLab2 } from './sims/lab2.js';
import { mountLab3 } from './sims/lab3.js';
import { mountLab3b } from './sims/lab3b.js';
import { mountLab4 } from './sims/lab4.js';
import { mountLab5 } from './sims/lab5.js';
import { mountLab6 } from './sims/lab6.js';
import { mountLab7 } from './sims/lab7.js';
import { mountLab8 } from './sims/lab8.js';
import { mountLab9 } from './sims/lab9.js';
import { mountLab10 } from './sims/lab10.js';
import { mountLab11 } from './sims/lab11.js';
import { mountLab12 } from './sims/lab12.js';
import { mountLab13 } from './sims/lab13.js';
import { mountLab14 } from './sims/lab14.js';
import { mountLab15 } from './sims/lab15.js';
import { mountLab16 } from './sims/lab16.js';
import { mountLab17 } from './sims/lab17.js';
import { mountLab1 } from './sims/lab1.js';
import { mountMetad } from './sims/metad.js';
import { mountMlpot } from './sims/mlpot.js';
import { mountRealcodes } from './sims/realcodes.js';
import { mountOnsager } from './sims/onsager.js';
import { mountConductivity } from './sims/conductivity.js';
import { mountVanhove } from './sims/vanhove.js';
import { mountVacf } from './sims/vacf.js';
import { mountMsd } from './sims/msd.js';
import { mountRdf } from './sims/rdf.js';
import { mountPolarizable } from './sims/polarizable.js';
import { mountPme } from './sims/pme.js';
import { mountEwald } from './sims/ewald.js';
import { mountCoulomb } from './sims/coulomb.js';
import { mountWater } from './sims/water.js';
import { mountShake } from './sims/shake.js';
import { mountForceField } from './sims/forcefield.js';
import { mountDihedrals } from './sims/dihedrals.js';
import { mountBonds } from './sims/bonds.js';
import { mountBarostats } from './sims/barostats.js';
import { mountLangevin } from './sims/langevin.js';
import { mountNoseHoover } from './sims/nosehoover.js';
import { mountBerendsen } from './sims/berendsen.js';
import { mountEnsembles } from './sims/ensembles.js';
import { mountTemperature } from './sims/temperature.js';
import { mountNeighbours } from './sims/neighbours.js';
import { mountCutoff } from './sims/cutoff.js';
import whatIsMd from './lessons/what-is-md.js';
import lennardJones from './lessons/lennard-jones.js';
import newtonInCode from './lessons/newton-in-code.js';
import timeStep from './lessons/time-step.js';
import verlet from './lessons/verlet.js';
import energyCheck from './lessons/energy-check.js';
import unitsLesson from './lessons/units.js';
import latticeStart from './lessons/lattices.js';
import maxwellBoltzmann from './lessons/maxwell-boltzmann.js';
import pbcLesson from './lessons/pbc.js';
import labIntegratorsLesson from './lessons/lab-integrators.js';
import labPackingLesson from './lessons/lab-packing.js';
import labMcLesson from './lessons/lab-mc.js';
import labDielectricLesson from './lessons/lab-dielectric.js';
import labStructureFactorLesson from './lessons/lab-structure-factor.js';
import labLifetimesLesson from './lessons/lab-lifetimes.js';
import labViscosityLesson from './lessons/lab-viscosity.js';
import labFiniteSizeLesson from './lessons/lab-finite-size.js';
import labFreeEnergyLesson from './lessons/lab-free-energy.js';
import labErrorsLesson from './lessons/lab-errors.js';
import labBeyondLesson from './lessons/lab-beyond.js';
import labAnalysisLesson from './lessons/lab-analysis.js';
import labElectrostaticsLesson from './lessons/lab-electrostatics.js';
import labMoleculesLesson from './lessons/lab-molecules.js';
import labBarostatsLesson from './lessons/lab-barostats.js';
import labThermostatsLesson from './lessons/lab-thermostats.js';
import labBuildingLesson from './lessons/lab-building.js';
import labBasicsLesson from './lessons/lab-basics.js';
import metadynamicsLesson from './lessons/metadynamics.js';
import mlPotentialsLesson from './lessons/ml-potentials.js';
import realCodesLesson from './lessons/real-codes.js';
import onsagerLesson from './lessons/onsager.js';
import conductivityLesson from './lessons/conductivity.js';
import vanHoveLesson from './lessons/van-hove.js';
import vacfLesson from './lessons/vacf.js';
import msdLesson from './lessons/msd.js';
import rdfLesson from './lessons/rdf.js';
import polarizableLesson from './lessons/polarizable.js';
import pmeLesson from './lessons/pme.js';
import ewaldLesson from './lessons/ewald.js';
import coulombPbcLesson from './lessons/coulomb-pbc.js';
import waterModelsLesson from './lessons/water-models.js';
import shakeLesson from './lessons/shake.js';
import forceFieldsLesson from './lessons/force-fields.js';
import dihedralsLesson from './lessons/dihedrals.js';
import bondsAnglesLesson from './lessons/bonds-angles.js';
import barostatsLesson from './lessons/barostats.js';
import langevinLesson from './lessons/langevin.js';
import noseHooverLesson from './lessons/nose-hoover.js';
import berendsenLesson from './lessons/berendsen.js';
import ensemblesLesson from './lessons/ensembles.js';
import temperatureLesson from './lessons/temperature.js';
import neighbourLists from './lessons/neighbour-lists.js';
import cutoffsLesson from './lessons/cutoffs.js';

const lessons = { 'what-is-md': whatIsMd, 'lennard-jones': lennardJones, 'newton-in-code': newtonInCode, 'time-step': timeStep, verlet, 'energy-check': energyCheck, units: unitsLesson, 'lattice-start': latticeStart, 'maxwell-boltzmann': maxwellBoltzmann, pbc: pbcLesson, cutoffs: cutoffsLesson, 'neighbour-lists': neighbourLists, temperature: temperatureLesson, ensembles: ensemblesLesson, berendsen: berendsenLesson, 'nose-hoover': noseHooverLesson, langevin: langevinLesson, barostats: barostatsLesson, 'bonds-angles': bondsAnglesLesson, dihedrals: dihedralsLesson, 'force-fields': forceFieldsLesson, shake: shakeLesson, 'water-models': waterModelsLesson, 'coulomb-pbc': coulombPbcLesson, ewald: ewaldLesson, pme: pmeLesson, polarizable: polarizableLesson, rdf: rdfLesson, msd: msdLesson, vacf: vacfLesson, 'van-hove': vanHoveLesson, conductivity: conductivityLesson, onsager: onsagerLesson, 'real-codes': realCodesLesson, 'ml-potentials': mlPotentialsLesson, metadynamics: metadynamicsLesson, 'lab-basics': labBasicsLesson, 'lab-building': labBuildingLesson, 'lab-thermostats': labThermostatsLesson, 'lab-barostats': labBarostatsLesson, 'lab-molecules': labMoleculesLesson, 'lab-electrostatics': labElectrostaticsLesson, 'lab-analysis': labAnalysisLesson, 'lab-beyond': labBeyondLesson, 'lab-errors': labErrorsLesson, 'lab-free-energy': labFreeEnergyLesson, 'lab-finite-size': labFiniteSizeLesson, 'lab-viscosity': labViscosityLesson, 'lab-lifetimes': labLifetimesLesson, 'lab-structure-factor': labStructureFactorLesson, 'lab-dielectric': labDielectricLesson, 'lab-mc': labMcLesson, 'lab-packing': labPackingLesson, 'lab-integrators': labIntegratorsLesson };   // add new lesson files here
const sims = { lj: mountLJ, pair: mountPair, newton: mountNewton, timestep: mountTimestep, integrators: mountIntegrators, energylab: mountEnergyLab, units: mountUnits, lattice: mountLattice, velocities: mountVelocities, pbc: mountPBC, cutoff: mountCutoff, neighbours: mountNeighbours, temperature: mountTemperature, ensembles: mountEnsembles, berendsen: mountBerendsen, nosehoover: mountNoseHoover, langevin: mountLangevin, barostats: mountBarostats, bonds: mountBonds, dihedrals: mountDihedrals, forcefield: mountForceField, shake: mountShake, water: mountWater, coulomb: mountCoulomb, ewald: mountEwald, pme: mountPme, polarizable: mountPolarizable, rdf: mountRdf, msd: mountMsd, vacf: mountVacf, vanhove: mountVanhove, conductivity: mountConductivity, onsager: mountOnsager, realcodes: mountRealcodes, mlpot: mountMlpot, metad: mountMetad, lab1: mountLab1, lab2: mountLab2, lab3: mountLab3, lab3b: mountLab3b, lab4: mountLab4, lab5: mountLab5, lab6: mountLab6, lab7: mountLab7, lab8: mountLab8, lab9: mountLab9, lab10: mountLab10, lab11: mountLab11, lab12: mountLab12, lab13: mountLab13, lab14: mountLab14, lab15: mountLab15, lab16: mountLab16, lab17: mountLab17 };
const researchPages = { 'lammps-start': lammpsStart, 'source-compile-run': sourceCompileRun, 'terminal-basics': terminalBasics, 'get-lammps': getLammps, 'build-lammps': buildLammps, 'read-cpp': readCpp, 'git-diff': gitDiff, 'input-script': inputScript, 'reading-a-line': readingALine, 'four-blocks': fourBlocks, units: unitsPage, 'building-the-system': buildingTheSystem, velocities: velocitiesPage, 'forces-in-script': forcesInScript, 'read-data': readData, 'log-and-errors': logAndErrors, 'source-tree': sourceTree, 'lammps-class': lammpsClass, styles: stylesPage, 'style-families': styleFamilies, packages: packagesPage, 'atom-storage': atomStorage };
const app = document.getElementById('app');
let cleanup = [];

const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const badge = lv => `<span class="badge ${lv}">${lv === 'adv' ? 'Advanced' : lv === 'lab' ? 'Lab' : 'Core'}</span>`;

// ---------- pages ----------
function home() {
  const total = allLessons.length, ready = allLessons.filter(l => l.ready).length;
  app.innerHTML = `
  <section class="hero">
    <div class="hero-text">
      <p class="eyebrow">A free, hands-on course</p>
      <h1>Molecular dynamics,<br><span class="c">from scratch.</span></h1>
      <p class="lead">Learn how MD really works by building it yourself, one small idea at a time. Each lesson has a 60-second video, a live simulation and code you can run.</p>
      <div class="row"><a class="btn primary" href="/lessons/what-is-md">Start lesson 1 →</a><a class="btn ghost" href="/lessons">See all lessons</a></div>
      <p class="meta">${chapters.length} chapters · ${total} lessons · ${ready} live so far</p>
    </div>
    <div class="hero-sim" id="hero-sim"></div>
  </section>

  <section class="how">
    ${[['▶', 'Watch', 'A 60-second animated video explains the idea.'], ['◉', 'Play', 'Change the live simulation and see what happens.'], ['{ }', 'Code', 'Read and run the Python, written from scratch.'], ['✓', 'Check', 'Three quick questions to test yourself.']]
      .map(([i, t, d]) => `<div class="card"><div class="ico">${i}</div><h3>${t}</h3><p>${d}</p></div>`).join('')}
  </section>

  <section class="paths">
    <div class="card"><span class="badge core">Core</span><h3>Undergraduate path</h3><p>Follow the Core lessons in order. You'll understand and write a working MD code.</p></div>
    <div class="card"><span class="badge adv">Advanced</span><h3>Graduate path</h3><p>Add the Advanced lessons for thermostats, Ewald sums, transport theory and ML potentials.</p></div>
  </section>

  <section id="lessons-list" class="chapters">
    <h2>Lessons</h2>
    ${chapters.map(c => `
      <div class="chapter">
        <div class="ch-head"><span class="ch-n">${String(c.n).padStart(2, '0')}</span><div><h3>${c.title}</h3><p>${c.blurb}</p></div></div>
        <ol class="lesson-list">
          ${c.lessons.map((l, i) => {
            const inner = `<span class="num">${c.n}.${i + 1}</span><span class="t">${l.title}</span>${badge(l.level)}${l.ready ? '<span class="go">→</span>' : '<span class="soon">Soon</span>'}`;
            return `<li>${l.ready ? `<a href="/lessons/${l.id}">${inner}</a>` : `<div class="off">${inner}</div>`}</li>`;
          }).join('')}
        </ol>
      </div>`).join('')}
  </section>

  <section class="research card">
    <h3>LAMMPS from scratch <span class="badge lab">New</span></h3>
    <p>Run the real LAMMPS in your browser, read its engine, change its code, and check every change with a hand calculation.</p>
    <a class="btn sm primary" href="/lammps">Open the LAMMPS section →</a>
  </section>`;
  cleanup.push(mountLJ(document.getElementById('hero-sim'), { mini: true }));
}

function lessonPage(id) {
  const meta = allLessons.find(l => l.id === id), L = lessons[id];
  if (!meta || !L) return notFound();
  const idx = allLessons.indexOf(meta), prev = allLessons[idx - 1], next = allLessons[idx + 1];
  if (L.lab) return labPage(meta, L, prev, next);
  const v = L.video, abs = p => p && !/^(\/|https?:)/.test(p) ? '/' + p : p;
  const video = v.youtube
    ? `<iframe src="https://www.youtube-nocookie.com/embed/${v.youtube}?rel=0&modestbranding=1&playsinline=1" title="${esc(meta.title)}" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>`
    : `<video src="${abs(v.mp4)}" poster="${abs(v.poster) || ''}" controls playsinline preload="metadata"></video>`;
  const navLink = (l, dir) => !l ? '<span></span>' : l.ready
    ? `<a class="pn ${dir}" href="/lessons/${l.id}"><small>${dir === 'prev' ? '← Previous' : 'Next →'}</small>${l.num} ${l.title}</a>`
    : `<div class="pn ${dir} off"><small>${dir === 'prev' ? '← Previous' : 'Next · coming soon'}</small>${l.num} ${l.title}</div>`;

  app.innerHTML = `
  <nav class="crumb"><a href="/">Home</a> / <a href="/lessons">Chapter ${meta.chapter}</a> / ${meta.num}</nav>
  <section class="lesson-top">
    <div class="video">${video}</div>
    <div class="intro">
      <p class="eyebrow">Lesson ${meta.num} · ${L.minutes} min ${badge(meta.level)}</p>
      <h1>${meta.title}</h1>
      <p class="lead">${L.summary}</p>
      <h4>In this lesson</h4>
      <ul class="learn">${L.learn.map(x => `<li>${x}</li>`).join('')}</ul>
    </div>
  </section>

  <section class="block">
    <h2>Try it</h2>
    <p class="muted">${L.simHint}</p>
    <div id="sim"></div>
  </section>

  <section class="block">
    <div class="tabs" role="tablist">
      ${['Idea', 'Algorithm', 'Code', 'Maths'].map((t, i) => `<button role="tab" class="${i ? '' : 'on'}" data-tab="${t}">${t}</button>`).join('')}
    </div>
    <div class="tab-body prose" data-panel="Idea">${L.idea}</div>
    <div class="tab-body prose" data-panel="Algorithm" hidden>${L.algorithm}</div>
    <div class="tab-body prose" data-panel="Code" hidden>
      <p>${L.code.note}</p>
      <div class="code"><button class="btn sm ghost copy">Copy</button><pre>${highlight(L.code.src)}</pre></div>
    </div>
    <div class="tab-body prose" data-panel="Maths" hidden>${L.maths()}</div>
  </section>

  <section class="block">
    <h2>Check yourself</h2>
    <div class="quiz">${L.quiz.map((q, qi) => `
      <div class="q" data-q="${qi}">
        <p class="qt"><b>${qi + 1}.</b> ${q.q}</p>
        <div class="opts">${q.options.map((o, oi) => `<button data-o="${oi}">${o}</button>`).join('')}</div>
        <p class="why" hidden></p>
      </div>`).join('')}
      <p class="score muted"></p>
    </div>
  </section>
  ${lessonBox(id, x => allResearch.find(p => p.id === x))}
  <nav class="prevnext">${navLink(prev, 'prev')}${navLink(next, 'next')}</nav>`;

  cleanup.push(sims[L.sim](document.getElementById('sim')));

  app.querySelectorAll('[data-tab]').forEach(b => b.onclick = () => {
    app.querySelectorAll('[data-tab]').forEach(x => x.classList.toggle('on', x === b));
    app.querySelectorAll('[data-panel]').forEach(p => p.hidden = p.dataset.panel !== b.dataset.tab);
  });
  const copy = app.querySelector('.copy');
  copy.onclick = async () => {
    try { await navigator.clipboard.writeText(L.code.src); copy.textContent = 'Copied ✓'; }
    catch { copy.textContent = 'Select the code to copy'; }
    setTimeout(() => copy.textContent = 'Copy', 2000);
  };
  let right = 0, done = 0;
  app.querySelectorAll('.q').forEach(qe => {
    const q = L.quiz[+qe.dataset.q];
    qe.querySelectorAll('[data-o]').forEach(b => b.onclick = () => {
      if (qe.classList.contains('done')) return;
      qe.classList.add('done'); done++;
      const ok = +b.dataset.o === q.answer; if (ok) right++;
      b.classList.add(ok ? 'right' : 'wrong');
      qe.querySelector(`[data-o="${q.answer}"]`).classList.add('right');
      const w = qe.querySelector('.why'); w.hidden = false; w.innerHTML = `<b>${ok ? 'Correct.' : 'Not quite.'}</b> ${q.why}`;
      if (done === L.quiz.length) app.querySelector('.score').textContent = `You got ${right} of ${L.quiz.length}. ${right === L.quiz.length ? 'Ready for the next lesson.' : 'Read the Idea tab again, then try the next lesson.'}`;
    });
  });
}

// ---------- lab pages: hand calculations on 1-10 particles (no video) ----------
function labPage(meta, L, prev, next) {
  const navLink = (l, dir) => !l ? '<span></span>' : l.ready
    ? `<a class="pn ${dir}" href="/lessons/${l.id}"><small>${dir === 'prev' ? '← Previous' : 'Next →'}</small>${l.num} ${l.title}</a>`
    : `<div class="pn ${dir} off"><small>${dir === 'prev' ? '← Previous' : 'Next · coming soon'}</small>${l.num} ${l.title}</div>`;
  app.innerHTML = `
  <nav class="crumb"><a href="/">Home</a> / <a href="/lessons">Chapter ${meta.chapter}</a> / ${meta.num}</nav>
  <section class="block lab-top">
    <p class="eyebrow">Lab ${meta.num} · ${L.minutes} min ${badge('lab')}</p>
    <h1>${meta.title}</h1>
    <p class="lead">${L.summary}</p>
    <h4>In this lab</h4>
    <ul class="learn">${L.learn.map(x => `<li>${x}</li>`).join('')}</ul>
    <p class="muted small">You need only a calculator (or the Python check at the end). Every number on this page is printed by that script.</p>
  </section>
  <section class="block lab-setup prose"><h2>The setup</h2>${L.setup}</section>
  ${L.steps.map((st, i) => `
  <section class="block lab-step prose" id="step-${i + 1}">
    <p class="eyebrow">Step ${i + 1}</p>
    <h2>${st.title}</h2>
    <p class="why"><b>Why:</b> ${st.why}</p>
    ${st.work}
    <div class="see"><b>What the numbers say:</b> ${st.see}</div>
  </section>`).join('')}
  <section class="block lab-summary prose"><h2>Summary</h2>${L.summaryTable}</section>
  <section class="block">
    <h2>Run it yourself</h2>
    <p class="muted">${L.simHint}</p>
    <div id="sim"></div>
  </section>
  <section class="block prose">
    <h2>Check every number in Python</h2>
    <p>${L.code.note}</p>
    <div class="code"><button class="btn sm ghost copy">Copy</button><pre>${highlight(L.code.src)}</pre></div>
  </section>
  <section class="block prose">
    <h2>Your turn</h2>
    ${L.exercises.map((e, i) => `<div class="ex"><p><b>${i + 1}.</b> ${e.q}</p><button class="btn sm ghost" data-ex="${i}">Show answer</button><div class="ans" hidden>${e.a}</div></div>`).join('')}
  </section>
  ${lessonBox(meta.id, x => allResearch.find(p => p.id === x))}
  <nav class="prevnext">${navLink(prev, 'prev')}${navLink(next, 'next')}</nav>`;
  cleanup.push(sims[L.sim](document.getElementById('sim')));
  app.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
  const copy = app.querySelector('.copy');
  copy.onclick = async () => {
    try { await navigator.clipboard.writeText(L.code.src); copy.textContent = 'Copied ✓'; }
    catch { copy.textContent = 'Select the code to copy'; }
    setTimeout(() => copy.textContent = 'Copy', 2000);
  };
}

// ---------- research: LAMMPS from scratch ----------
function researchHome() {
  const ready = allResearch.filter(p => p.ready).length;
  app.innerHTML = `
  <nav class="crumb"><a href="/">Home</a> / LAMMPS</nav>
  <section class="block lab-top prose">
    <p class="eyebrow">LAMMPS · ${ready} of ${allResearch.length} pages live</p>
    <h1>LAMMPS from scratch</h1>
    <p class="lead">Understand how LAMMPS, the most widely used molecular dynamics program, really works. Run it in your browser, read its engine, change its code, and check every change with a hand calculation. You need school maths and no programming.</p>
  </section>
  <section class="chapters">
    ${researchParts.map(p => `
      <div class="chapter research-part">
        <div class="ch-head"><span class="ch-n">${p.key}</span><div><h3>${p.title}</h3><p>${p.blurb}</p></div></div>
        <ol class="lesson-list">
          ${p.pages.map((pg, i) => {
            const inner = `<span class="num">${p.key === 'S' ? '▶' : `${p.key}.${i + 1}`}</span><span class="t">${pg.title}</span>${pg.ready ? '<span class="go">→</span>' : '<span class="soon">Soon</span>'}`;
            return `<li>${pg.ready ? `<a href="/lammps/${pg.id}">${inner}</a>` : `<div class="off">${inner}</div>`}</li>`;
          }).join('')}
        </ol>
      </div>`).join('')}
  </section>`;
}

function researchPage(id) {
  const meta = allResearch.find(p => p.id === id), P = researchPages[id];
  if (!meta || !P) return notFound();
  const ready = allResearch.filter(p => p.ready), i = ready.indexOf(meta), prev = ready[i - 1], next = ready[i + 1];
  const link = (p, dir) => p ? `<a class="pn ${dir}" href="/lammps/${p.id}"><small>${dir === 'prev' ? '← Previous' : 'Next →'}</small>${p.title}</a>` : `<a class="pn ${dir}" href="/lammps"><small>${dir === 'prev' ? '← Back' : 'More pages coming'}</small>All LAMMPS pages</a>`;
  app.innerHTML = `
  <nav class="crumb"><a href="/">Home</a> / <a href="/lammps">LAMMPS</a> / ${meta.num}</nav>
  <section class="block lab-top">
    <p class="eyebrow">LAMMPS from scratch · ${meta.num === 'Start' ? 'Start here' : meta.num} · ${P.minutes} min</p>
    <h1>${meta.title}</h1>
    <p class="lead">${P.summary}</p>
    <h4>On this page</h4>
    <ul class="learn">${P.learn.map(x => `<li>${x}</li>`).join('')}</ul>
  </section>
  ${P.html()}
  ${lammpsBox(id, x => allLessons.find(l => l.id === x))}
  <nav class="prevnext">${link(prev, 'prev')}${link(next, 'next')}</nav>`;
  cleanup.push(P.mount(app));
}

function notFound() { app.innerHTML = `<section class="block"><h1>Page not found</h1><p><a href="/">Back to home</a></p></section>`; }

function highlight(src) {
  const kw = /\b(import|as|def|return|for|in|if|range|round|print)\b/g;
  return esc(src).split('\n').map(line => {
    const h = line.indexOf('#');
    const code = h >= 0 ? line.slice(0, h) : line, com = h >= 0 ? line.slice(h) : '';
    const c = code.replace(/(f?"[^"]*")/g, '<span class="s">$1</span>').replace(kw, '<span class="k">$1</span>').replace(/\b(\d+\.?\d*)\b/g, '<span class="n">$1</span>');
    return c + (com ? `<span class="c">${com}</span>` : '');
  }).join('\n');
}

// ---------- page titles and descriptions (search engines and link previews) ----------
const plain = s => String(s).replace(/<[^>]+>/g, '').replace(/\s+/g, ' ').trim();
const short = s => { const p = plain(s); return p.length > 155 ? p.slice(0, 152).replace(/\s\S*$/, '') + '…' : p; };
function setMeta(title, desc, path) {
  document.title = title;
  const set = (sel, attr, val) => { const e = document.head.querySelector(sel); if (e) e.setAttribute(attr, val); };
  set('meta[name="description"]', 'content', desc); set('link[rel="canonical"]', 'href', SITE + path);
  set('meta[property="og:title"]', 'content', title); set('meta[property="og:description"]', 'content', desc); set('meta[property="og:url"]', 'content', SITE + path);
}

// ---------- router: /, /lessons, /lessons/<id>  (old #lesson-<id> links still work) ----------
function legacyHash() {
  const h = location.hash.slice(1);
  if (h.startsWith('lesson-')) history.replaceState(null, '', '/lessons/' + h.slice(7));
  else if (h === 'lessons') history.replaceState(null, '', '/lessons');
  else if (h === 'home') history.replaceState(null, '', '/');
}
function route() {
  cleanup.forEach(f => f()); cleanup = [];
  if (/^\/research(\/|$)/.test(location.pathname)) history.replaceState(null, '', location.pathname.replace(/^\/research/, '/lammps') + location.hash);
  const p = location.pathname.replace(/\/+$/, '') || '/';
  const m = p.match(/^\/lessons\/([\w-]+)$/), r = p.match(/^\/lammps(?:\/([\w-]+))?$/);
  if (r) {
    window.scrollTo(0, 0);
    if (!r[1]) { researchHome(); setMeta(`LAMMPS from scratch | ${NAME}`, 'Understand how LAMMPS works: run it in your browser, read its engine, change its code and check every change with a hand calculation.', '/lammps'); }
    else { const meta = allResearch.find(x => x.id === r[1]), P = researchPages[r[1]]; researchPage(r[1]); setMeta(meta && P ? `${meta.title} | LAMMPS from scratch` : `Page not found | ${NAME}`, P ? short(P.summary) : HOME_DESC, p); }
  } else if (m) {
    const meta = allLessons.find(l => l.id === m[1]), L = lessons[m[1]];
    lessonPage(m[1]); window.scrollTo(0, 0);
    if (meta && L) setMeta(`${meta.title} | ${NAME}`, short(L.summary), p); else setMeta(`Page not found | ${NAME}`, HOME_DESC, p);
  } else {
    home();
    if (p === '/lessons') { document.getElementById('lessons-list').scrollIntoView(); setMeta(`All lessons | ${NAME}`, HOME_DESC, '/lessons'); }
    else { window.scrollTo(0, 0); setMeta(HOME_TITLE, HOME_DESC, '/'); }
  }
}
document.addEventListener('click', e => {
  const a = e.target.closest('a[href^="/"]');
  if (!a || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey || a.target === '_blank' || /\.[a-z0-9]+$/i.test(a.getAttribute('href'))) return;
  e.preventDefault();
  if (a.getAttribute('href') !== location.pathname) history.pushState(null, '', a.getAttribute('href'));
  route();
});
addEventListener('popstate', route);
addEventListener('hashchange', () => { if (/^#(lesson-|lessons$|home$)/.test(location.hash)) { legacyHash(); route(); } });
legacyHash();
route();
