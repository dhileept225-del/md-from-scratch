// After `vite build`: write one real HTML page per lesson so search engines (and link previews)
// see a proper title, description and the lesson text. The app replaces the content on load.
// Output: dist/lessons/<id>.html, dist/lessons.html, dist/404.html, dist/sitemap.xml, dist/robots.txt
import { createServer } from 'vite';
import fs from 'node:fs';
import path from 'node:path';

import { SITE, NAME, HOME_TITLE, HOME_DESC } from '../src/site.js';
import { lessonBox, lammpsBox } from '../src/crosslinks.js';
const root = path.resolve(path.dirname(new URL(import.meta.url).pathname), '..');
const dist = path.join(root, 'dist');

const vite = await createServer({ root, logLevel: 'error', server: { middlewareMode: true }, appType: 'custom' });
const { chapters, allLessons } = await vite.ssrLoadModule('/src/data/chapters.js');
const lessons = {};
for (const f of fs.readdirSync(path.join(root, 'src/lessons')).filter(f => f.endsWith('.js'))) {
  const L = (await vite.ssrLoadModule(`/src/lessons/${f}`)).default;
  if (L && L.id) lessons[L.id] = L;
}
const { researchParts, allResearch } = await vite.ssrLoadModule('/src/research/index.js');
const researchPages = {};
for (const p of allResearch.filter(p => p.ready)) researchPages[p.id] = (await vite.ssrLoadModule(`/src/research/${p.id}.js`)).default;
await vite.close();

const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const plain = s => String(s).replace(/<[^>]+>/g, '').replace(/\s+/g, ' ').trim();
const short = s => { const p = plain(s); return p.length > 155 ? p.slice(0, 152).replace(/\s\S*$/, '') + '…' : p; };
const template = fs.readFileSync(path.join(dist, 'index.html'), 'utf8');

// share image: public/og/<key>.jpg made by scripts/make-og.py; falls back to the site card
const og = [];
const ogImg = (key, kicker, title) => { og.push({ key, kicker, title }); return fs.existsSync(path.join(root, 'public/og', key + '.jpg')) ? `${SITE}/og/${key}.jpg` : `${SITE}/og/home.jpg`; };
const crumbs = list => ({ '@context': 'https://schema.org', '@type': 'BreadcrumbList', itemListElement: list.map(([name, u], i) => ({ '@type': 'ListItem', position: i + 1, name, item: SITE + u })) });

function page({ title, desc, urlPath, body, type = 'website', jsonld, noindex, image, crumb }) {
  let h = template;
  const rep = (re, val) => { if (!re.test(h)) throw new Error(`template is missing ${re}`); h = h.replace(re, val); };
  rep(/<title>[^<]*<\/title>/, `<title>${esc(title)}</title>`);
  rep(/(<meta name="description" content=")[^"]*/, `$1${esc(desc)}`);
  rep(/(<link rel="canonical" href=")[^"]*/, `$1${SITE}${urlPath}`);
  rep(/(<meta property="og:type" content=")[^"]*/, `$1${type}`);
  rep(/(<meta property="og:title" content=")[^"]*/, `$1${esc(title)}`);
  rep(/(<meta property="og:description" content=")[^"]*/, `$1${esc(desc)}`);
  rep(/(<meta property="og:url" content=")[^"]*/, `$1${SITE}${urlPath}`);
  rep(/(<meta property="og:image" content=")[^"]*/, `$1${image || SITE + '/og/home.jpg'}`);
  const extra = (noindex ? '<meta name="robots" content="noindex">\n' : '') + (jsonld ? `<script type="application/ld+json">${JSON.stringify(jsonld)}</script>\n` : '')
    + (crumb ? `<script type="application/ld+json">${JSON.stringify(crumbs(crumb))}</script>\n` : '');
  h = h.replace('</head>', `${extra}</head>`);
  rep(/<main id="app"><\/main>/, `<main id="app">${body}</main>`);
  return h;
}

const listHtml = () => chapters.map(c => `
  <div class="chapter"><div class="ch-head"><span class="ch-n">${String(c.n).padStart(2, '0')}</span><div><h3>${c.title}</h3><p>${c.blurb}</p></div></div>
  <ol class="lesson-list">${c.lessons.map((l, i) => `<li>${l.ready ? `<a href="/lessons/${l.id}"><span class="num">${c.n}.${i + 1}</span><span class="t">${esc(l.title)}</span></a>` : `<div class="off"><span class="num">${c.n}.${i + 1}</span><span class="t">${esc(l.title)}</span></div>`}</li>`).join('')}</ol></div>`).join('');

const course = { '@type': 'Course', name: NAME, description: HOME_DESC, url: SITE + '/', inLanguage: 'en', isAccessibleForFree: true,
  provider: { '@type': 'Person', name: 'Thummuru Dhileep', affiliation: 'Okayama University' } };

// home (rewrite dist/index.html itself with static content) and /lessons
const homeBody = `<section class="hero"><div class="hero-text"><p class="eyebrow">A free, hands-on course</p><h1>Molecular dynamics, from scratch.</h1>
  <p class="lead">${HOME_DESC}</p></div></section><section id="lessons-list" class="chapters"><h2>Lessons</h2>${listHtml()}</section>`;
fs.writeFileSync(path.join(dist, 'index.html'), page({ title: HOME_TITLE, desc: HOME_DESC, urlPath: '/', body: homeBody, jsonld: { '@context': 'https://schema.org', ...course }, image: ogImg('home', 'A free, hands-on course', 'Molecular dynamics, from scratch') }));
fs.writeFileSync(path.join(dist, 'lessons.html'), page({ title: `All lessons | ${NAME}`, desc: HOME_DESC, urlPath: '/lessons', body: homeBody, jsonld: { '@context': 'https://schema.org', ...course, url: SITE + '/lessons' },
  image: ogImg('lessons', `${allLessons.filter(l => l.ready).length} lessons and labs`, 'Molecular dynamics, from scratch'), crumb: [['Home', '/'], ['Lessons', '/lessons']] }));
fs.writeFileSync(path.join(dist, '404.html'), page({ title: `Page not found | ${NAME}`, desc: HOME_DESC, urlPath: '/404', body: '', noindex: true }));

// one page per ready lesson
fs.mkdirSync(path.join(dist, 'lessons'), { recursive: true });
const ready = allLessons.filter(l => l.ready && lessons[l.id]);
const missing = allLessons.filter(l => l.ready && !lessons[l.id]).map(l => l.id);
if (missing.length) throw new Error('ready lessons without a lesson file: ' + missing.join(', '));
for (const meta of ready) {
  const L = lessons[meta.id], i = allLessons.indexOf(meta), prev = allLessons[i - 1], next = allLessons[i + 1];
  const kind = L.lab ? 'Lab' : 'Lesson';
  const nav = [prev, next].filter(l => l && l.ready).map(l => `<a href="/lessons/${l.id}">${l.num} ${esc(l.title)}</a>`).join(' · ');
  const body = `<nav class="crumb"><a href="/">Home</a> / <a href="/lessons">Chapter ${meta.chapter}</a> / ${meta.num}</nav>
  <section class="block"><p class="eyebrow">${kind} ${meta.num} · ${L.minutes} min</p><h1>${esc(meta.title)}</h1><p class="lead">${L.summary}</p>
  <ul class="learn">${L.learn.map(x => `<li>${x}</li>`).join('')}</ul></section>
  ${L.lab
    ? `<section class="block prose"><h2>The setup</h2>${L.setup}</section>${L.steps.map((s, k) => `<section class="block prose"><p class="eyebrow">Step ${k + 1}</p><h2>${s.title}</h2><p><b>Why:</b> ${s.why}</p>${s.work}<p><b>What the numbers say:</b> ${s.see}</p></section>`).join('')}<section class="block prose"><h2>Summary</h2>${L.summaryTable}</section>`
    : `<section class="block prose"><h2>The idea</h2>${L.idea}</section><section class="block prose"><h2>The algorithm</h2>${L.algorithm}</section>`}
  ${lessonBox(meta.id, x => allResearch.find(p => p.id === x))}
  <nav class="prevnext">${nav}</nav>`;
  const jsonld = { '@context': 'https://schema.org', '@type': 'LearningResource', name: meta.title, description: plain(L.summary), url: `${SITE}/lessons/${meta.id}`,
    learningResourceType: L.lab ? 'Lab' : 'Lesson', educationalLevel: meta.level === 'adv' ? 'Graduate' : 'Undergraduate', timeRequired: `PT${L.minutes}M`,
    inLanguage: 'en', isAccessibleForFree: true, teaches: L.learn.map(plain), isPartOf: course };
  fs.writeFileSync(path.join(dist, 'lessons', `${meta.id}.html`),
    page({ title: `${meta.title} | ${NAME}`, desc: short(L.summary), urlPath: `/lessons/${meta.id}`, body, type: 'article', jsonld,
      image: ogImg('lesson-' + meta.id, `${kind} ${meta.num}`, meta.title), crumb: [['Home', '/'], ['Lessons', '/lessons'], [meta.title, `/lessons/${meta.id}`]] }));
}

// research section
const RES_DESC = 'Understand how LAMMPS works: run it in your browser, read its engine, change its code and check every change with a hand calculation.';
const resList = researchParts.map(p => `<div class="chapter"><div class="ch-head"><span class="ch-n">${p.key}</span><div><h3>${p.title}</h3><p>${p.blurb}</p></div></div>
  <ol class="lesson-list">${p.pages.map(pg => `<li>${pg.ready ? `<a href="/lammps/${pg.id}"><span class="t">${esc(pg.title)}</span></a>` : `<div class="off"><span class="t">${esc(pg.title)}</span></div>`}</li>`).join('')}</ol></div>`).join('');
const lammpsCourse = { '@context': 'https://schema.org', '@type': 'Course', name: 'LAMMPS from scratch', description: RES_DESC, url: SITE + '/lammps', inLanguage: 'en', isAccessibleForFree: true, provider: course.provider };
fs.writeFileSync(path.join(dist, 'lammps.html'), page({ title: `LAMMPS from scratch | ${NAME}`, desc: RES_DESC, urlPath: '/lammps', jsonld: lammpsCourse,
  image: ogImg('lammps', 'Run it in your browser, check it by hand', 'LAMMPS from scratch'), crumb: [['Home', '/'], ['LAMMPS', '/lammps']],
  body: `<nav class="crumb"><a href="/">Home</a> / LAMMPS</nav><section class="block"><h1>LAMMPS from scratch</h1><p class="lead">${RES_DESC}</p></section><section class="chapters">${resList}</section>` }));
fs.mkdirSync(path.join(dist, 'lammps'), { recursive: true });
const readyRes = allResearch.filter(p => p.ready);
for (const meta of readyRes) {
  const P = researchPages[meta.id];
  const body = `<nav class="crumb"><a href="/">Home</a> / <a href="/lammps">LAMMPS</a> / ${meta.num}</nav>
  <section class="block"><p class="eyebrow">LAMMPS from scratch · ${P.minutes} min</p><h1>${esc(meta.title)}</h1><p class="lead">${P.summary}</p><ul class="learn">${P.learn.map(x => `<li>${x}</li>`).join('')}</ul></section>${P.html()}
  ${lammpsBox(meta.id, x => allLessons.find(l => l.id === x))}`;
  const jsonld = { '@context': 'https://schema.org', '@type': 'LearningResource', name: meta.title, description: plain(P.summary), url: `${SITE}/lammps/${meta.id}`,
    learningResourceType: 'Tutorial', educationalLevel: 'Graduate', timeRequired: `PT${P.minutes}M`, inLanguage: 'en', isAccessibleForFree: true, teaches: P.learn.map(plain), isPartOf: lammpsCourse };
  fs.writeFileSync(path.join(dist, 'lammps', `${meta.id}.html`), page({ title: `${meta.title} | LAMMPS from scratch`, desc: short(P.summary), urlPath: `/lammps/${meta.id}`, body, type: 'article', jsonld,
    image: ogImg('lammps-' + meta.id, `LAMMPS from scratch · ${meta.num}`, meta.title), crumb: [['Home', '/'], ['LAMMPS', '/lammps'], [meta.title, `/lammps/${meta.id}`]] }));
}

// sitemap and robots
const today = new Date().toISOString().slice(0, 10);
const urls = ['/', '/lessons', ...ready.map(l => `/lessons/${l.id}`), '/lammps', ...readyRes.map(p => `/lammps/${p.id}`)];
fs.writeFileSync(path.join(dist, 'sitemap.xml'), `<?xml version="1.0" encoding="UTF-8"?>\n<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n${urls.map(u => `  <url><loc>${SITE}${u}</loc><lastmod>${today}</lastmod></url>`).join('\n')}\n</urlset>\n`);
fs.writeFileSync(path.join(dist, 'robots.txt'), `User-agent: *\nAllow: /\n\nSitemap: ${SITE}/sitemap.xml\n`);
fs.writeFileSync(path.join(dist, 'og-list.json'), JSON.stringify(og, null, 1));
console.log(`prerendered ${ready.length} lesson pages + home, /lessons, 404; sitemap with ${urls.length} URLs`);
