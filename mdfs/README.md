# MD from Scratch

A free, hands-on course: learn molecular dynamics by building it.
Each lesson = 60-second video + live simulation + Idea / Algorithm / Code / Maths tabs + quiz.

## Run it on your computer
1. Install Node.js (LTS) from https://nodejs.org
2. In this folder: `npm install` then `npm run dev`
3. Open the link it prints (http://localhost:5173)

## Put it online (free, Vercel)
1. Create a GitHub repo and upload this folder (not `node_modules` or `dist`).
2. Go to https://vercel.com, sign in with GitHub, click **Add New → Project**, pick the repo.
3. Vercel detects Vite automatically. Click **Deploy**. Every push to GitHub updates the site.

## Add your YouTube video to a lesson
1. Upload the MP4 to YouTube (Public or Unlisted; embedding is on by default).
2. Open the video. The URL looks like `https://www.youtube.com/watch?v=AbC123xYz`.
3. Copy the part after `v=` (here `AbC123xYz`).
4. In `src/lessons/what-is-md.js`, set `video: { youtube: 'AbC123xYz', ... }`.
   Leave it empty to play the local MP4 in `public/videos/`.

## Add a new lesson
1. Copy `src/lessons/what-is-md.js` to `src/lessons/<id>.js` and edit the text.
2. In `src/main.js`, import it and add it to `lessons`.
3. In `src/data/chapters.js`, set `ready: true` on that lesson.

## Project layout
- `src/data/chapters.js`: the course map (chapters, lessons, Core/Advanced)
- `src/lessons/`: one file per lesson
- `src/sims/lj.js`: the live Lennard-Jones MD engine (velocity Verlet + thermostat)
- `public/videos/`: local video files and posters
