# Make the share images (1200x630 JPEG) that LinkedIn, X and others show for each page.
# Run locally after `npm run build` (needs Python Playwright + Chrome); images go to public/og/,
# then build again so the pages point at them. Only missing images are made unless --force.
#   python3 scripts/make-og.py [--force]
import asyncio, json, html, sys, pathlib
from playwright.async_api import async_playwright

root = pathlib.Path(__file__).resolve().parent.parent
items = json.loads((root / 'dist/og-list.json').read_text())
out = root / 'public/og'; out.mkdir(parents=True, exist_ok=True)
force = '--force' in sys.argv
logo = (root / 'public/favicon.svg').read_text()
font = root / 'node_modules/@fontsource/space-grotesk/files'

def card(kicker, title):
    size = 76 if len(title) < 34 else 64 if len(title) < 56 else 54
    return f'''<!doctype html><html><head><style>
@font-face {{ font-family: SG; font-weight: 400; src: url("file://{font}/space-grotesk-latin-400-normal.woff2"); }}
@font-face {{ font-family: SG; font-weight: 700; src: url("file://{font}/space-grotesk-latin-700-normal.woff2"); }}
html, body {{ margin: 0; width: 1200px; height: 630px; }}
body {{ font-family: SG, sans-serif; color: #F2F5FA; display: flex; flex-direction: column; justify-content: space-between; box-sizing: border-box; padding: 64px 72px;
  background: radial-gradient(900px 520px at 85% -10%, #14305C 0%, transparent 62%), #060C1A; }}
.top {{ display: flex; align-items: center; gap: 18px; font-size: 30px; font-weight: 700; }}
.top svg {{ width: 56px; height: 56px; }}
.kicker {{ color: #38E1FF; font-size: 30px; letter-spacing: .5px; margin-bottom: 18px; }}
h1 {{ font-size: {size}px; line-height: 1.08; margin: 0; font-weight: 700; max-width: 1050px; }}
.foot {{ display: flex; justify-content: space-between; color: #8FA3C7; font-size: 26px; border-top: 2px solid rgba(56,225,255,.25); padding-top: 22px; }}
.foot b {{ color: #FFC24B; font-weight: 400; }}
</style></head><body>
<div class="top">{logo}MD from Scratch</div>
<div><div class="kicker">{html.escape(kicker)}</div><h1>{html.escape(title)}</h1></div>
<div class="foot"><span>Free · hands-on · every number checked by hand</span><b>LAMMPS runs in your browser</b></div>
</body></html>'''

async def main():
    async with async_playwright() as p:
        br = await p.chromium.launch(executable_path='/usr/bin/google-chrome', args=['--allow-file-access-from-files'])
        pg = await br.new_page(viewport={'width': 1200, 'height': 630})
        made = 0
        for it in items:
            f = out / f"{it['key']}.jpg"
            if f.exists() and not force: continue
            tmp = root / 'dist/_og.html'; tmp.write_text(card(it['kicker'], it['title']))
            await pg.goto(tmp.as_uri()); await pg.wait_for_timeout(150)
            await pg.screenshot(path=str(f), type='jpeg', quality=82)
            made += 1
        (root / 'dist/_og.html').unlink(missing_ok=True)
        await br.close()
        print(f'made {made} images, {len(items)} pages')

asyncio.run(main())
