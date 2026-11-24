#!/usr/bin/env node
/*
 * Build reader-independent versions of the Markdown documents.
 *
 *   docs/<name>.html  single self-contained file: math pre-rendered by KaTeX,
 *                     KaTeX fonts and all figures embedded as data URIs.
 *                     Opens identically in any browser, offline.
 *   docs/<name>.pdf   printed from that HTML (A4, page numbers).
 *
 * Source of truth stays the Markdown in lecture/ and review/, written in
 * GitHub math syntax:  inline $`...`$   and   ```math fenced blocks.
 *
 * Usage:  cd tools && npm install && node build_docs.js [--no-pdf]
 */
const fs = require('fs');
const path = require('path');
const katex = require('katex');
const MarkdownIt = require('markdown-it');

const ROOT = path.resolve(__dirname, '..');
const OUT = path.join(ROOT, 'docs');
const REPO_URL = 'https://github.com/turkeyck/Computational-EM-from-MBP-to-topological-photonics-and-adjoint-methods/blob/main/';

const DOCS = [
  { src: 'lecture/lecture_notes.md', name: 'lecture_notes', title: '計算電磁學與拓撲光子學：從平面波展開到可微分 RCWA' },
  { src: 'review/lecture_review.md', name: 'lecture_review', title: '講義審閱報告' },
  { src: 'review/changelog.md', name: 'changelog', title: '修訂紀錄：原稿 → 修訂版' },
];

const MIME = { '.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg', '.svg': 'image/svg+xml', '.gif': 'image/gif' };

function dataUri(file) {
  const ext = path.extname(file).toLowerCase();
  return `data:${MIME[ext] || 'application/octet-stream'};base64,${fs.readFileSync(file).toString('base64')}`;
}

function katexCss() {
  const cssPath = require.resolve('katex/dist/katex.min.css');
  const fontDir = path.join(path.dirname(cssPath), 'fonts');
  let css = fs.readFileSync(cssPath, 'utf8');
  // keep only woff2 sources and inline them
  css = css.replace(/src:([^;}]*)/g, (m, srcs) => {
    const w = srcs.match(/url\(fonts\/([^)]+\.woff2)\)/);
    return w ? `src:url(${dataUri(path.join(fontDir, w[1])).replace('application/octet-stream', 'font/woff2')}) format("woff2")` : m;
  });
  return css;
}

/* Replace math with placeholders before Markdown parsing, so Markdown never touches TeX. */
function extractMath(md) {
  const store = [];
  // private-use characters as placeholders: markdown-it rewrites NUL to U+FFFD
  const keep = (html) => { store.push(html); return `\uE000${store.length - 1}\uE001`; };
  const render = (tex, display) => katex.renderToString(tex, { displayMode: display, throwOnError: true, strict: false, output: 'html' });

  // display blocks: ```math ... ```
  md = md.replace(/^```math\n([\s\S]*?)\n```[ \t]*$/gm, (m, tex) => '\n' + keep(`<div class="math-display">${render(tex, true)}</div>`) + '\n');

  // inline $`...`$, but not inside other code spans or fenced code
  const lines = md.split('\n');
  let inFence = false;
  for (let i = 0; i < lines.length; i++) {
    if (/^\s*```/.test(lines[i])) { inFence = !inFence; continue; }
    if (inFence) continue;
    let out = '', rest = lines[i];
    const span = /(?<!`)(`+)(?!`)[\s\S]*?(?<!`)\1(?!`)/g;
    let pos = 0, m;
    while ((m = span.exec(rest))) {
      const before = rest.slice(pos, m.index);
      // a $`...`$ math span is a code span immediately wrapped by $ on both sides
      if (before.endsWith('$') && rest[m.index + m[0].length] === '$' && m[1].length === 1) {
        out += before.slice(0, -1) + keep(render(m[0].slice(1, -1), false));
        pos = m.index + m[0].length + 1;
        span.lastIndex = pos;
      } else {
        out += before + m[0];
        pos = m.index + m[0].length;
      }
    }
    lines[i] = out + rest.slice(pos);
  }
  return { md: lines.join('\n'), store };
}

// Short ASCII ids: long CJK anchors become over-long PDF name tokens that strict readers reject.
function slugify(text, used) {
  const k = `s${used.size + 1}`;
  used.add(k);
  return k;
}

function build(doc) {
  const srcPath = path.join(ROOT, doc.src);
  const srcDir = path.dirname(srcPath);
  const { md, store } = extractMath(fs.readFileSync(srcPath, 'utf8'));

  const mdit = new MarkdownIt({ html: true, linkify: false, typographer: false });
  const defaultImage = mdit.renderer.rules.image;
  mdit.renderer.rules.image = (tokens, idx, opts, env, self) => {
    const t = tokens[idx];
    const src = t.attrGet('src');
    if (src && !/^(https?:|data:)/.test(src)) t.attrSet('src', dataUri(path.resolve(srcDir, src)));
    return defaultImage(tokens, idx, opts, env, self);
  };
  const defaultLink = mdit.renderer.rules.link_open || ((tokens, idx, opts, env, self) => self.renderToken(tokens, idx, opts));
  mdit.renderer.rules.link_open = (tokens, idx, opts, env, self) => {
    const t = tokens[idx];
    const href = t.attrGet('href');
    if (href && !/^(https?:|mailto:|#)/.test(href)) {
      const rel = path.relative(ROOT, path.resolve(srcDir, href)).split(path.sep).join('/');
      const local = DOCS.find(d => d.src === rel);
      t.attrSet('href', local ? `${local.name}.html` : REPO_URL + rel);   // sibling docs stay local
    }
    return defaultLink(tokens, idx, opts, env, self);
  };

  let html = mdit.render(md);
  html = html.replace(/\uE000(\d+)\uE001/g, (m, i) => store[+i]);
  if (/[\uE000\uE001\uFFFD]/.test(html)) throw new Error(`${doc.src}: math placeholder left in output`);

  // heading ids + table of contents (h2/h3)
  const used = new Set();
  const toc = [];
  html = html.replace(/<h([23])>([\s\S]*?)<\/h\1>/g, (m, lvl, inner) => {
    const id = slugify(inner.replace(/<[^>]+>/g, ''), used);
    toc.push({ lvl: +lvl, id, text: inner.replace(/<[^>]+>/g, '') });
    const major = lvl === '2' && /^(第\s*\d+\s*講|附錄|參考文獻)/.test(inner.replace(/<[^>]+>/g, ''));
    return `<h${lvl} id="${id}"${major ? ' class="major"' : ''}>${inner}</h${lvl}>`;
  });
  const tocHtml = toc.length > 4
    ? `<nav class="toc"><p class="toc-title">目錄</p><ul>${toc.map(h => `<li class="l${h.lvl}"><a href="#${h.id}">${h.text}</a></li>`).join('')}</ul></nav>`
    : '';
  // TOC goes right before the first section heading (after title, byline and any intro note)
  html = html.replace(/<h2 /, `${tocHtml}<h2 `);

  const page = `<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${doc.title}</title>
<style>${katexCss()}</style>
<style>
:root{--fg:#1f2328;--muted:#59636e;--line:#d1d9e0;--bg:#fff;--soft:#f6f8fa;--accent:#0969da}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font-family:"Noto Sans TC","PingFang TC","Microsoft JhengHei","WenQuanYi Zen Hei","Noto Sans CJK TC",system-ui,sans-serif;
  font-size:16px;line-height:1.75}
main{max-width:860px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:1.9em;line-height:1.3;border-bottom:1px solid var(--line);padding-bottom:.3em}
h2{font-size:1.45em;margin-top:2.2em;border-bottom:1px solid var(--line);padding-bottom:.25em}
h3{font-size:1.15em;margin-top:1.6em}
a{color:var(--accent);text-decoration:none} a:hover{text-decoration:underline}
blockquote{margin:1em 0;padding:.6em 1em;border-left:4px solid var(--accent);background:var(--soft);color:var(--fg)}
blockquote p{margin:.3em 0}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.88em;background:var(--soft);padding:.1em .35em;border-radius:4px}
pre{background:var(--soft);padding:12px;border-radius:6px;overflow-x:auto} pre code{background:none;padding:0}
table{border-collapse:collapse;margin:1em 0;display:block;overflow-x:auto;max-width:100%}
th,td{border:1px solid var(--line);padding:6px 10px;vertical-align:top}
th{background:var(--soft)}
img{max-width:100%;height:auto;display:block;margin:1em auto}
em{color:var(--muted)}
.math-display{overflow-x:auto;overflow-y:hidden;margin:1em 0;padding:2px 0}
.katex{font-size:1.05em}
.toc{background:var(--soft);border:1px solid var(--line);border-radius:6px;padding:12px 20px;margin:1.5em 0}
.toc-title{font-weight:bold;margin:0 0 .4em}
.toc ul{list-style:none;padding:0;margin:0;columns:2;column-gap:24px}
.toc li{break-inside:avoid;line-height:1.6} .toc li.l3{padding-left:1.2em;font-size:.92em}
@media (max-width:640px){.toc ul{columns:1} body{font-size:15px}}
@media print{
  main{max-width:none;padding:0}
  h2.major{break-before:page}
  h2,h3{break-after:avoid} img,table,.math-display,blockquote{break-inside:avoid}
  a{color:var(--fg)}
}
</style>
</head>
<body><main>
${html}
<hr><p style="color:var(--muted);font-size:.85em">本檔由 <code>tools/build_docs.js</code> 從 <code>${doc.src}</code> 自動產生，請勿直接編輯。</p>
</main></body></html>`;

  fs.mkdirSync(OUT, { recursive: true });
  const outHtml = path.join(OUT, `${doc.name}.html`);
  fs.writeFileSync(outHtml, page);
  return { outHtml, nMath: store.length };
}

async function toPdf(htmlFile, pdfFile, title) {
  let chromium;
  try { ({ chromium } = require('playwright')); } catch { console.warn('playwright not installed: skipping PDF'); return false; }
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + htmlFile, { waitUntil: 'load' });
  await page.evaluate(() => document.fonts.ready);
  await page.pdf({
    path: pdfFile, format: 'A4', printBackground: true,
    margin: { top: '18mm', bottom: '18mm', left: '16mm', right: '16mm' },
    displayHeaderFooter: true,
    headerTemplate: `<div style="font-size:8px;color:#888;width:100%;text-align:center">${title}</div>`,
    footerTemplate: '<div style="font-size:8px;color:#888;width:100%;text-align:center"><span class="pageNumber"></span> / <span class="totalPages"></span></div>',
  });
  await browser.close();
  return true;
}

(async () => {
  const noPdf = process.argv.includes('--no-pdf');
  for (const doc of DOCS) {
    const { outHtml, nMath } = build(doc);
    let msg = `${doc.src} -> docs/${doc.name}.html (${nMath} formulas)`;
    if (!noPdf && await toPdf(outHtml, path.join(OUT, `${doc.name}.pdf`), doc.title)) msg += `, docs/${doc.name}.pdf`;
    console.log(msg);
  }
})().catch(e => { console.error(e); process.exit(1); });
