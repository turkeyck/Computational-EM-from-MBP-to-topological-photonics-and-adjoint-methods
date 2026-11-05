# Computational EM: from MPB to topological photonics and adjoint methods

## 直接閱讀（任何閱讀器都能開）

| 文件 | PDF | 單檔 HTML（離線可開） |
| --- | --- | --- |
| 講義修訂版 | [`docs/lecture_notes.pdf`](docs/lecture_notes.pdf) | [`docs/lecture_notes.html`](docs/lecture_notes.html) |
| 審閱報告 | [`docs/lecture_review.pdf`](docs/lecture_review.pdf) | [`docs/lecture_review.html`](docs/lecture_review.html) |
| 修訂紀錄 | [`docs/changelog.pdf`](docs/changelog.pdf) | [`docs/changelog.html`](docs/changelog.html) |

PDF 與 HTML 中的數學式已預先算繪，圖片與字型都已內嵌，不依賴閱讀器的 Markdown 或數學外掛。HTML 檔請下載後用瀏覽器開啟（GitHub 網頁上只會顯示原始碼）。

## 檔案

- `lecture/lecture_notes.md` — **講義修訂版**（2026-11-04，已改正錯誤並補上圖與 demo 連結）
- `lecture/lecture_notes_original.md` — 講義原稿（2026-10-29 版）
- `review/changelog.md` — 原稿 → 修訂版的逐條修改紀錄
- `review/lecture_review.md` — 審閱報告：編排、正確性、引用、術語、排版、圖與 demo
- `demos/` — 最小可執行範例，兼作講義公式的數值驗證（PWE、FHS Chern 數、1D RCWA TE/TM、示意圖）
- `figures/` — demo 產生的圖

```bash
pip install numpy scipy matplotlib
cd demos && for f in 0*.py; do python "$f"; done
```

## 重新產生 PDF／HTML

原稿是 `lecture/`、`review/` 下的 Markdown（GitHub 數學語法：行內 ``$`…`$``、獨立式 ```` ```math ````）。修改後執行：

```bash
cd tools
npm install
npx playwright install chromium   # 首次執行需要
node build_docs.js                # 輸出到 docs/；加 --no-pdf 只產生 HTML
```
