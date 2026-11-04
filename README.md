# Computational EM: from MPB to topological photonics and adjoint methods

- `lecture/lecture_notes_original.md` — 講義原稿（2026-10-29 版）
- `review/lecture_review.md` — 審閱報告：編排、正確性、引用、術語、排版、圖與 demo
- `demos/` — 最小可執行範例，兼作講義公式的數值驗證（PWE、FHS Chern 數、1D RCWA TE/TM、示意圖）
- `figures/` — demo 產生的圖

```bash
pip install numpy scipy matplotlib
cd demos && for f in 0*.py; do python "$f"; done
```
