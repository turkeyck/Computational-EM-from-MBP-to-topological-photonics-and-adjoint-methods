"""第 2 講習題 2：方晶格介電柱（ε = 12、R = 0.2a）的 2D TM 平面波展開。

TM（E_z）：|k+G|^2 e_G = (ω/c)^2 Σ_G' ε_{G-G'} e_G'  →  K e = λ [[ε]] e
用 scipy.linalg.eigh(K, Eps) 解廣義 Hermitian 本徵問題。
對照 MPB 官方 Python 教學：TM 能隙 0.2826 → 0.4193（ωa/2πc），相對寬度約 39%。
"""
import numpy as np
from scipy.linalg import eigh
from scipy.special import j1

from common import plt, savefig

EPS_ROD, EPS_BG, R = 12.0, 1.0, 0.2   # a = 1
NBANDS = 6


def g_set(n_max):
    """|G| 排序截斷的倒晶格向量（方晶格 b1 = 2π x̂, b2 = 2π ŷ）。"""
    m = np.arange(-n_max, n_max + 1)
    mx, my = np.meshgrid(m, m, indexing="ij")
    idx = np.stack([mx.ravel(), my.ravel()], axis=1)
    r2 = (idx ** 2).sum(1)
    return idx[r2 <= n_max ** 2]


def eps_matrix(idx):
    """圓柱的解析傅立葉係數 ε_{G-G'}（第 2.3 節公式）。"""
    d = idx[:, None, :] - idx[None, :, :]
    G = 2 * np.pi * np.hypot(d[..., 0], d[..., 1])
    f = np.pi * R ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        form = np.where(G == 0, 1.0, 2 * j1(G * R) / (G * R))
    return EPS_BG * (G == 0) + (EPS_ROD - EPS_BG) * f * form


def bands_tm(kpts, idx, Eps):
    out = []
    Gvec = 2 * np.pi * idx
    for k in kpts:
        K = np.diag(((k + Gvec) ** 2).sum(1))
        w2 = eigh(K, Eps, eigvals_only=True, subset_by_index=[0, NBANDS - 1])
        out.append(np.sqrt(np.abs(w2)) / (2 * np.pi))  # ωa/2πc
    return np.array(out)


def kpath(n_seg=20):
    G, X, M = np.array([0, 0]), np.array([np.pi, 0]), np.array([np.pi, np.pi])
    pts, ticks = [], [0]
    for A, B in [(G, X), (X, M), (M, G)]:
        for t in np.linspace(0, 1, n_seg, endpoint=False):
            pts.append(A + t * (B - A))
        ticks.append(len(pts))
    pts.append(G)
    return np.array(pts), ticks


if __name__ == "__main__":
    kpts, ticks = kpath()
    print(" n_max   N_PW   band1 max   band2 min   gap %")
    for n_max in (5, 7, 9, 11):
        idx = g_set(n_max)
        b = bands_tm(kpts, idx, eps_matrix(idx))
        lo, hi = b[:, 0].max(), b[:, 1].min()
        print(f"{n_max:5d} {len(idx):6d}   {lo:.4f}      {hi:.4f}     {200*(hi-lo)/(hi+lo):.1f}")
    print("MPB tutorial:      0.2826      0.4193      39.0")

    fig, ax = plt.subplots(figsize=(5.2, 4))
    x = np.arange(len(kpts))
    ax.axhspan(lo, hi, color="tab:orange", alpha=0.2, label=f"TM gap {lo:.3f}–{hi:.3f}")
    ax.plot(x, b, color="tab:blue")
    ax.set_xticks(ticks, ["Γ", "X", "M", "Γ"])
    ax.set_xlim(0, x[-1])
    ax.set_ylabel("ωa / 2πc")
    ax.set_title(f"Square lattice, ε=12 rods, R=0.2a, TM ({len(idx)} plane waves)")
    ax.legend(loc="upper left")
    savefig(fig, "l02_pwe_tm_bands.png")
