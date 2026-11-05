"""第 1 講習題 3 + 第 2 講習題 1、3、4：一維 PWE 對 Rytov 解析解，比較 Laurent 規則與逆規則。

一維 Bragg 堆疊，ε 只隨 x 變化，介面法向為 x̂，正入射（k 沿 x）：
  (a) E_y 形式：-E'' = (ω/c)^2 ε E        → E_y 為切向、連續，用 Laurent [[ε]]（Li 定理 1）
  (b) H_z 形式：-(η H')' = (ω/c)^2 H      → η 作用在 D_y = ε E_y 上得到 E_y
       D_y 與 1/ε 同時跳躍、乘積 E_y 連續 → 逆規則 η = [[ε]]^{-1}（Li 定理 2）
       Laurent η = [[1/ε]] 收斂較慢。
注意：這裡關鍵是「切向」E_y 連續，而不是 D_x 連續（原講義第 2 講習題 4 的敘述需修正）。
"""
import numpy as np
from scipy.linalg import eigh, toeplitz
from scipy.optimize import brentq

from common import plt, savefig

n1, n2 = 3.5, 1.5           # a = Λ = 1
d1 = 0.3; d2 = 1 - d1
eps1, eps2 = n1 ** 2, n2 ** 2
K_BLOCH = 0.3 * 2 * np.pi    # 取一個一般的 Bloch 波數
NB = 4


def fourier_step(val_in, val_out, width, nmax):
    """週期 1、x∈[0,width) 取 val_in，其餘 val_out 的傅立葉係數 c_n（n = -nmax..nmax）。"""
    n = np.arange(-nmax, nmax + 1)
    c = np.empty(n.size, complex)
    c[n == 0] = val_in * width + val_out * (1 - width)
    nz = n != 0
    c[nz] = (val_in - val_out) * (1 - np.exp(-2j * np.pi * n[nz] * width)) / (2j * np.pi * n[nz])
    return c


def conv_matrix(c, M):
    """由係數 c_{-2M..2M} 組 (2M+1)×(2M+1) 的 Toeplitz 矩陣 [[f]]_{mn} = c_{m-n}。"""
    mid = len(c) // 2
    col = c[mid: mid + 2 * M + 1]
    row = c[mid - 2 * M: mid + 1][::-1]
    return toeplitz(col, row)


def rytov_freqs(k, nb=NB):
    """解 cos(kΛ) = cos(k1 d1)cos(k2 d2) - ½(n1/n2+n2/n1) sin sin，回傳前 nb 個 ωΛ/2πc。"""
    def f(w):  # w = ωΛ/2πc，k_i = n_i ω/c = 2π n_i w
        a, b = 2 * np.pi * n1 * w * d1, 2 * np.pi * n2 * w * d2
        return np.cos(a) * np.cos(b) - 0.5 * (n1 / n2 + n2 / n1) * np.sin(a) * np.sin(b) - np.cos(k)
    ws = np.linspace(1e-6, 1.2, 20001)
    fs = f(ws)
    roots = [brentq(f, ws[i], ws[i + 1]) for i in range(len(ws) - 1) if fs[i] * fs[i + 1] < 0]
    return np.array(roots[:nb])


def pwe(M, form):
    Gs = 2 * np.pi * np.arange(-M, M + 1)
    kG = np.diag(K_BLOCH + Gs)
    if form == "E":
        Eps = conv_matrix(fourier_step(eps1, eps2, d1, 2 * M), M)
        w2 = eigh(kG @ kG, Eps, eigvals_only=True)
    else:
        if form == "H-laurent":
            eta = conv_matrix(fourier_step(1 / eps1, 1 / eps2, d1, 2 * M), M)
        else:  # H-inverse
            eta = np.linalg.inv(conv_matrix(fourier_step(eps1, eps2, d1, 2 * M), M))
        w2 = np.linalg.eigvalsh(kG @ eta @ kG)
    return np.sqrt(np.abs(np.sort(w2)[:NB])) / (2 * np.pi)


if __name__ == "__main__":
    exact = rytov_freqs(K_BLOCH)
    print("Rytov exact bands:", np.round(exact, 6))
    Ms = np.unique(np.logspace(0.5, 2.6, 18).astype(int))
    errs = {f: [] for f in ("E", "H-laurent", "H-inverse")}
    for M in Ms:
        for f in errs:
            errs[f].append(abs(pwe(M, f)[1] - exact[1]) / exact[1])  # 第 2 條能帶
    labels = {"E": "E_y form, Laurent [[ε]]  (correct)",
              "H-laurent": "H_z form, Laurent [[1/ε]]",
              "H-inverse": "H_z form, inverse [[ε]]⁻¹  (correct)"}
    fig, ax = plt.subplots(figsize=(5.4, 4))
    for f, e in errs.items():
        e = np.array(e)
        p = np.polyfit(np.log(Ms[-8:]), np.log(e[-8:]), 1)[0]
        ax.loglog(2 * Ms + 1, e, "o-", ms=3, label=f"{labels[f]}, slope≈{p:.1f}")
        print(f"{f:10s} slope {p:.2f}  final rel.err {e[-1]:.2e}")
    ax.set_xlabel("number of plane waves N = 2M+1")
    ax.set_ylabel("relative error of band 2")
    ax.set_title(f"1D Bragg stack n=({n1},{n2}), k=0.3·2π/Λ")
    ax.legend(fontsize=7.5)
    savefig(fig, "l02_pwe_1d_laurent_vs_inverse.png")
