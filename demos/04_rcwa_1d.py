"""第 6、7 講：單層一維二元光柵的 RCWA（TE 與 TM），全部採用 e^{-iωt} 慣例。

慣例（講義應明確統一）：時諧 e^{-iωt}；橫向 e^{+i k_{x,m} x}；
沿 +z 傳播或衰減的模態為 e^{+i k_z z}，分支取 Im k_z > 0（傳播模取 Re k_z > 0）。

幾何：z<0 入射區 n_I；0<z<d 光柵層（寬 f·Λ 的 ε_r 脊，其餘 ε_g）；z>d 基板 n_II。
層內 k_z^2 為下列矩陣的本徵值（以 k0 無因次化）：
  TE：  E - Kx^2
  TM：  Li／Lalanne–Morris：[[1/ε]]^{-1} (I - Kx [[ε]]^{-1} Kx)
        舊式（Moharam 1995）：[[ε]] (I - Kx [[ε]]^{-1} Kx)
匹配 z=0、z=d 的切向場（TE：E_y、∂z E_y；TM：H_y、(1/ε)∂z H_y）。
層內場寫成 W[e^{iKz z'} c+ + e^{-iKz (z'-d')} c-]，只出現衰減指數。
"""
import numpy as np
from scipy.linalg import toeplitz

from common import plt, savefig


def fourier_coeffs(eps_ridge, eps_groove, fill, nmax):
    n = np.arange(-nmax, nmax + 1)
    c = np.empty(n.size, complex)
    c[n == 0] = eps_ridge * fill + eps_groove * (1 - fill)
    nz = n != 0
    c[nz] = (eps_ridge - eps_groove) * np.sin(np.pi * n[nz] * fill) / (np.pi * n[nz])  # 脊置中
    return c


def conv(c, M):
    mid = len(c) // 2
    return toeplitz(c[mid: mid + 2 * M + 1], c[mid - 2 * M: mid + 1][::-1])


def kz_branch(kz2):
    kz = np.sqrt(kz2.astype(complex))
    return np.where(kz.imag < 0, -kz, kz)


def rcwa_1d(wl, period, d, n_I, n_II, eps_ridge, eps_groove, fill, M,
            theta=0.0, pol="TE", li_rule=True):
    k0 = 2 * np.pi / wl
    m = np.arange(-M, M + 1)
    kx = n_I * np.sin(theta) - m * wl / period            # 已除以 k0
    Kx = np.diag(kx)
    I = np.eye(len(m))
    kzI = kz_branch(n_I ** 2 - kx ** 2)
    kzII = kz_branch(n_II ** 2 - kx ** 2)

    Ec = conv(fourier_coeffs(eps_ridge, eps_groove, fill, 2 * M), M)
    if pol == "TE":
        A = Ec - Kx @ Kx
        P = I                      # ∂z E_y 直接連續
        pI = pII = 1.0
    else:
        Einv = np.linalg.inv(Ec)
        B = I - Kx @ Einv @ Kx
        if li_rule:
            Ac = conv(fourier_coeffs(1 / eps_ridge, 1 / eps_groove, fill, 2 * M), M)
            A = np.linalg.solve(Ac, B)
            P = Ac                 # E_x = [[1/ε]] D_x （Li 定理 1）
        else:
            A = Ec @ B
            P = Einv
        pI, pII = 1 / n_I ** 2, 1 / n_II ** 2

    lam, W = np.linalg.eig(A)
    kz = kz_branch(lam)
    X = np.diag(np.exp(1j * kz * k0 * d))
    V = P @ W @ np.diag(kz)        # 層內「∂z 場」的係數
    N = len(m)
    delta = (m == 0).astype(complex)
    Z = np.zeros((N, N), complex)
    # 未知量：[R, T, c+, c-]
    Amat = np.block([
        [I, Z, -W, -W @ X],
        [pI * np.diag(kzI) * -1, Z, -V, V @ X],
        [Z, -I, W @ X, W],
        [Z, -pII * np.diag(kzII), V @ X, -V],
    ])
    rhs = np.concatenate([-delta, -pI * kzI * delta, np.zeros(N), np.zeros(N)])
    # 第二列：pI kzI (δ - R) = V (c+ - X c-)  →  -pI kzI R - V c+ + V X c- = -pI kzI δ
    sol = np.linalg.solve(Amat, rhs)
    R, T = sol[:N], sol[N:2 * N]
    c0 = n_I * np.cos(theta)
    DEr = np.abs(R) ** 2 * np.real(kzI) / c0
    DEt = np.abs(T) ** 2 * np.real(kzII * pII) / (c0 * pI)
    return m, DEr, DEt


def thin_film_R(wl, d, n1, n2, n3, theta=0.0, pol="TE"):
    """Airy 公式：均勻薄膜的反射率（解析對照）。"""
    s1 = n1 * np.sin(theta)
    c = [np.sqrt(complex(n ** 2 - s1 ** 2)) for n in (n1, n2, n3)]   # n cosθ
    if pol == "TE":
        r = lambda i, j: (c[i] - c[j]) / (c[i] + c[j])
    else:
        ns = (n1, n2, n3)
        r = lambda i, j: (ns[j] ** 2 * c[i] - ns[i] ** 2 * c[j]) / (ns[j] ** 2 * c[i] + ns[i] ** 2 * c[j])
    ph = np.exp(2j * 2 * np.pi / wl * c[1] * d)
    rr = (r(0, 1) + r(1, 2) * ph) / (1 + r(0, 1) * r(1, 2) * ph)
    return abs(rr) ** 2


if __name__ == "__main__":
    geo = dict(period=1.0, d=0.5, n_I=1.0, n_II=1.5, eps_ridge=4.0, eps_groove=1.0, fill=0.5)

    print("(a) energy conservation, lossless grating, θ = 20°")
    for pol in ("TE", "TM"):
        _, r, t = rcwa_1d(1.3, M=20, theta=np.radians(20), pol=pol, **geo)
        print(f"   {pol}: ΣDE_r + ΣDE_t - 1 = {r.sum() + t.sum() - 1:+.2e}")

    print("(b) thin-film limit (ridge = groove = 2.25), θ = 30°")
    for pol in ("TE", "TM"):
        _, r, _ = rcwa_1d(0.9, 1.0, 0.37, 1.0, 1.5, 2.25, 2.25, 0.5, M=5, theta=np.radians(30), pol=pol)
        print(f"   {pol}: RCWA R0 = {r[5]:.12f}   Airy = {thin_film_R(0.9, 0.37, 1.0, 1.5, 1.5, np.radians(30), pol):.12f}")

    # (c) 波長掃描與 Rayleigh anomaly
    wls = np.linspace(0.8, 2.0, 1201)
    res = np.array([rcwa_1d(w, M=15, pol="TE", **geo)[1:] for w in wls])  # (n, 2, N)
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.plot(wls, res[:, 0, 15], label="R₀")
    ax.plot(wls, res[:, 1, 15], label="T₀")
    ax.plot(wls, res[:, 1, 14] + res[:, 1, 16], label="T₊₁ + T₋₁")
    for wl_ra, lab in ((1.0, r"RA: $\lambda = n_I\Lambda$"), (1.5, r"RA: $\lambda = n_{II}\Lambda$")):
        ax.axvline(wl_ra, color="gray", ls="--", lw=0.8)
        ax.text(wl_ra + 0.01, 0.92, lab, fontsize=8)
    ax.set_xlabel("λ / Λ"); ax.set_ylabel("diffraction efficiency")
    ax.set_title("1D binary grating, TE, normal incidence")
    ax.legend(loc="center right")
    savefig(fig, "l06_rcwa_te_spectrum.png")

    # (d) 第 7 講習題 2：金屬光柵 TM 收斂，舊公式 vs Li 法則
    metal = dict(period=1.0, d=0.3, n_I=1.0, n_II=1.0, eps_ridge=(0.2 + 3j) ** 2, eps_groove=1.0, fill=0.5)
    Ms = np.arange(2, 81, 3)
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    for li, lab in ((False, "old: [[ε]] (Moharam 1995)"), (True, "Li 1996: [[1/ε]]⁻¹")):
        R0 = [rcwa_1d(0.8, M=M, theta=np.radians(10), pol="TM", li_rule=li, **metal)[1][M] for M in Ms]
        ax.plot(2 * Ms + 1, R0, "o-", ms=3, label=lab)
        print(f"   TM metal, {lab:28s}: R0(M=20)={R0[6]:.5f}  R0(M=80)={R0[-1]:.5f}")
    ax.set_xlabel("number of harmonics 2M+1"); ax.set_ylabel("R₀ (TM)")
    ax.set_title("Metallic lamellar grating, n = 0.2+3i, λ = 0.8Λ")
    ax.legend()
    savefig(fig, "l07_tm_metal_convergence.png")
