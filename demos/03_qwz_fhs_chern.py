"""第 3 講習題 2、4：用 Fukui–Hatsugai–Suzuki (FHS) 方法算 QWZ 模型下能帶的 Chern 數。

H(k) = sin kx σx + sin ky σy + (m + cos kx + cos ky) σz
FHS：U_μ(k) = <u(k)|u(k+μ)>/|...|，F12 = Im ln[U1(k) U2(k+1) U1(k+2)^-1 U2(k)^-1]，C = ΣF12/2π
本徵解法給的任意相位在乘積中抵消，所以不必做任何規範固定。
"""
import numpy as np

from common import plt, savefig

sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)


def lower_band_states(m, L):
    ks = 2 * np.pi * np.arange(L) / L
    kx, ky = np.meshgrid(ks, ks, indexing="ij")
    H = (np.sin(kx)[..., None, None] * sx + np.sin(ky)[..., None, None] * sy
         + (m + np.cos(kx) + np.cos(ky))[..., None, None] * sz)
    _, v = np.linalg.eigh(H)
    return v[..., :, 0], kx, ky           # 下能帶 (L, L, 2)


def fhs(m, L=24, return_field=False):
    u, kx, ky = lower_band_states(m, L)
    def link(a, b):
        z = np.sum(a.conj() * b, axis=-1)
        return z / np.abs(z)
    u1 = np.roll(u, -1, axis=0)           # k + x̂（週期邊界：BZ 為環面）
    u2 = np.roll(u, -1, axis=1)
    U1, U2 = link(u, u1), link(u, u2)
    F = np.angle(U1 * np.roll(U2, -1, axis=0) / np.roll(U1, -1, axis=1) / U2)
    C = F.sum() / (2 * np.pi)
    return (C, F, kx, ky) if return_field else C


if __name__ == "__main__":
    ms = np.linspace(-3, 3, 241)
    Cs = np.array([fhs(m) for m in ms])
    for m in (-2.5, -1, 1, 2.5):
        print(f"m = {m:+.1f}:  C = {fhs(m):+.6f}")
    print("grid-size test at m = 1:", {L: round(fhs(1.0, L), 10) for L in (4, 6, 12, 48)})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    ax1.plot(ms, Cs, lw=2)
    for x in (-2, 0, 2):
        ax1.axvline(x, color="gray", ls="--", lw=0.8)
    ax1.set_xlabel("m"); ax1.set_ylabel("Chern number (lower band)")
    ax1.set_title("QWZ phase diagram via FHS (24×24 grid)")
    ax1.set_yticks([-1, 0, 1])

    C, F, kx, ky = fhs(1.0, L=48, return_field=True)
    im = ax2.pcolormesh(kx / np.pi, ky / np.pi, F / (2 * np.pi / 48) ** 2, shading="auto", cmap="RdBu_r")
    fig.colorbar(im, ax=ax2, label="Berry curvature Ω (F12 / plaquette area)")
    ax2.set_xlabel("kx / π"); ax2.set_ylabel("ky / π")
    ax2.set_title(f"Ω(k) at m = 1, C = {C:+.0f}")
    ax2.set_aspect("equal"); ax2.grid(False)
    savefig(fig, "l03_qwz_fhs.png")
