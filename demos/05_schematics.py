"""講義缺少的示意圖：布里淵區與高對稱路徑、FHS 格子、RCWA 幾何與繞射階、S 矩陣串接、1D Yee 網格。"""
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle, Polygon

from common import plt, savefig


def arrow(ax, p, q, **kw):
    kw.setdefault("arrowstyle", "-|>")
    kw.setdefault("mutation_scale", 12)
    kw.setdefault("lw", 1.4)
    ax.add_patch(FancyArrowPatch(p, q, **kw))


# 1. 方晶格布里淵區 + Γ–X–M–Γ，以及 FHS 的 plaquette
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4))
a1.add_patch(Rectangle((-1, -1), 2, 2, fill=False, lw=1.5))
a1.add_patch(Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color="tab:orange", alpha=0.25,
                     label="irreducible BZ"))
for (p, q) in (((0, 0), (1, 0)), ((1, 0), (1, 1)), ((1, 1), (0, 0))):
    arrow(a1, p, q, color="tab:red")
for lab, (x, y) in {"Γ": (0, 0), "X": (1, 0), "M": (1, 1)}.items():
    a1.plot(x, y, "ko", ms=4)
    a1.annotate(lab, (x, y), textcoords="offset points", xytext=(-12, 6), fontsize=13)
a1.set_xlim(-1.3, 1.3); a1.set_ylim(-1.3, 1.3); a1.set_aspect("equal")
a1.set_xlabel("kx  (π/a)"); a1.set_ylabel("ky  (π/a)")
a1.set_title("1st Brillouin zone (square lattice)\npath Γ–X–M–Γ used in Lecture 2")
a1.legend(loc="lower left", fontsize=8)

L = 5
for i in range(L + 1):
    a2.plot([0, L], [i, i], color="gray", lw=0.6)
    a2.plot([i, i], [0, L], color="gray", lw=0.6)
x0, y0 = 2, 2
a2.add_patch(Rectangle((x0, y0), 1, 1, color="tab:blue", alpha=0.2))
arrow(a2, (x0, y0), (x0 + 1, y0), color="tab:red")
arrow(a2, (x0 + 1, y0), (x0 + 1, y0 + 1), color="tab:red")
arrow(a2, (x0 + 1, y0 + 1), (x0, y0 + 1), color="tab:red")
arrow(a2, (x0, y0 + 1), (x0, y0), color="tab:red")
a2.text(x0 + 0.5, y0 - 0.3, "U₁(k)", ha="center", fontsize=10)
a2.text(x0 + 1.1, y0 + 0.5, "U₂(k+1̂)", fontsize=10)
a2.text(x0 + 0.5, y0 + 1.12, "U₁(k+2̂)⁻¹", ha="center", fontsize=10)
a2.text(x0 - 0.1, y0 + 0.5, "U₂(k)⁻¹", ha="right", fontsize=10)
a2.text(x0 + 0.5, y0 + 0.5, "F₁₂(k)", ha="center", va="center", fontsize=10, color="tab:blue")
a2.text(0, -0.55, r"edges identified (torus): last-column links use $u_{k+G} = e^{-iG\cdot r}u_k$",
        fontsize=7.5)
a2.set_xlim(-0.5, L + 0.5); a2.set_ylim(-0.8, L + 0.5); a2.set_aspect("equal"); a2.axis("off")
a2.set_title("FHS lattice: Berry flux through one plaquette")
savefig(fig, "l01_l03_bz_and_fhs.png")

# 2. RCWA 幾何與繞射階
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.axhspan(0, 1, color="tab:blue", alpha=0.08)
ax.axhspan(-1.4, 0, color="white")
ax.axhspan(-2.6, -1.4, color="tab:green", alpha=0.10)
for x in np.arange(-3, 3, 1.0):
    ax.add_patch(Rectangle((x, -1.4), 0.5, 1.4, color="tab:gray", alpha=0.6))
ax.text(-3.4, 0.55, r"Region I  ($n_I$)"+"\nincident + reflected\nRayleigh expansion", fontsize=8)
ax.text(-3.4, -0.75, "grating layer, thickness d\nε(x) → [[ε]], modes W, Q", fontsize=8,
        bbox=dict(fc="white", ec="none", alpha=0.8))
ax.text(-3.4, -2.3, r"Region II  ($n_{II}$)"+"\ntransmitted orders", fontsize=8)
arrow(ax, (-0.9, 1.6), (-0.1, 0.05), color="k")
ax.text(-1.25, 1.5, "incident, θ", fontsize=9)
for mo, ang in ((-1, 35), (0, 0), (1, -35)):
    t = np.radians(ang)
    arrow(ax, (0, 0.05), (0.9 * np.sin(-t) * 1.4 + 0.6 * np.sin(np.radians(25)), 1.2), color="tab:blue")
    arrow(ax, (0, -1.4), (1.1 * np.sin(-t), -2.5), color="tab:green")
    ax.text(1.1 * np.sin(-t) + 0.05, -2.55, f"T{mo:+d}" if mo else "T0", fontsize=8, color="tab:green")
ax.text(1.2, 1.25, "R orders", fontsize=8, color="tab:blue")
ax.annotate("", xy=(2.5, -1.55), xytext=(3.5, -1.55), arrowprops=dict(arrowstyle="<->"))
ax.text(2.8, -1.8, "Λ", fontsize=10)
ax.text(1.6, 0.6, r"$k_{x,m} = k_0(n_I\sin\theta - m\lambda/\Lambda)$", fontsize=9)
ax.text(1.6, 0.3, r"$k_{z,m} = \sqrt{k_0^2n^2 - k_{x,m}^2}$", fontsize=9)
ax.set_xlim(-3.5, 3.6); ax.set_ylim(-2.7, 1.8); ax.axis("off")
ax.set_title("RCWA geometry (Lecture 6): Floquet orders in each half-space")
savefig(fig, "l06_rcwa_geometry.png")

# 3. S 矩陣串接（Redheffer 星積）
fig, ax = plt.subplots(figsize=(8, 2.8))
for i, (x, lab) in enumerate(((0, "S^A"), (3.2, "S^B"))):
    ax.add_patch(Rectangle((x, 0), 1.6, 1.6, fc="tab:blue", alpha=0.15, ec="tab:blue"))
    ax.text(x + 0.8, 0.8, f"${lab}$", ha="center", va="center", fontsize=14)
arrow(ax, (-1.2, 1.2), (0, 1.2)); ax.text(-1.2, 1.35, "$c_1^+$")
arrow(ax, (0, 0.4), (-1.2, 0.4)); ax.text(-1.2, 0.15, "$c_1^-$")
arrow(ax, (1.6, 1.2), (3.2, 1.2), color="tab:red"); arrow(ax, (3.2, 0.4), (1.6, 0.4), color="tab:red")
ax.text(1.75, 0.75, "multiple\nreflections\n$(I - S_{11}^B S_{22}^A)^{-1}$", fontsize=8, color="tab:red")
arrow(ax, (4.8, 1.2), (6, 1.2)); ax.text(5.6, 1.35, "$c_2^+$")
arrow(ax, (6, 0.4), (4.8, 0.4)); ax.text(5.6, 0.15, "$c_2^-$")
ax.text(6.4, 0.8, r"$S^{AB} = S^A \star S^B$" "\nonly $e^{-Qd}$ appears\n→ unconditionally stable",
        fontsize=10, va="center")
ax.set_xlim(-1.5, 9.5); ax.set_ylim(-0.1, 1.8); ax.axis("off")
ax.set_title("Lecture 8: cascading two blocks with the Redheffer star product")
savefig(fig, "l08_smatrix_cascade.png")

# 4. 1D Yee 交錯網格
fig, ax = plt.subplots(figsize=(8, 2.6))
for i in range(5):
    ax.plot(i, 1, "s", color="tab:blue", ms=9)
    ax.text(i, 1.18, f"$E_z^n({i})$" if i < 4 else "", ha="center", fontsize=8)
    if i < 4:
        ax.plot(i + 0.5, 0, "o", color="tab:red", ms=9)
        ax.text(i + 0.5, -0.3, f"$H_y^{{n+1/2}}({i}+\\frac{{1}}{{2}})$", ha="center", fontsize=8)
        arrow(ax, (i, 0.92), (i + 0.45, 0.08), color="gray", lw=0.8)
        arrow(ax, (i + 1, 0.92), (i + 0.55, 0.08), color="gray", lw=0.8)
ax.text(5.2, 1, "E at integer time steps n", color="tab:blue", va="center")
ax.text(5.2, 0, "H at half steps n+½\n(leapfrog)", color="tab:red", va="center")
ax.set_xlim(-0.4, 7.3); ax.set_ylim(-0.6, 1.5); ax.axis("off")
ax.set_title("Lecture 10: 1D Yee grid — E and H staggered by Δx/2 in space and Δt/2 in time")
savefig(fig, "l10_yee_1d.png")
