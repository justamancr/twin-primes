import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
G = "#9a9a9a"; B = "#1f5fbf"

# ---- Fig 1: Brun bounds
z = [7,11,13,17,19,23,29,31,37,41,47,53,61,71,83,97,113,131,151,181,211,251,307,367,439,523]
lo = [5238090,5844129,4891710,4363041,3902666,3562326,3312376,3087690,2905235,2743966,2442811,2304218,2020472,1695903,1041338,435333,-2131633,-3769472,-4926831,-6114841,-7036011,-8216825,-9631880,-11076450,-12713685,-14580026]
up = [7142866,6190495,4945136,4369763,3904720,3566630,3324207,3117475,2964063,2846900,2720469,2736704,2927265,2977532,3122832,3256124,3720191,3952288,4513568,5644482,6923916,9377684,14067812,16666668,16666668,16666668]
tr = [7142856,5844155,4945055,4363285,3903985,3564504,3318662,3104557,2936752,2793489,2550231,2454015,2293050,2161774,1999614,1913995,1740117,1686516,1594734,1486054,1413941,1334218,1253598,1184969,1119875,1057485]
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.axhline(0, color="#555", lw=0.8)
ax.plot(z, np.array(up)/1e6, "--", color=G, lw=1, marker="o", ms=2.5, label="Brun upper bound")
ax.plot(z, np.array(tr)/1e6, "-", color=G, lw=1.2, marker="o", ms=2.5, label="true survivor count")
ax.plot(z, np.array(lo)/1e6, "-", color=B, lw=1.8, marker="o", ms=2.5, label="Brun lower bound")
ax.set_xscale("log"); ax.set_xlabel("sieving depth $z$ (all primes $\\leq z$ applied)"); ax.set_ylabel("twin centers below $10^8$ (millions)")
ax.set_xticks([7,17,29,41,61,97,151,251,439]); ax.set_xticklabels([7,17,29,41,61,97,151,251,439])
ax.annotate("lower bound $<0$ from $z=113=x^{0.257}$", xy=(113, -2.13), xytext=(12, -8), color=B, fontsize=8, arrowprops=dict(arrowstyle="-", color=B, lw=0.6))
ax.legend(frameon=False, loc="upper left", fontsize=8); fig.tight_layout()
fig.savefig("fig1_brun.pdf"); fig.savefig("fig1_brun.png")

# ---- Fig 2: GPY ratio
a = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]; pair=[0.255,0.333,0.409,0.482,0.552,0.622,0.691,0.760]
a6=[0.15,0.2,0.25,0.3]; t6=[0.429,0.489,0.558,0.633]
fig, ax = plt.subplots(figsize=(6.2, 3.4))
ax.axhline(1, color="#555", lw=0.8, ls="--"); ax.text(0.555, 0.975, "$\\rho=1$: a second prime is certified", fontsize=8, ha="right", va="top", color="#555")
for xv, lab in [(0.25, "provable (Bombieri–Vinogradov)"), (0.5, "Elliott–Halberstam")]:
    ax.axvline(xv, color="#aaa", lw=0.7, ls=":"); ax.text(xv-0.004, 1.085, lab, rotation=0, fontsize=7.5, color="#666", va="top", ha="right")
ax.plot(a6, t6, "-o", color=G, ms=3, lw=1.2, label="6-tuple, diameter 16")
ax.plot(a, pair, "-o", color=B, ms=3, lw=1.8, label="pair $(n, n+2)$")
for xv, yv in [(0.25, 0.409), (0.5, 0.760)]: ax.text(xv-0.006, yv+0.03, f"{yv}", ha="right", fontsize=8)
ax.set_xlim(0.1, 0.56); ax.set_ylim(0.1, 1.12); ax.set_xlabel("exponent $a$ in the sieve level $R = x^a$"); ax.set_ylabel("$\\rho = S_1/(S_0 \\log x)$")
ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0,0.02), fontsize=8); fig.tight_layout()
fig.savefig("fig2_gpy.pdf"); fig.savefig("fig2_gpy.png")

# ---- Fig 3: rigidity primes vs twins (x=1e9) + HL predictions
q = [31,37,41,47,53,61,71,83,97,113,131,151,179,211,251,307,401,503,701,1009]
pr = [0.050,0.076,0.105,0.089,0.081,0.090,0.130,0.102,0.117,0.123,0.108,0.146,0.130,0.151,0.161,0.152,0.193,0.201,0.210,0.213]
tw = [0.694,0.522,0.812,0.571,0.446,0.721,0.493,0.482,0.513,0.485,0.542,0.649,0.572,0.659,0.586,0.653,0.658,0.647,0.610,0.629]
qp = [31,61,113,251,503,1009]; pred9=[0.486,0.523,0.554,0.593,0.625,0.655]
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.axhline(1, color="#555", lw=0.8, ls="--"); ax.text(32, 1.02, "Poisson (random sample)", fontsize=8, color="#555")
ax.axhline(0.933, color="#999", lw=0.8, ls="--"); ax.text(32, 0.945, "primes thinned at random to twin density", fontsize=8, color="#777")
ax.axhline(0.5, color=B, lw=0.8, ls=":"); ax.text(32, 0.36, "limit 1/2 as $x\\to\\infty$ (conditional on HL)", fontsize=8, color=B)
ax.plot(q, tw, "-o", color=B, ms=3, lw=1.4, label="twin centers, measured ($x=10^9$)")
ax.plot(qp, pred9, "s", color="k", ms=5, mfc="none", label="twin centers, HL second-order prediction")
ax.plot(q, pr, "-o", color=G, ms=3, lw=1.2, label="primes, measured ($x=10^9$)")
ax.set_xscale("log"); ax.set_xticks([31,61,113,251,503,1009]); ax.set_xticklabels([31,61,113,251,503,1009])
ax.set_xlabel("modulus $q$"); ax.set_ylabel("$\\chi^2$ per degree of freedom"); ax.set_ylim(0, 1.1)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5,-0.18), ncol=3, fontsize=7.5); fig.tight_layout()
fig.savefig("fig3_rigidity.pdf", bbox_inches="tight"); fig.savefig("fig3_rigidity.png", bbox_inches="tight")
print("figures written")
