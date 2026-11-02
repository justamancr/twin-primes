"""Phase 1 — Brun's pure sieve on the twin-center wheel: watch the lower bound die.

S(x,z) = #{ n<=x : 6|n, n != ±1 (mod q) for all 5<=q<=z }
Inclusion-exclusion over squarefree d | P(z):  S = sum_d mu(d) * N_d
where N_d = #{n<=x, 6|n, n in one of 2^w(d) forbidden classes mod d} = (x/6) * 2^w(d)/d + err, |err| <= 2^w(d).
Truncating at w(d) <= r gives (Bonferroni) an upper bound for even r, lower bound for odd r.
Main term M_r(z) = (x/6) * sum_{j<=r} (-1)^j e_j(2/q),   error E_r(z) = sum_{j<=r} 2^j C(m, j),  m = #primes in [5,z].
Lower bound: max over odd r of M_r - E_r.   Upper bound: min over even r of M_r + E_r.
"""
import numpy as np, sys, time
from math import log, comb, exp

x = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8

# primes to sqrt(x)+ via simple sieve
lim = int(x**0.5) + 10
s = np.ones(lim + 1, dtype=bool); s[:2] = False
for i in range(2, int(lim**0.5) + 1):
    if s[i]: s[i*i::i] = False
primes = np.nonzero(s)[0]

def elem_sym(vals, rmax):
    """elementary symmetric polynomials e_0..e_rmax of vals"""
    e = np.zeros(rmax + 1); e[0] = 1.0
    for v in vals:
        e[1:] = e[1:] + v * e[:-1]       # e_j += v*e_{j-1}, descending order handled by vectorized shift
    return e

def brun_bounds(x, z, rmax=12):
    qs = primes[(primes >= 5) & (primes <= z)]
    m = len(qs)
    e = elem_sym(2.0 / qs, min(rmax, m))
    M = np.cumsum([(x / 6) * ((-1) ** j) * e[j] for j in range(len(e))])
    E = np.cumsum([(2 ** j) * comb(m, j) for j in range(len(e))])
    lower = max((M[r] - E[r], r) for r in range(1, len(e), 2))
    upper = min((M[r] + E[r], r) for r in range(0, len(e), 2))
    full = (x / 6) * np.prod(1 - 2.0 / qs)      # untruncated product = "heuristic"
    return lower, upper, full, m

def true_survivors(x, z):
    mm = np.arange(6, x + 1, 6, dtype=np.int64)
    alive = np.ones(len(mm), dtype=bool)
    for q in primes[(primes >= 5) & (primes <= z)]:
        r = mm % q
        alive &= (r != 1) & (r != q - 1)
    return int(alive.sum())

print(f"x = {x:.0e},  sqrt(x) = {x**0.5:.0f}")
print(f"{'z':>6} {'#q':>4} {'LOWER bound':>14} {'(r)':>4} {'UPPER bound':>14} {'(r)':>4} {'true S(x,z)':>12} {'product':>12} {'z = x^a, a=':>12} {'max #factors':>12}")
zs = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 61, 71, 83, 97, 113, 131, 151, 181, 211, 251, 307, 367, 439, 523, 631, 751, 907, 1091, 1301, 1559, 1871, 2243, 2693, 3229, 3877, 4651, 5581, 6701, 8039, 9649]
died = None
for z in zs:
    if z > x**0.5: break
    (lo, rlo), (up, rup), full, m = brun_bounds(x, z)
    t = true_survivors(x, z) if z <= 1600 else -1
    a = log(z) / log(x)
    print(f"{z:>6} {m:>4} {lo:>14.0f} {rlo:>4} {up:>14.0f} {rup:>4} {t if t>=0 else '  (skipped)':>12} {full:>12.0f} {a:>12.3f} {int(log(x)/log(z)):>12}")
    if lo <= 0 and died is None: died = z

print()
print(f"Lower bound first dies at z = {died}  (z = x^{log(died)/log(x):.3f}); survivors there have all prime factors > {died},")
print(f"so each of n-1, n+1 has at most {int(log(x)/log(died))} prime factors.  Certification needs z = sqrt(x) = x^0.5 -> at most 1 factor = prime.")
