"""Phase 4a — hunt for structure in twin centers beyond the wheel's prediction. x = 1e8."""
import numpy as np, sys
from math import log, exp, sqrt
from scipy.stats import chi2

x = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
s = np.ones(x // 2 + 2, dtype=bool); s[0] = False
for i in range(1, int(x**0.5)//2 + 1):
    if s[i]:
        p = 2*i+1; s[p*p//2::p] = False
odd_primes = 2*np.nonzero(s)[0] + 1
primes = np.concatenate(([2], odd_primes))
isp = lambda n: s[n//2]   # odd n only

def pairs_at(h):
    p = odd_primes[odd_primes + h <= x]
    return p[isp(p + h)]
tw = pairs_at(2); centers = tw[tw > 1000] + 1
print(f"x={x:.0e}: twins={len(tw)}")

# ---------- A. residue equidistribution of twin centers mod q (allowed classes only) ----------
print("\nA. twin centers mod q: are the q-2 allowed residues equally used?  (chi2 p-value; random model => uniform p-values)")
print(f"{'q':>4} {'chi2':>9} {'dof':>4} {'p':>7}   top deviation")
for q in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]:
    r = centers % q
    counts = np.bincount(r, minlength=q).astype(float)
    assert counts[1] == 0 and counts[q-1] == 0
    allowed = np.array([c for c in range(q) if c not in (1, q-1)])
    obs = counts[allowed]; expc = obs.sum()/len(allowed)
    c2 = ((obs-expc)**2/expc).sum(); dof = len(allowed)-1
    pval = chi2.sf(c2, dof)
    dev = allowed[np.argmax(np.abs(obs-expc))]
    print(f"{q:>4} {c2:>9.1f} {dof:>4} {pval:>7.3f}   residue {dev}: {(obs.max()-expc)/expc*100:+.2f}%")

# ---------- B. joint independence of residues mod q1, q2 ----------
print("\nB. independence of (n mod q1, n mod q2) over allowed classes (chi2 p-value)")
qs = [5,7,11,13,17,19,23]
for i in range(len(qs)):
    for j in range(i+1, len(qs)):
        q1, q2 = qs[i], qs[j]
        a1 = [c for c in range(q1) if c not in (1,q1-1)]; a2 = [c for c in range(q2) if c not in (1,q2-1)]
        r1 = centers % q1; r2 = centers % q2
        T = np.zeros((q1, q2)); np.add.at(T, (r1, r2), 1)
        T = T[np.ix_(a1, a2)]
        E = np.outer(T.sum(1), T.sum(0)) / T.sum()
        c2 = ((T-E)**2/E).sum(); dof = (len(a1)-1)*(len(a2)-1)
        print(f"  ({q1:>2},{q2:>2}) chi2={c2:7.1f} dof={dof:3d} p={chi2.sf(c2,dof):.3f}", end="")
    print()

# ---------- C. how often is a twin CENTER divisible by q?  wheel predicts 1/(q-2), random numbers 1/q ----------
print("\nC. P(q | center): wheel prediction 1/(q-2) vs naive 1/q")
print(f"{'q':>4} {'observed':>9} {'1/(q-2)':>9} {'1/q':>7}")
for q in [5,7,11,13,17,19,23,29,31,37,41,43,47]:
    obs = (centers % q == 0).mean()
    print(f"{q:>4} {obs:>9.5f} {1/(q-2):>9.5f} {1/q:>7.5f}")

# ---------- D. 2-adic valuation of centers: random model P(v2 = j | even) = 2^-j ----------
print("\nD. 2-adic valuation v2(center): observed fraction vs 2^-j")
c = centers.copy(); v = np.zeros(len(c), dtype=int)
while True:
    m = (c % 2 == 0)
    if not m.any(): break
    v[m] += 1; c[m] //= 2
for j in range(1, 13):
    print(f"  v2={j:>2}: {(v==j).mean():.5f}   model {2.0**-j:.5f}")

# ---------- E. prime pairs at distance h: HL predicts count(h)/count(2) = prod_{p|h, p>2} (p-1)/(p-2) ----------
print("\nE. pairs (p, p+h) below x; ratio to h=2, vs Hardy–Littlewood singular series")
def ss(h):
    r = 1.0
    for p in primes[1:60]:
        if h % p == 0: r *= (p-1)/(p-2)
    return r
for h in [2,4,6,8,10,12,14,30,210,246]:
    n = len(pairs_at(h))
    print(f"  h={h:>3}: {n:>8}  ratio {n/len(tw):.4f}  HL {ss(h):.4f}")

# ---------- F. two-point Möbius sum ----------
print("\nF. sum_{n<=N} mu(n) mu(n+2)  (Chowla predicts o(N); random-walk scale sqrt(N))")
N = min(x, 10**8)
mu = np.ones(N + 3, dtype=np.int8)
for p in primes:
    if p*p > N+2: break
    mu[p*p::p*p] = 0
# count prime factors parity via sieve of omega
om = np.zeros(N + 3, dtype=np.int8)
for p in primes:
    om[p::p] += 1
# for squarefree n, mu = (-1)^omega
mu = np.where(mu == 0, 0, np.where(om % 2 == 0, 1, -1)).astype(np.int8)
mu[0] = 0
prod = mu[1:N+1].astype(np.int32) * mu[3:N+3].astype(np.int32)
cs = np.cumsum(prod, dtype=np.int64)
for e in range(4, int(log(N,10))+1):
    M = 10**e
    print(f"  N=1e{e}: sum={cs[M-1]:>8}   sum/sqrt(N)={cs[M-1]/sqrt(M):+.3f}   sum/N={cs[M-1]/M:+.2e}")
# logarithmic average
lw = (prod[1:] / np.arange(2, N+1)).cumsum()
print(f"  log-averaged: (1/log N) sum mu(n)mu(n+2)/n at N=1e8 = {lw[-1]/log(N):+.5f}  (Tao: -> 0)")
