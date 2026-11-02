"""Phase 4b-ii — bilinear / multiplicative-structure probe at x=1e8.
Control group: wheel survivors at depth z (centers n, 6|n, n != ±1 mod q for 5<=q<=z). Twins are a subset.
Null (Type I only): twin indicator is independent of multiplicative functions of n, given survival.
"""
import numpy as np, sys
from math import log, sqrt
x = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
s = np.ones(x//2+2, dtype=bool); s[0]=False
for i in range(1, int(x**0.5)//2+1):
    if s[i]:
        p=2*i+1; s[p*p//2::p]=False
odd_primes = 2*np.nonzero(s)[0]+1
primes = np.concatenate(([2], odd_primes))

# multiplicative functions on all n <= x : omega, Omega parity (Liouville), squarefree
M = x+2
Omega = np.zeros(M, dtype=np.int8); omega = np.zeros(M, dtype=np.int8); sqf = np.ones(M, dtype=bool)
for p in primes:
    omega[p::p] += 1
    pk = p
    while pk < M:
        Omega[pk::pk] += 1; pk *= p
    if p*p < M: sqf[p*p::p*p] = False
liou = np.where(Omega % 2 == 0, 1, -1).astype(np.int8)
mu = np.where(sqf, liou, 0).astype(np.int8)
print(f"x={x:.0e}: multiplicative tables built")

centers_all = np.arange(6, x+1, 6, dtype=np.int64)
isp = lambda n: s[n//2]
twin = isp(centers_all-1) & isp(centers_all+1)
twin &= centers_all > 1000

def survivors(z):
    alive = np.ones(len(centers_all), dtype=bool)
    for q in primes[(primes>=5)&(primes<=z)]:
        r = centers_all % q; alive &= (r!=1)&(r!=q-1)
    return alive

print("\n(i) mean of multiplicative functions: twins vs depth-z survivors (z-score of difference)")
print(f"{'z':>6} {'#surv':>9} {'#twin':>7} | {'<lambda> surv':>13} {'twin':>8} {'zsc':>6} | {'<mu> surv':>10} {'twin':>8} {'zsc':>6} | {'<omega> surv':>12} {'twin':>7} {'zsc':>6} | {'P(n/6 prime) surv':>17} {'twin':>7} {'zsc':>6}")
for z in [5, 31, 101, 331, 1009, 3163, 10007]:
    al = survivors(z); tw = al & twin
    S = centers_all[al]; T = centers_all[tw]
    out = [z, len(S), len(T)]
    for f in (liou, mu, omega):
        a = f[S].astype(float); b = f[T].astype(float)
        zsc = (b.mean()-a.mean())/sqrt(a.var()/len(b)) if len(b) else 0
        out += [a.mean(), b.mean(), zsc]
    a = isp(S//6 | 1) & (S//6 % 2 == 1); b = isp(T//6 | 1) & (T//6 % 2 == 1)   # n/6 odd prime
    zsc = (b.mean()-a.mean())/sqrt(a.mean()*(1-a.mean())/len(b))
    out += [a.mean(), b.mean(), zsc]
    print(f"{out[0]:>6} {out[1]:>9} {out[2]:>7} | {out[3]:>13.4f} {out[4]:>8.4f} {out[5]:>6.2f} | {out[6]:>10.4f} {out[7]:>8.4f} {out[8]:>6.2f} | {out[9]:>12.4f} {out[10]:>7.4f} {out[11]:>6.2f} | {out[12]:>17.5f} {out[13]:>7.5f} {out[14]:>6.2f}")

print("\n(i-b) KILL TEST: same comparison restricted to windows [x/2, x] and [x/4, x/2]  (size bias removed)")
for lo, hi in [(x//2, x), (x//4, x//2)]:
    win = (centers_all > lo) & (centers_all <= hi)
    for z in [331, 1009, 3163]:
        al = survivors(z) & win; tw = al & twin
        S = centers_all[al]; T = centers_all[tw]
        a = omega[S].astype(float); b = omega[T].astype(float)
        zo = (b.mean()-a.mean())/sqrt(a.var()/len(b))
        a2 = isp(S//6 | 1) & (S//6 % 2 == 1); b2 = isp(T//6 | 1) & (T//6 % 2 == 1)
        zp = (b2.mean()-a2.mean())/sqrt(a2.mean()*(1-a2.mean())/len(b2))
        print(f"  window ({lo:.0e},{hi:.0e}] z={z:>5}: #surv={len(S):>7} #twin={len(T):>6} | <omega> surv={a.mean():.4f} twin={b.mean():.4f} zsc={zo:+.2f} | P(n/6 prime) surv={a2.mean():.5f} twin={b2.mean():.5f} zsc={zp:+.2f}")

print("\n(ii) direct Type II sum: B(A) = sum_{A<=a<2A} lambda(a) * #{b : a*b is a twin center}  / sqrt(total pairs)   (null: O(1) noise)")
tc = centers_all[twin]
tcset = np.zeros(x+1, dtype=bool); tcset[tc] = True
print(f"{'A':>8} {'pairs':>9} {'B(A)':>8} {'B(A)/sqrt':>10}   {'same with mu(a)':>15}")
for A in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]:
    tot = 0; B = 0.0; Bm = 0.0
    for a in range(A, 2*A):
        cnt = int(tcset[a::a].sum())
        tot += cnt; B += int(liou[a])*cnt; Bm += int(mu[a])*cnt
    print(f"{A:>8} {tot:>9} {B:>8.0f} {B/sqrt(tot):>10.2f}   {Bm/sqrt(tot):>15.2f}")
