"""Phase 2 — GPY weighted sieve measured on real primes.

For an admissible tuple H = {h_1..h_k}, weights
    w(n) = sum_{d | prod(n+h_i), d <= R, d squarefree}  mu(d) * (log(R/d)/log R)^(k+l)
Quantities
    S0 = sum_{n<=x} w(n)^2
    S1 = sum_{n<=x} ( sum_i theta(n+h_i) ) w(n)^2       theta(p)=log p if p prime else 0
Criterion: if S1 > S0 * log x  then some n has >= 2 primes among n+h_i   (since each theta <= log x).
rho = S1 / (S0 log x).   The method proves rho>1 only for R <= x^(theta/2); theta=1/2 is what's proven (BV).
"""
import numpy as np, sys, time
from math import log, gcd
from itertools import product

x = int(sys.argv[1]) if len(sys.argv) > 1 else 10**7

# primes & theta
s = np.ones(x + 64, dtype=bool); s[:2] = False
for i in range(2, int((x + 64) ** 0.5) + 1):
    if s[i]: s[i*i::i] = False
theta = np.zeros(x + 64); idx = np.nonzero(s)[0]; theta[idx] = np.log(idx)
primes = idx

def admissible(H):
    for p in primes[primes <= len(H)]:
        if len(set(h % p for h in H)) == p: return False
    return True

def crt_classes(H, d_primes):
    """all residues r mod d=prod(d_primes) with p | prod(r+h) for every p in d_primes"""
    d = 1; classes = [0]
    for p in d_primes:
        res = sorted(set((-h) % p for h in H))
        new = []
        inv = pow(d, -1, p)
        for c in classes:
            for r in res:
                # solve y ≡ c (mod d), y ≡ r (mod p)
                t = ((r - c) * inv) % p
                new.append(c + d * t)
        d *= p; classes = new
    return d, classes

def gpy(H, a, l, verbose=False):
    k = len(H); R = x ** a; logR = log(R)
    w = np.zeros(x + 1)
    # enumerate squarefree d <= R built from primes, with lambda_d
    ps = [int(p) for p in primes[primes <= R]]
    def rec(start, d, dp):
        # add contribution of current d (if d>1)
        if d > 1:
            lam = ((-1) ** len(dp)) * (log(R / d) / logR) ** (k + l)
            dd, classes = crt_classes(H, dp)
            for c in classes:
                w[c::dd] += lam
        for i in range(start, len(ps)):
            p = ps[i]
            if d * p > R: break
            # skip primes where tuple occupies all classes (nu_p = p) -> then every n hit; fine, keep
            rec(i + 1, d * p, dp + [p])
    rec(0, 1, [])
    w[1:] += 1.0   # d = 1 term
    w[0] = 0
    n = np.arange(0, x + 1)
    w2 = w[1:] ** 2
    S0 = w2.sum()
    S1 = sum((theta[1 + h: x + 1 + h] * w2).sum() for h in H)
    return S1 / (S0 * log(x))

print(f"x = {x:.0e}")
tuples = {
    "pair {0,2}":        [0, 2],
    "6-tuple diam 16":   [0, 4, 6, 10, 12, 16],
}
for name, H in tuples.items():
    assert admissible(H), name
    print(f"\n{name}:  rho = S1/(S0 log x);  rho>1 <=> weights certify 2 primes in the tuple")
    print(f"{'R = x^a':>10} {'l=0':>8} {'l=1':>8} {'l=2':>8}   note")
    for a in ([0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5] if len(H) == 2 else [0.15, 0.2, 0.25, 0.3]):
        t = time.time()
        vals = [gpy(H, a, l) for l in (0, 1, 2)]
        note = "<- provable today (theta=1/2)" if abs(a - 0.25) < 1e-9 else ("<- needs Elliott-Halberstam (theta=1)" if abs(a-0.5)<1e-9 else "")
        print(f"{a:>10.2f} {vals[0]:>8.3f} {vals[1]:>8.3f} {vals[2]:>8.3f}   {note}   ({time.time()-t:.0f}s)")
