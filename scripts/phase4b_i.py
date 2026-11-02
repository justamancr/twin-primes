"""Phase 4b-i — scale the sub-Poisson uniformity. Segmented sieve to X, histograms of twin centers mod q."""
import numpy as np, sys, time
from math import log
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**10
SEG = 2 * 10**8                       # numbers per segment
lim = int(X**0.5) + 2
b = np.ones(lim+1, dtype=bool); b[:2] = False
for i in range(2, int(lim**0.5)+1):
    if b[i]: b[i*i::i] = False
base = np.nonzero(b)[0]; base_odd = base[base > 2]
QS = [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251]
hist = {q: np.zeros(q, dtype=np.int64) for q in QS}
checkpoints = [10**8, 10**9, 10**10]
ntw = 0; t0 = time.time()
def report(xnow):
    print(f"\n=== x = {xnow:.0e}: twins = {ntw}  ({time.time()-t0:.0f}s) ===")
    print(f"{'q':>4} {'chi2/dof':>9} {'log q/log x':>12}")
    for q in QS:
        h = hist[q].astype(float); allowed = [c for c in range(q) if c not in (1,q-1)]
        obs = h[allowed]; e = obs.sum()/len(allowed)
        print(f"{q:>4} {((obs-e)**2/e).sum()/(len(allowed)-1):>9.3f} {log(q)/log(xnow):>12.3f}")
lo = 1; ci = 0
while lo <= X:
    hi = min(lo + SEG - 1, X + 2)
    # odd numbers in [lo, hi]
    start = lo if lo % 2 else lo + 1
    n_odd = (hi - start)//2 + 1
    seg = np.ones(n_odd, dtype=bool)
    if start == 1: seg[0] = False
    for p in base_odd:
        if p*p > hi: break
        first = max(p*p, ((start + p - 1)//p)*p)
        if first % 2 == 0: first += p
        if first > hi: continue
        seg[(first - start)//2::p] = False
    # twins: odd p with p+2 prime, both in segment; centers = p+1
    idx = np.nonzero(seg[:-1] & seg[1:])[0]
    centers = start + 2*idx + 1
    centers = centers[(centers > 1000) & (centers <= X)]
    ntw += len(centers)
    for q in QS:
        hist[q] += np.bincount(centers % q, minlength=q)
    lo = hi + 1
    while ci < len(checkpoints) and hi >= checkpoints[ci]:
        report(checkpoints[ci]); ci += 1
