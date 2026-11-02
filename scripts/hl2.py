"""Second-order Hardy–Littlewood prediction for twin-center variance across residue classes mod q."""
import numpy as np, sys, time
from math import log
x = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**8
t0 = time.time()
# primes to x
b = np.ones(x//2+2, dtype=bool); b[0]=False
for i in range(1, int(x**0.5)//2+1):
    if b[i]:
        p=2*i+1; b[p*p//2::p]=False
op = 2*np.nonzero(b)[0]+1
# measured twins
tw = op[:-1][b[(op[:-1]+2)//2]]; tw = tw[tw>1000]; ce = tw+1; T = len(ce)
# rho(h) for h = 6k, k = 1..x/6
Kmax = x//6
rho = np.full(Kmax+1, 6.0)          # index k <-> h = 6k ; rho[0] unused
for p in op[op>=5]:
    p = int(p)
    g = (1-4/p)/(1-2/p)**2; m2 = (1-3/p)/(1-2/p)**2; z0 = 1/(1-2/p)
    rho *= g                         # generic factor for all h, then fix the special classes
    # h = 6k ≡ 0 mod p  <=> k ≡ 0 mod p
    rho[p::p] *= z0/g
    # h ≡ ±2 mod p <=> 6k ≡ ±2 mod p <=> k ≡ ±2*inv6 mod p
    inv6 = pow(6, -1, p)
    for r in ((2*inv6) % p, (-2*inv6) % p):
        rho[r::p] *= m2/g
    if p > 10**7 and p % 1000003 == 0: print(f"  p={p} {time.time()-t0:.0f}s", flush=True)
rho[0] = 0
print(f"x={x:.0e}: rho built in {time.time()-t0:.0f}s; mean rho over h (all h, incl zeros) = {rho[1:].sum()/x:.5f} (should be ~1)")
h = 6*np.arange(Kmax+1, dtype=np.float64)
w = rho * np.maximum(x - h, 0)       # rho(h)*(x-h)
A = 2*w[1:].sum()
print(f"T={T}, (T/x)^2 * A = {(T/x)**2*A:.0f}  vs T^2 = {T**2}  (ratio {(T/x)**2*A/T**2:.4f})")
print(f"\n{'q':>5} {'K':>5} {'HL predicted chi2/dof':>22} {'measured':>10}")
for q in [31,37,41,47,53,61,71,83,97,113,131,151,179,211,251,307,401,503,701,1009]:
    K = q-2
    # h ≡ 0 mod q <=> 6k ≡ 0 mod q <=> k ≡ 0 mod q
    Aq = 2*w[q::q].sum()
    EV = T*(1-1/K) + (T/x)**2*(Aq - A/K)
    pred = EV/((K-1)*T/K)
    al=[c for c in range(q) if c not in (1,q-1)]
    hh=np.bincount(ce%q, minlength=q)[al].astype(float); e=hh.mean(); meas=((hh-e)**2/e).sum()/(K-1)
    print(f"{q:>5} {K:>5} {pred:>22.3f} {meas:>10.3f}")
