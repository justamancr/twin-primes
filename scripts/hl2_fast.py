"""Exact-sum HL prediction of chi2/dof for twin centers mod q, at large x (vectorized rho builder)."""
import numpy as np, sys, time
from math import log
x = int(float(sys.argv[1])); C2 = 0.6601618158468696; t0=time.time()
b = np.ones(x//2+2, dtype=bool); b[0]=False
for i in range(1, int(x**0.5)//2+1):
    if b[i]:
        p=2*i+1; b[p*p//2::p]=False
op = 2*np.nonzero(b)[0]+1
tw = op[:-1][b[(op[:-1]+2)//2]]; tw = tw[tw>1000]; T = len(tw); ce = tw+1
del tw
Kmax = x//6
logrho = np.full(Kmax+1, log(6.0))
small = op[(op>=5)&(op<=Kmax)]; big = op[op>Kmax]
for p in small:
    p=int(p); g=(1-4/p)/(1-2/p)**2; m2=(1-3/p)/(1-2/p)**2; z0=1/(1-2/p)
    logrho += log(g); logrho[p::p] += log(z0/g)
    inv6 = pow(6,-1,p)
    for r in ((2*inv6)%p, (-2*inv6)%p): logrho[r::p] += log(m2/g)
print(f"small primes done ({time.time()-t0:.0f}s)", flush=True)
# big primes: global factor + at most one special class element each
pb = big.astype(np.float64)
logrho += np.sum(np.log((1-4/pb)/(1-2/pb)**2))
inv6 = np.where(big % 6 == 1, (5*big+1)//6, (big+1)//6)
r1 = (2*inv6) % big; r2 = big - r1
corr = np.log((1-3/pb)/(1-2/pb)**2) - np.log((1-4/pb)/(1-2/pb)**2)
for r in (r1, r2):
    m = r <= Kmax
    np.add.at(logrho, r[m], corr[m])
rho = np.exp(logrho); rho[0]=0; del logrho
print(f"x={x:.0e}: rho built ({time.time()-t0:.0f}s); mean rho = {rho.sum()/x:.6f}; T={T}", flush=True)
h = 6*np.arange(Kmax+1, dtype=np.float64); w = rho*(x-h); del rho, h
Aall = 2*w[1:].sum()
print(f"{'q':>5} {'prediction':>11} {'measured':>9}")
for q in [31,61,113,251,503,1009]:
    k=q-2; Aq = 2*w[q::q].sum()
    EV = T*(1-1/k) + (T/x)**2*(Aq - Aall/k)
    al=[c for c in range(q) if c not in (1,q-1)]
    hh=np.bincount(ce%q, minlength=q)[al].astype(float); e=hh.mean(); meas=((hh-e)**2/e).sum()/(k-1)
    print(f"{q:>5} {EV/((k-1)*T/k):>11.3f} {meas:>9.3f}")
