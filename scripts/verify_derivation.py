"""Numerical checks of the derivation:
 (1) S(Y) = sum_{d<=Y,(d,6)=1} mu^2(d) 2^omega(d)/phi2(d)  ~  (1/(24 C2)) log^2 Y
 (2) H - sum_{h<=H} rho(h)  ~  3 (S(H/6) - 1)
 (3) H/kappa - sum_{h<=H, q|h} rho(h)  ~  3 (q/kappa) (S_q(H/(6q)) - 1),  S_q excludes q|d
 (4) chi2/dof prediction: exact-sum vs limit 1/2
"""
import numpy as np, sys, time
from math import log
x = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**6
C2 = 0.6601618158468696
t0 = time.time()
b = np.ones(x//2+2, dtype=bool); b[0]=False
for i in range(1, int(x**0.5)//2+1):
    if b[i]:
        p=2*i+1; b[p*p//2::p]=False
op = 2*np.nonzero(b)[0]+1
# --- (1) S(Y): multiplicative sieve of a(d)=mu^2 2^omega/phi2 over d<=Y, (d,6)=1
Y = x//6
a = np.ones(Y+1); a[0]=0
sqf = np.ones(Y+1, dtype=bool)
for p in op[op>=5]:
    p=int(p)
    if p > Y: break
    a[p::p] *= 2.0/(p-2)
    if p*p <= Y: sqf[p*p::p*p] = False
# kill d with factor 2 or 3
a[2::2] = 0; a[3::3] = 0; a[~sqf] = 0
S = np.cumsum(a)                     # S[Y'] = sum_{d<=Y'} a(d)
A = 1/(24*C2)
print("(1) S(Y) vs A log^2 Y  (A = 1/(24 C2) = %.5f)" % A)
for Yp in [10**3, 10**4, 10**5, 10**6, Y]:
    if Yp > Y: break
    L = log(Yp)
    print(f"   Y={Yp:>9}: S={S[Yp]:9.4f}  A log^2 Y={A*L*L:9.4f}  S - A log^2 Y = {S[Yp]-A*L*L:8.4f}  (should be ~ B log Y + C)")
# fit B, C on the tail
Ys = np.array([2**k for k in range(8, int(np.log2(Y))+1)]); Ls = np.log(Ys)
resid = S[Ys] - A*Ls*Ls
Bc, Cc = np.polyfit(Ls, resid, 1)
print(f"   fit: S(Y) ≈ A log^2 Y + {Bc:.4f} log Y + {Cc:.4f}")

# --- rho(h) exact for h = 6k <= x
Kmax = x//6
rho = np.full(Kmax+1, 6.0)
for p in op[op>=5]:
    p=int(p)
    g=(1-4/p)/(1-2/p)**2; m2=(1-3/p)/(1-2/p)**2; z0=1/(1-2/p)
    rho *= g; rho[p::p] *= z0/g
    inv6 = pow(6,-1,p)
    for r in ((2*inv6)%p, (-2*inv6)%p): rho[r::p] *= m2/g
rho[0]=0
R = np.cumsum(rho)                   # R[k] = sum_{h<=6k} rho(h)
print(f"\nrho built ({time.time()-t0:.0f}s); mean rho = {R[-1]/x:.6f}")
print("(2) H - sum_{h<=H} rho(h)   vs   3 (S(H/6) - 1)")
for H in [10**4, 10**5, 10**6, 10**7]:
    if H > x: break
    K = H//6
    print(f"   H=1e{round(log(H,10))}: Delta={H - R[K]:9.3f}   3(S(K)-1)={3*(S[K]-1):9.3f}   (1/(8C2))log^2 H={log(H)**2/(8*C2):8.3f}")
print("(3) restricted: H/kappa - sum_{q|h} rho   vs   3 (q/kappa)(S_q(H/(6q)) - 1)")
for q in [31, 113, 1009]:
    kappa = q-2
    aq = a.copy(); aq[q::q] = 0; Sq = np.cumsum(aq)
    for H in [10**5, 10**6, 10**7]:
        if H > x: break
        K = H//6
        # sum over h<=H with q|h: h=6k, q|k
        s = rho[q:K+1:q].sum()
        J = K//q
        print(f"   q={q:>5} H=1e{round(log(H,10))}: Delta_q={H/kappa - s:9.3f}   3(q/k)(S_q(J)-1)={3*q/kappa*(Sq[J]-1):9.3f}")
# --- (4) chi2/dof: exact-sum prediction vs limit
print("\n(4) chi2/dof: exact-sum HL prediction (unweighted, (T/x)^2 approx) vs limit 1/2")
tw = op[:-1][b[(op[:-1]+2)//2]]; tw = tw[tw>1000]; T = len(tw)
h = 6*np.arange(Kmax+1, dtype=np.float64); w = rho*np.maximum(x-h,0); Aall = 2*w[1:].sum()
for q in [31,61,113,251,503,1009]:
    kappa=q-2; Aq = 2*w[q::q].sum()
    EV = T*(1-1/kappa) + (T/x)**2*(Aq - Aall/kappa)
    print(f"   q={q:>5}: prediction {EV/((kappa-1)*T/kappa):.3f}   leading-order 1 - 6A (T/x) log^2(x/6) = {1 - 6*A*(T/x)*log(x/6)**2:.3f}   limit 0.500")
