"""Type II probe done right: residual of #{b: a*b twin center} against the wheel model, lambda-weighted."""
import numpy as np, sys
from math import sqrt, log
x = 10**8
s = np.ones(x//2+2, dtype=bool); s[0]=False
for i in range(1, int(x**0.5)//2+1):
    if s[i]:
        p=2*i+1; s[p*p//2::p]=False
isp = lambda n: s[n//2]
c = np.arange(6, x+1, 6, dtype=np.int64)
tw = c[isp(c-1) & isp(c+1) & (c > 1000)]
N = len(tw)
tcset = np.zeros(x+1, dtype=bool); tcset[tw] = True
# small primes list and Liouville for a up to 2^15
AMAX = 2**15
spf = np.zeros(AMAX+1, dtype=np.int64)
for p in range(2, AMAX+1):
    if spf[p]==0:
        spf[p::p][spf[p::p]==0] = p
def factor(a):
    f = {}
    while a > 1:
        p = spf[a]; k = 0
        while a % p == 0: a//=p; k+=1
        f[p]=k
    return f
def model(a):
    """wheel-model probability that a random twin center is divisible by a"""
    pr = 1.0
    for p,k in factor(a).items():
        if p in (2,3): pr *= float(p)**(-(k-1))
        else: pr *= (1.0/(p-2)) * float(p)**(-(k-1))
    return pr
print(f"N twins = {N}")
print(f"{'A':>6} {'#a':>5} {'sum obs':>9} {'sum model':>10} {'chi2/dof':>9} {'lam-weighted resid / sd':>24} {'mu-weighted / sd':>17}")
for e in range(1, 15):
    A = 2**e
    obs=[]; exp_=[]; lam=[]; mus=[]
    for a in range(A, 2*A):
        cnt = int(tcset[a::a].sum()); m = N*model(a)
        f = factor(a); Om = sum(f.values()); sqf = all(k==1 for k in f.values())
        obs.append(cnt); exp_.append(m); lam.append((-1)**Om); mus.append((-1)**Om if sqf else 0)
    obs=np.array(obs,float); exp_=np.array(exp_); lam=np.array(lam,float); mus=np.array(mus,float)
    res = obs-exp_
    chi2 = ((res**2)/np.maximum(exp_,1)).sum()/len(obs)
    sd = sqrt((exp_).sum())          # Poisson-scale sd of the weighted sum
    print(f"{A:>6} {A:>5} {obs.sum():>9.0f} {exp_.sum():>10.0f} {chi2:>9.3f} {(lam*res).sum()/sd:>24.2f} {(mus*res).sum()/sd:>17.2f}")
