import numpy as np, sys
from math import log
X = int(float(sys.argv[1]))
s = np.ones(X//2+2, dtype=bool); s[0]=False
for i in range(1, int(X**0.5)//2+1):
    if s[i]:
        p=2*i+1; s[p*p//2::p]=False
op = 2*np.nonzero(s)[0]+1
pr = op[op > 1000]
tw = op[:-1][s[(op[:-1]+2)//2] ] ; tw = tw[tw>1000]; ce = tw+1
print(f"x={X:.0e}  primes={len(pr)}  twins={len(tw)}")
print(f"{'q':>4} {'PRIMES chi2/dof':>16} {'TWINS chi2/dof':>15} {'logq/logx':>10}")
for q in [31,37,41,47,53,61,71,83,97,113,131,151,179,211,251,307,401,503,701,1009]:
    h = np.bincount(pr % q, minlength=q)[1:].astype(float); e=h.mean()
    cp = ((h-e)**2/e).sum()/(q-2)
    al=[c for c in range(q) if c not in (1,q-1)]
    h2 = np.bincount(ce % q, minlength=q)[al].astype(float); e2=h2.mean()
    ct = ((h2-e2)**2/e2).sum()/(len(al)-1)
    print(f"{q:>4} {cp:>16.3f} {ct:>15.3f} {log(q)/log(X):>10.3f}")
