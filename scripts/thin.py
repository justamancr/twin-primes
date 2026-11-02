import numpy as np
from math import log
X = 10**9
s = np.ones(X//2+2, dtype=bool); s[0]=False
for i in range(1, int(X**0.5)//2+1):
    if s[i]:
        p=2*i+1; s[p*p//2::p]=False
op = 2*np.nonzero(s)[0]+1; pr = op[op>1000]
tw = op[:-1][s[(op[:-1]+2)//2]]; tw = tw[tw>1000]; ce = tw+1
rng = np.random.default_rng(1)
rho = len(tw)/len(pr)
print(f"thinning prob rho = {rho:.4f}; model ratio for independent thinning = 1-rho = {1-rho:.3f}")
print(f"{'q':>5} {'twins':>7} {'thinned primes (3 seeds)':>26} {'shifted twins p-1':>18}")
for q in [31,61,113,179,251,401,701,1009]:
    al=[c for c in range(q) if c not in (1,q-1)]
    def chi(v, classes):
        h=np.bincount(v%q, minlength=q)[classes].astype(float); e=h.mean(); return ((h-e)**2/e).sum()/(len(classes)-1)
    ct = chi(ce, al)
    th = [chi(pr[rng.random(len(pr))<rho], list(range(1,q))) for _ in range(3)]
    # control: use lower member p (classes exclude 0 and q-2)
    al2=[c for c in range(q) if c not in (0,q-2)]
    print(f"{q:>5} {ct:>7.3f} {th[0]:>8.3f} {th[1]:>8.3f} {th[2]:>8.3f} {chi(tw, al2):>18.3f}")
