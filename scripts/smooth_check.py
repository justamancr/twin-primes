import numpy as np
from math import log
exec(open('verify_derivation.py').read().split('# --- (4)')[0].replace('print(','(lambda *a,**k: None)(') )  # reuse arrays silently
# smoothed: W(H) = sum_{h<=H} (H-h) rho(h);  prediction  H^2/2 - 3*int_0^H (S(u/6)-1) du
print("smoothed:  H^2/2 - W(H)   vs   3 * int_0^H (S(u/6)-1) du")
Su = np.concatenate(([0.0], np.repeat(S[1:]-1, 6)))[:x+1]   # S(u/6)-1 as function of u (step)
I = np.cumsum(Su)
for H in [10**4, 10**5, 10**6]:
    K = H//6
    hh = 6*np.arange(K+1); W = (rho[:K+1]*(H-hh)).sum()
    print(f"   H=1e{round(log(H,10))}: {H*H/2 - W:12.1f}   {3*I[H]:12.1f}   ratio {(H*H/2-W)/(3*I[H]):.4f}")
print("restricted smoothed, q=113: H^2/(2k) - sum_{q|h}(H-h)rho  vs  3(q/k) int (S_q(u/(6q))-1) du")
q=113; k=q-2; aq=a.copy(); aq[q::q]=0; Sq=np.cumsum(aq)
Squ = np.concatenate(([0.0], np.repeat(Sq[1:]-1, 6*q)))[:x+1]; Iq=np.cumsum(Squ)
for H in [10**5, 10**6]:
    K=H//6; hh=6*np.arange(K+1); Wq=(rho[q:K+1:q]*(H-hh[q::q])).sum()
    print(f"   H=1e{round(log(H,10))}: {H*H/(2*k) - Wq:12.1f}   {3*q/k*Iq[H]:12.1f}   ratio {(H*H/(2*k)-Wq)/(3*q/k*Iq[H]):.4f}")
