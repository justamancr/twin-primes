"""Phase 0 — the instrument.
Twin-center wheel, prime sieve, Hardy–Littlewood check, Mertens, twin-paradox ratio.
"""
import numpy as np, math, sys, time
from math import log, exp, pi

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
GAMMA = 0.5772156649015329
C2 = 0.6601618158468696      # twin prime constant
TWO_C2 = 2 * C2

# ---------- 1. prime sieve (odd-only bitmap) ----------
t = time.time()
sieve = np.ones(N // 2 + 1, dtype=np.bool_)   # index i -> 2i+1
sieve[0] = False                               # 1 is not prime
for i in range(1, int(N**0.5) // 2 + 1):
    if sieve[i]:
        p = 2 * i + 1
        sieve[p * p // 2::p] = False
primes = 2 * np.nonzero(sieve)[0] + 1
primes = np.concatenate(([2], primes))
print(f"sieve to {N:.0e}: {len(primes)} primes in {time.time()-t:.1f}s")

# ---------- 2. twin counts vs. published table ----------
is_p = sieve  # odd n prime <=> sieve[n//2]
twin_lower = primes[1:][is_p[(primes[1:] + 2) // 2][:len(primes)-1] if False else None] if False else None
# simpler: twin pairs (p, p+2) with p odd prime
p_odd = primes[1:]
p_odd = p_odd[p_odd + 2 <= N]
mask = is_p[(p_odd + 2) // 2]
twins = p_odd[mask]              # lower member of each pair
published = {10**4:205, 10**5:1224, 10**6:8169, 10**7:58980, 10**8:440312, 10**9:3424506, 10**10:27412679}
print("\n--- twin pairs below x: computed vs published vs Hardy–Littlewood ---")
def li2(x, step=1000):
    # 2*C2 * integral_2^x dt/ln^2 t  (Simpson on log scale)
    us = np.linspace(log(2), log(x), step*2+1)
    f = np.exp(us) / us**2          # dt = e^u du, 1/ln^2 t = 1/u^2
    h = us[1]-us[0]
    return TWO_C2 * h/3 * (f[0] + f[-1] + 4*f[1:-1:2].sum() + 2*f[2:-1:2].sum())
x = 10**4
while x <= N:
    c = int(np.searchsorted(twins, x - 2, side='right'))
    hl = li2(x)
    print(f"x=1e{int(log10(x)) if False else round(log(x,10))}: computed={c:>10}  published={published.get(x,'?'):>10}  HL≈{hl:>12.0f}  ratio={c/hl:.4f}")
    x *= 10

# ---------- 3. Mertens and the 1.32 constant from the wheel ----------
print("\n--- Mertens: prod(1-1/p) * ln(p) * e^gamma -> 1 ---")
prod1 = 1.0
for i, p in enumerate(primes):
    prod1 *= (1 - 1/p)
    if p in (10**k for k in range(1, 8)) or (i < 10):
        pass
checks = [10**k for k in range(1,8)]
prod1 = 1.0; j = 0
for p in primes:
    prod1 *= (1 - 1/p)
    if j < len(checks) and p > checks[j]:
        print(f"p≤1e{j+1}: prod(1-1/p)·ln(p)·e^γ = {prod1*log(p)*exp(GAMMA):.5f}")
        j += 1

print("\n--- twin constant from the wheel: prod_{3≤q≤p} (1-2/q)/(1-1/q)^2 -> C2 = 0.66016 ---")
r = 1.0; j = 0
for p in primes[1:]:
    r *= (1 - 2/p) / (1 - 1/p)**2
    if j < len(checks) and p > checks[j]:
        print(f"q≤1e{j+1}: {r:.6f}   (x2 = {2*r:.5f}, target 1.32032)")
        j += 1

# ---------- 4. the twin-center wheel ----------
print("\n--- twin-center wheel: n ≡ 0 (mod 6), n ≢ ±1 (mod q) for all q ≤ depth ---")
def center_survivors(x, depth):
    """count centers n<=x surviving sieve by odd primes 3<q<=depth (n multiple of 6)."""
    m = np.arange(6, x + 1, 6, dtype=np.int64)      # candidate centers
    alive = np.ones(len(m), dtype=np.bool_)
    for q in primes[2:]:                            # start from 5
        if q > depth: break
        r = m % q
        alive &= (r != 1) & (r != q - 1)
    return int(alive.sum()), m[alive]

def heuristic(x, depth):
    """x * (1/2) * (1/3) * prod_{5≤q≤depth}(1-2/q)   -- 'each prime kills 2/q of the remainder'"""
    h = x / 6
    for q in primes[2:]:
        if q > depth: break
        h *= (1 - 2/q)
    return h

X = min(N, 10**7)
print(f"window x = {X:.0e}; true twin pairs with lower member < x: {int(np.searchsorted(twins, X-2, side='right'))}")
print(f"{'depth':>8} {'survivors':>10} {'heuristic':>12} {'surv/heur':>10} {'certified<depth²':>18}")
for depth in [5, 7, 11, 13, 17, 23, 31, 47, 71, 101, 151, 229, 347, 523, 787, 1187, 1789, 2699, 3163]:
    if depth > X**0.5 + 1: break
    s, surv = center_survivors(X, depth)
    nxt = primes[np.searchsorted(primes, depth, side='right')]
    cert = int((surv < nxt*nxt).sum())
    print(f"{depth:>8} {s:>10} {heuristic(X, depth):>12.0f} {s/heuristic(X, depth):>10.4f} {cert:>18}")

print("\n--- the twin '1.12 paradox': sieve to √x, compare survivors to truth and to heuristic ---")
print("prediction: heuristic overcounts; truth/heuristic -> 1/(4e^{-2γ}) =", f"{1/(4*exp(-2*GAMMA)):.4f}", " (primes: 1/(2e^{-γ}) = 0.8905)")
for X in [10**5, 10**6, 10**7]:
    depth = int(X**0.5)
    s, surv = center_survivors(X, depth)
    true = int(np.searchsorted(twins, X-2, side='right'))
    h = heuristic(X, depth)
    print(f"x=1e{round(log(X,10))}: survivors={s} (these ARE the twins + nothing else: true={true})  heuristic={h:.0f}  true/heur={true/h:.4f}")
