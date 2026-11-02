from sympy import factorint, primerange, isprime
k=2996863034895; m=1290000
print("k factors:", factorint(k))
print("m factors:", factorint(m))
# residues of center n = k*2^m mod small primes; twin needs n != ±1 mod q
rows=[]
zero=0; dodge=0
for q in primerange(3,200):
    r=(k%q)*pow(2,m,q)%q
    assert r not in (1,q-1)
    rows.append((q,r)); 
    if r==0: zero+=1
    else: dodge+=1
print(rows[:25])
print("primes<200 where center ≡0:",zero," where center lands on some other residue:",dodge)
# previous records, same question
for k,m in [(3756801695685,666669),(65516468355,333333),(2003663613,195000)]:
    print(k, factorint(k), [q for q in primerange(3,100) if (k%q)*pow(2,m,q)%q==0])
# the primorial-neighbor construction for comparison
from sympy import primorial
for i in range(1,16):
    P=primorial(i); print(i, P, isprime(P-1), isprime(P+1))
