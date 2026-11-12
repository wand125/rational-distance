import math
B=10**7; N=B; n=B//2
I64=2**63-1; I128=2**127-1; U64=2**64-1
def bits(x): return x.bit_length()
rows=[
 ("N^2 (nn,U,V,|Q|)",N**2),
 ("T<=N^3",N**3),
 ("R<=2N^4",2*N**4),
 ("W=2N^2-d^2<=2N^2",2*N**2),
 ("U*V<=N^4",N**4),
 ("K<=2N^6",2*N**6),
 ("sqrtK<=sqrt2 N^3",math.isqrt(2*N**6)+1),
 ("T+sqrtK",N**3+math.isqrt(2*N**6)+1),
 ("A<=N^3",N**3),
 ("2Ap<=2N^4",2*N**4),
 ("|X|,|Y|<=N^3",N**3),
 ("X^2+Y^2 / (A+p)^2 <= 2(N^3+N)^2",2*(N**3+N)**2),
 ("D^2=(N-x)^2+(N-y)^2<=2(N^3+N)^2",2*(N**3+N)**2),
 ("u*u (u32 side, n)",n*n),
 ("dense_even_pairs sum",sum((12*k+1)*(12*k+2)//2 for k in range(1,B//24+1))),
]
for name,v in rows:
    print(f"{name:40s} bits={bits(v):3d} i64={'ok' if v<=I64 else 'OVER'} i128={'ok' if v<=I128 else 'OVER'}")
# CRT distinct primes for n<=5e6
P=[2,3,5,7,11,13,17,19,23]; prod=1
for i,p in enumerate(P):
    prod*=p
    if prod>n: print("max distinct primes of n<=",n,":",i); break
print("critical primes p=1mod4 with 12p^2<=n: max p", max(p for p in range(5,700) if all(p%d for d in range(2,int(p**.5)+1)) and p%4==1 and 12*p*p<=n))
print("plans memory MB per worker ~", (n+1)*(7*3*4+4)/1e6)
