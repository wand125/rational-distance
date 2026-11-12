#!/usr/bin/env python3
"""Independent brute force for rpart --all --dump-three.
Every (Q, A, C) with 1<=Q<=N, signed A, C, B=Q-A, D=Q-C, max(|A|,|B|,|C|,|D|)<=M.
Row (Q, min(a,b), max(a,b), c, d) when c>=1, a,b>=1, a!=b, A²+C² and B²+C² squares,
and A²+D² or B²+D² square. Also prints full four-square points (any a,b,c,d, zeros allowed)."""
import sys, math
N, M = int(sys.argv[1]), int(sys.argv[2])
sq = lambda x: math.isqrt(x)**2 == x
rows, four = set(), []
for Q in range(1, N+1):
    for A in range(Q-M, M+1):
        B = Q-A
        if abs(A) > M or abs(B) > M: continue
        a, b = abs(A), abs(B)
        for C in range(Q-M, M+1):
            D = Q-C
            if abs(C) > M or abs(D) > M: continue
            c, d = abs(C), abs(D)
            s1, s2 = sq(a*a+c*c), sq(b*b+c*c)
            if not (s1 and s2): continue
            s3, s4 = sq(a*a+d*d), sq(b*b+d*d)
            if s3 and s4: four.append((Q, A, C))
            if c >= 1 and a >= 1 and b >= 1 and a != b and (s3 or s4):
                rows.add((Q, min(a, b), max(a, b), c, d))
out = sys.argv[3]
with open(out, 'w') as f:
    f.write('Q,a,b,c,d\n')
    for r in sorted(rows): f.write(','.join(map(str, r))+'\n')
print(len(rows), 'rows;', 'four-square points (Q,A,C):', four[:20], len(four))
