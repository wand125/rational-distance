#!/usr/bin/env python3
"""Local check of the filters for primitive tuples (signed integers).
For m in (48, 144): every residue triple (A, C, Q) mod m with B=Q-A, D=Q-C, not all of
A,B,C,D even, and A²+C², B²+C², A²+D², B²+D² all squares mod m must satisfy:
12 | Q, C mod 4 != 2, (C odd => A,B ≡ 0 mod 4), (4|C => A,B odd), (3|A <=> 3|B)."""
for m in (48, 144):
    sq = {(i*i) % m for i in range(m)}
    bad = 0; ok = 0
    for A in range(m):
        for C in range(m):
            for Q in range(m):
                B, D = (Q-A) % m, (Q-C) % m
                if A % 2 == 0 and B % 2 == 0 and C % 2 == 0 and D % 2 == 0: continue
                if not all((x*x+y*y) % m in sq for x, y in ((A, C), (B, C), (A, D), (B, D))): continue
                ok += 1
                cond = (Q % 12 == 0 and C % 4 != 2 and (C % 2 == 0 or (A % 4 == 0 and B % 4 == 0))
                        and (C % 4 != 0 or (A % 2 == 1 and B % 2 == 1)) and ((A % 3 == 0) == (B % 3 == 0)))
                if not cond: bad += 1
    print(m, 'residue tuples', ok, 'violations', bad)
