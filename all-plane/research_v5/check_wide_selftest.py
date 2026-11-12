#!/usr/bin/env python3
"""Independent Python big-integer check of the v5 256-bit product and integer square root."""
import csv, math, sys
def val(s):
    h, b = s.split('*2+'); return int(h)*2+int(b)
n = 0
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f):
        a, b, c = val(row['a']), val(row['b']), int(row['c'])
        assert a*b == int(row['product']), row
        assert math.isqrt(a*b+c) == val(row['isqrt_product_plus_c']), row
        n += 1
assert n == 20000
print('wide selftest OK', n)
