#!/usr/bin/env python3
"""Bounded exhaustive proof check of gamma-compatible local CRT residue sets."""
import json
import math


def valuation_capped(t,p,e):
    if t==0:
        return e
    k=0
    while k<e and t%p==0:
        k+=1
        t//=p
    return k


def direct_gamma(n,t):
    d=math.gcd(n,t)
    return math.gcd(d,n//d)


def branch_condition(p,e,h,u):
    k=valuation_capped(h,p,e)
    if k==0:
        return True
    v=u-h
    if k==e:
        return u%p!=0 or u%(p**e)==0
    if 2*k<e:
        modulus=p**(e-k)
        high=u if u%modulus==0 else v
        return high%(p*modulus)!=0
    if 2*k==e:
        modulus=p**k
        return (u//modulus)%p!=0 and (v//modulus)%p!=0
    low=p**(e-k)
    high=p**k
    return u%(p*low)!=0 or (
        u%high==0 and u%(p*high)!=0 and v%(p*high)!=0)


def main():
    results=[]
    for p in [3,5,7,13,17]:
        for e in range(1,7):
            modulus=p**e
            if modulus>500:
                break
            checked=0
            for h in range(modulus):
                k=valuation_capped(h,p,e)
                crt_modulus=p**(e-k)
                for u in range(modulus):
                    if u%crt_modulus not in (0,h%crt_modulus):
                        continue
                    expected=(direct_gamma(modulus,u)==direct_gamma(modulus,u-h)
                              ==direct_gamma(modulus,h))
                    actual=branch_condition(p,e,h,u)
                    if expected!=actual:
                        raise AssertionError((p,e,h,u,expected,actual))
                    checked+=1
            results.append({'prime':p,'exponent':e,'modulus':modulus,
                            'compatible_pairs_checked':checked,'exact_match':True})
    print(json.dumps({'checks':results,'all_exact_matches':True},indent=2))


if __name__=='__main__':
    main()
