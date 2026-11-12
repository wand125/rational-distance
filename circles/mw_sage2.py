# mw_sage2.py (mathematics-52): retry the rank proof for lines that mw_sage.py marked rank_unproved.
# usage: sage -python mw_sage2.py EXPORT.txt SAGE_OUT.txt OUT.txt
# Rank is an isogeny invariant: run mwrank (proof=True) on every curve in the isogeny class of E (2-isogenies exist
# since E has rational 2-torsion), each with second_limit 13 and 16; the first one that proves the rank decides.
# line: r j ok|BAD rank=<..> ngens=<..> via=<index in isogeny class>|-
import sys, os
try:
    from sage.all import EllipticCurve, QQ, sage_eval
except ImportError:
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_flint import *
    from sage.all__sagemath_pari import *
    from sage.misc.sage_eval import sage_eval
inp, sg, out = sys.argv[1:4]
want = set(tuple(l.split()[:2]) for l in open(sg) if "rank_unproved" in l)
done = set(tuple(l.split()[:2]) for l in open(out)) if os.path.exists(out) else set()
for l in open(inp):
    w = l.rstrip("\n").split("\t")
    r, j = w[0], w[1]
    if (r, j) not in want or (r, j) in done: continue
    gens = sage_eval(w[3]); ng = len(gens) if gens else 0
    E = EllipticCurve(QQ, sage_eval(w[2]))
    rk, via = -1, "-"
    for k, C in enumerate(E.isogeny_class().curves):
        for sl in (13, 16):
            try:
                mw = C.mwrank_curve(); mw.two_descent(verbose=False, second_limit=sl)
                if mw.certain(): rk, via = mw.rank(), "%d/%d" % (k, sl); break
            except RuntimeError:
                pass
        if rk >= 0: break
    with open(out, "a") as f:
        f.write("%s %s %s rank=%d ngens=%d via=%s\n" % (r, j, "ok" if rk == ng else "BAD", rk, ng, via))
