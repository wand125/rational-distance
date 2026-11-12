# mw_sage4.py (mathematics-52): mwrank (proof) on curve number i of the isogeny class of the exported curve.
# usage: sage -python mw_sage4.py EXPORT.txt r j i second_limit   -> "r j iso=i/n sl=.. rank=<r|-> certain=<0|1> ngens=.."
import sys
try:
    from sage.all import EllipticCurve, QQ, sage_eval
except ImportError:
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_flint import *
    from sage.all__sagemath_pari import *
    from sage.misc.sage_eval import sage_eval
inp, r, j, i, sl = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
for l in open(inp):
    w = l.rstrip("\n").split("\t")
    if (w[0], w[1]) != (r, j): continue
    gens = sage_eval(w[3]); ng = len(gens) if gens else 0
    cs = EllipticCurve(QQ, sage_eval(w[2])).minimal_model().isogeny_class().curves
    if i >= len(cs): print(r, j, "iso=%d/%d" % (i, len(cs)), "none", flush=True); break
    m = cs[i].mwrank_curve(); m.two_descent(verbose=False, second_limit=sl)
    print(r, j, "iso=%d/%d" % (i, len(cs)), "sl=%d" % sl, "rank=%d" % m.rank(), "certain=%d" % int(m.certain()), "ngens=%d" % ng, flush=True)
