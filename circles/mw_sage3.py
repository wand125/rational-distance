# mw_sage3.py (mathematics-52): rank proof for one curve of the export, by several methods; prints one line.
# usage: sage -python mw_sage3.py EXPORT.txt r j
#   an: root number + rigorously bounded L(E,1) (w=+1) or L'(E,1) (w=-1) with Sage's lseries().at1/deriv_at1 error bound;
#       |value| > err gives analytic rank 0/1, then rank = analytic rank by Gross-Zagier-Kolyvagin
#   mw: mwrank on each curve of the isogeny class (second_limit 8..12)
# line: r j N=<log10 conductor> ngens=<..> an=<rank|-> mw=<rank|-> via=<..>
import sys, math
try:
    from sage.all import EllipticCurve, QQ, sage_eval
except ImportError:
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_flint import *
    from sage.all__sagemath_pari import *
    from sage.misc.sage_eval import sage_eval
inp, r, j = sys.argv[1], sys.argv[2], sys.argv[3]
mode = sys.argv[4] if len(sys.argv) > 4 else "an"
for l in open(inp):
    w = l.rstrip("\n").split("\t")
    if (w[0], w[1]) != (r, j): continue
    gens = sage_eval(w[3]); ng = len(gens) if gens else 0
    E = EllipticCurve(QQ, sage_eval(w[2])).minimal_model()
    N = E.conductor(); res = "-"; via = "-"
    if mode == "an" and ng <= 1:
        k = int(4 * math.sqrt(float(N))) + 10
        if k < 3 * 10**7:
            eps = E.root_number()
            v, err = (E.lseries().at1(k) if eps == 1 else E.lseries().deriv_at1(k))
            if abs(v) > err: res, via = (0 if eps == 1 else 1), "eps=%d,|v|=%.4g,err=%.2g,k=%d" % (eps, abs(v), err, k)
            else: via = "eps=%d,|v|=%.4g,err=%.2g" % (eps, abs(v), err)
        else: via = "k=%d too big" % k
    if mode == "mw":
        for i, C in enumerate(E.isogeny_class().curves):
            for sl in (8, 10, 12):
                try:
                    m = C.mwrank_curve(); m.two_descent(verbose=False, second_limit=sl)
                    if m.certain(): res, via = m.rank(), "iso%d/sl%d" % (i, sl); break
                except RuntimeError: pass
            if res != "-": break
    print(r, j, "logN=%.1f" % math.log10(float(N)), "ngens=%d" % ng, "%s=%s" % (mode, res), "ok" if res == ng else "no", via, flush=True)
