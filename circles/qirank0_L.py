# qirank0_L.py (mathematics-52): second system (Sage) for qirank0.gp.
# usage: sage -python qirank0_L.py CURVES.txt   (lines "r j [a4,a6] [a4,a6] ntors_K" from qirank0.gp: J and its twist by -1)
# For J and J^(-1): root number and L(E,1) with Sage's rigorous truncation error (lseries().at1, 4*sqrt(N)+10 terms);
# |L| > err gives rank 0 (Kolyvagin). Also #J(Q(i))_tors in Sage, compared with PARI's count over Q(i).
# line: r j rank0_J=<0|1> rank0_Jt=<0|1> torsK=<n> parity=<ok|BAD> |L|,err for both
import sys, math
try:
    from sage.all import EllipticCurve, QQ, NumberField, PolynomialRing, sage_eval
except ImportError:
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_flint import *
    from sage.all__sagemath_pari import *
    from sage.misc.sage_eval import sage_eval
R = PolynomialRing(QQ, 'x'); K = NumberField(R.gen()**2 + 1, 'i')
for l in open(sys.argv[1]):
    w = l.rstrip("\n").split("\t")
    r, j, nt = w[0], w[1], int(w[4])
    out = []
    ok = True
    for ab in (w[2], w[3]):
        a4, a6 = sage_eval(ab)
        E = EllipticCurve(QQ, [a4, a6]).minimal_model()
        N = E.conductor(); k = int(4 * math.sqrt(float(N))) + 10
        if E.root_number() != 1: out.append("eps=-1"); ok = False; continue
        v, err = E.lseries().at1(k)
        if not abs(v) > err: ok = False
        out.append("|L|=%.4g,err=%.2g,logN=%.1f" % (abs(v), err, math.log10(float(N))))
    J = EllipticCurve(QQ, sage_eval(w[2]))
    tk = J.change_ring(K).torsion_order()
    print(r, j, "ok" if ok and tk == nt else "BAD", "torsK=%d/%d" % (tk, nt), " ".join(out), flush=True)
