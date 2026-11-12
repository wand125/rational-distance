# mw_sage.py (mathematics-52): second system for the Mordell-Weil sieve inputs (Sage, mwrank / Cremona-Siksek).
# usage: sage -python mw_sage.py EXPORT.txt OUT.txt   (resumable)
# For each tab-separated line "r j ainvs gens ntors quartic" of mw_export.gp, checks without PARI's ellrank/ellsaturation:
#   j(E) == j of the quartic (from its invariants I, J), gens on E, rank (mwrank, proof=True) == #gens, torsion order, saturation index 1 (max_prime=-1:
#   Sage computes the bound from the height pairing, so the index of <gens> + torsion in E(Q) is 1).
# line: r j ok|BAD rank=<..> tors=<..> sat=<index> [reason]
import sys, os
try:
    from sage.all import EllipticCurve, QQ, sage_eval
except ImportError:  # passagemath (modular distributions)
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_flint import *
    from sage.all__sagemath_pari import *
    from sage.misc.sage_eval import sage_eval
inp, out = sys.argv[1], sys.argv[2]
done = set(tuple(l.split()[:2]) for l in open(out)) if os.path.exists(out) else set()
for l in open(inp):
    w = l.rstrip("\n").split("\t")
    r, j = w[0], w[1]
    if (r, j) in done: continue
    a = sage_eval(w[2]); gens = sage_eval(w[3]); nt = int(w[4]); qa, qb, qc, qd, qe = [QQ(x) for x in sage_eval(w[5])]
    I = 12*qa*qe - 3*qb*qd + qc**2; J = 72*qa*qc*qe + 9*qb*qc*qd - 27*qa*qd**2 - 27*qe*qb**2 - 2*qc**3
    E = EllipticCurve(QQ, a); why = []
    if E.j_invariant() != EllipticCurve(QQ, [-27*I, -27*J]).j_invariant(): why.append("j")
    P = [E(g) for g in gens] if gens else []
    # mwrank's 2-descent cannot prove the rank when Sha[2] is nontrivial; then retry with a deeper second descent,
    # and if it still fails record rk = -1 ("rank" not confirmed by this system) and go on with the other checks
    try:
        rk = E.rank(only_use_mwrank=True, proof=True)
    except RuntimeError:
        try:
            E.two_descent(verbose=False, second_limit=13)
            rk = E.rank(only_use_mwrank=True, proof=True)
        except RuntimeError:
            rk = -1
    if rk != len(P): why.append("rank" if rk >= 0 else "rank_unproved")
    to = E.torsion_order()
    if to != nt: why.append("tors")
    idx = 1
    if P:
        sat, idx, reg = E.saturation(P, max_prime=-1)
        if idx != 1: why.append("sat")
    with open(out, "a") as f:
        f.write("%s %s %s rank=%d tors=%d sat=%d %s\n" % (r, j, "BAD" if why else "ok", rk, to, idx, ",".join(why)))
