# mw_one.py (mathematics-52): one Mordell-Weil sieve job (circle r, pair j) in a fresh PARI; prints one line.
# usage: python mw_one.py r j Lmax ncontrol
# line: r j ok rank ntor res=[status,M,#primes,#surv] controls=n wrong=k blowup=b ms=setup,primes,sieve,controls
#   status 1 = all classes killed, 0 = survivors remain, -1 = too many survivors (MW_CAP)
#   cross=[#points,#checks,#violations]: exact small points vs class tests (only after a kill)
#   a control is "wrong" if the sieve killed everything although a rational point satisfies the fake test
import sys, os, time, cypari2
p = cypari2.Pari(); p.allocatemem(400 * 10**6, silent=True)
p.default("parisizemax", 2 * 10**9); p.default("nbthreads", 1)
p('read("circles.gp"); read("mwsieve.gp")')
if os.environ.get("MW_ALLOW_AN") == "1": p('MW_ALLOW_AN = 1')
if os.environ.get("MW_HEEGNER") in ("1", "2", "3"): p('MW_HEEGNER = %s' % os.environ["MW_HEEGNER"])
if os.environ.get("MW_COVB"): p('MW_COVB = %s' % int(float(os.environ["MW_COVB"])))
r, j, Lmax, nc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t = [time.time()]
p('MWS = iferr(mw_setup(%s, %d), e, [0, Str(e)])' % (r, j)); t.append(time.time())
if p('MWS == 0 || MWS[1] == 0'):
    print(r, j, "skip", "nopt" if p('MWS == 0') else str(p('MWS[2..#MWS]')).replace(" ", ""), flush=True); sys.exit()
p('MWPD = mw_primes(MWS, %d)' % Lmax); t.append(time.time())
res = p('mw_sieve(MWS, MWPD)'); t.append(time.time())
pts = p('[x | x <- mw_points(cq(%s)[%d], 300), type(x) != "t_VEC"]' % (r, j))
wrong = blow = done = 0
for i in range(min(nc, len(pts))):
    c = int(p('mw_sieve(mw_control(MWS, %s, %d), MWPD)[1]' % (pts[i], i + 1))); done += 1
    wrong += (c == 1); blow += (c == -1)
cross = p('mw_cross(MWS, MWPD, if (#MWS[9] <= 1, 4, 2))') if int(res[0]) == 1 else None
t.append(time.time())
# record exactly what the sieve used (curve, generators, torsion, quartic) for the second system (mw_sage.py)
if int(res[0]) == 1 and os.environ.get("MW_DUMP"):
    p('my(qq = cq(%s)[%d]); write("%s", "%s", "\\t", %d, "\\t", Vec(MWS[7][1..5]), "\\t", MWS[9], "\\t", #MWS[10], "\\t", [polcoef(qq, k, \'t) | k <- [4,3,2,1,0]])'
      % (r, j, os.environ["MW_DUMP"], r, j))
ms = [int(1000 * (t[k + 1] - t[k])) for k in range(len(t) - 1)]
print(r, j, "ok", p('MWS[13]'), p('#MWS[10]'), "res=" + str(res).replace(" ", ""), "controls=%d" % done, "wrong=%d" % wrong,
      "blowup=%d" % blow, "cross=" + (str(cross).replace(" ", "") if cross is not None else "-"), "ms=" + ",".join(map(str, ms)), "an=%s" % p('MWS[14]'), flush=True)
