\\ mwsieve.gp (mathematics-52, 2026-11-09): Mordell-Weil sieve for the vertex circles |P| = r.
\\ A four-distance point on the circle gives t with f_a(t) f_b(t) a square for every pair, where
\\ f = [1+t^2, N2, N3, N4]. Pick one pair (a,b) whose quartic Y: v^2 = f_a f_b has a known rational
\\ point and whose elliptic curve E has rank lo = hi. Then E(Q) = <saturated generators> + torsion.
\\ For primes l, the class of a point in (Z/M)^k x torsion decides its reduction mod l, and the
\\ reduction decides t mod l (via Washington, Elliptic Curves, Thm 2.17). A class dies at l when the
\\ other two conditions f_a f_c (c not in {a,b}) are neither a square nor 0 mod l at that t.
\\ Unresolved reductions (y = 0, bad formulas) always pass, so a kill is sound given the generators.
\\ Needs circles.gp (cq).

MW_CAP = 200000;
MW_ALLOW_AN = 0;
\\ MW_HEEGNER = 1 (11/11): skip ellrank; rank = 1 from root number -1 and L'(E,1) != 0 (checked separately with an error bound,
\\ mw_sage3.py), generator = Heegner point (ellheegner), then saturated. Records MW_AN = 2.
MW_HEEGNER = 0;
MW_PAIRS = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];
mw_pairidx(a, b) = { my(p = vecsort([a, b])); for (j = 1, 6, if (MW_PAIRS[j] == p, return(j))); 0; }

\\ Washington Thm 2.17: v^2 = a u^4 + b u^3 + c u^2 + d u + q^2, q != 0.
mw_weier(F) =
{
  my(a = polcoef(F,4,'u), b = polcoef(F,3,'u), c = polcoef(F,2,'u), d = polcoef(F,1,'u), e = polcoef(F,0,'u), q);
  if (!issquare(e, &q) || q == 0, error("constant term not a nonzero square"));
  my(a2 = c - d^2/(4*q^2), a4 = -4*q^2*a);
  [[d/q, a2, 2*q*b, a4, a2*a4], [a, b, c, d, q]];
}
mw_fwd(W, u, v) = { my([a,b,c,d,q] = W); [(2*q*(v+q) + d*u)/u^2, (4*q^2*(v+q) + 2*q*(d*u + c*u^2) - d^2*u^2/(2*q))/u^3]; }
mw_inv(W, P) = { my([a,b,c,d,q] = W); (2*q*(P[1]+c) - d^2/(2*q)) / P[2]; }

\\ known rational t on the quartic qq (t = [1,0] means infinity), searched up to height B
mw_points(qq, B = 1000) =
{
  my(D = denominator(content(qq)), F = D^2 * subst(qq, 't, 'x), P = hyperellratpoints(F, B), L = List());
  for (k = 1, #P, if (P[k][2] != 0, listput(L, P[k][1])));
  if (issquare(polcoef(qq, 4, 't)), listput(L, [1,0]));
  Vec(L);
}

\\ setup for circle r, pair index j. Returns 0 if unusable, else a record.
mw_setup(r, j, B = 1000, satB = 100) =
{
  my(Q = cq(r), qq = Q[j], pts = mw_points(qq, B), t0, F, D = denominator(content(qq)), W, E, Em, ch, R, gens, tor, Tl);
  if (#pts == 0, return(0));
  \\ prefer a finite point
  t0 = pts[1]; for (k = 1, #pts, if (type(pts[k]) != "t_VEC", t0 = pts[k]; break));
  if (type(t0) == "t_VEC", F = D^2 * polrecip(subst(qq, 't, 'u)) , F = D^2 * subst(qq, 't, t0 + 'u));
  my(WW = mw_weier(F)); W = WW[2];
  E = ellinit(WW[1]); Em = ellminimalmodel(E, &ch);
  MW_AN = 0;
  if (MW_HEEGNER,
    if (ellrootno(Em) != -1, return([0, "rootno"]));
    R = [1, 1, 0, [ellheegner(Em)]]; MW_AN = 2,
    R = ellrank(Em);
    for (eff = 2, 4, if (R[1] == R[2], break); R = ellrank(Em, eff)));
  \\ optional (MW_ALLOW_AN = 1): accept rank = lower bound when the numerical analytic rank equals it and is <= 1
  \\ (then rank = analytic rank by Gross-Zagier-Kolyvagin, provided the numerical analytic rank is right)
  if (R[1] != R[2] && MW_ALLOW_AN && R[1] <= 1 && #R[4] == R[1] && ellanalyticrank(Em)[1] == R[1], MW_AN = 1; R[2] = R[1]);
  if (R[1] != R[2], return([0, "rank", R[1..2]]));
  gens = R[4];
  if (#gens < R[1], R = ellrank(Em, 3); if (R[1] == R[2], gens = R[4]));
  if (#gens < R[1] && R[1] == 1, gens = iferr([ellheegner(Em)], e, []));
  \\ the sieve needs a full set of independent generators: refuse otherwise
  if (#gens != R[1] || (#gens > 0 && matdet(ellheightmatrix(Em, gens)) == 0), return([0, "gens", [R[1], #gens]]));
  gens = if (R[1] == 0, [], ellsaturation(Em, gens, satB));
  tor = elltors(Em);
  Tl = List([[0]]);
  if (tor[1] > 1,
    my(G = tor[3], S = tor[2]); Tl = List();
    if (#S == 1, for (a = 0, S[1]-1, listput(Tl, ellmul(Em, G[1], a))),
      for (a = 0, S[1]-1, for (b = 0, S[2]-1, listput(Tl, elladd(Em, ellmul(Em, G[1], a), ellmul(Em, G[2], b)))))));
  my(ab = MW_PAIRS[j], oth = [c | c <- [1..4], c != ab[1] && c != ab[2]], tests = [Q[mw_pairidx(ab[1], c)] | c <- oth]);
  [1, r, j, t0, W, E, Em, ch, gens, Vec(Tl), tests, F, R[1], MW_AN];
}

\\ self-check of the map on small points of the quartic: forward lands on E, inverse returns u
mw_selfcheck(S) =
{
  my([ok, r, j, t0, W, E, Em, ch, gens, Tl, tests, F] = S, P = hyperellratpoints(denominator(content(F))^2*subst(F,'u,'x), 200), n = 0);
  for (k = 1, #P, my(u = P[k][1], v = P[k][2] / denominator(content(F)));
    if (u == 0, next);
    my(X = mw_fwd(W, u, v)); if (!ellisoncurve(E, X), error("fwd off curve"));
    if (X[2] != 0 && mw_inv(W, X) != u, error("inv mismatch")); n++);
  n;
}

\\ reduction mod l of a rational point on Em
mw_red(P, l) = { if (P == [0], return([0])); if (valuation(denominator(P[1]), l) > 0, return([0])); P * Mod(1, l); }

\\ does t-value pass the tests mod l (t = [1,0] or integral-in-l value as Mod)
mw_tpass(tests, tl, l) =
{
  for (k = 1, #tests, my(g = tests[k], v);
    \\ a coefficient with l in the denominator (possible only for the fake tests of mw_control): unresolved, pass
    v = iferr(if (tl === "inf", polcoef(g, 4, 't) * Mod(1, l), subst(g, 't, tl)), e, 0);
    if (v != 0 && !issquare(v), return(0)));
  1;
}

\\ per-prime data: [l, El, gens mod l, orders, torsion mod l, cache]
mw_prime(S, l) =
{
  my([ok, r, j, t0, W, E, Em, ch, gens, Tl, tests, F] = S);
  \\ l must be good for Em and the map constants and the tests
  if (Em.disc % l == 0, return(0));
  my(cs = concat([Vec(F), W, ch, [if (type(t0) == "t_VEC", 1, t0)], concat([Vec(tests[k]) | k <- [1..#tests]])]));
  for (k = 1, #cs, if (type(cs[k]) == "t_INT" || type(cs[k]) == "t_FRAC", if (valuation(denominator(cs[k]), l) > 0, return(0))));
  if (valuation(W[5], l) > 0 || valuation(ch[1], l) != 0, return(0));
  my(El = ellinit(Em, l), g = [mw_red(P, l) | P <- gens], o = [ellorder(El, P) | P <- g], T = [mw_red(P, l) | P <- Tl]);
  [l, El, g, o, T, Map()];
}

\\ pass test of a point R of El
mw_rpass(S, D, R) =
{
  my([ok, r, j, t0, W, E, Em, ch, gens, Tl, tests, F] = S, l = D[1], key = if (R == [0], [0], lift(R)), v);
  if (mapisdefined(D[6], key, &v), return(v));
  if (R == [0],
    v = mw_tpass(tests, if (type(t0) == "t_VEC", "inf", t0 * Mod(1, l)), l),
    my(X = ellchangepointinv(R, ch * Mod(1, l)));
    if (X[2] == 0,
      \\ y = 0: u is a root of (-q + u(ux-d)/(2q))^2 - F(u), or infinity; pass if any candidate passes
      \\ (-q + u(ux-d)/(2q))^2 - F(u) = u^2 * h(u) identically; h has degree 2 and its roots are the u of
      \\ the two points of E with this x. Infinity is a candidate when the degree drops.
      my([a,b,c,d,q] = W * Mod(1, l), x = X[1], h = (x^2/(4*q^2) - a)*'u^2 + (-d*x/(2*q^2) - b)*'u + (-x - c + d^2/(4*q^2)), cand = List());
      if (h == 0, v = 1,
        if (poldegree(h, 'u) < 2, listput(cand, "inf"));
        if (poldegree(h, 'u) > 0, my(rt = polrootsmod(lift(h), l)); for (k = 1, #rt, listput(cand, rt[k])));
        v = 0;
        for (k = 1, #cand, my(uu = cand[k], tl);
          if (type(t0) == "t_VEC", tl = if (uu === "inf" || uu == 0, if (uu === "inf", 0 * Mod(1, l), "inf"), 1/uu),
            tl = if (uu === "inf", "inf", t0 * Mod(1, l) + uu));
          if (mw_tpass(tests, tl, l), v = 1; break))),
      my(uu = mw_inv(W * Mod(1, l), X));
      if (type(t0) == "t_VEC",
        v = if (uu == 0, mw_tpass(tests, "inf", l), mw_tpass(tests, 1/uu, l)),
        v = mw_tpass(tests, t0 + uu, l))));
  mapput(D[6], key, v); v;
}

\\ point of El for class (n vector, torsion index)
mw_img(D, n, ti) = { my(El = D[2], R = D[5][ti]); for (i = 1, #n, R = elladd(El, R, ellmul(El, D[3][i], n[i]))); R; }

mw_primes(S, Lmax) = { my(PD = List()); forprime (l = 3, Lmax, my(D = mw_prime(S, l)); if (D != 0, listput(PD, D))); Vec(PD); }

\\ the sieve (Lmax: bound for l, or a vector from mw_primes). Mchain: list of primes p; M grows by each p. L: primes l to consider.
mw_sieve(S, Lmax = 20000, Mchain = [2,2,2,3,2,5,3,7,2,11,13,5,17,19,3,23,2,29,31,7,37,41,43,47], verb = 0) =
{
  my(k = #S[9], nT = #S[10], PD = List(), used = Map(), surv, M = 1, nprime = 0);
  if (type(Lmax) == "t_INT", forprime (l = 3, Lmax, my(D = mw_prime(S, l)); if (D != 0, listput(PD, D))),
    \\ precomputed prime data: copy with fresh caches (the cache depends on the tests)
    for (i = 1, #Lmax, my(D = Lmax[i]); D[6] = Map(); listput(PD, D)));
  surv = List(); for (ti = 1, nT, listput(surv, [vector(k), ti]));
  \\ stage 0: classes of torsion only (M = 1) with primes where all generator orders are 1 (rare); go straight to lifting
  for (st = 1, #Mchain,
    my(p = Mchain[st], Mn = M * p, ns = List());
    if (k == 0 && st > 1, break);
    if (k > 0,
      \\ lift: n -> n + M*e, e in [0,p-1]^k
      forvec (e = vector(k, i, [0, p-1]),
        for (s = 1, #surv, listput(ns, [surv[s][1] + M * e, surv[s][2]])));
      surv = ns; M = Mn, M = 1);
    if (#surv > MW_CAP, MW_SURV = []; return([-1, M, nprime, #surv]));
    \\ apply primes whose orders divide M
    for (i = 1, #PD, my(D = PD[i]);
      if (mapisdefined(used, D[1]), next);
      if (#D[4] && vecmax([M % o | o <- D[4]]) != 0, next);
      mapput(used, D[1], 1); nprime++;
      my(ns2 = List()); for (s = 1, #surv, if (mw_rpass(S, D, mw_img(D, surv[s][1] % M, surv[s][2])), listput(ns2, surv[s])));
      surv = ns2;
      if (#surv == 0, return([1, M, nprime, 0])));
    if (verb, printf("  M=%d primes=%d surv=%d\n", M, nprime, #surv)));
  MW_SURV = Vec(surv);
  [0, M, nprime, #surv];
}

\\ positive control: replace the tests by one random quartic g with g(t1) = s^2 at a known rational
\\ point t1 of the quartic (finite, t1 != t0 allowed). A sound sieve must keep at least one class.
mw_control(S, t1, seed = 1) =
{
  setrand(seed);
  my(S2 = S, G = vector(2, i, my(k = sum(i = 0, 3, (random(21) - 10) * 't^i), s2 = (random(9) + 1)^2, g = s2 + ('t - t1) * k);
     if (poldegree(g, 't) < 4, g += ('t - t1)^4);
     \\ clear denominators (a square factor): the real tests are l-integral at the primes used, and so must the fake ones be
     g * denominator(content(g))^2));
  S2[11] = G; S2;
}

\\ driver for one circle: try pairs in order of increasing rank
mw_circle(r, Lmax = 20000, verb = 0) =
{
  my(Q = cq(r), info = List(), res);
  for (j = 1, 6,
    my(S = iferr(mw_setup(r, j), e, [0, "err"]));
    if (S == 0, listput(info, [j, "nopt"]); next);
    if (S[1] == 0, listput(info, [j, S[2]]); next);
    listput(info, [j, S[13]]);
    if (S[13] > 3, next);
    res = mw_sieve(S, Lmax, , verb);
    if (res[1] == 1, return([1, j, S[13], res, Vec(info)])));
  [0, Vec(info)];
}

\\ rational t of a rational point P of Em (list of 1 or 2 candidates when y = 0 in the Washington model)
mw_tof(S, P) =
{
  my([ok, r, j, t0, W, E, Em, ch, gens, Tl, tests, F] = S, us);
  if (P == [0], return([t0]));
  my(X = ellchangepointinv(P, ch), [a,b,c,d,q] = W);
  if (X[2] != 0, us = [mw_inv(W, X)],
    my(h = (X[1]^2/(4*q^2) - a)*'u^2 + (-d*X[1]/(2*q^2) - b)*'u + (-X[1] - c + d^2/(4*q^2)));
    us = [rt | rt <- nfroots(, h)]; if (poldegree(h,'u) < 2, us = concat(us, ["inf"])));
  [if (type(t0) == "t_VEC", if (u === "inf", 0, if (u == 0, [1,0], 1/u)), if (u === "inf", [1,0], t0 + u)) | u <- us];
}
\\ cross-check: for all class vectors with |n_i| <= nb and every torsion point, compute the rational point,
\\ its exact t, and test the implication (t passes the tests mod l) => (class passes at l) for each prime of PD.
\\ Also checks that each t lies on the quartic. Returns [#points, #prime checks, #violations].
mw_cross(S, PD, nb = 3) =
{
  my(k = #S[9], Em = S[7], qq = cq(S[2])[S[3]], np = 0, nc = 0, nv = 0);
  forvec (n = vector(k, i, [-nb, nb]), for (ti = 1, #S[10],
    my(P = S[10][ti]); for (i = 1, k, P = elladd(Em, P, ellmul(Em, S[9][i], n[i])));
    my(T = mw_tof(S, P), good = 0);
    for (m = 1, #T, my(tt = T[m]); if (type(tt) == "t_VEC" || issquare(subst(qq, 't, tt)), good = 1));
    if (!good, nv++; next);
    np++;
    for (i = 1, #PD, my(D = PD[i], l = D[1], cls = mw_rpass(S, D, mw_img(D, n, ti)), dir = 0);
      for (m = 1, #T, my(tt = T[m], tl);
        if (type(tt) == "t_VEC", tl = "inf", if (valuation(denominator(tt), l) > 0, tl = "inf", tl = tt * Mod(1, l)));
        if (mw_tpass(S[11], tl, l), dir = 1));
      nc++; if (dir && !cls, nv++))));
  [np, nc, nv];
}
