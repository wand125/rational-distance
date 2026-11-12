\\ circles.gp  (mathematics-52, 2026-11-07)
\\ Points P with rational distances to the four vertices of the unit square, on the circle |P| = r (r rational):
\\ every such point lies on the circle of radius r = |P| about (0,0).
\\ P = r ((1 - t^2), 2t) / (1 + t^2). With N2 = (r-1)^2 + (r+1)^2 t^2, N3 = (r^2+1)(1+t^2) - 4 r t,
\\ N4 = (r^2+2)(1+t^2) - 2r(1-t^2) - 4 r t, the squared distances to (1,0), (0,1), (1,1) are N_j / (1+t^2), so the
\\ conditions are: A_j = (1+t^2) N_j is a square, j = 2, 3, 4. Every root lies in Q(i).
\\ Genus-1 quotients: A2, A3, A4 and B23 = N2 N3, B24 = N2 N4, B34 = N3 N4 (quartics in t).
\\ Jacobian of y^2 = a t^4 + b t^3 + c t^2 + d t + e: E : Y^2 = X^3 - 27 I X - 27 J with
\\ I = 12ae - 3bd + c^2, J = 72ace + 9bcd - 27ad^2 - 27eb^2 - 2c^3.

cq(r) = { my(N2 = (r-1)^2 + (r+1)^2 * t^2, N3 = (r^2+1)*(1+t^2) - 4*r*t, N4 = (r^2+2)*(1+t^2) - 2*r*(1-t^2) - 4*r*t);
  [(1+t^2)*N2, (1+t^2)*N3, (1+t^2)*N4, N2*N3, N2*N4, N3*N4]; }
qjac(q) = { my(a = polcoef(q,4,t), b = polcoef(q,3,t), c = polcoef(q,2,t), d = polcoef(q,1,t), e = polcoef(q,0,t),
  I = 12*a*e - 3*b*d + c^2, J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3);
  ellinit([0, 0, 0, -27*I, -27*J]); }
\\ per circle: for each of the 6 quartics [rank lo, rank hi, #small points (t of height <= B, both infinity points counted
\\ when the leading coefficient is a square), #E_tors]
circ_info(r, B = 200) =
{
  my(Q = cq(r));
  vector(#Q, j, my(q = Q[j], E = qjac(q), Em = ellminimalmodel(E), R = ellrank(Em), D = denominator(content(q)), F = D^2 * subst(q, t, 'x),
     P = hyperellratpoints(F, B), ninf = if (issquare(polcoef(q, 4, t)), 2, 0));
    [R[1], R[2], #P + ninf, elltors(Em)[1]]);
}

\\ exact test of a parameter t (t = oo allowed as [1,0]) on circle r: are all 4 distances rational?
circ_pt(r, t) = { my(P = if (type(t) == "t_VEC", [-r, 0], r * [1 - t^2, 2*t] / (1 + t^2)), x = P[1], y = P[2]);
  issquare(x^2 + y^2) && issquare((x-1)^2 + y^2) && issquare(x^2 + (y-1)^2) && issquare((x-1)^2 + (y-1)^2); }
\\ verdict for circle r: [status, j, ranks]; status 1 = excluded by a rank-0 genus-1 quotient whose point list is complete
\\ (#affine points + #infinity points = #E(Q)_tors > 0) and no listed point is a 4-distance point; 0 = open; -1 = SOLUTION
circ_verdict(r, B = 200) =
{
  my(Q = cq(r), rk = vector(6), best = 0);
  for (j = 1, 6,
    my(q = Q[j], Em = ellminimalmodel(qjac(q)), R = ellrank(Em)); rk[j] = R[2];
    if (R[2] == 0 && !best,
      my(D = denominator(content(q)), F = D^2 * subst(q, t, 'x), P = hyperellratpoints(F, B),
         ninf = if (issquare(polcoef(q, 4, t)), 2, 0), nt = elltors(Em)[1]);
      if (#P + ninf == nt && nt > 0,
        for (k = 1, #P, if (circ_pt(r, P[k][1]), return([-1, j, rk, P[k][1]])));
        if (ninf && circ_pt(r, [1, 0]), return([-1, j, rk, "oo"]));
        best = j)));
  [best > 0, best, rk];
}

\\ ---- local solubility of y^2 = f(x), f in Z[x] squarefree -------------------------------------------
\\ ls_disc(f, p, a, k, d): is there x in a + p^k Z_p with f(x) a nonzero square in Q_p (or f(x) = 0)?
\\ returns 1 (yes), 0 (no), -1 (undecided at depth limit d)
ls_disc(f, p, a, k, d) =
{
  my(h = subst(f, 'x, a + p^k * 'x), v0 = if (polcoef(h, 0) == 0, oo, valuation(polcoef(h, 0), p)), m = oo, need = if (p == 2, 3, 1));
  for (i = 1, poldegree(h), my(c = polcoef(h, i)); if (c, m = min(m, valuation(c, p))));
  if (v0 == oo, return(1));
  if (v0 < m,
    if (v0 % 2, return(0));
    if (m - v0 >= need, my(u = polcoef(h, 0) / p^v0); return(if (p == 2, u % 8 == 1, kronecker(u, p) == 1))));
  \\ Hensel: v(f(a)) > 2 v(f'(a)) gives a root in the disc (y = 0)
  my(a1 = polcoef(h, 1)); if (a1 && v0 > 2 * valuation(a1, p), return(1));
  if (d == 0, return(-1));
  my(und = 0); for (r = 0, p - 1, my(s = ls_disc(f, p, a + r * p^k, k + 1, d - 1)); if (s == 1, return(1)); if (s == -1, und = 1));
  -und;
}
\\ y^2 = q(t) (degree 4, rational coefficients) is soluble over Q_p, including t = infinity?
ls_quartic(q, p, d = 14) =
{
  my(F = subst(q, t, 'x), D = denominator(content(F)), G); F = F * D^2; G = polrecip(F);
  my(s1 = ls_disc(F, p, 0, 0, d)); if (s1 == 1, return(1));
  my(s2 = ls_disc(G, p, 0, 1, d)); if (s2 == 1, return(1));
  if (s1 == -1 || s2 == -1, -1, 0);
}
\\ primes where an integral model of q can fail to have Q_p points
ls_badprimes(q) = { my(F = subst(q, t, 'x)); F = F * denominator(content(F))^2; F = F / content(F);
  select(p -> p > 1, factor(2 * poldisc(F) * pollead(F))[, 1]~); }
\\ 1 = everywhere locally soluble, 0 = not (with the prime), -1 = undecided
ls_els(q) = { my(ps = ls_badprimes(q), und = 0); for (k = 1, #ps, my(s = ls_quartic(q, ps[k])); if (s == 0, return([0, ps[k]])); if (s == -1, und = 1)); [if (und, -1, 1), 0]; }

\\ verdict v2: [status, how, detail, ranks]; how = "nols" (quotient j has no Q_p points), "rank0" (quotient j rank 0, list complete)
circ_verdict2(r, B = 200, pmax = 2000) =
{
  my(Q = cq(r));
  for (j = 1, 6, my(ps = select(p -> p <= pmax, ls_badprimes(Q[j])));
    for (k = 1, #ps, if (ls_quartic(Q[j], ps[k]) == 0, return([1, "nols", [j, ps[k]], []]))));
  my(v = circ_verdict(r, B)); if (v[1] == -1, return([-1, "SOL", v[4], v[3]]));
  if (v[1] == 1, [1, "rank0", [v[2]], v[3]], [0, "open", [], v[3]]);
}

\\ ---- descent by the common square class --------------------------------------------------------------
\\ At a 4-distance point (t rational, or t = oo), 1+t^2, N2, N3, N4 are all d * (rational square) for one squarefree
\\ d > 0 (all four are positive definite). So for each such d there must be t in P^1(Q_p), for every p, with
\\ d f(t) a square in Q_p for all f in fs. Only "no" answers (0) are used as proofs; "undecided" counts as soluble.
\\ lm_disc(fs, d, p, a, k, dep): 1 = all four determined and square somewhere, 0 = impossible on a + p^k Z_p, -1 = undecided
lm_disc(fs, d, p, a, k, dep) =
{
  my(und = 0, need = if (p == 2, 3, 1));
  for (i = 1, #fs,
    my(h = subst(fs[i], 'x, a + p^k * 'x), c0 = polcoef(h, 0), m = oo);
    for (e = 1, poldegree(h), my(c = polcoef(h, e)); if (c, m = min(m, valuation(c, p))));
    if (c0 == 0, und = 1; next);
    my(v0 = valuation(c0, p));
    if (v0 < m && m - v0 >= need,
      my(z = d * c0, w = valuation(z, p), u = z / p^w);
      if (w % 2, return(0));
      if (if (p == 2, u % 8 != 1, kronecker(u, p) != 1), return(0)),
    if (v0 < m && (v0 + valuation(d, p)) % 2, return(0), und = 1)));
  if (!und, return(1));
  if (dep == 0, return(-1));
  my(r1 = 0);
  for (r = 0, p - 1, my(s = lm_disc(fs, d, p, a + r * p^k, k + 1, dep - 1)); if (s == 1, return(1)); if (s == -1, r1 = 1));
  -r1;
}
\\ fs: polynomials in t of even degree (2 here); t = oo uses the reversed polynomials on p Z_p
lm_p(fs, d, p, dep = 12) =
{
  my(F = vector(#fs, i, my(g = subst(fs[i], t, 'x)); g * denominator(content(g))^2),
     G = vector(#fs, i, polrecip(F[i]) * if (poldegree(F[i]) < 2, 'x^(2 - poldegree(F[i])), 1)));
  my(s1 = lm_disc(F, d, p, 0, 0, dep)); if (s1 == 1, return(1));
  my(s2 = lm_disc(G, d, p, 0, 1, dep)); if (s2 == 1, return(1));
  if (s1 == -1 || s2 == -1, -1, 0);
}
circ_fs(r) = [1 + t^2, (r-1)^2 + (r+1)^2 * t^2, (r^2+1)*(1+t^2) - 4*r*t, (r^2+2)*(1+t^2) - 2*r*(1-t^2) - 4*r*t];
\\ returns [1, list of [d, p]] if every d is killed at some prime (circle excluded), else [0, surviving d's]
circ_cdesc(r, pmax = 2000, idx = [1, 2, 3, 4]) =
{
  my(fs = vecextract(circ_fs(r), idx), S = [2], n = #fs);
  for (i = 1, n, for (j = i + 1, n, my(a = subst(fs[i], t, 'x), b = subst(fs[j], t, 'x));
    a *= denominator(content(a))^2; b *= denominator(content(b))^2; S = setunion(S, select(p -> p > 1, factor(polresultant(a, b))[, 1]~))));
  for (i = 1, n, my(a = subst(fs[i], t, 'x)); a *= denominator(content(a))^2; S = setunion(S, select(p -> p > 1, factor(pollead(a) * polcoef(a, 0))[, 1]~)));
  \\ d runs over all squarefree products of S; only the local tests are limited to p <= pmax
  my(kill = List(), live = List(), St = select(p -> p <= pmax, S));
  forvec(e = vector(#S, i, [0, 1]), my(d = prod(i = 1, #S, S[i]^e[i]), dead = 0);
    for (k = 1, #St, if (lm_p(fs, d, St[k]) == 0, listput(kill, [d, St[k]]); dead = 1; break));
    if (!dead, listput(live, d)));
  if (#live, [0, Vec(live), S], [1, Vec(kill), S]);
}

\\ ---- fast version: lm_p depends on d only through its class in Q_p^* / Q_p^{*2} ----------------------
lp_cls(d, p) = { my(v = valuation(d, p), u = d / p^v); if (p == 2, [v % 2, u % 8], [v % 2, kronecker(u, p)]); }
lp_reps(p) = { if (p == 2, return([1, 3, 5, 7, 2, 6, 10, 14])); my(n = 2); while (kronecker(n, p) != -1, n++); [1, n, p, p * n]; }
circ_S(fs) = { my(S = [2], n = #fs);
  for (i = 1, n, for (j = i + 1, n, my(a = subst(fs[i], t, 'x), b = subst(fs[j], t, 'x));
    a *= denominator(content(a))^2; b *= denominator(content(b))^2; S = setunion(S, select(p -> p > 1, factor(polresultant(a, b))[, 1]~))));
  for (i = 1, n, my(a = subst(fs[i], t, 'x)); a *= denominator(content(a))^2; S = setunion(S, select(p -> p > 1, factor(pollead(a) * polcoef(a, 0))[, 1]~)));
  S; }
\\ [status (1 = every d killed), #S, #alive d, alive d's (at most 20), kill counts per prime as [p, n]]
circ_cdesc2(r, pmax = 2000, idx = [1, 2, 3, 4]) =
{
  my(fs = vecextract(circ_fs(r), idx), S = circ_S(fs), St = select(p -> p <= pmax, S), bad = vector(#St), kc = vector(#St), nal = 0, al = List());
  for (k = 1, #St, my(p = St[k], R = lp_reps(p)); bad[k] = Set(apply(c -> lp_cls(c, p), select(c -> lm_p(fs, c, p) == 0, R)));
    if (#bad[k] == #R, kc[k] = 2^#S; return([1, #S, 0, [], vector(#St, j, [St[j], kc[j]])])));
  forvec(e = vector(#S, i, [0, 1]), my(d = prod(i = 1, #S, S[i]^e[i]), dead = 0);
    for (k = 1, #St, if (#bad[k] && setsearch(bad[k], lp_cls(d, St[k])), kc[k]++; dead = 1; break));
    if (!dead, nal++; if (#al < 20, listput(al, d))));
  [nal == 0, #S, nal, Vec(al), vector(#St, k, [St[k], kc[k]])];
}

\\ ---- linear-algebra version ----------------------------------------------------------------------------
\\ The class of d = prod S[i]^e[i] in Q_p^*/Q_p^{*2} is an F_2-linear function of e. Bits per prime:
\\ odd p: (v, s) with s = [Legendre of the unit part is -1]; p = 2: (v, a, b) with a = [u = 3 mod 4], b = [u = 3, 5 mod 8].
cl_bits(c, p) = { my(v = valuation(c, p), u = c / p^v);
  if (p == 2, [v % 2, u % 4 == 3, (u % 8 == 3) || (u % 8 == 5)], [v % 2, kronecker(u, p) == -1]); }
\\ echelon basis over F_2 as a Map pivot -> [row bitmask, rhs]; add an equation, return 0 if inconsistent
ech_add(~B, row, rhs) =
{
  while (row, my(pv = #binary(row) - 1);
    if (mapisdefined(B, pv), my(z = mapget(B, pv)); row = bitxor(row, z[1]); rhs = bitxor(rhs, z[2]),
      mapput(B, pv, [row, rhs]); return(1)));
  rhs == 0;
}
\\ DFS over primes, choosing a good class at each; returns the number of alive d
ld_count(rows, good, k, B, nS) =
{
  if (k > #rows, return(2^(nS - #B)));
  my(tot = 0);
  for (g = 1, #good[k], my(B2 = B, ok = 1, cls = good[k][g]);
    for (b = 1, #cls, if (!ech_add(~B2, rows[k][b], cls[b]), ok = 0; break));
    if (ok, tot += ld_count(rows, good, k + 1, B2, nS)));
  tot;
}
\\ same output shape as circ_cdesc2 except the per-prime kill counts (empty)
circ_cdesc3(r, pmax = 2000, idx = [1, 2, 3, 4]) =
{
  my(fs = vecextract(circ_fs(r), idx), S = circ_S(fs), St = select(p -> p <= pmax, S), rows = List(), good = List(), nb);
  for (k = 1, #St, my(p = St[k], R = lp_reps(p), gd = List());
    for (c = 1, #R, if (lm_p(fs, R[c], p) != 0, listput(gd, cl_bits(R[c], p))));
    if (#gd == 0, return([1, #S, 0, [], [[p, 2^#S]]]));
    nb = if (p == 2, 3, 2);
    \\ row b = bitmask over S of the contributions to bit b
    my(rw = vector(nb, b, sum(i = 1, #S, (cl_bits(S[i], p)[b] % 2) * 2^(i - 1))));
    listput(rows, rw); listput(good, Vec(gd)));
  my(n = ld_count(Vec(rows), Vec(good), 1, Map(), #S));
  [n == 0, #S, n, [], []];
}
