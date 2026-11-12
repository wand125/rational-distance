\\ qirank0.gp (mathematics-52, 2026-11-11): close a circle with a genus-1 quotient C: v^2 = f(t) whose Jacobian J has
\\ rank 0 over Q(i), i.e. rank J(Q) = rank J^(-1)(Q) = 0 (J^(-1) = quadratic twist by -1; the ranks are proved elsewhere,
\\ by L(J,1) != 0 and L(J^(-1),1) != 0 with error bounds, see qirank0_L.py).
\\ Every quartic of cq(r) has a root t0 in Q(i) (all roots of 1+t^2, N2, N3, N4 lie in Q(i)), so C is isomorphic over K = Q(i)
\\ to the cubic model E below, and C(K) <-> E(K) = E(K)_tors (finite). We list E(K)_tors (PARI elltors over a number field),
\\ map every point back to t, check that it lies on C, and keep the rational t: these are all of C(Q).
\\ Requires: read("circles.gp") first (defines cq, circ_pt and the variable t).
't; 'u;   \\ create t and u before the Q(i) variable
QIW = varlower("i_");   \\ the variable of Q(i) must have lower priority than t and u
QI_K = nfinit(QIW^2 + 1);
\\ returns [status, r, j, t0, #E(K)_tors, cyc, rational t list, solutions]; status 1 = no 4-distance point on this circle
qi_rank0(r, j) =
{
  my(f = cq(r)[j], D = denominator(content(f)), F = D^2 * f, rts = nfroots(QI_K, F), t0, G, c, E, T, pts, rat = List(), sol = List());
  if (#rts == 0, error("no root in Q(i)"));
  t0 = Mod(rts[1], QIW^2 + 1);
  \\ u^4 F(t0 + 1/u) = c3 u^3 + c2 u^2 + c1 u + c0 since F(t0) = 0
  G = 'u^4 * subst(F, t, t0 + 1/'u); G = simplify(lift(G)); G = Pol(G, 'u);
  c = vector(5, k, Mod(polcoef(G, k-1, 'u), QIW^2 + 1));   \\ c[k] = coefficient of u^(k-1)
  if (c[5] != 0, error("F(t0) != 0"));
  \\ w'^2 = c3 u^3 + c2 u^2 + c1 u + c0 with w' = v u^2; X = c3 u, Y = c3 w': Y^2 = X^3 + c2 X^2 + c1 c3 X + c0 c3^2
  E = ellinit([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2], QI_K);
  if (E.j != ellinit(qi_curves(r, j)[1]).j, error("j mismatch"));
  T = elltors(E);
  pts = List([[0]]);
  if (T[1] > 1,
    my(S = T[2], g = T[3]); pts = List();
    if (#S == 1, for (a = 0, S[1]-1, listput(pts, ellmul(E, g[1], a))),
      for (a = 0, S[1]-1, for (b = 0, S[2]-1, listput(pts, elladd(E, ellmul(E, g[1], a), ellmul(E, g[2], b)))))));
  for (k = 1, #pts, my(P = pts[k]);
    if (P == [0], next);                      \\ u = oo <-> t = t0, not rational (t0 is not real)
    my(X = Mod(lift(P[1]), QIW^2 + 1), Y = Mod(lift(P[2]), QIW^2 + 1), u = X / c[4]);
    if (u == 0,
      \\ t = oo: C has points at infinity iff the leading coefficient of F is a square in K; they are rational t = oo
      listput(rat, "oo"); if (circ_pt(r, [1, 0]), listput(sol, "oo")); next);
    my(tt = t0 + 1/u, ww = Y / c[4], vv = ww / u^2);
    if (vv^2 != subst(F, t, tt), error("point not on C"));
    my(tl = lift(tt));
    if (poldegree(tl, QIW) <= 0, my(tq = polcoef(tl, 0, QIW)); listput(rat, tq); if (circ_pt(r, tq), listput(sol, tq))));
  [#sol == 0, r, j, lift(t0), T[1], T[2], Vec(rat), Vec(sol)];
}
\\ the two curves whose ranks must be 0: J and its twist by -1, as [a4, a6] of short Weierstrass models
qi_curves(r, j) = { my(q = cq(r)[j], a = polcoef(q,4,t), b = polcoef(q,3,t), c = polcoef(q,2,t), d = polcoef(q,1,t), e = polcoef(q,0,t),
  I = 12*a*e - 3*b*d + c^2, J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3); [[-27*I, -27*J], [-27*I, 27*J]]; }

\\ ---- J(Q) of rank 0 only (the twist J^(-1) may have positive rank) ----
\\ C, t0, E, the iso phi: C -> E over K (phi(P0) = O, P0 = (t0, 0)) as above. Galois conjugation acts on E(K) through C:
\\ sig(Q) = phi(conj(phi^-1(Q))) = tau(Q) + c with tau a group automorphism (tau^2 = 1) and c = sig(O) = phi(conj P0).
\\ For P in C(Q) and Q = phi(P): sig(Q) = Q, so S = Q + tau(Q) is tau-fixed and 2Q = S + c (since Q - tau(Q) = c).
\\ The tau-fixed subgroup of E(K) is isomorphic to J(Q) (finite when rank J(Q) = 0, so inside E(K)_tors).
\\ So phi(C(Q)) is contained in { Q : 2Q = S + c, S in E(K)_tors with tau(S) = S }, a finite set that we list
\\ (ellisdivisible gives one half, the others differ by E(K)[2]); then keep the Q with sig(Q) = Q and map back to t.
qi_cj(z) = Mod(subst(lift(z), QIW, -QIW), QIW^2 + 1);
\\ point of E -> [t, v] on C (t = "oo" at infinity: then v stands for v/t^2)
qi_toC(Q, t0, c) = { if (Q == [0], return([t0, 0]));
  my(X = Mod(lift(Q[1]), QIW^2+1), Y = Mod(lift(Q[2]), QIW^2+1), u = X / c[4], ww = Y / c[4]);
  if (u == 0, ["oo", ww], [t0 + 1/u, ww / u^2]); }
qi_toE(P, t0, c) = { if (P[1] != "oo" && P[1] == t0 && P[2] == 0, return([0]));
  if (P[1] == "oo", return([0, c[4] * P[2]]));
  my(u = 1 / (P[1] - t0)); [c[4] * u, c[4] * P[2] * u^2]; }
qi_sig(Q, t0, c) = { my(P = qi_toC(Q, t0, c)); qi_toE([if (P[1] == "oo", "oo", qi_cj(P[1])), qi_cj(P[2])], t0, c); }
qi_peq(P, Q) = { if (P == [0] || Q == [0], return(P == Q)); Mod(lift(P[1]), QIW^2+1) == Mod(lift(Q[1]), QIW^2+1) && Mod(lift(P[2]), QIW^2+1) == Mod(lift(Q[2]), QIW^2+1); }
qi_rank0b(r, j) =
{
  my(f = cq(r)[j], D = denominator(content(f)), F = D^2 * f, rts = nfroots(QI_K, F), t0, G, c, E, T, tors, two, cc, rat = List(), sol = List(), nS = 0, nQ = 0);
  t0 = Mod(rts[1], QIW^2 + 1);
  G = 'u^4 * subst(F, t, t0 + 1/'u); G = simplify(lift(G)); G = Pol(G, 'u);
  c = vector(5, k, Mod(polcoef(G, k-1, 'u), QIW^2 + 1));
  if (c[5] != 0, error("F(t0) != 0"));
  E = ellinit([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2], QI_K);
  if (E.j != ellinit(qi_curves(r, j)[1]).j, error("j mismatch"));
  T = elltors(E); tors = List([[0]]);
  if (T[1] > 1, my(S = T[2], g = T[3]); tors = List();
    if (#S == 1, for (a = 0, S[1]-1, listput(tors, ellmul(E, g[1], a))),
      for (a = 0, S[1]-1, for (b = 0, S[2]-1, listput(tors, elladd(E, ellmul(E, g[1], a), ellmul(E, g[2], b)))))));
  two = [P | P <- Vec(tors), P == [0] || ellmul(E, P, 2) == [0]];
  cc = qi_sig([0], t0, c);
  if (!ellisoncurve(E, cc), error("c not on E"));
  for (k = 1, #tors, my(S = tors[k], tS = elladd(E, qi_sig(S, t0, c), ellneg(E, cc)));
    if (!qi_peq(tS, S), next); nS++;
    my(H, Z = elladd(E, S, cc));
    if (!ellisdivisible(E, Z, 2, &H), next);
    for (m = 1, #two, my(Q = elladd(E, H, two[m])); nQ++;
      if (!qi_peq(qi_sig(Q, t0, c), Q), next);
      my(P = qi_toC(Q, t0, c));
      if (P[1] == "oo", listput(rat, "oo"); if (circ_pt(r, [1, 0]), listput(sol, "oo")); next);
      my(tl = lift(P[1])); if (poldegree(tl, QIW) > 0, error("fixed point with t not rational"));
      my(tq = polcoef(tl, 0, QIW), vq = polcoef(lift(P[2]), 0, QIW));
      if (vq^2 != subst(F, t, tq), error("rational point not on C"));
      listput(rat, tq); if (circ_pt(r, tq), listput(sol, tq))));
  [#sol == 0, r, j, T[1], nS, nQ, Vec(rat), Vec(sol)];
}
