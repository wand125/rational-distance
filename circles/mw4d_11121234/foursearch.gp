\\ foursearch.gp (mathematics-52, 2026-11-12): second 2-descent on a 2-covering y^2 = g(x) whose quartic splits over Q
\\ as g = q1 q2 (happens for curves with full rational 2-torsion). A rational point gives q1(x) = d w1^2, q2(x) = d w2^2
\\ for a squarefree d built from the primes of 2 lc(q1) lc(q2) Res(q1,q2). For each d whose conic q1 = d w^2 has a
\\ rational point (qfsolve), parametrize it by s (qfparam) and search the 4-covering v^2 = d q2(x(s), z(s)) with
\\ hyperellratpoints. Every point found is mapped back to x = x(s)/z(s) on the 2-covering; the caller checks it.
\\ This is a search (it only has to find a point of infinite order); nothing here is used as a proof of non-existence.
FS_DBG = 0; FS_MIN = 1; FS_DMAX = 0; FS_DSET = []; FS_DREL = 0;
fs_hom(q, X, Z) = { polcoef(q, 2) * X^2 + polcoef(q, 1) * X * Z + polcoef(q, 0) * Z^2; }
\\ GL2(Z) reduction of a binary quartic by a covariant positive definite quadratic form (sum over the roots of
\\ |x - a z|^2 / |Q'(a)|^(1/2), as in Stoll-Cremona), Gauss-reduced. Returns [Q(U.(x,1)) as a polynomial, U].
fs_reduce(Q) =
{
  my(R = polroots(Q), dQ = Q', A = 0., B = 0., C = 0.);
  for (i = 1, #R, my(a = R[i], w = abs(subst(dQ, variable(Q), a))^(-1/2)); A += w; B += -2 * w * real(a); C += w * abs(a)^2);
  \\ Gauss reduction of the real form A x^2 + B x z + C z^2 (x = U11 x' + U12 z', z = U21 x' + U22 z')
  my(U = matid(2), it = 0);
  while (it++ < 10000,
    my(m = round(-B / (2*A)));
    if (m != 0, U = U * [1, m; 0, 1]; C = A*m^2 + B*m + C; B = B + 2*A*m);
    if (A > C, U = U * [0, -1; 1, 0]; [A, C] = [C, A]; B = -B, break));
  my(v = variable(Q), Xs = U[1,1]*'x + U[1,2], Zs = U[2,1]*'x + U[2,2]);
  my(Qr = sum(k = 0, 4, polcoef(Q, k, v) * Xs^k * Zs^(4-k)));
  [Qr, U];
}
\\ returns the list of x (rational, or [1,0] for infinity) on y^2 = q1 q2 found from 4-covering points of height <= B
fs_points(q1, q2, B) =
{
  my(P = select(p -> p > 0, factor(2 * pollead(q1) * pollead(q2) * polresultant(q1, q2))[,1]), L = List(), ds = List(), cand = List());
  forsubset(#P, S, my(d = prod(i = 1, #S, P[S[i]])); listput(ds, d); listput(ds, -d));
  if (#FS_DSET, ds = List(FS_DSET));
  for (k = 1, #ds, my(d = ds[k], G = [polcoef(q1,2), polcoef(q1,1)/2, 0; polcoef(q1,1)/2, polcoef(q1,0), 0; 0, 0, -d], sol = qfsolve(G));
    if (type(sol) != "t_COL", next);
    my(M = qfparam(G, sol), X = M[1,1]*'s^2 + M[1,2]*'s + M[1,3], Z = M[2,1]*'s^2 + M[2,2]*'s + M[2,3]);
    my(Q4 = d * fs_hom(q2, X, Z));
    if (poldegree(Q4, 's) < 3 || poldisc(Q4) == 0, next);   \\ (a square factor: not a genus-1 covering)
    \\ integral, square-free content, then GL2(Z)-reduce (fs_reduce) so small points stay small
    my(Qi = Q4 * denominator(content(Q4))^2, c = content(Qi), Qr, M2, pts);
    Qi = Qi / (c / core(c));
    \\ minimize (hyperellminimalmodel: x = (M11 x' + M12)/(M21 x' + M22), new form y'^2 + H y' = P' -> v^2 = 4P' + H^2),
    \\ then reduce; M2 is the composite map from the searched variable back to s
    my(Mm = matid(2));
    if (FS_MIN, my(mm, PH = iferr(hyperellminimalmodel(Qi, &mm), E, 0));
      if (PH != 0, my(Qm = if (type(PH) == "t_VEC", 4*PH[1] + PH[2]^2, PH)); Qm = Qm / (content(Qm) / core(content(Qm)));
        if (poldegree(Qm) >= 3, Qi = Qm; Mm = mm[2])));
    my(RU = iferr(fs_reduce(Qi), E, 0)); if (RU == 0, next); Qr = RU[1]; M2 = Mm * RU[2];
    if (FS_DBG, print("   d=", d, " log10|disc| ", round(log(abs(poldisc(Qr)))/log(10)), " log10 max|coef| ", round(log(vecmax(abs(Vec(Qr))))/log(10))));
    if (FS_DMAX && abs(poldisc(Qr)) > FS_DMAX, next);
    listput(cand, [abs(poldisc(Qr)), d, Qr, M2, M, X, Z, Q4]));
  \\ second pass: smallest discriminants first; with FS_DREL, only those within that factor of the smallest
  \\ (a locally soluble 4-covering minimizes to the level of the curve; the much larger ones are not locally soluble)
  cand = vecsort(Vec(cand), 1);
  for (k = 1, #cand, my([D4, d, Qr, M2, M, X, Z, Q4] = cand[k], pts);
    if (FS_DREL && D4 > FS_DREL * cand[1][1], break);
    my(T1 = getabstime()); pts = hyperellratpoints(Qr, B);
    if (FS_DBG, print("   d=", d, " points ", #pts, " ms ", getabstime() - T1, " Qr ", Qr));
    \\ (no closure here: a GP closure captures L by value, so listput inside it would be lost)
    my(svs = List());
    for (n = 1, #pts, my(sp = pts[n][1], den = M2[2,1] * sp + M2[2,2]);
      listput(svs, if (den == 0, [1, 0], (M2[1,1] * sp + M2[1,2]) / den)));
    if (issquare(pollead(Qr)), listput(svs, if (M2[2,1] == 0, [1, 0], M2[1,1] / M2[2,1])));
    for (n = 1, #svs, my(sv = svs[n]);
      if (type(sv) == "t_VEC", listput(L, if (M[2,1] != 0, M[1,1] / M[2,1], [1, 0])); next);
      if (!issquare(subst(Q4, 's, sv)), error("fs_points: mapped s not on the 4-cover"));
      my(xv = subst(X, 's, sv), zv = subst(Z, 's, sv)); listput(L, if (zv == 0, [1, 0], xv / zv))));
  Vec(L);
}
