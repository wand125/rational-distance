\\ qidesc.gp  (mathematics-52, 2026-11-07)
\\ Refined descent on the circle |P| = r over K = Q(i). Read circles.gp first.
\\ fs = [f_1 .. f_m] (here 1+t^2, N2, N3, N4, or a subset), f_j = c_j (t - al_j)(t - conj al_j), al_j in K.
\\ At a rational point t = a/b, the class of t - al_j in K^*/K^*2 is b * delta_j with delta_j in K(S,2) and b in Q^*.
\\ Unknowns: the exponent vectors of delta_1..delta_m on the basis [i, generators of the primes of K above S] (e-vars),
\\ and, for each rational prime p in S, the class of the common rational factor in Q_p^*/Q_p^*2 (lambda-vars).
\\ Conditions:
\\   (global) c_j N(delta_j) = c_1 N(delta_1) in Q^*/Q^*2 (both are d, the common class of the values f_j(t));
\\   (local, p in S) (bits of delta_j at each prime P | p) + lambda_p = bits of (t - al_j) for some t in P^1(Q_p).
\\ Local images are over-approximated (unknown or root discs are "free"), so "no solution" is a proof.
\\ Local classes in K_P^*/K_P^*2 are coordinatised by Hilbert symbols against a local basis W_P (nondegenerate pairing).

QI_K = nfinit(y^2 + 1); QI_B = bnfinit(y^2 + 1);

\\ local basis of K_P^*/K_P^*2 by greedy rank of the Hilbert pairing Gram matrix
\\ local basis of K_P^*/K_P^*2: columns of the Hilbert Gram matrix of small elements (they span), pick a column basis.
\\ Since the pairing is nondegenerate, the picked elements span K_P^*/K_P^*2 and z -> ((z, w))_w is injective.
QI_LB = Map();
qi_locbasis(P) =
{
  my(key = [P.p, P.gen]);
  if (mapisdefined(QI_LB, key), return(mapget(QI_LB, key)));
  my(p = P.p, dim = if (p == 2, 4, 2), C = List());
  forvec(v = [[-3, 3], [-3, 3]], if (v != [0, 0], listput(C, v[1] + v[2] * y)));
  listput(C, p); listput(C, p * y); listput(C, p + y); listput(C, 1 + p * y);
  \\ make sure a unit of non-square norm (inert p) and a rational non-residue (split p) are among the candidates
  if (p > 2, listput(C, qi_nonres(p)); my(x0 = 1); while (kronecker(x0^2 + 1, p) != -1 && x0 < 10 * p, x0++); listput(C, x0 + y));
  C = Vec(C);
  my(G = matrix(#C, #C, a, b, nfhilbert(QI_K, C[a], C[b], P) == -1), Gm = Mod(G, 2), sel = List(), cur = matrix(#C, 0));
  if (matrank(Gm) != dim, error("qi_locbasis: rank ", matrank(Gm), " at ", P));
  for (b = 1, #C, my(c2 = concat(cur, Gm[, b])); if (matrank(c2) > matrank(cur), cur = c2; listput(sel, C[b])); if (#sel == dim, break));
  my(W = Vec(sel)); mapput(QI_LB, key, W); W;
}
qi_nonres(p) = { my(u = 2); while (kronecker(u, p) != -1, u++); u; }
qi_bits(z, P, W) = vector(#W, k, nfhilbert(QI_K, z, W[k], P) == -1);

\\ local pieces at p: vectors over components (j, P) of bit vectors, or "F" (free: any class)
\\ comps: list of [j, P, W]; al: roots; mode 0: t in a + p^k Z_p; mode 1: t = 1/s, s in a + p^k Z_p (k >= 1)
\\ large odd p: only the children containing a root are explored; all other children are summarised by one "generic" piece
\\ in which a root component is ["U", e] (z = e * unit, unit class arbitrary) and a non-constant component is free.
QI_PBIG = 400;
QI_RHO = Map();
\\ image of i in Z_p for a prime P of degree 1 above p (p-adic, precision QI_PREC)
QI_PREC = 30;
qi_rho(P) =
{
  my(key = [P.p, P.gen]); if (mapisdefined(QI_RHO, key), return(mapget(QI_RHO, key)));
  my(p = P.p, R = polrootspadic('x^2 + 1, p, QI_PREC), rho = 0);
  for (k = 1, #R, if (nfeltval(QI_K, y - truncate(R[k]), P) >= QI_PREC - 2, rho = R[k]));
  if (rho == 0, error("qi_rho ", P)); mapput(QI_RHO, key, rho); rho;
}
qi_img(z, P) = subst(lift(z), y, qi_rho(P));
qi_disc(al, comps, p, a, k, mode, dep, ~out) =
{
  my(nc = #comps, fx = vector(nc), unk = 0, nroot = 0, isr = vector(nc));
  for (c = 1, nc, my(j = comps[c][1], P = comps[c][2], W = comps[c][3], z, vz, vv, isroot = 0);
    if (mode == 0, z = a - al[j], z = 1 - al[j] * a);
    \\ root containment (only for P of degree 1, unramified: K_P = Q_p)
    if (P.f == 1 && P.e == 1,
      my(rt = if (mode == 0, al[j], 1 / al[j]));
      if (nfeltval(QI_K, a - rt, P) >= k, isroot = 1));
    if (isroot, nroot++; isr[c] = 1; fx[c] = "F"; next);
    vz = nfeltval(QI_K, z, P);
    vv = if (mode == 0, k * P.e, k * P.e + nfeltval(QI_K, al[j], P));
    if (vv > vz + 2 * nfeltval(QI_K, 2, P), fx[c] = qi_bits(z, P, W), unk = 1; fx[c] = -1));
  if ((unk || nroot > 1) && dep > 0,
    if (p <= QI_PBIG || p == 2,
      for (r = 0, p - 1, qi_disc(al, comps, p, a + r * p^k, k + 1, mode, dep - 1, ~out)); return);
    \\ large p: children holding a root, then one generic piece for the rest
    my(rs = Set());
    for (c = 1, nc, if (isr[c], my(j = comps[c][1], P = comps[c][2], rt = if (mode == 0, al[j], 1 / al[j]), im = qi_img(rt, P));
      rs = setunion(rs, [lift(Mod(truncate((im - a) / p^k), p))])));
    foreach(rs, r, qi_disc(al, comps, p, a + r * p^k, k + 1, mode, dep - 1, ~out));
    my(g = vector(nc));
    for (c = 1, nc, my(j = comps[c][1]);
      if (isr[c], g[c] = ["U", if (mode == 0, p^k, -al[j] * p^k)],
        g[c] = if (type(fx[c]) == "t_INT", "F", fx[c])));
    mapput(out, g, 1); return);
  for (c = 1, nc, if (type(fx[c]) == "t_INT", fx[c] = "F"));
  mapput(out, fx, 1);
}
qi_local(al, comps, p, dep = 10) =
{
  my(out = Map());
  qi_disc(al, comps, p, 0, 0, 0, dep, ~out);
  qi_disc(al, comps, p, 0, 1, 1, dep, ~out);
  Mat(out)[, 1];
}

\\ F_2 echelon (same as ech_add in circles.gp, repeated so this file is self-contained after circles.gp)
qi_add(~B, row, rhs) =
{
  while (row, my(pv = #binary(row) - 1);
    if (mapisdefined(B, pv), my(z = mapget(B, pv)); row = bitxor(row, z[1]); rhs = bitxor(rhs, z[2]),
      mapput(B, pv, [row, rhs]); return(1)));
  rhs == 0;
}
qi_dfs(eqs, k, B) =
{
  if (k > #eqs, return(1));
  for (g = 1, #eqs[k], my(B2 = B, ok = 1);
    for (b = 1, #eqs[k][g], if (!qi_add(~B2, eqs[k][g][b][1], eqs[k][g][b][2]), ok = 0; break));
    if (ok && qi_dfs(eqs, k + 1, B2), return(1)));
  0;
}

\\ main: returns [status (1 = no rational t), #S, #pieces per prime]
qi_desc(r, idx = [1, 2, 3, 4], pmax = 400, dep = 10) =
{
  my(fs = vecextract(circ_fs(r), idx), m = #fs, S = circ_S(fs), al = vector(m), cj = vector(m));
  for (j = 1, m, my(F = subst(fs[j], t, 'x), R = nfroots(QI_K, F)); al[j] = Mod(lift(R[1]), y^2 + 1); cj[j] = pollead(F);
    S = setunion(S, select(p -> p > 1, factor(poldisc(F * denominator(content(F))^2))[, 1]~));
    S = setunion(S, select(p -> p > 1, factor(denominator(content(lift(al[j]))) * numerator(cj[j]) * denominator(cj[j]))[, 1]~)));
  \\ global basis of K(S,2): i and a generator of each prime above S
  my(gb = List([y]), gN = List([1]));
  for (k = 1, #S, foreach(idealprimedec(QI_K, S[k]), P,
    \\ K has class number 1: a generator of P
    my(g = nfbasistoalg(QI_K, bnfisprincipal(QI_B, P)[2]));
    listput(gb, lift(g)); listput(gN, norm(Mod(lift(g), y^2 + 1)))));
  my(n = #gb, St = select(p -> p <= pmax, S), nE = m * n, eqs = List(), npc = List());
  \\ global norm equations: c_j N(delta_j) / (c_1 N(delta_1)) is a square: one row per prime of S
  my(geq = List());
  for (j = 2, m, for (k = 1, #S, my(p = S[k], row = 0);
    for (b = 1, n, if (valuation(gN[b], p) % 2, row = bitxor(row, 2^((j - 1) * n + b - 1)); row = bitxor(row, 2^(b - 1))));
    listput(geq, [row, valuation(cj[j] / cj[1], p) % 2])));
  \\ (sign: norms and c_j are positive, nothing to add)
  my(B = Map()); for (k = 1, #geq, if (!qi_add(~B, geq[k][1], geq[k][2]), return([1, #S, "global"])));
  my(off = nE);
  for (kp = 1, #St, my(p = St[kp], Ps = idealprimedec(QI_K, p), comps = List(), lam = if (p == 2, [2, -1, 5], [p, qi_nonres(p)]));
    foreach(Ps, P, my(W = qi_locbasis(P)); for (j = 1, m, listput(comps, [j, P, W])));
    comps = Vec(comps);
    my(pcs = qi_local(al, comps, p, dep), opts = List());
    \\ rows for each component bit: e-vars of delta_j and the lambda-vars of p
    my(rowsc = vector(#comps, c, my(j = comps[c][1], P = comps[c][2], W = comps[c][3], gbits = vector(n, b, qi_bits(Mod(gb[b], y^2 + 1), P, W)),
         lbits = vector(#lam, l, qi_bits(lam[l], P, W)));
       vector(#W, w, sum(b = 1, n, gbits[b][w] * 2^((j - 1) * n + b - 1)) + sum(l = 1, #lam, lbits[l][w] * 2^(off + l - 1)))));
    for (g = 1, #pcs, my(fx = pcs[g], E = List());
      for (c = 1, #comps,
        if (type(fx[c]) == "t_STR", next);
        if (type(fx[c][1]) == "t_STR",
          \\ ["U", e]: (z, u) = (e, u) for a non-square unit u, i.e. sum_w uc[w] * bit_w(z) = [(e, u) = -1]
          my(P = comps[c][2], W = comps[c][3], u = qi_nonres(p), G = matrix(#W, #W, a, b, nfhilbert(QI_K, W[a], W[b], P) == -1),
             uc = lift(matsolve(Mod(G, 2), Mod(qi_bits(u, P, W)~, 2))), row = 0);
          for (w = 1, #W, if (uc[w], row = bitxor(row, rowsc[c][w])));
          listput(E, [row, nfhilbert(QI_K, fx[c][2], u, P) == -1]),
          for (w = 1, #fx[c], listput(E, [rowsc[c][w], fx[c][w]]))));
      listput(opts, Vec(E)));
    off += #lam; listput(eqs, Vec(opts)); listput(npc, [p, #pcs]));
  \\ eliminate lambda-vars and dedupe per prime, then search with forward checking
  my(C = vector(#eqs, k, my(U = Map()); foreach(eqs[k], E, my(c = qi_canon(E, 2^nE)); if (c[1] == 1, mapput(U, c[2], 1))); if (#U, Mat(U)[, 1]~, [])));
  for (k = 1, #C, if (#C[k] == 0, return([1, #S, Vec(npc), "local"])));
  [!qi_dfs2(C, B), #S, Vec(npc), vector(#C, k, #C[k])];
}

\\ ---- faster search: eliminate the local lambda-vars, dedupe, forward checking with fewest-candidates-first ----------
\\ a piece is a list of [row, rhs]; return [1, reduced equations in the e-vars only (bits < nE), canonical], or [0] if inconsistent.
\\ (An empty equation list means "no condition"; it must not be confused with inconsistency: in gp, [] == 0 is true.)
qi_canon(E, nE) =
{
  my(B = Map());
  for (k = 1, #E, if (!qi_add(~B, E[k][1], E[k][2]), return([0])));
  if (#B == 0, return([1, []]));
  my(M = Mat(B), rows = List());
  \\ keep rows whose pivot is an e-var; they contain no lambda bits (pivot = highest bit)
  for (k = 1, #M~, if (M[k, 1] < nE, listput(rows, M[k, 2])));
  \\ reduced echelon among themselves, sorted by pivot
  my(R = vecsort(Vec(rows), (a, b) -> cmp(#binary(b[1]), #binary(a[1]))));
  for (a = 1, #R, my(pa = #binary(R[a][1]) - 1); for (b = 1, #R, if (b != a && bittest(R[b][1], pa), R[b] = [bitxor(R[b][1], R[a][1]), bitxor(R[b][2], R[a][2])])));
  [1, vecsort(R)];
}
qi_ok(B, E) = { my(B2 = B); for (k = 1, #E, if (!qi_add(~B2, E[k][1], E[k][2]), return(0))); 1; }
qi_dfs2(C, B) =
{
  \\ C: list of candidate lists (one per prime not yet chosen); filter, pick the smallest
  if (#C == 0, return(1));
  my(F = vector(#C), best = 0, bn = oo);
  for (k = 1, #C, F[k] = select(E -> qi_ok(B, E), C[k]); if (#F[k] == 0, return(0)); if (#F[k] < bn, bn = #F[k]; best = k));
  my(rest = vector(#C - 1, k, F[if (k < best, k, k + 1)]));
  for (g = 1, #F[best], my(B2 = B, E = F[best][g]); for (b = 1, #E, qi_add(~B2, E[b][1], E[b][2])); if (qi_dfs2(rest, B2), return(1)));
  0;
}
