\\ circ_ptsearch.gp (mathematics-52, 2026-11-09): point search on open circles (a search, not a proof).
\\ All six genus-1 quotients share the parameter t, so a 4-distance point with t = a/b, max(|a|, |b|) <= B appears among the
\\ rational points of any one of them. For each circle take the quotient with the fewest points up to 10^4, list its points
\\ with hyperellratpoints(F, B), and test every t exactly (circ_pt). t = oo is P = (-r, 0), checked too.
\\ Set LIST, NSH, SH, OUT, B before reading. Resumable. Line: r j #points SOL-or-none ms
read("circles.gp");
{
  my(seen = Set());
  iferr(my(L0 = readstr(OUT)); seen = Set(apply(s -> strsplit(s, " ")[1], select(s -> s != "DONE" && #s, L0))), E, );
  for (k = 1, #LIST, if ((k - 1) % NSH == SH && !setsearch(seen, Str(LIST[k])),
    my(r = LIST[k], T = getwalltime(), Q = cq(r), best = 0, bn = oo, F, P, sol = List());
    for (j = 1, 6, my(G = subst(Q[j], t, 'x)); G *= denominator(content(G))^2; my(n = #hyperellratpoints(G, 10^4)); if (n < bn, bn = n; best = j));
    F = subst(Q[best], t, 'x); F *= denominator(content(F))^2;
    P = hyperellratpoints(F, B);
    for (i = 1, #P, if (circ_pt(r, P[i][1]), listput(sol, P[i][1])));
    if (circ_pt(r, [1, 0]), listput(sol, "oo"));
    write(OUT, Str(r, " ", best, " ", #P, " ", if (#sol, Str("SOL ", Vec(sol)), "none"), " ", getwalltime() - T));
    if (#sol, print("SOLUTION ", r, " ", Vec(sol)))));
  write(OUT, "DONE");
}
