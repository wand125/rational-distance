\\ circ_stage2.gp (mathematics-52, 2026-11-07): second stage for the circles left by circ_cdesc3.
\\ Set LIST (vector of r), NSH, SH, OUT before reading. Resumable. One line per circle:
\\   r rank0 j          (quotient j has rank 0 and its complete point list has no 4-distance point)
\\   r qi [a,b] ms      (the refined descent over Q(i) on the factor pair [a,b] has no solution, so that quotient has no point)
\\   r open ranks ms    (nothing worked)
\\   r SOL ...          (a 4-distance point; should never happen)
\\ Pairs: [1,j] is (1+t^2)N_j (quotient j-1 in cq), [a,b] with a,b >= 2 is N_a N_b. Rank-0 quotients are tried first.
read("circles.gp"); read("qidesc.gp");
ST2_SEC = 120;   \\ seconds per pair; a pair that runs out of time is skipped (recorded as timeouts=)
ST2_PAIRS = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];   \\ in the order of cq(r): A2, A3, A4, B23, B24, B34
{
  my(seen = Set());
  iferr(my(L0 = readstr(OUT)); seen = Set(apply(s -> strsplit(s, " ")[1], select(s -> s != "DONE" && #s, L0))), E, );
  for (k = 1, #LIST, if ((k - 1) % NSH == SH && !setsearch(seen, Str(LIST[k])),
    my(r = LIST[k], T = getwalltime(), v = circ_verdict(r, 300), line);
    if (v[1] == -1, line = Str(r, " SOL ", v),
      if (v[1] == 1, line = Str(r, " rank0 ", v[2]),
        my(rk = v[3], ord = vecsort(vector(6, j, [rk[j], j])), done = 0);
        my(nto = 0);
        for (o = 1, 6, my(j = ord[o][2], a = alarm(ST2_SEC, qi_desc(r, ST2_PAIRS[j])));
          if (type(a) == "t_ERROR", nto++; next);
          if (a[1], line = Str(r, " qi ", ST2_PAIRS[j], " ", getwalltime() - T); done = 1; break));
        if (!done, line = Str(r, " open ", rk, " ", getwalltime() - T, if (nto, Str(" timeouts=", nto), "")))));
    write(OUT, line); if (#line && strsplit(line, " ")[2] != "rank0" && strsplit(line, " ")[2] != "qi", print(line))));
  write(OUT, "DONE");
}
