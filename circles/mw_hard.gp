\\ mw_hard.gp (mathematics-52): deeper search for generators on pairs 1, 2 of the circles in MWH_LIST.
read("circles.gp"); read("mwsieve.gp");
{foreach(MWH_LIST, r, for (j = 1, 2, my(E = ellminimalmodel(qjac(cq(r)[j])), T = getabstime(), R, ar);
  ar = iferr(ellanalyticrank(E), e, [-1]);
  write(MWH_OUT, r, " ", j, " analytic ", ar[1], " ", getabstime() - T, " ms");
  for (eff = 2, 6, R = ellrank(E, eff); write(MWH_OUT, r, " ", j, " effort ", eff, " ", R[1..3], " #gens ", #R[4], " ", getabstime() - T, " ms"); if (R[1] == R[2], break))))}
