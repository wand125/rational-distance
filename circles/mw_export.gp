\\ mw_export.gp (mathematics-52): for each killed (r, j) in MWX_LIST, write one line
\\   r <tab> j <tab> [a1,a2,a3,a4,a6] <tab> gens <tab> ntors <tab> quartic coefficients [a,b,c,d,e]
\\ with the minimal model Em, the saturated generators used by the sieve, the torsion order, and the quartic
\\ (the second system recomputes its j-invariant from the invariants I, J and compares).
\\ Set MWX_LIST (vector of [r, j]) and MWX_OUT before reading.
read("circles.gp"); read("mwsieve.gp");
{for (i = 1, #MWX_LIST, my(r = MWX_LIST[i][1], j = MWX_LIST[i][2], S = mw_setup(r, j), Em = S[7], qq = cq(r)[j]);
  write(MWX_OUT, r, "\t", j, "\t", Vec(Em[1..5]), "\t", S[9], "\t", #S[10], "\t", [polcoef(qq, k, 't) | k <- [4,3,2,1,0]]));}
