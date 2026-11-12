\\ positive controls for qidesc.gp: real points with rational distance to (0,0) and to the vertices of the factor set;
\\ the refined descent on those factors must never say "no point" (pmax 400 and pmax 1e6)
read("circles.gp"); read("qidesc.gp");
brute(D, V) = { my(R = Map()); for (d = 1, D, for (a = -3*d, 4*d, for (b = -3*d, 4*d, if (gcd([a, b, d]) == 1,
  my(x = a/d, y0 = b/d, r2 = x^2 + y0^2, ok = r2 != 0 && issquare(r2));
  if (ok, for (k = 1, #V, if (!issquare((x - V[k][1])^2 + (y0 - V[k][2])^2), ok = 0; break)));
  if (ok, my(r = sqrtint(numerator(r2)) / sqrtint(denominator(r2))); if (r != 1, mapput(R, r, [x, y0]))))))); Mat(R); }
VX = [[0,0], [1,0], [0,1], [1,1]];   \\ factor j <-> vertex VX[j] (factor 1 = 1+t^2 is automatic)
test(idx, D) = { my(V = vector(#idx - (idx[1] == 1), k, VX[idx[k + (idx[1] == 1)]]), L = brute(D, V), bad = [0, 0]);
  for (k = 1, #L~, for (m = 1, 2, if (qi_desc(L[k,1], idx, if (m == 1, 400, 10^6))[1], bad[m]++; print("  WRONG pmax=", if (m == 1, 400, "1e6"), " r=", L[k,1], " P=", L[k,2], " idx=", idx))));
  print(idx, ": ", #L~, " radii, wrongly killed (pmax 400, 1e6) ", bad); }
test([1,2], 30); test([1,3], 30); test([1,4], 30); test([2,3], 40); test([2,4], 40); test([3,4], 40);
test([1,2,3], 45); test([1,2,4], 40); test([1,3,4], 40);
print("DONE");
