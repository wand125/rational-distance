read("circles.gp"); read("qidesc.gp");
brute(D, V) = { my(R = Map()); for (d = 1, D, for (a = -3*d, 4*d, for (b = -3*d, 4*d, if (gcd([a, b, d]) == 1,
  my(x = a/d, y0 = b/d, r2 = x^2 + y0^2, ok = r2 != 0 && issquare(r2));
  if (ok, for (k = 1, #V, if (!issquare((x - V[k][1])^2 + (y0 - V[k][2])^2), ok = 0; break)));
  if (ok, my(r = sqrtint(numerator(r2)) / sqrtint(denominator(r2))); if (r != 1, mapput(R, r, [x, y0]))))))); Mat(R); }
test(V, idx, D, lim) = { my(L = brute(D, V), bad = 0, n = min(#L~, lim)); for (k = 1, n, if (qi_desc(L[k,1], idx)[1], bad++; print("  WRONG r=", L[k,1], " P=", L[k,2])));
  print(V, " idx ", idx, ": ", n, " radii, wrongly killed ", bad); }
T = getwalltime();
test([[1,0],[0,1]], [1, 2, 3], 45, 99);
test([[1,0],[1,1]], [1, 2, 4], 35, 99);
test([[1,1]], [1, 4], 25, 99);
print(getwalltime() - T, " ms"); print("DONE");
