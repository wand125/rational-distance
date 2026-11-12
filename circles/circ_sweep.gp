\\ circ_sweep.gp (mathematics-52, 2026-11-07): common-square-class descent on the circles |P| = r, r = p/q, max(p, q) <= H,
\\ p % NSH == SH. Set H, NSH, SH, OUT before reading. One line per circle:
\\ r status #S #alive  kills-per-prime   (status 1 = every d killed: no 4-distance point at distance r from (0,0))
read("circles.gp");
{
  for (q = 1, H, for (p = 1, H, if (p % NSH == SH && gcd(p, q) == 1 && p != q,
    my(r = p/q, a = circ_cdesc2(r));
    write(OUT, Str(r, " ", a[1], " ", a[2], " ", a[3], " ", select(z -> z[2], a[5]), if (a[3], Str(" alive=", a[4]), "")));
    if (!a[1], print("ALIVE ", r, " ", a[3])))));
  write(OUT, "DONE");
}
