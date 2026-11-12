\\ circ_sweep3.gp (mathematics-52, 2026-11-07): like circ_sweep.gp but with circ_cdesc3 (F_2 linear algebra) and resumable.
\\ Circles r = p/q, max(p, q) <= H, r != 1, p % NSH == SH. Set H, NSH, SH, OUT before reading.
\\ One line per circle: r status #S #alive  (status 1 = every d killed); "DONE" at the end.
read("circles.gp");
{
  my(seen = Set());
  iferr(my(L = readstr(OUT)); seen = Set(apply(s -> strsplit(s, " ")[1], select(s -> s != "DONE" && #s, L))), E, );
  if (#seen, print("resume: ", #seen, " circles already done"));
  for (q = 1, H, for (p = 1, H, if (p % NSH == SH && gcd(p, q) == 1 && p != q,
    my(r = p/q);
    if (!setsearch(seen, Str(r)),
      my(a = circ_cdesc3(r));
      write(OUT, Str(r, " ", a[1], " ", a[2], " ", a[3]));
      if (!a[1], print("ALIVE ", r, " ", a[3]))))));
  write(OUT, "DONE");
}
