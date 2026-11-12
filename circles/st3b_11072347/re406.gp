\\ recheck the stage-2 "qi" kills with the corrected qidesc.gp (same pmax 400, pairs in the stage-2 order)
read("circles.gp"); read("qidesc.gp"); read("st3b_11072347/list406.gp");
ID = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];
{for (k = 1, #LIST, my(r = LIST[k], rk = vector(6, j, ellrank(ellminimalmodel(qjac(cq(r)[j])))[2]), ord = vecsort(vector(6, j, [rk[j], j])), hit = 0);
  for (o = 1, 6, my(j = ord[o][2], a = alarm(120, qi_desc(r, ID[j]))); if (type(a) != "t_ERROR" && a[1], hit = j; break));
  write("st3b_11072347/re406.txt", Str(r, " ", if (hit, Str("qi ", ID[hit]), "open"))));}
write("st3b_11072347/re406.txt", "DONE");
