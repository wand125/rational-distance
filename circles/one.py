# one.py (mathematics-52): one refined descent in a fresh PARI; prints "1" (no point) or "0". usage: python one.py r pair pmax
import sys, os, cypari2
p = cypari2.Pari(); p.allocatemem(400 * 10**6, silent=True); p.default("parisizemax", 2 * 10**9)
p.default("nbthreads", 1)   # the threaded libpari deadlocked in ellrank (mt_queue / pthread_join) on 11/08
here = os.path.dirname(os.path.abspath(__file__))
p('read("%s/circles.gp"); read("%s/qidesc.gp");' % (here, here))
print(int(p('qi_desc(%s, %s, %s)' % (sys.argv[1], sys.argv[2], sys.argv[3]))[0]))
