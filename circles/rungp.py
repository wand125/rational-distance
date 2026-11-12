# rungp.py (mathematics-52): run a gp script through cypari2 (no gp binary on this host). usage: python rungp.py file.gp
import sys, cypari2
p = cypari2.Pari()
p.allocatemem(400 * 10**6, silent=True)
p.default("parisizemax", 2 * 10**9); p.default("nbthreads", 1)
p('read("%s")' % sys.argv[1])
