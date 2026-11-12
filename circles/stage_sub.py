# stage_sub.py (mathematics-52, 2026-11-08): circ_stage2.gp / circ_stage3.gp driven from Python through cypari2.
# Each refined descent runs in its own child process (one.py) killed after SEC seconds: interrupting PARI inside
# cypari2 (cysignals alarm) left PARI in a broken state in stage_py.py, so the parent PARI is never interrupted.
# usage: python stage_sub.py LISTFILE NSH SH OUT PMAX FULL SEC
#   LISTFILE: one r per line; PMAX: prime bound of the refined descent; FULL: 1 = also try [1,2,3,4] last; SEC: seconds per descent
# Output lines are the same as circ_stage2.gp / circ_stage3.gp:  r rank0 j | r qi [a, b] ms | r open ranks ms [timeouts=n] | r SOL ...
import sys, os, time, cypari2
import subprocess
lst, nsh, sh, out, pmax, full, sec = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6] == "1", float(sys.argv[7])
p = cypari2.Pari(); p.allocatemem(400 * 10**6, silent=True); p.default("parisizemax", 2 * 10**9)
p.default("nbthreads", 1)   # the threaded libpari deadlocked in ellrank (mt_queue / pthread_join) on 11/08
here = os.path.dirname(os.path.abspath(__file__))
p('read("%s/circles.gp"); read("%s/qidesc.gp");' % (here, here))
PAIRS = ["[1,2]", "[1,3]", "[1,4]", "[2,3]", "[2,4]", "[3,4]"]
R = [l.strip() for l in open(lst) if l.strip()]
seen = set()
if os.path.exists(out):
    seen = {l.split()[0] for l in open(out) if l.strip() and l.strip() != "DONE"}
f = open(out, "a")
for k, r in enumerate(R):
    if k % nsh != sh or r in seen: continue
    T = time.time()
    v = p('circ_verdict(%s, 300)' % r)
    if int(v[0]) == -1: line = "%s SOL %s" % (r, v)
    elif int(v[0]) == 1: line = "%s rank0 %s" % (r, v[1])
    else:
        rk = [int(z) for z in v[2]]
        order = sorted(range(6), key=lambda j: (rk[j], j))
        trial = [PAIRS[j] for j in order] + (["[1,2,3,4]"] if full else [])
        line, nto = None, 0
        for pr in trial:
            try:
                cp = subprocess.run([sys.executable, os.path.join(here, "one.py"), r, pr, pmax], capture_output=True, text=True, timeout=sec)
            except subprocess.TimeoutExpired:
                nto += 1; continue
            res = cp.stdout.strip().splitlines()
            if cp.returncode != 0 or not res or res[-1] not in ("0", "1"):
                raise RuntimeError("one.py failed for %s %s: rc=%s out=%r err=%r" % (r, pr, cp.returncode, cp.stdout[-200:], cp.stderr[-300:]))
            if res[-1] == "1":
                line = "%s qi %s %d" % (r, pr.replace(",", ", "), int((time.time() - T) * 1000)); break
        if line is None:
            line = "%s open %s %d%s" % (r, str(rk).replace(" ", "").replace(",", ", "), int((time.time() - T) * 1000), " timeouts=%d" % nto if nto else "")
    f.write(line + "\n"); f.flush()
f.write("DONE\n"); f.close()
