# mw_trial.py (mathematics-52): run mw_one.py for every (r, j) of a list, each with a timeout; resumable.
# usage: python mw_trial.py LIST OUT Lmax ncontrol sec [NSH SH]
import sys, subprocess, os
lst, out, Lmax, nc, sec = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
nsh, sh = (int(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (1, 0)
rs = [l.split()[0] for l in open(lst) if l.strip()][sh::nsh]
done = set()
if os.path.exists(out):
    done = {tuple(l.split()[:2]) for l in open(out)}
here = os.path.dirname(os.path.abspath(__file__))
killed = set()
if os.path.exists(out):
    for l in open(out):
        w = l.split()
        if len(w) > 3 and w[2] == "ok" and w[5].startswith("res=[1,") and "wrong=0" in l and l.split("cross=")[1].split()[0].endswith(",0]"):
            killed.add(w[0])
for r in rs:
    for j in range(1, 7):
        if r in killed: break
        if (r, str(j)) in done: continue
        try:
            q = subprocess.run([sys.executable, os.path.join(here, "mw_one.py"), r, str(j), Lmax, nc], cwd=here,
                               capture_output=True, text=True, timeout=sec)
            line = q.stdout.strip().splitlines()[-1] if q.stdout.strip() else "%s %d fail %s" % (r, j, q.stderr.strip().splitlines()[-1][:80].replace(" ", "_") if q.stderr.strip() else "noout")
        except subprocess.TimeoutExpired:
            line = "%s %d timeout %d" % (r, j, sec)
        with open(out, "a") as f: f.write(line + "\n")
        w = line.split()
        if len(w) > 3 and w[2] == "ok" and w[5].startswith("res=[1,") and "wrong=0" in line and line.split("cross=")[1].split()[0].endswith(",0]"):
            killed.add(r)
open(out + ".DONE", "w").close()
