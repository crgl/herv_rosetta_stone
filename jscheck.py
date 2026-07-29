import re, subprocess, sys, os
h = open("dash/index.html").read()
blocks = re.findall(r'<script(?![^>]*src=)[^>]*>(.*?)</script>', h, re.S)
print("inline script blocks in index.html:", len(blocks))
targets = []
for i, b in enumerate(blocks):
    p = f"/tmp/inline{i}.js"
    open(p, "w").write(b)
    targets.append(p)
targets.append("dash/detail.js")
for p in targets:
    r = subprocess.run(["node", "--check", p], capture_output=True, text=True)
    n = sum(1 for _ in open(p))
    if r.returncode == 0:
        print(f"OK    {p} ({n} lines)")
    else:
        print(f"FAIL  {p} ({n} lines)")
        print(r.stderr.strip()[:1200])
