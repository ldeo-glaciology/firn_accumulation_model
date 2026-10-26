"""Concatenate rendered scenes into one video with merged subtitles."""
import os, re, subprocess, sys

SCENES = [
    ("s00_intro", "Intro"), ("s01_firn", "WhatIsFirn"), ("s02_puzzle", "Puzzle"),
    ("s03a_frame", "ModelFrame"), ("s03b_mass", "StressAndMass"), ("s03c_grains", "GrainsAndThickness"),
    ("s04_scaling", "Scaling"), ("s05_numerics", "Numerics"), ("s06_ode", "SteadyODE"),
    ("s07_constgrain", "ConstantGrain"), ("s08a_growth_sims", "GrowthSims"), ("s08b_growth_theory", "GrowthTheory"),
    ("s09_proof", "FullModelProof"), ("s10_robust", "Robustness"), ("s11_conclusion", "Conclusion"),
]
qdir = sys.argv[1] if len(sys.argv) > 1 else "480p15"
out = sys.argv[2] if len(sys.argv) > 2 else f"output/firn_explainer_{qdir}.mp4"
os.makedirs(os.path.dirname(out), exist_ok=True)


def dur(p):
    return float(subprocess.check_output(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0", p]).decode())


def ts(t):
    h, r = divmod(t, 3600)
    m, s = divmod(r, 60)
    return f"{int(h):02d}:{int(m):02d}:{int(s):02d},{int(round((s - int(s)) * 1000)):03d}".replace(",1000", ",999")


def parse_srt(p):
    blocks = open(p).read().strip().split("\n\n")
    res = []
    for b in blocks:
        lines = b.strip().split("\n")
        if len(lines) < 3:
            continue
        a, c = lines[1].split(" --> ")
        f = lambda x: (lambda h, m, s: int(h) * 3600 + int(m) * 60 + float(s.replace(",", ".")))(*x.strip().split(":"))
        res.append((f(a), f(c), " ".join(lines[2:])))
    return res


lst, subs, offset = [], [], 0.0
chapters = []
for f, c in SCENES:
    p = f"media/videos/{f}/{qdir}/{c}.mp4"
    if not os.path.exists(p):
        print("missing", p)
        continue
    lst.append(f"file '{os.path.abspath(p)}'")
    srt = p[:-4] + ".srt"
    if os.path.exists(srt):
        for a, b, t in parse_srt(srt):
            subs.append((a + offset, b + offset, t))
    chapters.append((offset, c))
    offset += dur(p)
open("output/list.txt", "w").write("\n".join(lst) + "\n")
srt_out = out[:-4] + ".srt"
with open(srt_out, "w") as fh:
    for i, (a, b, t) in enumerate(subs, 1):
        # wrap long captions into two lines
        words, lines, cur = t.split(), [], ""
        for w in words:
            if len(cur) + len(w) > 60 and cur:
                lines.append(cur)
                cur = w
            else:
                cur = (cur + " " + w).strip()
        lines.append(cur)
        fh.write(f"{i}\n{ts(a)} --> {ts(b)}\n" + "\n".join(lines) + "\n\n")
tmp = out[:-4] + "_nosub.mp4"
subprocess.check_call(["ffmpeg", "-v", "error", "-y", "-f", "concat", "-safe", "0", "-i", "output/list.txt",
                       "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-ar", "48000", tmp])
subprocess.check_call(["ffmpeg", "-v", "error", "-y", "-i", tmp, "-i", srt_out, "-map", "0", "-map", "1",
                       "-c", "copy", "-c:s", "mov_text", "-metadata:s:s:0", "language=eng", "-movflags", "+faststart", out])
os.remove(tmp)
print(f"{out}: {offset/60:.1f} min, {len(subs)} captions")
for t, c in chapters:
    print(f"  {int(t//60):02d}:{int(t%60):02d}  {c}")
