"""Transcribe every narration clip with Whisper and flag ones that don't match the script."""
import glob, os, re, sys, difflib, json
import numpy as np, soundfile as sf
from scipy.signal import resample_poly
from faster_whisper import WhisperModel

m = WhisperModel("small.en", device="cpu", compute_type="int8", cpu_threads=4)
done = {}
if os.path.exists("test/audio_check.json"):
    done = json.load(open("test/audio_check.json"))


def norm(t):
    t = t.lower().replace("feern", "firn").replace("fern", "firn")
    t = re.sub(r"[^a-z0-9 ]", " ", t)
    return t.split()


for wav in sorted(glob.glob("media/tts/*.wav")):
    key = os.path.basename(wav)
    if key in done:
        continue
    ref = open(wav[:-4] + ".txt").read()
    a, sr = sf.read(wav)
    a = resample_poly(a, 16000, sr).astype(np.float32)
    segs, _ = m.transcribe(a, beam_size=5)
    hyp = " ".join(s.text for s in segs)
    r = difflib.SequenceMatcher(None, norm(ref), norm(hyp)).ratio()
    done[key] = dict(ratio=r, ref=ref, hyp=hyp, dur=len(a) / 16000)
    json.dump(done, open("test/audio_check.json", "w"), indent=1)

bad = sorted(done.items(), key=lambda kv: kv[1]["ratio"])
tot = sum(v["dur"] for v in done.values())
print(f"{len(done)} clips, {tot/60:.1f} min of narration")
for k, v in bad[:12]:
    print(f"{v['ratio']:.2f} {k}\n  REF: {v['ref']}\n  HYP: {v['hyp']}")
