# Explainer video: grain size and the accumulation dependence of firn thickness

A 3Blue1Brown-style explainer (Manim Community + a British text-to-speech voice-over) of
Kingslake, Skarbek, Case & McCarthy (2022), *Grain-size evolution controls the accumulation
dependence of modelled firn thickness*, The Cryosphere 16, 3413–3430.

## Pieces

| file | what it does |
|---|---|
| `firn_sim.py` | Python port of the nondimensional model (method of lines on the stretched grid, BDF time stepping) plus the steady ODE models and the closed-form solutions of Eqs. 25 and 27. Reproduces the saved MATLAB results to ~1e-9. |
| `make_data.py` | Runs every simulation shown in the video (and reads the paper's saved MATLAB results in `../savedResults`) into `data/*.npz`. |
| `make_images.py` | Renders the contour-map colour fields (Figs 1, 6, 7, 8) as PNGs. |
| `tau_creature.py` | The tau creature (blinks, looks, moods, speech/thought bubbles). |
| `common.py` | Colours, LaTeX macros, plotting helpers and the voice-over machinery (`with self.vo("..."):`). |
| `scenes/s*.py` | One file per chapter. |
| `render.sh`, `queue.sh` | Render a scene / a list of scenes (`ql` preview or `final` = 1080p30). |
| `assemble.py` | Concatenates the chapters and merges subtitles into `output/`. |
| `check_audio.py` | Transcribes every narration clip with Whisper and flags mismatches. |

## Requirements

`pip install manim numpy scipy matplotlib soundfile kokoro-onnx` (+ `faster-whisper` for the audio check),
a LaTeX install with `dvisvgm`, and the Kokoro-82M ONNX model plus voice files in `$KOKORO_DIR`
(default `/home/user/tts_models`: `model.onnx`, `voices.npz`). Narration uses the British voice `bm_george`.

```
python make_data.py && python make_images.py
./queue.sh final s00_intro.py:Intro s01_firn.py:WhatIsFirn ...   # see assemble.py for the full list
python assemble.py 1080p30 output/firn_explainer_1080p.mp4
```
