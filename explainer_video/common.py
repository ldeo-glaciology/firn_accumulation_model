"""Shared style, LaTeX macros and the voice-over machinery for the explainer."""
import hashlib
import os
import re
import sys
from contextlib import contextmanager

import numpy as np
from manim import *

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tau_creature import TauCreature, Snowflake, TAU_BLUE  # noqa: E402,F401

config.background_color = "#0E1014"

# manim's camera calls Image.fromarray(zeros((h, w)), mode="RGBA"), which newer Pillow rejects
import PIL.Image as _PILImage  # noqa: E402
_orig_fromarray = _PILImage.fromarray


def _fromarray(obj, mode=None):
    try:
        return _orig_fromarray(obj, mode) if mode is not None else _orig_fromarray(obj)
    except (ValueError, TypeError):
        a = np.asarray(obj)
        if mode == "RGBA" and a.ndim == 2:
            return _PILImage.new("RGBA", (a.shape[1], a.shape[0]))
        if mode is not None:
            return _orig_fromarray(a).convert(mode)
        raise


_PILImage.fromarray = _fromarray

# ------------------------------------------------------------------ colours
C_PHI = "#58C4DD"     # porosity
C_SIG = "#FC6255"     # stress
C_W = "#FFD84D"       # velocity
C_R = "#C08BFF"       # grain size
C_A = "#83C167"       # age
C_B = "#FF9A3C"       # accumulation
C_AL = "#FF7AB6"      # compaction number
C_DL = "#9CDBFF"      # saturation ratio
C_Z = "#5CD0B3"       # firn thickness z830
C_H = "#E8E8E8"
ICE = "#BFE6FF"
SNOW = "#F4F8FF"

TEX = TexTemplate()
TEX.add_to_preamble(r"""
\usepackage{xcolor}
\usepackage{amssymb}
\definecolor{cphi}{HTML}{58C4DD}
\definecolor{csig}{HTML}{FC6255}
\definecolor{cw}{HTML}{FFD84D}
\definecolor{cr}{HTML}{C08BFF}
\definecolor{cA}{HTML}{83C167}
\definecolor{cb}{HTML}{FF9A3C}
\definecolor{cal}{HTML}{FF7AB6}
\definecolor{cdl}{HTML}{9CDBFF}
\definecolor{cz}{HTML}{5CD0B3}
\newcommand{\ph}{{\color{cphi}\phi}}
\newcommand{\sg}{{\color{csig}\sigma}}
\newcommand{\ww}{{\color{cw}w}}
\newcommand{\rsq}{{\color{cr}r^2}}
\newcommand{\rr}{{\color{cr}r}}
\newcommand{\rs}{{\color{cr}r_s^2}}
\newcommand{\rf}{{\color{cr}r_f^2}}
\newcommand{\rz}{{\color{cr}r_0^2}}
\newcommand{\Ag}{{\color{cA}A}}
\newcommand{\bb}{{\color{cb}b}}
\newcommand{\bt}{{\color{cb}\beta}}
\newcommand{\al}{{\color{cal}\alpha}}
\newcommand{\dl}{{\color{cdl}\delta}}
\newcommand{\zc}{{\color{cz}z_{830}}}
\newcommand{\DDt}[1]{\frac{\mathrm{D}#1}{\mathrm{D}t}}
\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\dd}[2]{\frac{\mathrm{d} #1}{\mathrm{d} #2}}
""")


def M(*parts, size=40, **kw):
    """MathTex using the coloured macro template."""
    return MathTex(*parts, tex_template=TEX, font_size=size, **kw)


def T(*parts, size=36, **kw):
    return Tex(*parts, tex_template=TEX, font_size=size, **kw)


def label_box(mob, color=YELLOW, buff=0.15):
    return SurroundingRectangle(mob, color=color, buff=buff, corner_radius=0.08)


def eq_tag(n, size=26):
    return Tex(f"({n})", font_size=size, color=GREY_B)


def with_tag(eq, n, buff=0.4):
    tag = eq_tag(n).next_to(eq, RIGHT, buff=buff)
    return VGroup(eq, tag)


def chapter_card(title, number=None):
    t = Text(title, font="CMU Serif", font_size=54, weight=BOLD)
    g = VGroup(t)
    if number is not None:
        n = Text(f"Part {number}", font="CMU Serif", font_size=30, color=GREY_B)
        n.next_to(t, UP, buff=0.3)
        g.add(n)
    line = Line(LEFT * 3, RIGHT * 3, color=C_PHI, stroke_width=2).next_to(t, DOWN, buff=0.3)
    g.add(line)
    return g.move_to(ORIGIN)


# ------------------------------------------------------------- voice over
TTS_DIR = os.path.join(HERE, "media", "tts")
VOICE = os.environ.get("FIRN_VOICE", "bm_george")
SPEED = float(os.environ.get("FIRN_SPEED", "0.95"))
MODEL_DIR = os.environ.get("KOKORO_DIR", "/home/user/tts_models")

# spelling tweaks so the British espeak front-end says technical words properly
PRONOUNCE = [
    (r"\bEulerian\b", "Oilerian"),
    (r"\bLagrangian\b", "Lagrahnjian"),
    (r"\bArrhenius\b", "Arreenius"),
    (r"\bArthern\b", "Arthurn"),
    (r"\bLeibniz\b", "Libenitz"),
    (r"\bODEs\b", "O D Es"),
    (r"\bODE\b", "O D E"),
    (r"\bPDEs\b", "P D Es"),
]

_kokoro = None


def _engine():
    global _kokoro
    if _kokoro is None:
        from kokoro_onnx import Kokoro
        _kokoro = Kokoro(os.path.join(MODEL_DIR, "model.onnx"), os.path.join(MODEL_DIR, "voices.npz"))
    return _kokoro


def split_markup(text):
    """'[shown](spoken)' markup -> (caption text, spoken text)."""
    caption = re.sub(r"\[([^\]]*)\]\(([^)]*)\)", r"\1", text)
    spoken = re.sub(r"\[([^\]]*)\]\(([^)]*)\)", r"\2", text)
    for pat, rep in PRONOUNCE:
        spoken = re.sub(pat, rep, spoken)
    caption = re.sub(r"\s+", " ", caption).strip()
    spoken = re.sub(r"\s+", " ", spoken).strip()
    return caption, spoken


def tts(text):
    import soundfile as sf
    caption, spoken = split_markup(text)
    key = hashlib.sha1(f"{VOICE}|{SPEED}|{spoken}".encode()).hexdigest()[:16]
    os.makedirs(TTS_DIR, exist_ok=True)
    path = os.path.join(TTS_DIR, key + ".wav")
    if not os.path.exists(path):
        samples, sr = _engine().create(spoken, voice=VOICE, speed=SPEED, lang="en-gb")
        # trim leading/trailing near-silence, keep a short natural tail
        thr = 0.01 * np.max(np.abs(samples))
        idx = np.where(np.abs(samples) > thr)[0]
        if len(idx):
            a = max(idx[0] - int(0.03 * sr), 0)
            b = min(idx[-1] + int(0.12 * sr), len(samples))
            samples = samples[a:b]
        sf.write(path, samples, sr)
        with open(path[:-4] + ".txt", "w") as f:
            f.write(spoken)
    info = sf.info(path)
    return path, info.frames / info.samplerate, caption


class _Tracker:
    def __init__(self, scene, start, duration):
        self.scene, self.start, self.duration = scene, start, duration

    @property
    def remaining(self):
        return max(self.start + self.duration - self.scene.renderer.time, 0.0)

    def until(self, frac):
        """Seconds from now until `frac` of the clip has been spoken."""
        return max(self.start + frac * self.duration - self.scene.renderer.time, 0.0)


class VScene(Scene):
    """Scene with `with self.vo("..."):` voice-over blocks.

    The audio clip starts when the block is entered; animations inside run
    concurrently; on exit the scene waits for the clip to finish (+ gap).
    """
    gap = 0.35

    def wait(self, duration=1.0, *args, **kw):
        if duration < 0.04:
            return
        super().wait(duration, *args, **kw)

    @contextmanager
    def vo(self, text, gap=None):
        path, dur, caption = tts(text)
        start = self.renderer.time
        self.add_sound(path)
        self.add_subcaption(caption, duration=dur, offset=0)
        tr = _Tracker(self, start, dur)
        yield tr
        rest = start + dur + (self.gap if gap is None else gap) - self.renderer.time
        if rest > 0.04:
            self.wait(rest)

    def say(self, text, *anims, run_time=None, gap=None, **kw):
        """Speak `text` while playing `anims` (spread over the clip if run_time is None)."""
        with self.vo(text, gap=gap) as tr:
            if anims:
                rt = run_time if run_time is not None else min(tr.duration, max(1.0, tr.duration * 0.6))
                self.play(*anims, run_time=rt, **kw)

    def fill(self, tr, frac=1.0):
        t = tr.until(frac)
        if t > 0.04:
            self.wait(t)

    def clear_all(self, run_time=0.8):
        if self.mobjects:
            self.play(*[FadeOut(m) for m in self.mobjects], run_time=run_time)


# --------------------------------------------------------------- data/plots
DATA = os.path.join(HERE, "data")


def load(name):
    return dict(np.load(os.path.join(DATA, name + ".npz"), allow_pickle=True))


def depth_axes(x_range, y_range, x_len=5, y_len=4.2, x_label=None, y_label=None,
               x_step=None, y_step=None, label_size=28, num_size=22, flip_y=False):
    ax = Axes(x_range=[x_range[0], x_range[1], x_step or (x_range[1] - x_range[0]) / 4],
              y_range=[y_range[0], y_range[1], y_step or (y_range[1] - y_range[0]) / 4],
              x_length=x_len, y_length=y_len, tips=False,
              axis_config=dict(color=GREY_B, stroke_width=2, include_numbers=True,
                               font_size=num_size, decimal_number_config=dict(num_decimal_places=1)))
    labels = VGroup()
    if x_label is not None:
        xl = x_label if isinstance(x_label, Mobject) else M(x_label, size=label_size)
        xl.next_to(ax.x_axis, DOWN, buff=0.45)
        labels.add(xl)
    if y_label is not None:
        yl = y_label if isinstance(y_label, Mobject) else M(y_label, size=label_size)
        yl.next_to(ax.y_axis, LEFT, buff=0.45)
        labels.add(yl)
    return ax, labels


def curve(ax, x, y, color=WHITE, width=3, **kw):
    m = np.isfinite(x) & np.isfinite(y)
    pts = [ax.c2p(a, b) for a, b in zip(x[m], y[m])]
    v = VMobject(stroke_color=color, stroke_width=width, **kw)
    if len(pts) > 1:
        v.set_points_as_corners(pts)
    return v


# ------------------------------------------------------------ shared visuals
DEEP_ICE = "#2F7FB8"


def layer_color(k, scale=5.0):
    return interpolate_color(ManimColor(SNOW), ManimColor(DEEP_ICE), 1 - np.exp(-k / scale))


def layer_thick(k):
    return 0.13 + 0.55 * np.exp(-k / 3.5)


def build_stack(n_layers, top_y, width=3.6, x=0.0, offset=0):
    layers = VGroup()
    y = top_y
    for k in range(n_layers):
        t = layer_thick(k + offset)
        r = Rectangle(width=width, height=t)
        r.set_fill(layer_color(k + offset), 1).set_stroke(config.background_color, 1.2)
        r.move_to([x, y - t / 2, 0])
        y -= t
        layers.add(r)
    return layers


def snowfall(n=36, x_range=(-7, 7), y_bottom=1.5, y_top=4.4, speed=0.9, seed=1):
    rng = np.random.default_rng(seed)
    flakes = VGroup()
    for _ in range(n):
        f = Snowflake(size=rng.uniform(0.05, 0.11)).move_to(
            [rng.uniform(*x_range), rng.uniform(y_bottom, y_top), 0])
        f.v = speed * rng.uniform(0.7, 1.3)
        f.phase = rng.uniform(0, TAU)
        flakes.add(f)

    def upd(group, dt):
        for f in group:
            f.shift(DOWN * f.v * dt + RIGHT * 0.15 * np.sin(f.phase + 2 * f.get_y()) * dt)
            f.rotate(0.6 * dt)
            if f.get_y() < y_bottom:
                f.set_y(y_top)
    flakes.add_updater(upd)
    return flakes


def grain_cluster(radius=1.3, n=70, rmin=0.06, rmax=0.13, seed=3, color=SNOW, porosity=0.5):
    """Random non-overlapping grains inside a disc (a magnified view of firn)."""
    rng = np.random.default_rng(seed)
    grains = VGroup()
    placed = []
    tries = 0
    target_area = (1 - porosity) * PI * radius ** 2
    area = 0
    while area < target_area and tries < 20000:
        tries += 1
        r = rng.uniform(rmin, rmax)
        p = rng.uniform(-radius, radius, 2)
        if np.linalg.norm(p) > radius - r:
            continue
        if any(np.linalg.norm(p - q) < r + s + 0.01 for q, s in placed):
            continue
        placed.append((p, r))
        area += PI * r ** 2
        grains.add(Circle(radius=r).move_to([p[0], p[1], 0]).set_fill(color, 1).set_stroke(GREY_B, 1))
    return grains


class DepthPlot(VGroup):
    """Profile plot with depth increasing downwards and the x-axis along the top."""

    def __init__(self, x_range, z_range, width=4.0, height=5.0, x_ticks=None, z_ticks=None,
                 x_label=None, z_label=None, x_dec=1, z_dec=1, num_size=22, label_size=30,
                 color=GREY_B, **kw):
        super().__init__(**kw)
        self.x_range, self.z_range, self.w, self.h = x_range, z_range, width, height
        x0, x1 = x_range
        z0, z1 = z_range
        self.top = Line(ORIGIN, RIGHT * width, color=color, stroke_width=2)
        self.left = Line(ORIGIN, DOWN * height, color=color, stroke_width=2)
        self.add(self.top, self.left)
        self.xnums = VGroup()
        self.znums = VGroup()
        for xt in (x_ticks if x_ticks is not None else np.linspace(x0, x1, 5)):
            p = self.p(xt, z0)
            self.add(Line(p, p + DOWN * 0.1, color=color, stroke_width=2))
            n = DecimalNumber(xt, num_decimal_places=x_dec, font_size=num_size, color=GREY_A).next_to(p, UP, buff=0.12)
            self.xnums.add(n)
        for zt in (z_ticks if z_ticks is not None else np.linspace(z0, z1, 5)):
            p = self.p(x0, zt)
            self.add(Line(p, p + RIGHT * 0.1, color=color, stroke_width=2))
            n = DecimalNumber(zt, num_decimal_places=z_dec, font_size=num_size, color=GREY_A).next_to(p, LEFT, buff=0.12)
            self.znums.add(n)
        self.add(self.xnums, self.znums)
        self.labels = VGroup()
        if x_label is not None:
            xl = x_label if isinstance(x_label, Mobject) else M(x_label, size=label_size)
            xl.next_to(self.xnums, UP, buff=0.15)
            self.labels.add(xl)
            self.x_lab = xl
        if z_label is not None:
            zl = z_label if isinstance(z_label, Mobject) else M(z_label, size=label_size)
            zl.next_to(self.znums, LEFT, buff=0.15)
            self.labels.add(zl)
            self.z_lab = zl
        self.add(self.labels)

    def p(self, x, z):
        x0, x1 = self.x_range
        z0, z1 = self.z_range
        o = self.top.get_start()
        return o + RIGHT * (x - x0) / (x1 - x0) * self.w + DOWN * (z - z0) / (z1 - z0) * self.h

    def curve(self, x, z, color=WHITE, width=3, clip=True, **kw):
        x = np.asarray(x, float)
        z = np.asarray(z, float)
        m = np.isfinite(x) & np.isfinite(z)
        if clip:
            m &= (z >= self.z_range[0] - 1e-9) & (z <= self.z_range[1] + 1e-9)
        pts = [self.p(a, b) for a, b in zip(x[m], z[m])]
        v = VMobject(stroke_color=color, stroke_width=width, **kw)
        if len(pts) > 1:
            v.set_points_as_corners(pts)
        return v

    def hline(self, z, color=GREY_B, dashed=True, width=2):
        l = Line(self.p(self.x_range[0], z), self.p(self.x_range[1], z), color=color, stroke_width=width)
        return DashedLine(l.get_start(), l.get_end(), color=color, stroke_width=width, dash_length=0.08) if dashed else l

    def vline(self, x, color=GREY_B, dashed=True, width=2):
        a, b = self.p(x, self.z_range[0]), self.p(x, self.z_range[1])
        return DashedLine(a, b, color=color, stroke_width=width, dash_length=0.08) if dashed else Line(a, b, color=color, stroke_width=width)


def firn_column(width=1.6, height=5.0, phi_profile=None, top=ORIGIN, seed=5, n_grains=220, grain_growth=True):
    """A column whose colour and grain packing follow a porosity profile phi(depth fraction)."""
    if phi_profile is None:
        phi_profile = lambda s: 0.5 * np.exp(-((s / 0.35) ** 2))
    g = VGroup()
    nb = 40
    for i in range(nb):
        s = (i + 0.5) / nb
        phi = phi_profile(s)
        c = interpolate_color(ManimColor(DEEP_ICE), ManimColor(SNOW), min(phi / 0.5, 1) ** 0.8)
        r = Rectangle(width=width, height=height / nb + 0.005, stroke_width=0).set_fill(c, 1)
        r.move_to(top + DOWN * (s * height))
        g.add(r)
    col = VGroup(*g)
    col.move_to(top + DOWN * height / 2)
    outline = Rectangle(width=width, height=height, stroke_color=WHITE, stroke_width=2).move_to(col)
    rng = np.random.default_rng(seed)
    grains = VGroup()
    for _ in range(n_grains):
        s = rng.uniform(0.0, 1.0)
        phi = phi_profile(s)
        if rng.uniform() > (phi / 0.5) ** 0.7 + 0.05:
            continue
        rad = 0.025 + (0.05 * s if grain_growth else 0.0)
        x = rng.uniform(-width / 2 + rad, width / 2 - rad)
        grains.add(Circle(radius=rad, stroke_width=0.8, stroke_color=GREY_C).set_fill(WHITE, 0.55).move_to(
            col.get_top() + DOWN * s * height + RIGHT * x))
    return VGroup(col, grains, outline)


IMAGES = os.path.join(HERE, "images")


class MapPanel(Group):
    """A colour-field image (from make_images.py) inside manim-drawn axes."""

    def __init__(self, name, x_range, y_range, x_len=3.6, y_len=3.0, x_ticks=(), y_ticks=(),
                 x_tick_labels=None, y_tick_labels=None, x_label=None, y_label=None, title=None,
                 num_size=20, label_size=26, x_dec=0, y_dec=2, **kw):
        super().__init__(**kw)
        self.x_range, self.y_range, self.x_len, self.y_len = x_range, y_range, x_len, y_len
        img = ImageMobject(os.path.join(IMAGES, name + ".png"))
        img.stretch_to_fit_width(x_len)
        img.stretch_to_fit_height(y_len)
        self.img = img
        frame = Rectangle(width=x_len, height=y_len, stroke_color=GREY_B, stroke_width=2)
        self.frame = frame
        self.add(img, frame)
        nums = VGroup()
        for i, xt in enumerate(x_ticks):
            p = self.c2p(xt, y_range[0])
            nums.add(Line(p, p + UP * 0.08, color=GREY_B, stroke_width=2))
            lab = (x_tick_labels[i] if x_tick_labels else DecimalNumber(xt, num_decimal_places=x_dec, font_size=num_size))
            if isinstance(lab, str):
                lab = M(lab, size=num_size)
            nums.add(lab.next_to(p, DOWN, buff=0.1))
        for i, yt in enumerate(y_ticks):
            p = self.c2p(x_range[0], yt)
            nums.add(Line(p, p + RIGHT * 0.08, color=GREY_B, stroke_width=2))
            lab = (y_tick_labels[i] if y_tick_labels else DecimalNumber(yt, num_decimal_places=y_dec, font_size=num_size))
            if isinstance(lab, str):
                lab = M(lab, size=num_size)
            nums.add(lab.next_to(p, LEFT, buff=0.1))
        self.nums = nums
        self.add(nums)
        if x_label:
            xl = M(x_label, size=label_size) if isinstance(x_label, str) else x_label
            xl.next_to(frame, DOWN, buff=0.45)
            self.add(xl)
        if y_label:
            yl = M(y_label, size=label_size) if isinstance(y_label, str) else y_label
            if yl.width > 1.2:
                yl.rotate(PI / 2)
            yl.next_to(self.nums, LEFT, buff=0.15)
            self.add(yl)
        if title:
            tl = M(title, size=label_size) if isinstance(title, str) else title
            tl.next_to(frame, UP, buff=0.15)
            self.title = tl
            self.add(tl)

    def c2p(self, x, y):
        f = self.frame
        x0, x1 = self.x_range
        y0, y1 = self.y_range
        return f.get_corner(DL) + RIGHT * (x - x0) / (x1 - x0) * self.x_len + UP * (y - y0) / (y1 - y0) * self.y_len
