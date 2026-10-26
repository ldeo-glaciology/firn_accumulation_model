"""A 3Blue1Brown-style 'tau creature' built from manim primitives."""
import numpy as np
from manim import *
from manim.mobject.geometry.boolean_ops import Union

TAU_BLUE = "#2A7FA8"
TAU_BLUE_DARK = "#1D5E7E"


def _body_shape():
    from manim.mobject.geometry.boolean_ops import Intersection
    # top bar: a chunky rounded bar with a drooping, curled left tip
    bar = RoundedRectangle(width=2.3, height=0.62, corner_radius=0.3).move_to([0.05, 0.78, 0])
    droop = Ellipse(width=0.62, height=0.66).move_to([-0.95, 0.6, 0])
    # single stem that turns to the right at the bottom, like a hand-written tau
    th, inner = 0.56, 0.2
    cx, cy = 0.12 - th / 2 + inner + th, -0.45
    stem = Rectangle(width=th, height=0.95 - cy).move_to([0.12, (0.95 + cy) / 2, 0])
    ring = Annulus(inner_radius=inner, outer_radius=inner + th).move_to([cx, cy, 0])
    quarter = Intersection(ring, Rectangle(width=inner + th + 0.02, height=inner + th + 0.02).move_to(
        [cx - (inner + th) / 2 + 0.01, cy - (inner + th) / 2 + 0.01, 0]))
    cap = Circle(radius=th / 2).move_to([cx, cy - inner - th / 2, 0])
    return Union(bar, droop, stem, quarter, cap)


class TauCreature(VGroup):
    def __init__(self, color=TAU_BLUE, mood="plain", height=2.0, **kwargs):
        super().__init__(**kwargs)
        self.body = _body_shape().set_fill(color, 1).set_stroke(width=0)
        self.eye_centers = [np.array([-0.12, 1.08, 0]), np.array([0.42, 1.08, 0])]
        self.eyes = VGroup()
        self.pupils = VGroup()
        for c in self.eye_centers:
            eye = Ellipse(width=0.44, height=0.5).move_to(c)
            eye.set_fill(WHITE, 1).set_stroke(BLACK, 1.5)
            pupil = VGroup(
                Circle(radius=0.13).set_fill(BLACK, 1).set_stroke(width=0),
                Circle(radius=0.04).set_fill(WHITE, 1).set_stroke(width=0).shift(0.05 * UR),
            ).move_to(c + 0.04 * DOWN)
            self.eyes.add(eye)
            self.pupils.add(pupil)
        self.mouth = self._mouth(mood)
        self.add(self.body, self.eyes, self.pupils, self.mouth)
        self.mood = mood
        self.scale(height / 2.4)

    # ---------------------------------------------------------------- mouths
    def _mouth_anchor(self):
        return np.array([0.16, 0.66, 0])

    def _mouth(self, mood):
        c = self._mouth_anchor()
        if mood in ("plain", "happy", "speaking"):
            m = ArcBetweenPoints(c + 0.22 * LEFT, c + 0.22 * RIGHT, angle=PI / 2.2 if mood == "happy" else PI / 3.5)
            m.set_stroke(BLACK, 4)
            if mood == "happy":
                m = VMobject().set_points_as_corners([c + 0.24 * LEFT + 0.06 * UP, c + 0.24 * RIGHT + 0.06 * UP])
                m.add_points_as_corners([c + 0.24 * RIGHT + 0.06 * UP])
                arc = Arc(radius=0.24, start_angle=PI, angle=PI).move_to(c + 0.05 * DOWN)
                m = VGroup(arc).set_fill(BLACK, 1).set_stroke(BLACK, 2)
                m = Sector(radius=0.24, start_angle=PI, angle=PI).move_to(c + 0.02 * DOWN)
                m.set_fill(BLACK, 1).set_stroke(BLACK, 1)
            if mood == "speaking":
                m = Ellipse(width=0.34, height=0.2).move_to(c).set_fill(BLACK, 1).set_stroke(BLACK, 1)
            else:
                m = m if isinstance(m, Sector) else ArcBetweenPoints(c + 0.22 * LEFT, c + 0.22 * RIGHT, angle=PI / 3).set_stroke(BLACK, 4)
        elif mood == "surprised":
            m = Ellipse(width=0.2, height=0.26).move_to(c).set_fill(BLACK, 1).set_stroke(BLACK, 1)
        elif mood == "confused":
            m = VMobject().set_points_smoothly([c + 0.22 * LEFT + 0.03 * DOWN, c + 0.07 * LEFT + 0.03 * UP,
                                                c + 0.07 * RIGHT + 0.03 * DOWN, c + 0.22 * RIGHT + 0.04 * UP])
            m.set_stroke(BLACK, 4)
        elif mood == "pondering":
            m = Line(c + 0.18 * LEFT + 0.02 * UP, c + 0.14 * RIGHT + 0.02 * DOWN).set_stroke(BLACK, 4)
        elif mood == "sad":
            m = ArcBetweenPoints(c + 0.2 * LEFT, c + 0.2 * RIGHT, angle=-PI / 3).set_stroke(BLACK, 4).shift(0.04 * DOWN)
        else:
            m = ArcBetweenPoints(c + 0.22 * LEFT, c + 0.22 * RIGHT, angle=PI / 3).set_stroke(BLACK, 4)
        return m

    def _mouth_transformed_to_current(self, mood):
        """Build a new mouth for `mood` in the creature's current frame."""
        ref = TauCreature(mood=mood)
        new = ref.mouth
        # map from reference frame to this creature: use the body bounding boxes
        new.match_height(self.mouth, stretch=False) if False else None
        s = self.body.width / ref.body.width
        new.scale(s, about_point=ORIGIN)
        new.shift(self.body.get_center() - ref.body.get_center() * s)
        if self.body.get_width() and self._is_flipped():
            new.flip(UP, about_point=self.body.get_center())
        return new

    def _is_flipped(self):
        # left droop sits left of the eyes in the default orientation
        return self.eyes[0].get_center()[0] > self.eyes[1].get_center()[0]

    # ------------------------------------------------------------ animations
    def change(self, mood, look_at=None):
        """Return an animation changing the mouth (and optionally gaze)."""
        target = self.generate_target()
        new_mouth = self._mouth_transformed_to_current(mood)
        target.mouth.become(new_mouth)
        if mood == "surprised":
            for e in target.eyes:
                e.stretch(1.12, 1)
        if look_at is not None:
            target.look_at_point(look_at)
        self.mood = mood
        return MoveToTarget(self)

    def look_at_point(self, point):
        if isinstance(point, Mobject):
            point = point.get_center()
        for eye, pupil in zip(self.eyes, self.pupils):
            d = np.array(point) - eye.get_center()
            n = np.linalg.norm(d)
            if n > 0:
                d = d / n
            r = 0.32 * eye.width / 2 * 1.0
            pupil.move_to(eye.get_center() + r * d * np.array([1.0, 0.9, 0]))
        return self

    def look_at(self, point):
        return self.animate.look_at_point(point)

    def blink(self):
        target = self.copy()
        for eye, pupil in zip(target.eyes, target.pupils):
            y = eye.get_center()[1]
            eye.stretch(0.08, 1).set_y(y)
            pupil.stretch(0.08, 1).set_y(y)
        return Transform(self, target, rate_func=there_and_back, run_time=0.3)

    def bubble(self, content, kind="speech", width=None, height=None, direction=None):
        """Create (not animate) a speech/thought bubble above the creature."""
        if isinstance(content, str):
            content = Tex(content, font_size=34)
        content.set_color(BLACK)
        pad = 0.35
        w = width or max(content.width + 2 * pad, 2.0)
        h = height or max(content.height + 2 * pad, 1.0)
        flip = (direction is RIGHT) if direction is not None else (self.get_center()[0] < 0)
        side = 1 if flip else -1
        top = self.get_top()
        center = top + np.array([side * (w / 2 + 0.2), h / 2 + 0.55, 0])
        if kind == "speech":
            box = RoundedRectangle(width=w, height=h, corner_radius=0.35).move_to(center)
            tail = Polygon(center + np.array([-side * w * 0.28, -h / 2 + 0.3, 0]),
                           center + np.array([-side * w * 0.08, -h / 2 + 0.3, 0]),
                           top + np.array([side * 0.25, 0.12, 0]))
            shape = Union(box, tail)
            shape.set_fill(WHITE, 1).set_stroke(GREY_B, 2)
            bubble = VGroup(shape)
        else:
            box = Ellipse(width=w * 1.15, height=h * 1.25).move_to(center)
            box.set_fill(WHITE, 1).set_stroke(GREY_B, 2)
            dots = VGroup(*[Circle(radius=r).set_fill(WHITE, 1).set_stroke(GREY_B, 2).move_to(
                top + np.array([side * (0.15 + 0.22 * i), 0.1 + 0.2 * i, 0])) for i, r in enumerate([0.05, 0.08, 0.11])])
            bubble = VGroup(box, dots)
        content.move_to(center)
        bubble.add(content)
        bubble.content = content
        return bubble

    def says(self, content, mood="speaking", **kw):
        b = self.bubble(content, "speech", **kw)
        return b, AnimationGroup(self.change(mood, b), FadeIn(b[0], scale=0.8, shift=0.2 * DOWN), Write(b.content), lag_ratio=0.1)

    def thinks(self, content, mood="pondering", **kw):
        b = self.bubble(content, "thought", **kw)
        return b, AnimationGroup(self.change(mood, b), FadeIn(b[:-1], scale=0.8), Write(b.content), lag_ratio=0.1)


class Snowflake(VGroup):
    def __init__(self, size=0.15, color=WHITE, **kw):
        super().__init__(**kw)
        for k in range(3):
            self.add(Line(size * LEFT, size * RIGHT, stroke_width=1.6, color=color).rotate(k * PI / 3))
