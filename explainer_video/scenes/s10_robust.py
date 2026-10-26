from common import *


class Robustness(VScene):
    def construct(self):
        title = Text("How robust is this?", font="CMU Serif", font_size=40).to_edge(UP, buff=0.25)
        ap = M(r"\al' = \al\left(\frac{\sigma_0}{4}\right)^{1-n}", size=40).move_to([0, 1.8, 0])
        with self.vo("Two more checks. First: what if firn isn't linearly viscous? The authors re-ran the full model with stress exponents, n, "
                     "of two, three and four, adjusting alpha so that the strain rate matches the linear case at a typical stress, "
                     "a quarter of sigma zero.") as tr:
            self.play(Write(title), run_time=1)
            law = M(r"\DDt{\ph} \;\propto\; |\sg|^{n}", size=48).move_to([0, 0.2, 0])
            ns = [M(r"n = %d" % k, size=40, color=YELLOW).next_to(law, DOWN, buff=0.4) for k in [1, 2, 3, 4]]
            self.play(Write(law), FadeIn(ns[0]), run_time=1.2)
            for a, b in zip(ns[:-1], ns[1:]):
                self.wait(0.5)
                self.play(Transform(ns[0], b), run_time=0.5)
            self.wait(tr.until(0.6))
            self.play(Write(ap), FadeOut(law), FadeOut(ns[0]), run_time=1.5)
        kw = dict(x_range=(0.1, 10), y_range=(0.001, 0.1), x_len=3.3, y_len=3.0, x_ticks=[2, 4, 6, 8, 10],
                  y_ticks=[0.02, 0.04, 0.06, 0.08, 0.1], x_dec=0, y_dec=2, num_size=18, label_size=26)
        panels = Group(*[MapPanel(f"fig7_n{n}", title=r"n = %d" % n, x_label=r"\bt", y_label=(r"\rs" if n == 2 else None), **kw)
                         for n in [2, 3, 4]]).arrange(RIGHT, buff=0.55).move_to([0.2, -1.0, 0])
        with self.vo("The numbers shift, but the pattern is identical. Larger surface grains, stronger dependence on accumulation. "
                     "Just as the derivation predicts, the competition between advection of porosity and advection of grain size "
                     "doesn't care about n.") as tr:
            self.play(ap.animate.scale(0.75).move_to([0, 2.7, 0]).set_opacity(0.8), title.animate.set_opacity(0), run_time=0.6)
            self.play(LaggedStart(*[FadeIn(p) for p in panels], lag_ratio=0.3), run_time=2.5)
        self.play(FadeOut(panels), FadeOut(ap), title.animate.set_opacity(1), run_time=0.7)

        # surface height change
        base = -2.6
        col = Rectangle(width=1.6, height=3.6).set_fill(DEEP_ICE, 0.5).set_stroke(WHITE, 2)
        col.move_to([-4.6, base + 1.8, 0])
        h_tr = ValueTracker(3.6)
        col.add_updater(lambda m: m.stretch_to_fit_height(h_tr.get_value()).move_to([-4.6, base + h_tr.get_value() / 2, 0]))
        hdot = M(r"\dot h \ne 0", size=36).move_to([-4.6, 2.2, 0])
        mult = VGroup(M(r"\bb \times 1.1", r"\;\Rightarrow\; \dot h > 0", size=34), M(r"\bb \times 0.9", r"\;\Rightarrow\; \dot h < 0", size=34))
        mult.arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([-0.2, 0.6, 0])
        with self.vo("Second: so far, every result has been a steady state, with the ice-sheet surface at a fixed height. What if the surface is rising "
                     "or falling, because snowfall and ice flow are out of balance? The authors mimic this by multiplying the accumulation by one point one, "
                     "or by point nine, so the column steadily thickens or thins.") as tr:
            self.play(FadeIn(col), Write(hdot), run_time=1)
            self.wait(tr.until(0.5))
            self.play(FadeIn(mult[0]), h_tr.animate.set_value(4.4), run_time=1.5)
            self.play(FadeIn(mult[1]), h_tr.animate.set_value(2.8), run_time=1.5)
            self.play(h_tr.animate.set_value(3.6), run_time=1)
        col.clear_updaters()
        kw8 = dict(x_range=(0.1, 7.03), y_range=(0.001, 0.0703), x_len=3.9, y_len=3.3, x_ticks=[1, 3, 5, 7],
                   y_ticks=[0.01, 0.03, 0.05, 0.07], x_dec=0, y_dec=2, num_size=18, label_size=26)
        p_up = MapPanel("fig8_thicken", title=r"\dot h > 0 \text{ (thickening)}", x_label=r"\bt", y_label=r"\rs", **kw8)
        p_dn = MapPanel("fig8_thin", title=r"\dot h < 0 \text{ (thinning)}", x_label=r"\bt", **kw8)
        maps = Group(p_up, p_dn).arrange(RIGHT, buff=0.8).move_to([0.3, -0.7, 0])
        tau = TauCreature(height=1.3).move_to([6.2, -2.9, 0])
        with self.vo("A rising surface gives somewhat thicker firn, because a parcel ends up deeper by the time it has compacted. A sinking surface "
                     "gives somewhat thinner firn. But, once again, the same pattern.") as tr:
            self.play(FadeOut(VGroup(col, hdot, mult)), run_time=0.6)
            self.play(FadeIn(p_up), run_time=1.2)
            self.wait(tr.until(0.45))
            self.play(FadeIn(p_dn), run_time=1.2)
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            self.play(tau.change("happy", maps), run_time=0.6)
        self.wait(0.5)
        self.clear_all()
