from common import *
from s03a_frame import setup_column


class StressAndMass(VScene):
    def construct(self):
        col, zaxis = setup_column()
        self.add(col, zaxis)
        top, bot = col.get_top()[1], col.get_bottom()[1]
        H = top - bot

        # ---- stress = weight of overlying firn
        above = Rectangle(width=col.width, height=2.0, stroke_width=0).set_fill(C_SIG, 0.28).align_to(col, UP).set_x(col.get_x())
        slab = Rectangle(width=col.width, height=0.25, stroke_color=YELLOW, stroke_width=3).next_to(above, DOWN, buff=0)
        dz = M(r"\mathrm{d}z", size=28).next_to(slab, RIGHT, buff=0.15)
        weight = M(r"|\sg| = \text{weight above}", size=34).move_to([2.0, 2.4, 0])
        s1 = M(r"\mathrm{d}\sg", "=", r"-\,\rho_i\, g\,(1-\ph)\,\mathrm{d}z", size=44).move_to([2.0, 1.2, 0])
        s2 = M(r"\pd{\sg}{z}", "=", r"-\,\rho_i\, g\,(1-\ph)", r",\qquad \sg(0)=0", size=44).move_to([2.0, 1.2, 0])
        tag5 = eq_tag(5).next_to(s2, RIGHT, buff=0.4)
        with self.vo("First, we need the stress. That's simply the weight of everything above. Go down a small distance, d zed, and you add a slab "
                     "whose mass, per unit area, is the ice density times one minus phi times d zed. So d sigma by d zed equals minus rho i g "
                     "times one minus phi, with zero stress at the surface.") as tr:
            self.play(FadeIn(above), Write(weight), run_time=1.2)
            self.wait(tr.until(0.2))
            self.play(Create(slab), FadeIn(dz), run_time=0.8)
            self.wait(tr.until(0.35))
            self.play(Write(s1), run_time=1.5)
            self.wait(tr.until(0.62))
            self.play(TransformMatchingShapes(s1, s2[:3]), run_time=1.2)
            self.wait(tr.until(0.88))
            self.play(Write(s2[3]), FadeIn(tag5))
        stress_eq = VGroup(s2, tag5)
        self.play(FadeOut(above), FadeOut(slab), FadeOut(dz), FadeOut(weight),
                  stress_eq.animate.scale(0.65).to_corner(UR).shift(LEFT * 0.1), run_time=0.8)

        # ---- mass conservation in a slab
        ctr = np.array([-0.9, 0.4, 0])
        rng = np.random.default_rng(4)
        box0 = Rectangle(width=2.4, height=1.6).move_to(ctr).set_stroke(WHITE, 2)
        grains = VGroup(*[Circle(0.12).set_fill(SNOW, 0.9).set_stroke(GREY_B, 1).move_to(
            ctr + [rng.uniform(-1.05, 1.05), rng.uniform(-0.65, 0.65), 0]) for _ in range(22)])
        slab_g = VGroup(box0, grains)
        lab_ice = M(r"\text{ice content} = (1-\ph)\,\Delta z", size=32).next_to(box0, DOWN, buff=0.3)
        dzl = BraceLabel(box0, r"\Delta z", brace_direction=RIGHT, font_size=32)
        with self.vo("Next, mass. Take a thin slab of firn. The amount of ice in it, one minus phi times its thickness, can't change as it moves, "
                     "because the firn isn't flowing sideways. So if the slab loses pore space, it has to get thinner.") as tr:
            self.play(Create(box0), LaggedStart(*[GrowFromCenter(g) for g in grains], lag_ratio=0.03), run_time=1.5)
            self.play(FadeIn(dzl), Write(lab_ice), run_time=1.2)
            self.wait(tr.until(0.7))
            tgt = VGroup(box0.copy().stretch(0.55, 1, about_edge=DOWN))
            new_grains = grains.copy()
            for g in new_grains:
                g.move_to([g.get_x(), box0.get_bottom()[1] + 0.13 + (g.get_y() - box0.get_bottom()[1] - 0.13) * 0.5, 0])
            self.play(box0.animate.stretch(0.55, 1, about_edge=DOWN), Transform(grains, new_grains),
                      dzl.animate.stretch(0.55, 1, about_edge=DOWN), run_time=1.5)

        mc1 = M(r"\DDt{}\big[(1-\ph)\,\Delta z\big] = 0", size=38).move_to([3.7, 1.3, 0])
        mc2 = M(r"\DDt{\Delta z} = \Delta z\,\pd{\ww}{z}", size=38).next_to(mc1, DOWN, buff=0.35)
        mc3 = M(r"\DDt{\ph}", "=", r"(1-\ph)\,\pd{\ww}{z}", size=48).move_to([1.4, -1.6, 0])
        tag6 = eq_tag(6).next_to(mc3, RIGHT, buff=0.4)
        with self.vo("The slab's thickness changes at a rate set by the velocity gradient, partial w partial zed. Expand the derivative, "
                     "and out pops a beautifully simple statement of mass conservation: D phi by D t equals one minus phi, "
                     "times partial w partial zed.") as tr:
            self.play(FadeOut(lab_ice), Write(mc1), run_time=1.3)
            self.play(Write(mc2), run_time=1.3)
            self.wait(tr.until(0.55))
            self.play(FadeOut(VGroup(box0, grains, dzl)), run_time=0.5)
            self.play(TransformMatchingShapes(VGroup(mc1.copy(), mc2.copy()), mc3), FadeIn(tag6), run_time=2)
        self.play(FadeOut(mc1), FadeOut(mc2), VGroup(mc3, tag6).animate.scale(0.75).move_to([1.2, 2.3, 0]), run_time=0.8)

        # ---- combine with compaction law -> velocity gradient
        law = M(r"\DDt{\ph}", "=", r"k_c\,\mathrm{sign}(\sg)|\sg|^n\ph^m", r"(1-\ph)", r"\,\frac{e^{-E_c/RT_s}}{\rsq}", size=40).next_to(mc3, DOWN, buff=0.5)
        comb = M(r"(1-\ph)", r"\pd{\ww}{z}", "=", r"k_c\,\mathrm{sign}(\sg)|\sg|^n\ph^m", r"(1-\ph)", r"\,\frac{e^{-E_c/RT_s}}{\rsq}", size=40).next_to(law, DOWN, buff=0.5)
        vel = M(r"\pd{\ww}{z}", "=", r"k_c\,\mathrm{sign}(\sg)|\sg|^n\ph^m", r"\,\frac{e^{-E_c/RT_s}}{\rsq}", size=44).move_to(comb)
        with self.vo("Now combine this with the compaction law. Both say something about D phi by D t, so set them equal. "
                     "The factors of one minus phi cancel, and we're left with an equation for how the velocity changes with depth.") as tr:
            self.play(FadeIn(law, shift=UP * 0.2), run_time=1)
            self.play(TransformMatchingShapes(VGroup(mc3[2].copy(), law[2:].copy()), comb), run_time=1.6)
            self.wait(tr.until(0.6))
            x1 = Cross(comb[0], stroke_color=RED, stroke_width=4)
            x2 = Cross(comb[4], stroke_color=RED, stroke_width=4)
            self.play(Create(x1), Create(x2), run_time=0.8)
            self.play(FadeOut(VGroup(x1, x2, comb[0], comb[4])),
                      *[ReplacementTransform(comb[i], vel[j]) for i, j in [(1, 0), (2, 1), (3, 2), (5, 3)]], run_time=1.5)

        self.play(FadeOut(VGroup(mc3, tag6, law)), vel.animate.move_to([1.2, 1.9, 0]), run_time=0.8)
        # velocity arrows on the column
        d = load("fig2")
        z, w = d["z"][-1], d["w"][-1]
        arrows = VGroup()
        for f in np.linspace(0.05, 0.95, 8):
            ww = np.interp(f * 0.9, z, w)
            p = np.array([col.get_right()[0] + 0.45, top - f * H, 0])
            arrows.add(Arrow(p + UP * 0.3 * ww / 2, p + DOWN * 0.3 * ww / 2, buff=0, color=C_W, stroke_width=7,
                             max_tip_length_to_length_ratio=0.3))
        with self.vo("Wherever firn is compacting, the downward velocity drops. Fluffy snow near the surface moves down quickly; "
                     "dense ice at depth moves more slowly.") as tr:
            self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.12), run_time=2.5)

        bc = M(r"\ww(0) = \frac{\bb}{1-\ph_s}", size=46).next_to(vel, DOWN, buff=0.7)
        tag7 = eq_tag(7).next_to(bc, RIGHT, buff=0.6)
        box = label_box(bc, color=C_B)
        b_arrow = Arrow([col.get_x(), top + 1.0, 0], [col.get_x(), top + 0.05, 0], buff=0, color=C_B, stroke_width=6)
        b_lab = M(r"\bb", size=40).next_to(b_arrow, LEFT, buff=0.1)
        tau = TauCreature(height=1.5).to_corner(DR)
        with self.vo("And here, at last, is where the accumulation rate enters. The velocity at the surface is set by the snowfall: "
                     "w of zero equals b over one minus phi s. Here b is measured in metres of ice equivalent per year, "
                     "and dividing by one minus phi s converts that into the thickness of fluffy snow that arrives each year.") as tr:
            self.play(Write(bc), FadeIn(tag7), GrowArrow(b_arrow), FadeIn(b_lab), run_time=1.5)
            self.play(Create(box), run_time=0.8)
            self.wait(tr.until(0.6))
            self.play(Indicate(bc, color=C_B, scale_factor=1.1), run_time=1.2)
        with self.vo("This one boundary condition is the only place that b appears, anywhere in the model.") as tr:
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.7)
            b, anim = tau.says("That's it?!", mood="surprised")
            self.play(anim, run_time=1)
            self.play(tau.blink())
        self.wait(0.6)
        self.clear_all()
