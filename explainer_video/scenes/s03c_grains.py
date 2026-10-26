from common import *


class GrainsAndThickness(VScene):
    def construct(self):
        # ---- normal grain growth: coarsening cluster
        rng = np.random.default_rng(11)
        ctr = np.array([-3.8, 0.3, 0])
        grains = VGroup()
        pos = []
        while len(pos) < 26:
            p = ctr + np.array([rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5), 0])
            if all(np.linalg.norm(p - q) > 0.5 for q in pos):
                pos.append(p)
        sizes = rng.uniform(0.12, 0.26, len(pos))
        for p, r in zip(pos, sizes):
            grains.add(Circle(radius=r).move_to(p).set_fill(SNOW, 0.9).set_stroke(GREY_B, 1.5))
        title = Text("Normal grain growth", font="CMU Serif", font_size=36).move_to([-3.8, 2.8, 0])
        ax = Axes(x_range=[0, 1, 0.5], y_range=[0, 1, 0.5], x_length=3.6, y_length=2.6, tips=False,
                  axis_config=dict(color=GREY_B, stroke_width=2)).move_to([2.8, 0.8, 0])
        axl = VGroup(M(r"t", size=30).next_to(ax.x_axis, DOWN, buff=0.15),
                     M(r"\rsq", size=30).next_to(ax.y_axis, LEFT, buff=0.15))
        lin = ax.plot(lambda t: 0.1 + 0.8 * t, color=C_R, stroke_width=4)
        with self.vo("Next, grain size. Snow grains grow. In what's called normal grain growth, small grains are gradually consumed by their larger "
                     "neighbours, driven by the energy stored in the boundaries between grains, and the square of the mean grain radius "
                     "increases linearly with time, at a rate that, once again, follows an Arrhenius law.") as tr:
            self.play(Write(title), LaggedStart(*[GrowFromCenter(g) for g in grains], lag_ratio=0.03), run_time=1.6)
            order = np.argsort(sizes)
            for stage in range(3):
                kill = order[stage * 5:(stage + 1) * 5]
                anims = []
                for i, g in enumerate(grains):
                    if i in kill:
                        anims.append(g.animate.scale(0.01))
                    elif i in order[(stage + 1) * 5:]:
                        anims.append(g.animate.scale(1.12))
                self.play(*anims, run_time=1.4)
                if stage == 0:
                    self.play(Create(ax), FadeIn(axl), run_time=0.8)
                    self.play(Create(lin), run_time=1.2)
            self.fill(tr)

        sat = ax.plot(lambda t: 0.85 - 0.75 * np.exp(-2.2 * t), color=C_R, stroke_width=4)
        rf_line = DashedLine(ax.c2p(0, 0.85), ax.c2p(1, 0.85), color=GREY_B)
        rf_lab = M(r"\rf", size=30).next_to(rf_line, RIGHT, buff=0.1)
        eq8 = M(r"\pd{\rsq}{t}", "=", r"k_g\,e^{-E_g/RT_s}\,\big(\rf - \rsq\big)", r"\;-\;\ww\pd{\rsq}{z}", size=44).move_to([0.5, -2.2, 0])
        bc8 = M(r"\rsq(0) = \rs", size=36).next_to(eq8, DOWN, buff=0.25)
        tag8 = eq_tag(8).next_to(eq8, RIGHT, buff=0.3)
        bg = Brace(eq8[2], UP, buff=0.08, color=C_R)
        bgt = Text("growth, with saturation", font="CMU Serif", font_size=22, color=C_R).next_to(bg, UP, buff=0.08)
        ba = Brace(eq8[3], UP, buff=0.08, color=C_W)
        bat = Text("advection", font="CMU Serif", font_size=22, color=C_W).next_to(ba, UP, buff=0.08)
        with self.vo("Grains can't keep growing for ever, so the model adds a saturation size, r f, with growth slowing down as r squared "
                     "approaches r f squared. And since grains are carried downward with the firn, there's an advection term too. "
                     "At the surface, the grains have a prescribed size, r s squared. Keep an eye on that number. "
                     "It's the second key input, alongside the accumulation rate.") as tr:
            self.play(Transform(lin, sat), Create(rf_line), FadeIn(rf_lab), run_time=1.5)
            self.wait(tr.until(0.22))
            self.play(FadeOut(grains), FadeOut(title), VGroup(ax, axl, lin, rf_line, rf_lab).animate.scale(0.8).move_to([4.6, 1.8, 0]),
                      run_time=1)
            self.play(Write(eq8[:3]), FadeIn(tag8), run_time=1.5)
            self.play(GrowFromCenter(bg), FadeIn(bgt), run_time=0.6)
            self.wait(tr.until(0.45))
            self.play(Write(eq8[3]), GrowFromCenter(ba), FadeIn(bat), run_time=1.2)
            self.wait(tr.until(0.65))
            self.play(Write(bc8), run_time=1)
            self.play(Circumscribe(bc8, color=C_R), run_time=1.2)
        grain_block = VGroup(eq8, bc8, tag8)
        self.play(FadeOut(VGroup(bg, bgt, ba, bat, ax, axl, lin, rf_line, rf_lab)),
                  grain_block.animate.scale(0.7).to_corner(UL).shift(DOWN * 0.1), run_time=0.8)

        eq9 = M(r"\pd{\Ag}{t} = 1 - \ww\pd{\Ag}{z}, \qquad \Ag(0) = 0", size=36).next_to(grain_block, DOWN, buff=0.35).align_to(grain_block, LEFT)
        tag9 = eq_tag(9).next_to(eq9, RIGHT, buff=0.3)
        with self.vo("There's also an equation for age, which says that, as it's carried downward, firn gets one year older every year.") as tr:
            self.play(Write(eq9), FadeIn(tag9), run_time=2)

        # ---- Appendix A: kinematic condition for h
        hd = Text("How thick is the domain?", font="CMU Serif", font_size=30, color=GREY_A).move_to([2.2, 2.7, 0])
        steps = [
            r"\pd{}{t}\left[\int_0^{z_b}(1-\ph)\,\mathrm{d}z\right] = 0",
            r"(1-\ph_b)\,\dot h \;-\; \int_0^{z_b}\pd{\ph}{t}\,\mathrm{d}z = 0",
            r"\pd{\ph}{t} = \pd{}{z}\big[(1-\ph)\,\ww\big]",
            r"\int_0^{z_b}\pd{\ph}{t}\,\mathrm{d}z = (1-\ph_b)\,\ww(z_b) - \underbrace{(1-\ph_s)\,\ww(0)}_{\bb}",
            r"(1-\ph_b)\,\dot h - (1-\ph_b)\,\ww(z_b) + \bb = 0",
            r"\pd{h}{t} = \ww(z_b) - \frac{\bb}{1-\ph_b}",
        ]
        eqs = [M(s, size=38) for s in steps]
        x0 = 2.2
        eqs[0].move_to([x0, 1.6, 0])
        eqs[1].move_to([x0, 0.4, 0])
        eqs[2].move_to([x0, -0.7, 0])
        eqs[3].move_to([x0, -1.9, 0])
        with self.vo("Finally, the base of the domain is free to move. Take the total amount of ice in the column, the integral of one minus phi, "
                     "and demand that it stays fixed: whatever snow arrives at the top must leave through the base.") as tr:
            self.play(FadeIn(hd), run_time=0.6)
            self.play(Write(eqs[0]), run_time=2)
        with self.vo("Differentiate under the integral using the [Leibniz](Leibniz) rule. One term comes from the moving base, at speed h dot. "
                     "The other is the integral of partial phi by partial t.") as tr:
            self.play(TransformMatchingShapes(eqs[0].copy(), eqs[1]), run_time=2)
        with self.vo("Mass conservation tells us that partial phi by partial t is the depth derivative of the ice flux, one minus phi times w. "
                     "So that integral is just the flux at the base, minus the flux at the surface. And the flux at the surface is exactly b.") as tr:
            self.play(Write(eqs[2]), run_time=1.5)
            self.wait(tr.until(0.45))
            self.play(Write(eqs[3]), run_time=2.5)
        eqs[4].move_to([x0, 0.4, 0])
        eqs[5].move_to([x0, -0.9, 0])
        tag10 = eq_tag(10).next_to(eqs[5], RIGHT, buff=0.4)
        box10 = label_box(eqs[5], color=WHITE)
        with self.vo("Substitute, rearrange, and we arrive at equation ten: h dot equals the velocity at the base, minus the accumulation rate "
                     "converted to the porosity at the base. In a steady state these two balance. If they don't, the column thickens or thins, "
                     "which we'll use later to mimic a rising or falling ice-sheet surface.") as tr:
            self.play(FadeOut(eqs[0]), FadeOut(eqs[2]), ReplacementTransform(VGroup(eqs[1], eqs[3]), eqs[4]), run_time=2)
            self.play(TransformMatchingShapes(eqs[4].copy(), eqs[5]), FadeIn(tag10), run_time=1.8)
            self.play(Create(box10))
            self.fill(tr)

        # ---- summary of the full model
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=0.8)
        rows = [
            (r"\pd{\ph}{t} = k_c\frac{\mathrm{sign}(\sg)|\sg|^n\ph^m(1-\ph)}{\rsq}e^{-E_c/RT_s} - \ww\pd{\ph}{z}", 4),
            (r"\pd{\sg}{z} = -\rho_i g(1-\ph)", 5),
            (r"\pd{\ww}{z} = k_c\frac{\mathrm{sign}(\sg)|\sg|^n\ph^m}{\rsq}e^{-E_c/RT_s}", 7),
            (r"\pd{\rsq}{t} = k_g e^{-E_g/RT_s}(\rf-\rsq) - \ww\pd{\rsq}{z}", 8),
            (r"\pd{\Ag}{t} = 1-\ww\pd{\Ag}{z}", 9),
            (r"\pd{h}{t} = \ww(z_b)-\frac{\bb}{1-\ph_b}", 10),
        ]
        sys_ = VGroup()
        for s, n in rows:
            e = M(s, size=32)
            t = eq_tag(n, size=22).next_to(e, RIGHT, buff=0.3)
            sys_.add(VGroup(e, t))
        sys_.arrange(DOWN, aligned_edge=LEFT, buff=0.28).move_to([-1.6, 0.0, 0])
        inputs = VGroup(
            Text("Inputs", font="CMU Serif", font_size=32),
            M(r"\ph_s", r"\text{ surface porosity}", size=32),
            M(r"\rs", r"\text{ surface grain size}", size=32),
            M(r"T_s", r"\text{ temperature}", size=32),
            M(r"\bb", r"\text{ accumulation rate}", size=32),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([4.6, 0.8, 0])
        frame = SurroundingRectangle(inputs, color=GREY_B, buff=0.3, corner_radius=0.1)
        tau = TauCreature(height=1.4).to_corner(DR).shift(UP * 0.1)
        with self.vo("And that's the full model. Equations four to ten, describing six variables, which respond to just four inputs: "
                     "the surface porosity, the surface grain size, the temperature, and the accumulation rate.") as tr:
            self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.2) for r in sys_], lag_ratio=0.2), run_time=3)
            self.wait(tr.until(0.6))
            self.play(Create(frame), LaggedStart(*[FadeIn(i) for i in inputs], lag_ratio=0.25), FadeIn(tau, shift=UP * 0.2), run_time=2)
            self.play(Indicate(inputs[2], color=C_R), Indicate(inputs[4], color=C_B), tau.change("pondering", inputs), run_time=1.2)
        self.wait(0.5)
        self.clear_all()
