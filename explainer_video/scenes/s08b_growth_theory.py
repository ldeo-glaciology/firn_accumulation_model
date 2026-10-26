from common import *


class GrowthTheory(VScene):
    def construct(self):
        title = Text("Two effects of accumulation", font="CMU Serif", font_size=38).to_edge(UP, buff=0.25)
        e26 = VGroup(
            M(r"\dd{\ph}{z} = -\frac{z^n\ph^m(1-\ph)}{\al\,\ww\,\rsq}", size=38),
            M(r"\dd{\ww}{z} = -\frac{z^n\ph^m}{\al\,\rsq}", size=38),
            M(r"\dd{\rsq}{z} = \frac{1}{\ww}", size=38),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to([-3.4, 0.3, 0])
        t26 = eq_tag(26).next_to(e26, RIGHT, buff=0.3)
        assume = M(r"\dl = 0,\qquad \sg = -z", size=34, color=GREY_A).next_to(e26, UP, buff=0.45)
        with self.vo("Back to the ODE model. Set delta to zero, since grains rarely saturate, and replace sigma with minus zed, as before. "
                     "That leaves three equations, which already reproduce the pattern we just saw.") as tr:
            self.play(FadeIn(title), run_time=0.8)
            self.play(Write(assume), run_time=1)
            self.play(LaggedStart(*[FadeIn(e, shift=RIGHT * 0.2) for e in e26], lag_ratio=0.2), FadeIn(t26), run_time=2)

        wb = M(r"\ww = \bt", size=40, color=C_W).move_to([2.6, 2.0, 0])
        g1 = M(r"\dd{\rsq}{z} = \frac{1}{\bt}", size=42).move_to([2.6, 0.9, 0])
        g2 = M(r"\rsq = \rs + \frac{z}{\bt}", size=46).move_to([2.6, -0.4, 0])
        box = label_box(g2, color=C_R)
        with self.vo("Now set w equal to beta. The grain-size equation becomes trivial: d r squared by d zed equals one over beta. "
                     "So r squared equals r s squared, plus zed over beta.") as tr:
            self.play(Write(wb), run_time=0.8)
            self.play(TransformFromCopy(e26[2], g1), run_time=1.5)
            self.wait(tr.until(0.6))
            self.play(TransformMatchingShapes(g1.copy(), g2), Create(box), run_time=1.5)

        # two columns with grains growing with depth
        self.play(FadeOut(VGroup(e26, t26, assume, wb, g1)), VGroup(g2, box).animate.move_to([-4.2, 2.0, 0]), run_time=0.8)
        cols = VGroup()
        for k, (x, beta) in enumerate([(-1.0, 1.0), (2.6, 4.0)]):
            frame = Rectangle(width=2.2, height=4.6, stroke_color=WHITE, stroke_width=2).move_to([x, -0.6, 0]).set_fill(DEEP_ICE, 0.25)
            gr = VGroup()
            rng = np.random.default_rng(7)
            for row in range(9):
                s = (row + 0.5) / 9
                zz = s * 1.0
                rad = 0.07 + 0.16 * np.sqrt(0.02 + zz / beta) / np.sqrt(1.02)
                for j in range(3):
                    gr.add(Circle(radius=rad).set_fill(SNOW, 0.85).set_stroke(GREY_B, 1).move_to(
                        [x - 0.7 + 0.7 * j + rng.uniform(-0.08, 0.08), frame.get_top()[1] - s * 4.6, 0]))
            lab = M(r"\bt = %d" % beta, size=34).next_to(frame, UP, buff=0.2)
            cols.add(VGroup(frame, gr, lab))
        same = DashedLine([-2.4, -0.6, 0], [4.0, -0.6, 0], color=YELLOW, stroke_width=2)
        same_l = Text("same depth", font="CMU Serif", font_size=22, color=YELLOW).next_to(same, RIGHT, buff=0.1)
        with self.vo("Think about what that says. With higher accumulation, a parcel reaches any given depth sooner, so its grains have had "
                     "less time to grow. Smaller grains. And smaller grains, remember, compact faster.") as tr:
            self.play(FadeIn(cols[0][0]), FadeIn(cols[0][2]), LaggedStart(*[GrowFromCenter(g) for g in cols[0][1]], lag_ratio=0.02), run_time=1.5)
            self.play(FadeIn(cols[1][0]), FadeIn(cols[1][2]), LaggedStart(*[GrowFromCenter(g) for g in cols[1][1]], lag_ratio=0.02), run_time=1.5)
            self.play(Create(same), FadeIn(same_l), run_time=0.8)
            self.wait(tr.until(0.8))
            self.play(Indicate(cols[1][1][12:15], color=C_R), run_time=1.2)
        self.play(FadeOut(VGroup(cols, same, same_l, g2, box)), run_time=0.7)

        # Eq 27 with overbraces
        e27 = M(r"\dd{\ph}{z} = -\frac{z^n\ph^m(1-\ph)}{\al}\;", r"\overbrace{\frac{1}{\bt}}^{\text{from } w}", r"\;",
                r"\overbrace{\frac{1}{\left(\frac{z}{\bt}+\rs\right)}}^{\text{from } r^2}", size=46).move_to([0, 1.4, 0])
        e27[1].set_color(C_W)
        e27[3].set_color(C_R)
        bw = SurroundingRectangle(e27[1], color=C_W, buff=0.08)
        br = SurroundingRectangle(e27[3], color=C_R, buff=0.08)
        lw = Text("porosity advection:\nthicker firn", font="CMU Serif", font_size=24, color=C_W, line_spacing=0.9).next_to(bw, DOWN, buff=0.35).shift(LEFT * 1.0)
        lr = Text("grain-size advection:\nsmaller grains, thinner firn", font="CMU Serif", font_size=24, color=C_R, line_spacing=0.9).next_to(br, DOWN, buff=0.35).shift(RIGHT * 0.8)
        t27 = eq_tag(27).next_to(e27, RIGHT, buff=0.3)
        with self.vo("Substitute both into the porosity equation, and look at what we get. Beta appears twice. Once from the advection of porosity: "
                     "this one over beta, which on its own would make the firn thicker. And once from the advection of grain size, inside this "
                     "second fraction, which makes the grains smaller, and the firn thinner. Two effects of accumulation, pulling in opposite "
                     "directions.") as tr:
            self.play(Write(e27), FadeIn(t27), run_time=2.5)
            self.wait(tr.until(0.32))
            self.play(Create(bw), FadeIn(lw), run_time=1)
            self.wait(tr.until(0.6))
            self.play(Create(br), FadeIn(lr), run_time=1)
        e27b = M(r"\dd{\ph}{z} = -\frac{z^{n-1}\ph^m(1-\ph)}{\al\left(1+", r"\dfrac{\bt\,\rs}{z}", r"\right)}", size=48).move_to([0, -1.4, 0])
        with self.vo("Multiply out, and beta now only ever appears multiplied by r s squared, and compared with zed.") as tr:
            self.play(FadeOut(lw), FadeOut(lr), TransformMatchingShapes(e27.copy(), e27b), run_time=2)
            self.play(Circumscribe(e27b[1], color=C_B), run_time=1.2)

        # punchline: rs2 = 0
        tau = TauCreature(height=1.4).move_to([5.9, -2.9, 0])
        zero = M(r"\rs = 0", size=40, color=YELLOW).next_to(e27b, DOWN, buff=0.5).shift(LEFT * 2.5)
        e_final = M(r"\dd{\ph}{z} = -\frac{z^{n-1}\ph^m(1-\ph)}{\al}", size=48).move_to(e27b)
        with self.vo("And now, the punchline. If the surface grain size is zero, beta disappears completely. The porosity becomes a function "
                     "of depth alone, and accumulation has no effect on firn thickness whatsoever.") as tr:
            self.play(FadeOut(VGroup(bw, br)), Write(zero), run_time=1)
            cross = Cross(e27b[1], stroke_width=5)
            self.play(Create(cross), run_time=0.8)
            self.play(FadeOut(cross), TransformMatchingShapes(e27b, e_final), run_time=1.5)
            self.play(FadeIn(tau, shift=UP * 0.3), Circumscribe(e_final, color=YELLOW), run_time=1.2)
        with self.vo("None whatsoever. Not a weaker dependence: no dependence at all.") as tr:
            b, anim = tau.says("None at all?!", mood="surprised")
            self.play(anim, run_time=1)
            self.play(tau.blink())
        self.play(FadeOut(b), FadeOut(VGroup(e27, t27, zero)), e_final.animate.scale(0.75).to_corner(UL).shift(DOWN * 0.8), tau.change("plain"), run_time=0.8)

        # numerical check from the full model
        dc = load("appC")
        zg = dc["zgrid"]
        pl = []
        for k, (tag, lab) in enumerate([("small", r"\rs = 10^{-4}"), ("large", r"\rs = 0.1")]):
            p = DepthPlot([0, 0.5], [0, 1.0], width=3.0, height=3.6, x_ticks=[0, 0.25, 0.5], z_ticks=[0, 0.5, 1.0], x_dec=2, z_dec=1,
                          x_label=r"\ph", z_label=(r"z" if k == 0 else None), num_size=20, label_size=30)
            p.shift(np.array([-1.0 + 4.3 * k, 1.5, 0]) - p.top.get_start())
            curves = VGroup(*[p.curve(dc["phi_" + tag][i], zg, color=interpolate_color(ManimColor("#1E4F7A"), ManimColor(C_PHI), i / 4), width=3)
                              for i in range(len(dc["beta"]))])
            t_ = M(lab, size=30).next_to(p, DOWN, buff=0.2)
            pl.append((p, curves, t_))
        bl = M(r"\bt = 0.5,\ 1,\ 2,\ 5,\ 10", size=28, color=C_PHI).move_to([1.2, -3.4, 0])
        hd = Text("full model", font="CMU Serif", font_size=28, color=GREY_A).move_to([1.2, 2.9, 0])
        with self.vo("The full model agrees. With almost no grain size at the surface, five very different accumulation rates give almost identical "
                     "porosity profiles. With large surface grains, they spread right out.") as tr:
            self.play(FadeIn(hd), FadeIn(pl[0][0]), FadeIn(pl[0][2]), LaggedStart(*[Create(c) for c in pl[0][1]], lag_ratio=0.2), FadeIn(bl), run_time=2.5)
            self.wait(tr.until(0.6))
            self.play(FadeIn(pl[1][0]), FadeIn(pl[1][2]), LaggedStart(*[Create(c) for c in pl[1][1]], lag_ratio=0.2), run_time=2)
        self.play(FadeOut(VGroup(*[m for trip in pl for m in trip], bl, hd)), run_time=0.7)

        # intuition with one parcel
        steps = VGroup(
            M(r"|\sg| \approx z = \bt\, t", size=40),
            M(r"\rsq = t \qquad (\rs = 0)", size=40),
            M(r"\text{compaction rate} \;\propto\; \frac{|\sg|}{\rsq} = \frac{\bt\, t}{t} = \bt", size=40),
            M(r"\text{total compaction} \;\propto\; \bt\, t = z", size=40),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.45).move_to([1.0, -0.2, 0])
        with self.vo("Here's the picture for a single parcel. Its stress grows in proportion to its depth, beta t. But with no initial size, its grain "
                     "size squared also grows in proportion to time: t. So its compaction rate, stress over r squared, is just beta: constant. "
                     "Multiply by the time it's been travelling, and the total compaction is beta t, which is simply its depth. "
                     "Whatever the accumulation rate, a parcel has compacted by exactly the same amount by the time it reaches a given depth.") as tr:
            for i, f in enumerate([0.0, 0.2, 0.36, 0.6]):
                self.wait(tr.until(f))
                self.play(FadeIn(steps[i], shift=RIGHT * 0.2), run_time=1)
            self.wait(tr.until(0.8))
            self.play(Circumscribe(steps[3], color=YELLOW), tau.change("happy", steps), run_time=1.2)
        self.play(FadeOut(steps), FadeOut(e_final), run_time=0.6)

        # finite rs2: crossover depth beta*rs2
        ax = Axes(x_range=[0, 1, 0.25], y_range=[0, 1.1, 0.25], x_length=5.0, y_length=3.6, tips=False,
                  axis_config=dict(color=GREY_B, stroke_width=2)).move_to([-2.6, 0.0, 0])
        axl = VGroup(M(r"z", size=30).next_to(ax.x_axis, DOWN, buff=0.2), M(r"\rsq", size=30).next_to(ax.y_axis, LEFT, buff=0.2))
        rs = 0.15
        lines = VGroup()
        for beta, c in [(1.0, C_R), (4.0, "#E6D2FF")]:
            lines.add(ax.plot(lambda zz, b=beta: rs + zz / b, x_range=[0, 0.95], color=c, stroke_width=4))
        lbl = VGroup(M(r"\bt = 1", size=28, color=C_R).next_to(ax.c2p(0.95, rs + 0.95), RIGHT, buff=0.1),
                     M(r"\bt = 4", size=28, color="#E6D2FF").next_to(ax.c2p(0.95, rs + 0.95 / 4), RIGHT, buff=0.1))
        rsl = DashedLine(ax.c2p(0, rs), ax.c2p(1, rs), color=GREY_B)
        rsl_l = M(r"\rs", size=28, color=GREY_A).next_to(ax.c2p(0, rs), LEFT, buff=0.1)
        xover = Dot(ax.c2p(4 * rs, 2 * rs), color=YELLOW)
        xover_l = M(r"z \sim \bt\,\rs", size=30, color=YELLOW).next_to(xover, UP + LEFT, buff=0.1)
        reg = Text("surface grains dominate:\nconstant-grain-size behaviour,\nsquare-root dependence on β", font="CMU Serif", font_size=24,
                   line_spacing=0.9, color=GREY_A).move_to([3.6, 0.6, 0])
        with self.vo("With grains of finite size at the surface, the balance breaks. Near the surface, the grain size is dominated by r s squared, "
                     "which knows nothing about accumulation. So there we're back in the constant-grain-size world, with its square-root dependence "
                     "on beta. That region extends down to a depth of about beta times r s squared. So the larger beta, or r s squared, "
                     "the more of the firn lives in that regime, and the more its thickness depends on accumulation.") as tr:
            self.play(Create(ax), FadeIn(axl), Create(rsl), FadeIn(rsl_l), run_time=1.2)
            self.play(Create(lines), FadeIn(lbl), run_time=1.5)
            self.wait(tr.until(0.3))
            self.play(FadeIn(reg), run_time=1)
            self.wait(tr.until(0.55))
            self.play(FadeIn(xover, scale=2), Write(xover_l), run_time=1)
        self.play(FadeOut(VGroup(ax, axl, lines, lbl, rsl, rsl_l, xover, xover_l, reg)), run_time=0.6)

        cf = M(r"\zc - \bt\rs\,\ln\!\left(1+\frac{\zc}{\bt\rs}\right) = \al L", size=46).move_to([0, 1.6, 0])
        cf_note = Text("(n = m = 1; follows directly from Eq. 27)", font="CMU Serif", font_size=20, color=GREY_B).next_to(cf, DOWN, buff=0.2)
        lim1 = M(r"\bt\rs \ll \zc:", r"\quad \zc \approx \al L", r"\quad\text{(no } \bt \text{ at all)}", size=36)
        lim2 = M(r"\bt\rs \gg \zc:", r"\quad \zc \approx \sqrt{2\al\bt\rs L}", r"\quad\text{(square-root law)}", size=36)
        lims = VGroup(lim1, lim2).arrange(DOWN, aligned_edge=LEFT, buff=0.45).move_to([0, -0.9, 0])
        with self.vo("With n equal to one, this simplest model can even be solved exactly, giving this equation for [z830](zed eight-thirty). "
                     "When beta r s squared is small, it gives alpha L: no beta at all. When it's large, we recover the square-root law from before. "
                     "A smooth hand-over between the two regimes.") as tr:
            self.play(Write(cf), FadeIn(cf_note), run_time=2.5)
            self.wait(tr.until(0.45))
            self.play(FadeIn(lim1, shift=UP * 0.2), run_time=1)
            self.wait(tr.until(0.68))
            self.play(FadeIn(lim2, shift=UP * 0.2), run_time=1)
        self.play(FadeOut(VGroup(cf, cf_note, lims, title)), FadeOut(tau), run_time=0.7)

        # Fig 6
        common_kw = dict(x_range=(0.1, 10), y_range=(0.001, 0.1), x_len=3.3, y_len=3.2, x_ticks=[2, 4, 6, 8, 10],
                         y_ticks=[0.02, 0.04, 0.06, 0.08, 0.1], x_dec=0, y_dec=2, num_size=18, label_size=24)
        p1 = MapPanel("fig6_full", title=r"\text{full model}", x_label=r"\bt", y_label=r"\rs", **common_kw)
        p2 = MapPanel("fig6_eq26", title=r"\text{ODE model (Eq. 26)}", x_label=r"\bt", **common_kw)
        p3 = MapPanel("fig6_eq27", title=r"\text{closed form (Eq. 27)}", x_label=r"\bt", **common_kw)
        maps = Group(p1, p2, p3).arrange(RIGHT, buff=0.55).move_to([0.2, 0.0, 0])
        hd6 = M(r"\text{firn thickness } \zc(\bt, \rs)", size=36).to_edge(UP, buff=0.25)
        with self.vo("Here's firn thickness across the whole plane of accumulation and surface grain size: from the full model, from the ODE model, "
                     "and from that closed form. The numbers differ, but the pattern is the same in all three.") as tr:
            self.play(FadeIn(hd6), FadeIn(p1), run_time=1.2)
            self.wait(tr.until(0.35))
            self.play(FadeIn(p2), run_time=1)
            self.wait(tr.until(0.55))
            self.play(FadeIn(p3), run_time=1)
        arr_lo = Arrow(p1.c2p(1, 0.008), p1.c2p(9.5, 0.008), buff=0, color=WHITE, stroke_width=4)
        arr_hi = Arrow(p1.c2p(1, 0.09), p1.c2p(9.5, 0.09), buff=0, color=WHITE, stroke_width=4)
        with self.vo("Along the bottom, where surface grains are small, you can move a long way in beta without crossing a single contour. "
                     "Near the top, the contours come thick and fast. As the surface grain size increases, firn thickness depends more "
                     "and more strongly on accumulation. And the picture is roughly symmetric, because beta and r s squared work together, "
                     "through their product.") as tr:
            self.play(GrowArrow(arr_lo), run_time=1.2)
            self.wait(tr.until(0.3))
            self.play(GrowArrow(arr_hi), run_time=1.2)
            self.wait(tr.until(0.75))
            self.play(Indicate(p3.frame), run_time=1)
        self.wait(0.5)
        self.clear_all()
