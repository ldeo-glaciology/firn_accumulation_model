from common import *


class Scaling(VScene):
    def construct(self):
        card = chapter_card("Scaling", 2)
        self.play(FadeIn(card, shift=UP * 0.3))
        self.wait(0.8)
        self.play(FadeOut(card))

        sc = M(r"\ph=\ph^*\ph_0,\quad \sg=\sg^*\sigma_0,\quad \ww=\ww^*w_0,\quad \rsq=r^{2*}r_0^2,\quad \Ag=\Ag^*A_0,",
               size=38).move_to([0, 2.7, 0])
        sc2 = M(r"h=h^*h_0,\quad \bb=\bt\, b_0,\quad z=z^*z_0,\quad t=t^*t_0", size=38).next_to(sc, DOWN, buff=0.25)
        tag11 = eq_tag(11).next_to(sc2, RIGHT, buff=0.4)
        with self.vo("Now for one of the most powerful tricks in applied mathematics: scaling. Write every variable as a dimensionless number, "
                     "marked with a star, times a characteristic scale, marked with a zero.") as tr:
            self.play(Write(sc), run_time=2.5)
            self.play(Write(sc2), FadeIn(tag11), run_time=2)

        scales = VGroup(
            M(r"z_0 = h_0 = 100\ \mathrm{m}, \qquad \ph_0 = 1", size=36),
            M(r"w_0 = b_0", size=36),
            M(r"t_0 = \frac{z_0}{w_0}", r"\quad\text{(transit time)}", size=36),
            M(r"\sigma_0 = \rho_i\, g\, z_0", r"\quad\text{(weight of 100 m of ice)}", size=36),
            M(r"\bt = \bb / b_0", r"\quad\text{(dimensionless accumulation)}", size=36),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([0, -0.6, 0])
        with self.vo("Depth is scaled by zed zero, a hundred metres, and porosity by one. Velocity is scaled by the accumulation scale, b zero. "
                     "That makes the time scale, t zero, equal to zed zero over b zero: the time it takes firn to travel through the domain. "
                     "The stress scale is the weight of a hundred metres of ice, rho i g zed zero. And the dimensionless accumulation rate gets "
                     "its own symbol, beta. Beta equals one means the accumulation is exactly at the reference value, b zero.") as tr:
            for i, f in enumerate([0.0, 0.17, 0.3, 0.55, 0.72]):
                self.wait(tr.until(f))
                self.play(FadeIn(scales[i], shift=UP * 0.2), run_time=0.7)
        self.play(FadeOut(scales), VGroup(sc, sc2, tag11).animate.scale(0.6).to_edge(UP, buff=0.2), run_time=0.8)

        e4 = M(r"\pd{\ph}{t} = k_c\frac{\mathrm{sign}(\sg)|\sg|^n\ph^m(1-\ph)}{\rsq}e^{-E_c/RT_s} - \ww\pd{\ph}{z}", size=40).move_to([0, 1.3, 0])
        e4b = M(r"\frac{1}{t_0}\pd{\ph^*}{t^*} = -\frac{k_c\,\sigma_0^n\, e^{-E_c/RT_s}}{r_0^2}\,\frac{|\sg^*|^n\ph^{*m}(1-\ph^*)}{r^{2*}}"
                r" - \frac{w_0}{z_0}\,\ww^*\pd{\ph^*}{z^*}", size=40).move_to([0, -0.1, 0])
        e12 = M(r"\pd{\ph}{t}", "=", r"-\frac{1}{\al}", r"\frac{|\sg|^n\ph^m(1-\ph)}{\rsq}", r"-\ww\pd{\ph}{z}", size=48).move_to([0, -0.1, 0])
        e13 = M(r"\al = \frac{r_0^2}{k_c\, t_0\, \sigma_0^n\, e^{-E_c/RT_s}}", size=46).move_to([0, -2.0, 0])
        t12 = eq_tag(12).next_to(e12, RIGHT, buff=0.4)
        t13 = eq_tag(13).next_to(e13, RIGHT, buff=0.4)
        with self.vo("Substitute these into the porosity equation. Multiply through by t zero, and remember that w zero t zero is just zed zero. "
                     "Every dimensional constant in the compaction term gathers into one single dimensionless group.") as tr:
            self.play(Write(e4), run_time=1.5)
            self.play(TransformMatchingShapes(e4.copy(), e4b), run_time=2)
            self.wait(tr.until(0.75))
            self.play(FadeOut(e4), TransformMatchingShapes(e4b, e12), FadeIn(t12), run_time=2)
            self.play(Write(e13), FadeIn(t13), run_time=1.5)
        box = label_box(e12[2], color=C_AL)
        name = Text("compaction number", font="CMU Serif", font_size=30, color=C_AL).next_to(e13, DOWN, buff=0.3)
        with self.vo("That group is alpha, which the authors call the compaction number. It compares how fast firn is carried through the domain "
                     "with how fast it compacts. Small alpha means compaction wins, and the porosity drops to nearly zero well within the domain.") as tr:
            self.play(Create(box), FadeIn(name), run_time=1)
            self.play(Indicate(e13, color=C_AL, scale_factor=1.05), run_time=1.2)
        self.play(FadeOut(VGroup(e12, t12, box)), VGroup(e13, t13, name).animate.scale(0.75).to_corner(UL).shift(DOWN * 0.9), run_time=0.8)

        e14 = M(r"\pd{\sg}{z} = -(1-\ph), \qquad \sg(0)=0", size=38)
        e15 = M(r"\pd{\ww}{z} = -\frac{1}{\al}\frac{|\sg|^n\ph^m}{\rsq}, \qquad \ww(0) = \frac{\bt}{1-\ph_s}", size=38)
        e16 = M(r"r_0^2 = \frac{k_g\, z_0\, r_f^2}{b_0}\, e^{-E_g/RT_s}", size=38)
        e17 = M(r"\pd{\rsq}{t} = (1-\dl\,\rsq) - \ww\pd{\rsq}{z}, \qquad \dl = \frac{r_0^2}{r_f^2}", size=38)
        grp = VGroup(VGroup(e14, eq_tag(14)), VGroup(e15, eq_tag(15)), VGroup(e16, eq_tag(16)), VGroup(e17, eq_tag(17)))
        for g in grp:
            g[1].next_to(g[0], RIGHT, buff=0.4)
        grp.arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to([1.2, -0.6, 0])
        with self.vo("The other equations scale neatly too. The stress gradient becomes just minus one minus phi. "
                     "The velocity equation picks up the same factor of one over alpha, and its surface value becomes beta over one minus phi s.") as tr:
            self.play(FadeIn(grp[0], shift=UP * 0.2), run_time=1)
            self.wait(tr.until(0.4))
            self.play(FadeIn(grp[1], shift=UP * 0.2), run_time=1)
        with self.vo("But what's r zero, the grain-size scale? It's chosen by balancing the terms in the grain-growth equation, "
                     "so it's roughly how big a grain gets during one transit time. That leaves a tidy grain equation with a single parameter, "
                     "delta: the ratio of r zero squared to the saturation size, which measures how close grains get to saturating.") as tr:
            self.play(FadeIn(grp[2], shift=UP * 0.2), run_time=1)
            self.wait(tr.until(0.5))
            self.play(FadeIn(grp[3], shift=UP * 0.2), run_time=1)
        self.play(FadeOut(grp), FadeOut(VGroup(sc, sc2, tag11)), VGroup(e13, t13, name).animate.move_to([0, 2.3, 0]).scale(1.2), run_time=0.8)

        # --- the b0 cancellation
        num = M(r"k_g", r"z_0", r"r_f^2\, e^{-E_g/RT_s}", r"\;/\;", r"\bb_0", size=44)
        den = M(r"k_c", r"z_0", r"\;/\;", r"\bb_0", r"\;\sigma_0^n\, e^{-E_c/RT_s}", size=44)
        bar = Line(LEFT, RIGHT, stroke_width=2.5).set_width(max(num.width, den.width) + 0.2)
        frac = VGroup(num, bar, den).arrange(DOWN, buff=0.15)
        lhs = M(r"\al =", size=44).next_to(frac, LEFT, buff=0.2)
        sub = VGroup(lhs, frac).move_to([0, 0.0, 0])
        subs_note = M(r"r_0^2 = \frac{k_g z_0 r_f^2}{\bb_0}e^{-E_g/RT_s}, \qquad t_0 = \frac{z_0}{\bb_0}", size=32, color=GREY_A).next_to(sub, DOWN, buff=0.5)
        e20 = M(r"\al = \frac{k_g\, r_f^2}{k_c\,\sigma_0^n}\, \exp\!\left[\frac{E_c - E_g}{R\,T_s}\right]", size=48).move_to([0, -0.2, 0])
        t20 = eq_tag(20).next_to(e20, RIGHT, buff=0.4)
        tau = TauCreature(height=1.5).to_corner(DL)
        with self.vo("Now something lovely happens. Substitute the expressions for r zero and t zero into alpha. The zed zeros cancel. "
                     "And so does b zero, the accumulation scale! Alpha depends only on temperature.") as tr:
            self.play(FadeIn(subs_note), run_time=0.8)
            self.play(Write(sub), run_time=2)
            self.wait(tr.until(0.45))
            self.play(Create(Cross(num[1], stroke_width=4)), Create(Cross(den[1], stroke_width=4)), run_time=0.7)
            self.wait(tr.until(0.62))
            self.play(Create(Cross(num[4], stroke_width=4, stroke_color=C_B)), Create(Cross(den[3], stroke_width=4, stroke_color=C_B)), run_time=0.7)
            crosses = [m for m in self.mobjects if isinstance(m, Cross)]
            self.wait(tr.until(0.8))
            self.play(FadeOut(VGroup(sub, subs_note, *crosses)), FadeIn(e20), FadeIn(t20), FadeIn(tau, shift=UP * 0.3), run_time=1.2)
        with self.vo("Physically: higher accumulation gives each parcel less time to compact. But it also gives the grains less time to grow, "
                     "and smaller grains compact faster. In the scales, those two effects cancel exactly. "
                     "Hold on to that thought, because it's a preview of the main result.") as tr:
            b, anim = tau.thinks(r"less time, but\\smaller grains...", mood="pondering")
            self.play(anim, run_time=1.2)
            self.wait(tr.until(0.8))
            self.play(tau.change("happy"), FadeOut(b), run_time=0.8)
        self.play(FadeOut(tau), FadeOut(VGroup(e13, t13, name)), VGroup(e20, t20).animate.scale(0.7).to_edge(UP, buff=0.3), run_time=0.8)

        # --- Fig 1
        ylab = [r"10^{-2}", r"10^{-1}", r"10^{0}"]
        pa = MapPanel("fig1_alpha", (233, 273), (-2, 0), x_len=5.0, y_len=3.8, x_ticks=[240, 250, 260, 270], y_ticks=[-2, -1, 0],
                      y_tick_labels=ylab, x_label=r"\text{surface temperature } T_s\ [\mathrm{K}]",
                      y_label=r"b_0\ [\mathrm{m\,a^{-1}}]", title=r"\text{compaction number } \al")
        pd_ = MapPanel("fig1_delta", (233, 273), (-2, 0), x_len=5.0, y_len=3.8, x_ticks=[240, 250, 260, 270], y_ticks=[-2, -1, 0],
                       y_tick_labels=ylab, x_label=r"\text{surface temperature } T_s\ [\mathrm{K}]",
                       y_label=r"b_0\ [\mathrm{m\,a^{-1}}]", title=r"\log_{10}\dl")
        pa.move_to([0.3, -0.9, 0])
        with self.vo("Temperature enters through the difference between the two activation energies. Warmer firn compacts faster, "
                     "but it also grows bigger grains, which compact more slowly. Because the activation energy for compaction is larger, "
                     "compaction wins, and alpha falls as the temperature rises. And notice: no dependence on accumulation at all. "
                     "Just vertical stripes.") as tr:
            self.play(FadeIn(pa), run_time=1.5)
            self.wait(tr.until(0.8))
            self.play(Indicate(pa.img, scale_factor=1.03), run_time=1)

        rows = [
            ["", r"\text{maritime glacier}", r"\text{coastal Antarctica}", r"\text{East Antarctic interior}"],
            [r"b_0\ [\mathrm{m\,a^{-1}}]", "1", "0.1", "0.01"],
            [r"T_s", r"0^\circ\mathrm{C}", r"-20^\circ\mathrm{C}", r"-40^\circ\mathrm{C}"],
            [r"t_0", r"100\ \mathrm{yr}", r"1000\ \mathrm{yr}", r"10\,000\ \mathrm{yr}"],
            [r"\al", "0.044", "0.082", "0.17"],
            [r"\dl", "0.038", "0.088", "0.16"],
        ]
        tbl = VGroup()
        for r in rows:
            tbl.add(VGroup(*[M(c, size=30) for c in r]))
        col_x = [-4.6, -2.0, 1.0, 4.3]
        for i, r in enumerate(tbl):
            for j, c in enumerate(r):
                c.move_to([col_x[j], 0.6 - 0.62 * i, 0])
        hline = Line([-5.6, 0.3, 0], [5.9, 0.3, 0], color=GREY_B, stroke_width=1.5)
        with self.vo("Here are three example climates from the paper: a warm, snowy maritime glacier; a site near the coast of Antarctica; "
                     "and the cold, dry interior of East Antarctica. Alpha ranges from about point zero four to point one seven: "
                     "always much less than one, so compaction dominates. And the time scale ranges from a century to ten thousand years.") as tr:
            self.play(FadeOut(pa), run_time=0.6)
            self.play(LaggedStart(*[FadeIn(r) for r in tbl[:4]], lag_ratio=0.2), Create(hline), run_time=2)
            self.wait(tr.until(0.55))
            self.play(FadeIn(tbl[4]), run_time=0.8)
            self.play(Indicate(tbl[4][1:], color=C_AL), run_time=1)
            self.wait(tr.until(0.85))
            self.play(Indicate(tbl[3][1:]), run_time=1)
        with self.vo("And delta? It's small for almost all realistic conditions, and only approaches one where it's both warm and the accumulation is low, "
                     "an unusual combination on real ice sheets. So grains rarely get anywhere near saturation within the firn, "
                     "and later on we'll often just set delta to zero.") as tr:
            self.play(FadeIn(tbl[5]), run_time=0.8)
            self.wait(tr.until(0.25))
            self.play(FadeOut(VGroup(tbl, hline)), run_time=0.6)
            pd_.move_to([0.3, -0.9, 0])
            self.play(FadeIn(pd_), run_time=1.2)
            corner = Circle(radius=0.45, color=YELLOW, stroke_width=4).move_to(pd_.c2p(268, -1.75))
            self.wait(tr.until(0.45))
            self.play(Create(corner), run_time=0.8)
        self.wait(0.4)
        self.clear_all()
