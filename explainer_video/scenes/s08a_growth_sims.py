from common import *


class GrowthSims(VScene):
    def construct(self):
        card = chapter_card("Letting the grains grow", 4)
        self.play(FadeIn(card, shift=UP * 0.3))
        self.wait(0.8)
        self.play(FadeOut(card))

        d = load("fig4")
        betas, rs2, zg = d["beta"], d["rs2"], d["zgrid"]
        nb = len(betas)
        title = Text("Evolving grain size", font="CMU Serif", font_size=38).to_edge(UP, buff=0.25)
        pphi = DepthPlot([0, 0.5], [0, 1.2], width=3.2, height=4.4, x_ticks=[0, 0.25, 0.5], z_ticks=[0, 0.4, 0.8, 1.2], x_dec=2, z_dec=1,
                         x_label=r"\ph", z_label=r"z", num_size=20, label_size=32)
        pphi.shift(np.array([-5.6, 2.0, 0]) - pphi.top.get_start())
        pw = DepthPlot([0, 20], [0, 1.2], width=3.2, height=4.4, x_ticks=[0, 10, 20], z_ticks=[0, 0.4, 0.8, 1.2], x_dec=0, z_dec=1,
                       x_label=r"\ww", num_size=20, label_size=32)
        pw.shift(np.array([-1.7, 2.0, 0]) - pw.top.get_start())
        ax = Axes(x_range=[0, 10, 2], y_range=[0, 1.2, 0.4], x_length=3.6, y_length=3.0, tips=False,
                  axis_config=dict(color=GREY_B, stroke_width=2, include_numbers=True, font_size=20,
                                   decimal_number_config=dict(num_decimal_places=1))).move_to([4.6, 0.0, 0])
        axl = VGroup(M(r"\bt", size=32).next_to(ax.x_axis, DOWN, buff=0.4), M(r"\zc", size=30).next_to(ax.y_axis, LEFT, buff=0.35).shift(UP * 0.4))
        lev = pphi.vline(0.096, color=GREY_C)
        all_lines = VGroup(*[ax.plot_line_graph(betas, d["z830"][:, i], add_vertex_dots=False, line_color=GREY_D, stroke_width=1.5)
                             for i in range(len(rs2))])

        def fam(i_rs):
            g1, g2 = VGroup(), VGroup()
            for j in range(nb):
                c1 = interpolate_color(ManimColor("#1E4F7A"), ManimColor(C_PHI), j / (nb - 1))
                c2 = interpolate_color(ManimColor("#6B5B1E"), ManimColor(C_W), j / (nb - 1))
                g1.add(pphi.curve(d["phi"][j, i_rs], zg, color=c1, width=2.2))
                g2.add(pw.curve(d["w"][j, i_rs], zg, color=c2, width=2.2))
            dots = VGroup(*[Dot(pphi.p(0.096, d["z830"][j, i_rs]), radius=0.035, color=C_Z) for j in range(nb)])
            line = ax.plot_line_graph(betas, d["z830"][:, i_rs], add_vertex_dots=False, line_color=C_Z, stroke_width=4)
            return g1, g2, dots, line

        def rs_label(i_rs, val):
            return M(r"\rs = %s" % val, size=36).move_to([4.6, 2.6, 0])

        with self.vo("Now let's switch grain growth back on, and repeat the experiment: accumulation from a tenth up to ten. "
                     "This time, we'll also vary the size of the grains at the surface, r s squared.") as tr:
            self.play(Write(title), FadeIn(pphi), FadeIn(pw), Create(lev), Create(ax), FadeIn(axl), run_time=1.5)
            self.play(Create(all_lines), run_time=1.5)

        cases = [(0, "0.001"), (2, "0.01"), (20, "0.1")]
        texts = [
            "Start with tiny surface grains, r s squared equal to point zero zero one. As the accumulation rate goes up a hundred-fold, "
            "the porosity profiles hardly move. Firn thickness barely depends on accumulation at all.",
            "With r s squared equal to point zero one, there's a clear dependence.",
            "And with large surface grains, point one, firn thickness is very sensitive to accumulation: it nearly quadruples "
            "across this range of beta.",
        ]
        prev = None
        for (i_rs, val), txt in zip(cases, texts):
            g1, g2, dots, line = fam(i_rs)
            lab = rs_label(i_rs, val)
            with self.vo(txt) as tr:
                if prev is None:
                    self.play(FadeIn(lab), run_time=0.5)
                    self.play(LaggedStart(*[AnimationGroup(Create(a), Create(b)) for a, b in zip(g1, g2)], lag_ratio=0.15), run_time=3)
                    self.play(FadeIn(dots), Create(line), run_time=1)
                else:
                    pg1, pg2, pdots, pline, plab = prev
                    self.play(Transform(pg1, g1), Transform(pg2, g2), Transform(pdots, dots), Transform(pline, line),
                              Transform(plab, lab), run_time=2)
                    g1, g2, dots, line, lab = pg1, pg2, pdots, pline, plab
            prev = (g1, g2, dots, line, lab)
        g1, g2, dots, line, lab = prev

        # Fig 5: slope vs surface grain size
        ax5 = Axes(x_range=[0, 0.1, 0.02], y_range=[0, 0.08, 0.02], x_length=6.0, y_length=4.0, tips=False,
                   axis_config=dict(color=GREY_B, stroke_width=2, include_numbers=True, font_size=22,
                                    decimal_number_config=dict(num_decimal_places=2))).move_to([0.6, -0.3, 0])
        ax5l = VGroup(M(r"\text{surface grain size } \rs", size=30).next_to(ax5.x_axis, DOWN, buff=0.45),
                      M(r"\text{sensitivity } \frac{\mathrm{d}\zc}{\mathrm{d}\bt}", size=30).rotate(PI / 2).next_to(ax5.y_axis, LEFT, buff=0.5))
        pts = VGroup(*[Dot(ax5.c2p(r, s), radius=0.06, color=C_SIG) for r, s in zip(rs2, d["slope"])])
        ln5 = ax5.plot_line_graph(rs2, d["slope"], add_vertex_dots=False, line_color=C_SIG, stroke_width=3)
        lo = M(r"0.005", size=28, color=C_SIG).next_to(pts[0], RIGHT, buff=0.15)
        hi = M(r"0.075", size=28, color=C_SIG).next_to(pts[-1], UP, buff=0.15)
        with self.vo("Do this for a whole range of surface grain sizes, and fit a straight line to each of these curves. The slope tells us how "
                     "strongly firn thickness depends on accumulation. It rises from about point zero zero five to point zero seven five: "
                     "a fifteen-fold increase. The bigger the grains at the surface, the more firn thickness cares about accumulation.") as tr:
            self.play(FadeOut(VGroup(pphi, pw, g1, g2, dots, lev, lab)), Group(ax, axl, all_lines, line).animate.scale(0.6).to_corner(UR).shift(DOWN * 0.6),
                      run_time=1)
            self.play(Create(ax5), FadeIn(ax5l), run_time=1)
            self.play(LaggedStart(*[GrowFromCenter(p) for p in pts], lag_ratio=0.08), Create(ln5), run_time=2.5)
            self.wait(tr.until(0.5))
            self.play(FadeIn(lo), run_time=0.6)
            self.play(FadeIn(hi), run_time=0.6)
        tau = TauCreature(height=1.4).move_to([-5.8, -2.6, 0])
        with self.vo("Why on earth would the size of the snow grains at the surface control how firn responds to snowfall?") as tr:
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            b, anim = tau.thinks(r"Why would grain\\size matter?", mood="confused")
            self.play(anim, run_time=1.2)
        self.wait(0.5)
        self.clear_all()
