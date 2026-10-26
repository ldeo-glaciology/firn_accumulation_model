from common import *


class WhatIsFirn(VScene):
    def construct(self):
        d = load("fig2")
        z = d["z"][-1] * 100.0          # steady profile, metres
        phi = d["phi"][-1]
        rho = 918 * (1 - phi)
        z830 = d["z830"][-1] * 100

        def phi_of_s(s):
            return np.interp(s * 80.0, z, phi)

        col = firn_column(width=1.7, height=5.6, phi_profile=phi_of_s, top=np.array([-4.1, 2.8, 0]))
        plot = DepthPlot([400, 950], [0, 80], width=4.6, height=5.6, x_ticks=[400, 550, 700, 850],
                         z_ticks=[0, 20, 40, 60, 80], x_dec=0, z_dec=0,
                         x_label=r"\text{density } \rho \;[\mathrm{kg\,m^{-3}}]", z_label=r"\text{depth [m]}")
        plot.shift(np.array([-0.9, 2.8, 0]) - plot.top.get_start())
        prof = plot.curve(rho, z, color=WHITE, width=4)
        ice_line = plot.vline(918, color=ICE)
        ice_lab = M(r"\rho_i = 918", size=28, color=ICE).next_to(plot.p(918, 80), UP, buff=0.1).shift(LEFT * 0.75)

        with self.vo("Let's picture a column of firn. At the surface, fresh snow is mostly air. "
                     "Its density might be around half that of solid ice. As you go deeper, the grains pack together and bond, "
                     "and the density rises towards that of ice: about nine hundred and eighteen kilograms per cubic metre.") as tr:
            self.play(FadeIn(col[0]), Create(col[2]), run_time=1)
            self.play(LaggedStart(*[FadeIn(g) for g in col[1]], lag_ratio=0.005), run_time=1.5)
            self.play(FadeIn(plot), run_time=1)
            self.play(Create(prof), run_time=3, rate_func=linear)
            self.wait(tr.until(0.85))
            self.play(Create(ice_line), FadeIn(ice_lab))

        l830 = plot.vline(830, color=C_Z)
        h830 = plot.hline(z830, color=C_Z)
        dot = Dot(plot.p(830, z830), color=C_Z)
        cz_top = col.get_top()[1]
        y830 = plot.p(0, z830)[1]
        brace = BraceBetweenPoints(col.get_left() * [1, 0, 0] + [0, cz_top, 0], [col.get_left()[0], y830, 0], LEFT, color=C_Z)
        zlab = M(r"\zc", size=40).next_to(brace, LEFT)
        col_line = DashedLine([col.get_left()[0], y830, 0], [col.get_right()[0], y830, 0], color=C_Z, stroke_width=3)
        rng = np.random.default_rng(2)
        bubbles = VGroup(*[Ellipse(width=rng.uniform(0.05, 0.1), height=rng.uniform(0.03, 0.06)).set_fill(BLACK, 0.55).set_stroke(width=0)
                           .move_to([col.get_x() + rng.uniform(-0.7, 0.7), rng.uniform(col.get_bottom()[1] + 0.1, y830 - 0.15), 0])
                           for _ in range(40)])
        with self.vo("A handy measure of how thick the firn is, is the depth at which the density reaches eight hundred and thirty "
                     "kilograms per cubic metre. We'll call it [z830](zed eight-thirty). It's roughly where the air spaces between the grains "
                     "pinch off into isolated bubbles.") as tr:
            self.play(Create(l830), run_time=1)
            self.wait(tr.until(0.3))
            self.play(Create(h830), FadeIn(dot, scale=2), run_time=1)
            self.play(Create(col_line), GrowFromEdge(brace, UP), Write(zlab), run_time=1.5)
            self.wait(tr.until(0.75))
            self.play(LaggedStart(*[GrowFromCenter(b) for b in bubbles], lag_ratio=0.03), run_time=1.8)

        # why it matters: ice core + satellite
        core = RoundedRectangle(width=0.6, height=3.2, corner_radius=0.25).set_fill("#9ED3F0", 0.85).set_stroke(WHITE, 2)
        cb = VGroup(*[Circle(radius=rng.uniform(0.03, 0.06)).set_fill(WHITE, 0.9).set_stroke(width=0)
                      .move_to(core.get_center() + [rng.uniform(-0.2, 0.2), rng.uniform(-1.4, 1.4), 0]) for _ in range(25)])
        core_g = VGroup(core, cb).move_to([5.6, 0.6, 0])
        core_lab = Text("ice core", font="CMU Serif", font_size=26).next_to(core_g, DOWN)
        air = Text("ancient air", font="CMU Serif", font_size=22, color=GREY_A).next_to(core_g, LEFT, buff=0.15).shift(UP * 0.9)
        sat_body = Rectangle(width=0.5, height=0.35).set_fill(GREY_B, 1).set_stroke(WHITE, 1)
        panels = VGroup(*[Rectangle(width=0.7, height=0.25).set_fill("#2E5C8A", 1).set_stroke(WHITE, 1).next_to(sat_body, s, buff=0.05)
                          for s in [LEFT, RIGHT]])
        sat = VGroup(sat_body, panels).move_to([5.4, 3.1, 0])
        beam = DashedLine(sat.get_bottom(), [5.4, 0.4, 0], color=YELLOW, stroke_width=2)
        alt_lab = Text("satellite altimetry", font="CMU Serif", font_size=24).next_to(sat, UP, buff=0.12)
        with self.vo("Those bubbles are a big part of why firn matters. They're tiny samples of ancient atmosphere, "
                     "and to interpret them in an ice core, you need to know how deep the firn was when they were sealed in. "
                     "And when satellites measure the height of an ice sheet, to work out whether it's gaining or losing mass, "
                     "changes in the thickness of the firn can masquerade as changes in the amount of ice.") as tr:
            self.play(FadeIn(core_g, shift=UP * 0.3), FadeIn(core_lab), run_time=1.2)
            self.play(LaggedStart(*[Indicate(b, scale_factor=1.6) for b in cb[:10]], lag_ratio=0.1), FadeIn(air), run_time=1.6)
            self.wait(tr.until(0.5))
            self.play(FadeOut(core_g), FadeOut(core_lab), FadeOut(air), run_time=0.6)
            self.play(FadeIn(sat, shift=DOWN * 0.3), FadeIn(alt_lab), run_time=1)
            
            self.play(Create(beam), run_time=1)
            self.play(beam.animate.shift(DOWN * 0.4), rate_func=there_and_back, run_time=1.5)
        self.play(FadeOut(VGroup(sat, beam, alt_lab)), run_time=0.6)

        # porosity
        eq = M(r"\ph", "=", r"1 - \frac{\rho}{\rho_i}", size=48).move_to([5.4, 1.8, 0])
        eqlab = Text("porosity", font="CMU Serif", font_size=28, color=C_PHI).next_to(eq, UP)
        plot2 = DepthPlot([0, 0.5], [0, 80], width=4.6, height=5.6, x_ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5],
                          z_ticks=[0, 20, 40, 60, 80], x_dec=1, z_dec=0,
                          x_label=r"\text{porosity } \ph", z_label=r"\text{depth [m]}")
        plot2.shift(plot.top.get_start() - plot2.top.get_start())
        prof2 = plot2.curve(phi, z, color=C_PHI, width=4)
        l096 = plot2.vline(0.096, color=C_Z)
        lab096 = M(r"\ph = 0.096", size=28, color=C_Z).next_to(plot2.p(0.096, 75), RIGHT, buff=0.1)
        with self.vo("Rather than density, the paper mostly works with porosity, phi: the fraction of the volume that's air. "
                     "Phi is one minus rho over rho i. So at the surface, phi is about a half. Deep down, it's zero. "
                     "And [z830](zed eight-thirty) is just the depth where phi falls to about zero point zero nine six.") as tr:
            self.play(FadeIn(eqlab), Write(eq), run_time=1.5)
            self.wait(tr.until(0.35))
            self.play(ReplacementTransform(plot, plot2), ReplacementTransform(prof, prof2),
                      FadeOut(ice_line), FadeOut(ice_lab), FadeOut(l830),
                      dot.animate.move_to(plot2.p(0.096, z830)), run_time=2)
            self.wait(tr.until(0.8))
            self.play(Create(l096), FadeIn(lab096), Flash(dot, color=C_Z), run_time=1.2)
        self.wait(0.5)
        self.clear_all()
