from common import *


class Numerics(VScene):
    def construct(self):
        card = chapter_card("Simulations", 3)
        self.play(FadeIn(card, shift=UP * 0.3))
        self.wait(0.8)
        self.play(FadeOut(card))

        # ---- method of lines
        top, bot, x = 2.8, -2.8, -4.6
        col = Rectangle(width=1.0, height=top - bot).move_to([x, (top + bot) / 2, 0]).set_fill(DEEP_ICE, 0.35).set_stroke(WHITE, 2)
        nodes = VGroup(*[Dot([x, y, 0], radius=0.045, color=WHITE) for y in np.linspace(top, bot, 26)])
        lab = M(r"N = 101 \text{ nodes}", size=30).next_to(col, DOWN, buff=0.2)
        zoom = VGroup()
        ys = [1.6, 0.4, -0.8]
        names = [r"\ph_{i-1}", r"\ph_{i}", r"\ph_{i+1}"]
        for yy, nm in zip(ys, names):
            zoom.add(VGroup(Dot([-0.6, yy, 0], radius=0.09, color=C_PHI), M(nm, size=34).next_to([-0.6, yy, 0], LEFT, buff=0.2)))
        spine = Line([-0.6, 2.3, 0], [-0.6, -1.5, 0], color=GREY_B)
        dz_br = BraceBetweenPoints([-0.4, ys[0], 0], [-0.4, ys[1], 0], RIGHT)
        dz_lab = M(r"\Delta z", size=30).next_to(dz_br, RIGHT, buff=0.1)
        flow = Arrow([-2.5, 2.2, 0], [-2.5, -1.2, 0], color=C_W, buff=0)
        flow_l = M(r"\ww>0", size=28).next_to(flow, LEFT, buff=0.1)
        up = M(r"\pd{\ph}{z}\Big|_i \approx \frac{\ph_i - \ph_{i-1}}{\Delta z}", size=40).move_to([3.6, 1.6, 0])
        up_l = Text("upwind difference", font="CMU Serif", font_size=26, color=GREY_A).next_to(up, UP, buff=0.25)
        odes = M(r"\dd{\ph_i}{t} = F_i\big(\ph_1,\dots,\ph_N,\ \rsq_1,\dots,\ h\big)", size=34).move_to([3.7, -0.3, 0])
        odes_l = M(r"3N+1 \text{ coupled ODEs in time}", size=30, color=GREY_A).next_to(odes, DOWN, buff=0.3)
        conn = DashedLine(nodes[10].get_center(), zoom[1][0].get_center(), color=GREY_C)
        with self.vo("To solve the full model, the authors use the method of lines. Chop the column into a hundred and one grid points. "
                     "Approximate each depth derivative with an upwind finite difference: firn only ever flows downward, so information "
                     "arrives from above, and we difference using the point above. What's left is a big system of ordinary differential "
                     "equations in time, one for each variable at each grid point, which a standard stiff solver can march forward.") as tr:
            self.play(Create(col), LaggedStart(*[GrowFromCenter(n) for n in nodes], lag_ratio=0.03), FadeIn(lab), run_time=1.8)
            self.play(Create(conn), Create(spine), LaggedStart(*[FadeIn(z) for z in zoom], lag_ratio=0.2), run_time=1.2)
            self.play(GrowFromCenter(dz_br), FadeIn(dz_lab), run_time=0.6)
            self.wait(tr.until(0.3))
            self.play(GrowArrow(flow), FadeIn(flow_l), run_time=0.8)
            self.play(FadeIn(up_l), Write(up), run_time=1.5)
            self.play(Indicate(zoom[0]), Indicate(zoom[1]), run_time=1)
            self.wait(tr.until(0.68))
            self.play(Write(odes), run_time=1.5)
            self.play(FadeIn(odes_l), run_time=0.6)
        self.play(FadeOut(VGroup(spine, zoom, dz_br, dz_lab, flow, flow_l, up, up_l, odes, odes_l, conn)), run_time=0.6)

        b1 = M(r"z = h(t)\,\hat z", size=40)
        b3 = M(r"\pd{}{z} = \frac{1}{h}\pd{}{\hat z}", size=40)
        b6 = M(r"\pd{}{t} \;\to\; \pd{}{t} - \frac{\hat z\,\dot h}{h}\pd{}{\hat z}", size=40)
        bs = VGroup(b1, b3, b6).arrange(DOWN, buff=0.4).move_to([2.2, 0.4, 0])
        stretch_l = Text("stretched coordinate (Appendix B)", font="CMU Serif", font_size=26, color=GREY_A).next_to(bs, UP, buff=0.4)
        with self.vo("There's one twist. The thickness of the domain, h, changes with time, so they stretch the coordinate: zed hat equals zed over h, "
                     "which always runs from zero to one. The chain rule then adds extra terms, proportional to h dot, "
                     "which account for the grid itself stretching or shrinking.") as tr:
            self.play(FadeIn(stretch_l), Write(b1), run_time=1.2)
            self.play(VGroup(col, nodes).animate.stretch(0.8, 1, about_edge=UP), rate_func=there_and_back, run_time=2)
            self.play(Write(b3), run_time=1)
            self.wait(tr.until(0.6))
            self.play(Write(b6), run_time=1.5)
        self.play(FadeOut(VGroup(bs, stretch_l, col, nodes, lab)), run_time=0.6)

        # ---- the time-dependent simulation
        d = load("fig2")
        t = d["t"]
        panels = []
        specs = [("phi", [0, 0.5], [0, 0.25, 0.5], r"\ph", C_PHI, 2),
                 ("w", [0, 2.2], [0, 1, 2], r"\ww", C_W, 0),
                 ("r2", [0, 1.0], [0, 0.5, 1.0], r"\rsq", C_R, 1),
                 ("A", [0, 1.0], [0, 0.5, 1.0], r"\Ag", C_A, 1)]
        for i, (k, xr, xt, lab_, c, dec) in enumerate(specs):
            p = DepthPlot(xr, [0, 1.0], width=2.6, height=4.6, x_ticks=xt, z_ticks=[0, 0.5, 1.0], x_dec=dec, z_dec=1,
                          x_label=lab_, z_label=(r"z" if i == 0 else None), num_size=20, label_size=34)
            p.shift(np.array([-5.6 + 3.25 * i, 1.9, 0]) - p.top.get_start())
            panels.append(p)
        idx = ValueTracker(0)

        def cur(k, c, p):
            def f():
                i = int(round(idx.get_value()))
                return p.curve(d[k][i], d["z"][i], color=c, width=4)
            return f
        curves = [always_redraw(cur(k, c, p)) for (k, *_r, c, _d), p in zip(specs, panels)]
        clock = always_redraw(lambda: M(r"t = %.2f \;\;(\approx %d\ \text{years})" % (t[int(round(idx.get_value()))],
                                                                                   round(t[int(round(idx.get_value()))] * 1000, -1)),
                                        size=34).to_edge(DOWN, buff=0.35))
        hline = always_redraw(lambda: DashedLine(panels[0].p(0, d["H"][int(round(idx.get_value()))]),
                                                 panels[-1].p(1.0, d["H"][int(round(idx.get_value()))]),
                                                 color=GREY_C, stroke_width=1.5))
        hlab = always_redraw(lambda: M(r"h", size=26, color=GREY_B).next_to(panels[-1].p(1.0, d["H"][int(round(idx.get_value()))]), RIGHT, buff=0.1))
        setup = M(r"b_0 = 0.1\ \mathrm{m\,a^{-1}},\;\; T_s = -20^\circ\mathrm{C},\;\; \bt = 1,\;\; r_s = 0.5\ \mathrm{mm}", size=30).to_edge(UP, buff=0.15)
        n = len(t)
        with self.vo("Let's run it. The intermediate climate: ten centimetres of ice-equivalent accumulation per year at minus twenty degrees, "
                     "beta equal to one, and surface grains with a radius of half a millimetre. We start from a straight-line porosity profile, "
                     "and let it go. After a dimensionless time of around point eight, which here is about eight hundred years, "
                     "everything settles into a steady state.") as tr:
            self.play(FadeIn(setup), *[FadeIn(p) for p in panels], run_time=1.2)
            self.add(*curves, clock, hline, hlab)
            self.wait(tr.until(0.5))
            self.play(idx.animate.set_value(n - 1), run_time=max(tr.until(0.97), 6), rate_func=linear)
        for m in curves + [clock, hline, hlab]:
            m.clear_updaters()

        P = panels
        a1 = Arrow(P[0].p(0.5, 0.0) + RIGHT * 0.5 + DOWN * 0.3, P[0].p(0.5, 0.0) + DOWN * 0.05, buff=0.05, color=WHITE, stroke_width=3)
        z830 = d["z830"][-1]
        lz = P[0].hline(z830, color=C_Z)
        lzl = M(r"\zc \approx 36\ \mathrm{m}", size=26).next_to(P[0].p(0.25, z830), DOWN, buff=0.1)
        w2 = M(r"\ww(0) = \tfrac{\bt}{1-\ph_s} = 2", size=24).move_to(P[1].p(0.75, 0.2))
        w1 = M(r"\ww \to \bt = 1", size=26).next_to(P[1].p(1.0, 0.85), RIGHT, buff=0.1)
        with self.vo("The porosity falls from a half at the surface to almost nothing about forty metres down, with [z830](zed eight-thirty) at about "
                     "thirty-six metres. The velocity drops from two at the surface to one at depth: exactly beta over one minus phi. "
                     "And the age and grain size both increase roughly linearly with depth.") as tr:
            self.play(Create(lz), FadeIn(lzl), run_time=1)
            self.wait(tr.until(0.42))
            self.play(FadeIn(w2), run_time=0.8)
            self.play(FadeIn(w1), run_time=0.8)
            self.wait(tr.until(0.78))
            self.play(Indicate(curves[2], scale_factor=1.0), Indicate(curves[3], scale_factor=1.0), run_time=1.5)

        circ = Circle(radius=0.45, color=YELLOW, stroke_width=4).move_to(P[0].p(0.49, 0.03))
        tau = TauCreature(height=1.3).move_to([5.9, -2.9, 0])
        with self.vo("Notice one thing. Right at the surface, the porosity profile is flat. That's because the stress there is zero, "
                     "so nothing is compacting yet. Real firn doesn't behave like this. It's a limitation shared by any model in which "
                     "compaction depends only on the overburden stress, and the paper flags it clearly.") as tr:
            self.play(Create(circ), run_time=1)
            self.wait(tr.until(0.5))
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            self.play(tau.change("confused", circ), run_time=0.8)
            self.play(tau.blink())
        self.wait(0.4)
        self.clear_all()
