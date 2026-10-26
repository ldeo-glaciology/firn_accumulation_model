from common import *


def beta_color(i, n):
    return interpolate_color(ManimColor("#1E4F7A"), ManimColor(C_PHI), i / (n - 1))


class ConstantGrain(VScene):
    def construct(self):
        d = load("fig3")
        betas, zg = d["beta"], d["zgrid"]
        nb = len(betas)
        title = Text("Uniform, constant grain size", font="CMU Serif", font_size=40).to_edge(UP, buff=0.25)
        rc = M(r"\rsq(z,t) = \rs", size=40).next_to(title, DOWN, buff=0.3)
        with self.vo("Let's start with the case that most firn models assume: a grain size that's uniform and never changes. "
                     "Freeze r squared at its surface value, everywhere.") as tr:
            self.play(Write(title), run_time=1.2)
            self.play(Write(rc), run_time=1.2)
        self.play(FadeOut(rc), run_time=0.4)

        pphi = DepthPlot([0, 0.5], [0, 1], width=3.2, height=4.2, x_ticks=[0, 0.25, 0.5], z_ticks=[0, 0.5, 1], x_dec=2, z_dec=1,
                         x_label=r"\ph", z_label=r"z", num_size=20, label_size=32)
        pphi.shift(np.array([-5.6, 2.1, 0]) - pphi.top.get_start())
        pw = DepthPlot([0, 20], [0, 1], width=3.2, height=4.2, x_ticks=[0, 10, 20], z_ticks=[0, 0.5, 1], x_dec=0, z_dec=1,
                       x_label=r"\ww", num_size=20, label_size=32)
        pw.shift(np.array([-1.7, 2.1, 0]) - pw.top.get_start())
        ax = Axes(x_range=[0, 10, 2], y_range=[0, 0.6, 0.2], x_length=3.6, y_length=3.0, tips=False,
                  axis_config=dict(color=GREY_B, stroke_width=2, include_numbers=True, font_size=20,
                                   decimal_number_config=dict(num_decimal_places=1))).move_to([4.6, 0.3, 0])
        axl = VGroup(M(r"\bt", size=32).next_to(ax.x_axis, DOWN, buff=0.4), M(r"\zc", size=30).next_to(ax.y_axis, LEFT, buff=0.35).shift(UP * 0.4))
        lev = pphi.vline(0.096, color=GREY_C)

        def fam(phis, ws, zz):
            g1, g2, dots1, dots2 = VGroup(), VGroup(), VGroup(), VGroup()
            for i in range(nb):
                c = beta_color(i, nb)
                g1.add(pphi.curve(phis[i], zz, color=c, width=2.5))
                g2.add(pw.curve(ws[i], zz, color=interpolate_color(ManimColor("#6B5B1E"), ManimColor(C_W), i / (nb - 1)), width=2.5))
            return g1, g2

        f_phi, f_w = fam(d["full_phi"], d["full_w"], zg)
        zdots = VGroup(*[Dot(ax.c2p(b, z), radius=0.05, color=C_Z) for b, z in zip(betas, d["z830_full"])])
        pdots = VGroup(*[Dot(pphi.p(0.096, z), radius=0.04, color=C_Z) for z in d["z830_full"]])
        full_line = ax.plot_line_graph(betas, d["z830_full"], add_vertex_dots=False, line_color=C_Z, stroke_width=3)
        with self.vo("Now turn up the accumulation rate, beta, from a half all the way up to ten, and find the steady state each time.") as tr:
            self.play(FadeIn(pphi), FadeIn(pw), Create(lev), run_time=1)
            self.play(LaggedStart(*[AnimationGroup(Create(a), Create(b)) for a, b in zip(f_phi, f_w)], lag_ratio=0.25),
                      run_time=max(tr.remaining - 0.3, 3))
        arrow = Arrow(pphi.p(0.15, 0.18), pphi.p(0.3, 0.6), buff=0, color=WHITE, stroke_width=3)
        with self.vo("More accumulation, thicker firn. But the relationship is sub-linear. Increase beta twenty-fold, and the firn thickness grows "
                     "by only about a factor of four and a half.") as tr:
            self.play(GrowArrow(arrow), LaggedStart(*[FadeIn(p) for p in pdots], lag_ratio=0.05), run_time=1.5)
            self.play(Create(ax), FadeIn(axl), run_time=1)
            self.play(LaggedStart(*[TransformFromCopy(p, z) for p, z in zip(pdots, zdots)], lag_ratio=0.05), run_time=2)
            self.play(Create(full_line), run_time=1)
        key = VGroup()

        def key_row(style, text):
            if style == "solid":
                l = Line(ORIGIN, RIGHT * 0.5, stroke_width=3, color=C_Z)
            elif style == "dash":
                l = DashedLine(ORIGIN, RIGHT * 0.5, stroke_width=3, color=C_Z, dash_length=0.08)
            else:
                l = DashedLine(ORIGIN, RIGHT * 0.5, stroke_width=3, color=C_Z, dash_length=0.02, dashed_ratio=0.4)
            return VGroup(l, M(text, size=24).next_to(l, RIGHT, buff=0.12))
        k1 = key_row("solid", r"\text{full model}")
        k2 = key_row("dash", r"\text{Eq. 24}")
        k3 = key_row("dot", r"\text{Eq. 25}")
        VGroup(k1, k2, k3).arrange(DOWN, aligned_edge=LEFT, buff=0.12).next_to(ax, UP, buff=0.2).align_to(ax, LEFT).shift(RIGHT * 0.3)

        e24 = VGroup(M(r"\dd{\ph}{z} = -\frac{z^n\ph^m(1-\ph)}{\al\,\ww\,\rs}", size=36),
                     M(r"\dd{\ww}{z} = -\frac{z^n\ph^m}{\al\,\rs}", size=36)).arrange(RIGHT, buff=0.6).move_to([-1.7, -2.85, 0])
        t24 = eq_tag(24).next_to(e24, RIGHT, buff=0.3)
        o_phi, o_w = fam(d["ode_phi"], d["ode_w"], zg)
        dash24 = DashedVMobject(ax.plot_line_graph(betas, d["z830_ode"], add_vertex_dots=False, line_color=C_Z, stroke_width=3)["line_graph"], num_dashes=30)
        with self.vo("Why sub-linear? Let's simplify. The stress turns out to be close to linear in depth, so replace sigma with minus zed. "
                     "Here's what that model gives. Same qualitative picture.") as tr:
            self.play(FadeIn(k1), FadeOut(arrow), run_time=0.5)
            self.play(Write(e24), FadeIn(t24), run_time=2)
            self.wait(tr.until(0.6))
            self.play(Transform(f_phi, o_phi), Transform(f_w, o_w), FadeOut(pdots), run_time=1.5)
            self.play(Create(dash24), FadeIn(k2), run_time=1.2)
        e25 = M(r"\dd{\ph}{z}", "=", r"-\frac{z^n\ph^m(1-\ph)}{\al\,", r"\bt", r"\,\rs}", size=44).move_to([-1.7, -2.85, 0])
        t25 = eq_tag(25).next_to(e25, RIGHT, buff=0.3)
        a_phi, _ = fam(d["an_phi"], d["ode_w"], zg)
        a_w = VGroup(*[pw.curve(np.full_like(zg, b), zg, color=interpolate_color(ManimColor("#6B5B1E"), ManimColor(C_W), i / (nb - 1)), width=2.5)
                       for i, b in enumerate(betas)])
        dot25 = DashedVMobject(ax.plot_line_graph(betas, d["z830_an"], add_vertex_dots=False, line_color=C_Z, stroke_width=3)["line_graph"], num_dashes=60, dashed_ratio=0.35)
        with self.vo("Simplify again. Ignore the way compaction slows the flow, and just set w equal to beta. Still the same picture.") as tr:
            self.play(FadeOut(e24), FadeOut(t24), FadeIn(e25), FadeIn(t25), run_time=1.2)
            self.play(Transform(f_phi, a_phi), Transform(f_w, a_w), run_time=1.5)
            self.play(Create(dot25), FadeIn(k3), run_time=1)
        with self.vo("And now look at the equation. Beta sits in the denominator. Higher accumulation means a gentler porosity gradient, "
                     "so the firn is thicker.") as tr:
            self.wait(tr.until(0.3))
            self.play(Circumscribe(e25[3], color=C_B), e25[3].animate.scale(1.3), run_time=1.5)

        # exact solution
        self.play(FadeOut(VGroup(pphi, pw, f_phi, f_w, lev)), VGroup(e25, t25).animate.move_to([-3.2, 2.2, 0]), run_time=1)
        sol = M(r"\frac{\ph}{1-\ph} = \frac{\ph_s}{1-\ph_s}\,\exp\!\left(-\frac{z^2}{2\,\al\,\bt\,\rs}\right)", size=40).move_to([-3.2, 0.7, 0])
        z830 = M(r"\zc = \sqrt{2\,\al\,\bt\,\rs\,L}\;\propto\;\sqrt{\bt}", size=44).move_to([-3.2, -0.8, 0])
        Ldef = M(r"L = \ln\!\left[\frac{\ph_s/(1-\ph_s)}{\ph_{830}/(1-\ph_{830})}\right]\approx 2.24", size=28, color=GREY_A).next_to(z830, DOWN, buff=0.3)
        note = Text("(n = m = 1; follows directly from Eq. 25)", font="CMU Serif", font_size=20, color=GREY_B).next_to(Ldef, DOWN, buff=0.2)
        with self.vo("In fact, this one we can solve exactly. With n and m equal to one, the equation separates. Phi over one minus phi decays "
                     "like a Gaussian in depth, with a width set by the square root of alpha, beta, r s squared. "
                     "So [z830](zed eight-thirty) grows like the square root of beta. And that's exactly the dotted curve.") as tr:
            self.play(Write(sol), run_time=2.5)
            self.wait(tr.until(0.55))
            self.play(Write(z830), FadeIn(Ldef), FadeIn(note), run_time=2)
            self.wait(tr.until(0.88))
            self.play(Indicate(dot25, color=WHITE), run_time=1.2)

        # intuition with a parcel
        self.play(FadeOut(VGroup(sol, Ldef, note, e25, t25)), z830.animate.scale(0.8).move_to([-3.2, -2.9, 0]), run_time=0.8)
        steps = VGroup(
            M(r"\text{depth: } z = \bt\, t", size=36),
            M(r"\text{compaction rate} \;\propto\; |\sg| \;\approx\; z = \bt\, t", size=36),
            M(r"\text{total compaction} \;\propto\; \int_0^t \bt\, t'\,\mathrm{d}t' = \frac{\bt t^2}{2} = \frac{z^2}{2\bt}", size=36),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.45).move_to([-1.6, 0.8, 0])
        colm = firn_column(width=0.9, height=4.6, top=np.array([-6.4, 2.6, 0]), n_grains=60)
        parcel = Square(0.25).set_fill(YELLOW, 1).set_stroke(BLACK, 1).move_to(colm.get_top() + DOWN * 0.2)
        with self.vo("Here's the intuition. A parcel buried at speed beta reaches depth zed after a time zed over beta. Its compaction rate is "
                     "proportional to the stress it feels, which is proportional to its depth. Add up that compaction over its journey, "
                     "and you get zed squared over two beta. To reach the density of bubble close-off, that total has to hit a fixed amount. "
                     "So zed scales like the square root of beta.") as tr:
            self.play(FadeIn(colm), FadeIn(parcel), run_time=0.8)
            self.play(parcel.animate.move_to(colm.get_bottom() + UP * 0.3), Write(steps[0]), run_time=2, rate_func=linear)
            self.wait(tr.until(0.25))
            self.play(Write(steps[1]), run_time=2)
            self.wait(tr.until(0.48))
            self.play(Write(steps[2]), run_time=2.5)
            self.wait(tr.until(0.85))
            self.play(Circumscribe(z830, color=C_Z), run_time=1.2)
        tau = TauCreature(height=1.3).move_to([6.1, -2.9, 0])
        with self.vo("Faster burial does carry low-density firn deeper. But it also means each parcel feels a large stress sooner, "
                     "which partly makes up for it. Hence: sub-linear.") as tr:
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            b, anim = tau.says(r"$\sqrt{\beta}$, not $\beta$!", mood="happy")
            self.wait(tr.until(0.6))
            self.play(anim, run_time=1)
        self.wait(0.5)
        self.clear_all()
