from common import *


def setup_column(scene=None):
    d = load("fig2")
    z, phi = d["z"][-1], d["phi"][-1]
    col = firn_column(width=1.5, height=5.4, phi_profile=lambda s: np.interp(s * 0.9, z, phi),
                      top=np.array([-4.4, 2.6, 0]))
    x = col.get_left()[0] - 0.45
    zax = Arrow([x, 2.6, 0], [x, -3.0, 0], buff=0, color=GREY_A, stroke_width=3, max_tip_length_to_length_ratio=0.04)
    z0 = M(r"z=0", size=26).next_to([x, 2.6, 0], LEFT, buff=0.1)
    zb = M(r"z=z_b", size=26).next_to([x, -2.75, 0], LEFT, buff=0.1)
    zl = M(r"z", size=32).next_to(zax.get_end(), DOWN, buff=0.1)
    return col, VGroup(zax, z0, zb, zl)


class ModelFrame(VScene):
    def construct(self):
        card = chapter_card("Building the model", 1)
        self.play(FadeIn(card, shift=UP * 0.3))
        self.wait(0.8)
        self.play(FadeOut(card))

        col, zaxis = setup_column()
        temp = M(r"T = T_s \text{ everywhere}", size=28).next_to(col, DOWN, buff=0.3).shift(RIGHT * 0.5)
        with self.vo("Here's the setup. A single, one-dimensional column of firn and ice, at one uniform temperature: the surface temperature. "
                     "Depth, zed, points downward, with zed equals zero at the surface, and zed b at the base of the domain.") as tr:
            self.play(FadeIn(col), run_time=1.2)
            self.play(GrowArrow(zaxis[0]), FadeIn(zaxis[1:]), run_time=1.2)
            self.wait(tr.until(0.45))
            self.play(Write(temp))

        # Eulerian frame: fixed grid, material flows through it
        top, bot = col.get_top()[1], col.get_bottom()[1]
        H = top - bot
        xl, xr = col.get_left()[0], col.get_right()[0]
        grid = VGroup(*[Line([xr + 0.02, top - f * H, 0], [xr + 0.32, top - f * H, 0], color=WHITE, stroke_width=3)
                        for f in np.linspace(0, 1, 11)])
        grid_lab = Text("fixed grid", font="CMU Serif", font_size=24).next_to(grid, RIGHT, buff=0.15).shift(UP * 2.2)
        stripes = VGroup(*[Line([xl, top - f * H, 0], [xr, top - f * H, 0], color="#0E2A40", stroke_width=2.5)
                           for f in np.linspace(0.02, 0.98, 14)])

        def wprof(s):
            return 0.55 * (1.0 + np.exp(-(s / 0.3) ** 2))

        def flow(g, dt):
            for ln in g:
                s = (top - ln.get_y()) / H
                ln.shift(DOWN * wprof(s) * dt)
                if ln.get_y() < bot:
                    ln.set_y(top)
        wlab = M(r"\ww", size=40).next_to(col, RIGHT, buff=1.6).shift(UP * 0.4)
        warr = Arrow(wlab.get_bottom() + UP * 0.05, wlab.get_bottom() + DOWN * 1.0, color=C_W, buff=0)
        lag = Text("Lagrangian: follow each layer", font="CMU Serif", font_size=28, color=GREY_B)
        eul = Text("Eulerian: fixed grid, firn flows through", font="CMU Serif", font_size=28, color=WHITE)
        VGroup(lag, eul).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to([2.8, 2.2, 0])
        with self.vo("Most firn models are [Lagrangian](Lagrangian): they follow each layer of snow as it's buried. "
                     "This one is Eulerian. The grid stays fixed relative to the surface, and the firn flows down through it, "
                     "with a velocity, w. That's what will let us see advection and compaction side by side, "
                     "as separate terms in the same equation.") as tr:
            self.play(FadeIn(lag), run_time=0.8)
            self.wait(tr.until(0.25))
            self.play(FadeIn(eul), lag.animate.set_opacity(0.4), Create(grid), FadeIn(grid_lab), FadeIn(stripes), run_time=1.2)
            stripes.add_updater(flow)
            self.play(Write(wlab), GrowArrow(warr), run_time=1)
            self.fill(tr)

        rows = [
            (r"\ph", "porosity"), (r"\sg", "vertical stress"), (r"\ww", "downward velocity"),
            (r"\rsq", "grain size (radius squared)"), (r"\Ag", "age"), (r"h", "domain thickness"),
        ]
        table = VGroup()
        for sym, desc in rows:
            r = VGroup(M(sym, size=40), Text(desc, font="CMU Serif", font_size=28))
            r[1].next_to(r[0], RIGHT, buff=0.4)
            table.add(r)
        table.arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([2.6, -0.2, 0])
        for r in table:
            r[1].align_to(table[3][1], LEFT)
        with self.vo("The model tracks six quantities: the porosity, phi; the vertical stress, sigma; the downward velocity, w; "
                     "the grain size, measured by the square of the mean grain radius, r squared; the age of the firn, A; "
                     "and the thickness of the whole domain, h.") as tr:
            self.play(FadeOut(lag), FadeOut(eul), FadeOut(wlab), FadeOut(warr), FadeOut(grid_lab), run_time=0.6)
            for i, r in enumerate(table):
                self.wait(tr.until([0.0, 0.2, 0.36, 0.52, 0.77, 0.88][i]))
                self.play(FadeIn(r, shift=RIGHT * 0.2), run_time=0.5)
        self.wait(0.3)
        self.play(FadeOut(table), run_time=0.6)

        # material derivative
        md = M(r"\DDt{\ph}", "=", r"\pd{\ph}{t}", "+", r"\ww\pd{\ph}{z}", size=50).move_to([2.3, 2.3, 0])
        b1 = Brace(md[2], DOWN, buff=0.1)
        b1t = Text("change at a\nfixed depth", font="CMU Serif", font_size=22, line_spacing=0.9).next_to(b1, DOWN, buff=0.1).shift(LEFT * 0.45)
        b2 = Brace(md[4], DOWN, buff=0.1)
        b2t = Text("carried down\ninto new firn", font="CMU Serif", font_size=22, line_spacing=0.9).next_to(b2, DOWN, buff=0.1).shift(RIGHT * 0.35)
        with self.vo("Now, one key idea. Imagine riding along with a parcel of firn as it's buried. The rate at which its porosity changes, "
                     "which we write capital D phi by D t, has two parts. Partial phi by partial t is the change you'd see sitting still "
                     "at one depth. And w times partial phi by partial zed accounts for the parcel being carried down into firn with a "
                     "different porosity.") as tr:
            self.play(Write(md[0]), run_time=1)
            self.wait(tr.until(0.35))
            self.play(Write(md[1:]), run_time=1.5)
            self.wait(tr.until(0.55))
            self.play(GrowFromCenter(b1), FadeIn(b1t), md[2].animate.set_color(YELLOW), run_time=0.8)
            self.wait(tr.until(0.75))
            self.play(GrowFromCenter(b2), FadeIn(b2t), md[2].animate.set_color(WHITE), run_time=0.8)

        d = load("fig2")
        z, phi = d["z"][-1], d["phi"][-1]
        plot = DepthPlot([0, 0.5], [0, 1], width=3.2, height=3.0, x_ticks=[0, 0.25, 0.5], z_ticks=[0, 0.5, 1],
                         x_dec=2, z_dec=1, x_label=r"\ph", z_label=r"z", num_size=20)
        plot.shift(np.array([1.6, -0.6, 0]) - plot.top.get_start())
        prof = plot.curve(phi, z, color=C_PHI, width=4)
        parcel = Square(0.3).set_fill(YELLOW, 0.9).set_stroke(BLACK, 1).move_to([col.get_x(), top - 0.15, 0])
        pdot = Dot(plot.p(phi[0], 0), color=YELLOW)
        tracker = ValueTracker(0.0)
        parcel.add_updater(lambda m: m.move_to([col.get_x(), top - 0.15 - tracker.get_value() * (H - 0.3), 0]))
        pdot.add_updater(lambda m: m.move_to(plot.p(np.interp(tracker.get_value(), z, phi), tracker.get_value())))
        steady = Text("steady profile", font="CMU Serif", font_size=24, color=GREY_B).next_to(plot, RIGHT, buff=0.2).shift(UP * 0.5)
        with self.vo("So even if the profile isn't changing at all, a parcel moving down through it still sees its porosity fall.") as tr:
            self.play(FadeIn(plot), Create(prof), FadeIn(steady), FadeIn(parcel), FadeIn(pdot), run_time=1.2)
            self.play(tracker.animate.set_value(0.95), run_time=max(tr.remaining - 0.2, 1.5), rate_func=lambda t: 1 - (1 - t) ** 1.6)
        parcel.clear_updaters()
        pdot.clear_updaters()
        self.play(FadeOut(VGroup(plot, prof, steady, parcel, pdot, b1, b1t, b2, b2t)),
                  md.animate.scale(0.6).to_corner(UR).shift(LEFT * 0.2), run_time=0.8)

        # constitutive law
        law = M(r"\DDt{\ph}", "=", r"k_c\,", r"\mathrm{sign}(\sg)|\sg|^n", r"\,\ph^m", r"(1-\ph)", r"\,f_T(T_s)", r"\,f_r(\rr)",
                size=46).move_to([2.0, 1.4, 0])
        tag2 = eq_tag(2).next_to(law, DOWN, buff=0.15).align_to(law, RIGHT)
        ref = Text("viscous compaction (Arthern et al., 2010)", font="CMU Serif", font_size=24, color=GREY_B).next_to(law, UP, buff=0.35)
        with self.vo("So how fast does a parcel compact? The model uses a viscous law, following Arthern and colleagues, in 2010.") as tr:
            self.play(FadeIn(ref), Write(law), FadeIn(tag2), run_time=2.5)
        notes = [
            (3, "stress, to a power n\n(negative: compressive)", C_SIG),
            (4, "more pore space,\neasier to squash", C_PHI),
            (5, "solid fraction", C_PHI),
            (6, "temperature", WHITE),
            (7, "grain size", C_R),
        ]
        note_mobs = []
        with self.vo("The compaction rate is proportional to the stress, sigma, raised to a power n. The stress is negative because it's compressive, "
                     "hence the sign function, so porosity goes down. Then there's phi to the m: firn with more open space is easier to squash. "
                     "A factor of one minus phi. And two functions, of temperature and of grain size.") as tr:
            fr = [0.0, 0.43, 0.68, 0.8, 0.9]
            for (i, txt, c), f in zip(notes, fr):
                self.wait(tr.until(f))
                box = SurroundingRectangle(law[i], color=c, buff=0.06)
                t = Text(txt, font="CMU Serif", font_size=20, color=c, line_spacing=0.85).next_to(box, DOWN, buff=0.15)
                if note_mobs:
                    self.play(FadeOut(note_mobs[-1]), Create(box), FadeIn(t), run_time=0.6)
                else:
                    self.play(Create(box), FadeIn(t), run_time=0.6)
                note_mobs.append(VGroup(box, t))
        self.play(FadeOut(note_mobs[-1]), run_time=0.4)

        fT = M(r"f_T(T_s) = \exp\!\left[-\frac{E_c}{R\,T_s}\right]", size=40).move_to([2.0, -0.4, 0])
        warm = Text("warmer  ⟹  faster compaction", font="CMU Serif", font_size=24, color=GREY_A).next_to(fT, DOWN, buff=0.2)
        with self.vo("The temperature dependence is an Arrhenius law: e to the minus E c over R T s. Warmer firn compacts faster.") as tr:
            self.play(Indicate(law[6]), TransformFromCopy(law[6], fT), run_time=1.5)
            self.play(FadeIn(warm), run_time=0.8)
        self.play(FadeOut(fT), FadeOut(warm), run_time=0.5)

        fr_eq = M(r"f_r(\rr) = \frac{1}{\rsq}", size=46).move_to([0.6, -0.5, 0])
        small = Circle(0.35).set_fill(SNOW, 0.85).set_stroke(GREY_B, 2).move_to([3.4, -0.6, 0])
        big = Circle(0.8).set_fill(SNOW, 0.85).set_stroke(GREY_B, 2).move_to([5.3, -0.6, 0])
        pa = Arrow(small.get_left(), small.get_right(), buff=0.05, color=C_R, stroke_width=4)
        pb = Arrow(big.get_left(), big.get_right(), buff=0.05, color=C_R, stroke_width=4)
        fast = Text("fast", font="CMU Serif", font_size=24, color=C_R).next_to(small, DOWN)
        slow = Text("slow", font="CMU Serif", font_size=24, color=C_R).next_to(big, DOWN)
        nh = Text("Nabarro–Herring creep:\ndiffusion through grains", font="CMU Serif", font_size=22, color=GREY_B, line_spacing=0.9).next_to(VGroup(small, big), DOWN, buff=0.7)
        tau = TauCreature(height=1.5).to_corner(DL).shift(RIGHT * 2.3)
        with self.vo("And the grain-size dependence is one over r squared. It comes from Nabarro–Herring creep, in which ice deforms by molecules "
                     "diffusing through the crystals. Bigger grains mean longer diffusion paths, so bigger grains compact more slowly. "
                     "Remember this one over r squared. It's going to be the hero of our story.") as tr:
            self.play(Indicate(law[7], color=C_R), TransformFromCopy(law[7], fr_eq), run_time=1.5)
            self.play(FadeIn(small), FadeIn(big), FadeIn(nh), run_time=0.8)
            self.play(GrowArrow(pa), GrowArrow(pb), FadeIn(fast), FadeIn(slow), run_time=1)
            self.wait(tr.until(0.78))
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            self.play(Circumscribe(fr_eq, color=C_R), tau.change("happy", fr_eq), run_time=1.2)
        self.play(FadeOut(VGroup(small, big, pa, pb, fast, slow, nh)), FadeOut(tau), run_time=0.6)

        eq4 = M(r"\pd{\ph}{t}", "=", r"k_c\,\frac{\mathrm{sign}(\sg)\,|\sg|^n\,\ph^m(1-\ph)}{\rsq}\,e^{-E_c/RT_s}", r"\;-\;\ww\pd{\ph}{z}",
                size=38).move_to([1.8, -0.9, 0])
        bc4 = M(r"\ph(0) = \ph_s", size=34).next_to(eq4, DOWN, buff=0.3).align_to(eq4, LEFT)
        tag4 = eq_tag(4).next_to(eq4, UP, buff=0.15).align_to(eq4, RIGHT)
        bc_ = Brace(eq4[2], DOWN, buff=0.1, color=C_SIG)
        ba_ = Brace(eq4[3], DOWN, buff=0.1, color=C_W)
        bct = Text("compaction", font="CMU Serif", font_size=26, color=C_SIG).next_to(bc_, DOWN, buff=0.1)
        bat = Text("advection", font="CMU Serif", font_size=26, color=C_W).next_to(ba_, DOWN, buff=0.1)
        with self.vo("For most of the paper, n and m are both one: a linear viscous law. Put the pieces together with the material derivative, "
                     "and we get an equation for how the porosity evolves at each depth.") as tr:
            self.play(FadeOut(fr_eq), run_time=0.5)
            self.play(TransformMatchingShapes(VGroup(law.copy(), md.copy()), eq4), run_time=2.5)
            self.play(FadeIn(bc4), FadeIn(tag4))
        with self.vo("Two terms. Compaction, and advection. And the competition between those two is the story of this whole video.") as tr:
            self.play(GrowFromCenter(bc_), FadeIn(bct), run_time=0.8)
            self.play(GrowFromCenter(ba_), FadeIn(bat), run_time=0.8)
        self.wait(0.5)
        stripes.clear_updaters()
        self.clear_all()
