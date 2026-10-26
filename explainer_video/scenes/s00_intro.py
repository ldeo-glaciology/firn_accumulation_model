from common import *


class Intro(VScene):
    def construct(self):
        surface_y = 1.3
        ground = Rectangle(width=16, height=6, stroke_width=0).set_fill("#16324A", 1)
        ground.move_to([0, surface_y - 3, 0])
        surf = Line([-8, surface_y, 0], [8, surface_y, 0], color=WHITE, stroke_width=3)
        flakes = snowfall(y_bottom=surface_y)
        stack = build_stack(15, surface_y, width=3.4)
        window = VGroup(ground, surf)

        self.add(window, stack, flakes)
        with self.vo("Every year, snow falls on the great ice sheets of Antarctica and Greenland. "
                     "And in the coldest places it never melts. It just gets buried, by next year's snow, "
                     "and the year after that, and the year after that.") as tr:
            for i in range(5):
                new = build_stack(15, surface_y, width=3.4, offset=0)
                fresh = new[0].copy().shift(UP * 0.6).set_opacity(0)
                old = stack
                self.play(*[Transform(old[k], new[k + 1]) for k in range(14)],
                          FadeOut(old[14]), new[0].animate(rate_func=linear).set_opacity(1) if False else FadeIn(new[0], shift=DOWN * 0.3),
                          run_time=1.5)
                self.remove(old[14])
                stack = VGroup(new[0], *old[:14])
                self.add(stack)
            self.fill(tr)

        # labels
        top = stack[0].get_top()[1]
        bottom = stack[-1].get_bottom()[1]
        x_r = stack.get_right()[0] + 0.3
        snow_br = BraceBetweenPoints([x_r, top, 0], [x_r, top - 0.95, 0], RIGHT)
        firn_br = BraceBetweenPoints([x_r, top - 0.95, 0], [x_r, bottom + 0.6, 0], RIGHT)
        ice_br = BraceBetweenPoints([x_r, bottom + 0.6, 0], [x_r, bottom, 0], RIGHT)
        snow_t = Text("snow", font="CMU Serif", font_size=30).next_to(snow_br, RIGHT)
        firn_t = Text("firn", font="CMU Serif", font_size=40, color=C_PHI).next_to(firn_br, RIGHT)
        ice_t = Text("ice", font="CMU Serif", font_size=30, color=ICE).next_to(ice_br, RIGHT)
        with self.vo("Squeezed under the weight of everything above it, that snow slowly turns into solid glacial ice. "
                     "The in-between stuff, snow that has survived at least one whole year but isn't yet ice, has a name. "
                     "It's called firn.") as tr:
            self.play(GrowFromCenter(snow_br), FadeIn(snow_t), run_time=1)
            self.wait(tr.until(0.35))
            self.play(GrowFromCenter(ice_br), FadeIn(ice_t), run_time=1)
            self.wait(tr.until(0.82))
            self.play(GrowFromCenter(firn_br), Write(firn_t), run_time=1.2)

        q = Text("More snowfall...", font="CMU Serif", font_size=36).to_corner(UL).shift(DOWN * 0.2)
        arrow = Arrow(q.get_bottom(), q.get_bottom() + DOWN * 1.1, color=C_B)
        thick = Text("thicker firn?", font="CMU Serif", font_size=34, color=C_PHI)
        thin = Text("thinner firn?", font="CMU Serif", font_size=34, color=C_SIG)
        opts = VGroup(thick, thin).arrange(DOWN, buff=0.25).next_to(firn_t, RIGHT, buff=0.6)
        with self.vo("Now, the climate is warming, and over large parts of Antarctica and Greenland, "
                     "a warmer atmosphere is expected to bring more snowfall. So here's a natural question. "
                     "If it snows more, does the firn layer get thicker, or thinner?") as tr:
            self.play(Write(q), GrowArrow(arrow), flakes.animate.set_opacity(1), run_time=1.5)
            self.wait(tr.until(0.7))
            self.play(FadeIn(thick, shift=RIGHT * 0.3), run_time=0.7)
            self.play(FadeIn(thin, shift=RIGHT * 0.3), run_time=0.7)

        tau = TauCreature(height=1.9).to_corner(DL).shift(UP * 0.1 + RIGHT * 0.3)
        with self.vo("Our friend here, the tau creature, thinks the answer is obvious.") as tr:
            self.play(FadeIn(tau, shift=UP * 0.5), run_time=1)
            self.play(tau.look_at(stack), run_time=0.5)
            self.play(tau.blink())
        bub, anim = tau.says(r"More snow $\Rightarrow$ thicker firn.\\Obviously!", mood="happy")
        with self.vo("Pile more snow on top, and surely you get a thicker layer of firn. Right?") as tr:
            self.play(anim, run_time=1.2)
            self.play(Indicate(thick, color=C_PHI))

        mag_center = np.array([5.5, -2.3, 0])
        disc = Circle(radius=1.35).move_to(mag_center).set_fill("#20384F", 1).set_stroke(WHITE, 3)
        grains = grain_cluster(radius=1.25, porosity=0.55).move_to(mag_center)
        src = stack[9].get_right() + LEFT * 0.3
        connector = Line(src, disc.get_left(), color=WHITE, stroke_width=2)
        dot = Dot(src, radius=0.06)
        with self.vo("Well, it turns out the answer depends, rather surprisingly, on something much smaller: "
                     "the size of the individual snow grains.") as tr:
            self.play(FadeOut(bub), FadeOut(opts), FadeOut(q), FadeOut(arrow), tau.change("surprised", disc), run_time=0.8)
            self.play(FadeIn(dot), Create(connector), FadeIn(disc), run_time=0.8)
            self.play(LaggedStart(*[GrowFromCenter(g) for g in grains], lag_ratio=0.02), run_time=2.2)
            self.play(tau.blink())

        title = Text("Grain-size evolution controls the accumulation\ndependence of modelled firn thickness",
                     font="CMU Serif", font_size=38, line_spacing=1.1)
        authors = Text("Jonathan Kingslake, Robert Skarbek, Elizabeth Case & Christine McCarthy",
                       font="CMU Serif", font_size=26, color=GREY_A)
        journal = Text("The Cryosphere, 16, 3413–3430 (2022)", font="CMU Serif", font_size=24, color=GREY_B, slant=ITALIC)
        card = VGroup(title, authors, journal).arrange(DOWN, buff=0.4).move_to(UP * 0.6)
        with self.vo("In this video we'll work through a paper by Kingslake, Skarbek, Case and McCarthy, "
                     "published in The Cryosphere in 2022. It builds a deliberately simple model of firn, "
                     "then scales it, simplifies it and simulates it, to untangle exactly how the rate of snow accumulation "
                     "controls firn thickness.") as tr:
            flakes.clear_updaters()
            self.play(*[FadeOut(m) for m in [window, stack, flakes, snow_br, snow_t, firn_br, firn_t, ice_br, ice_t,
                                              disc, grains, connector, dot, thick, thin]],
                      tau.change("plain"), run_time=1.2)
            self.play(Write(title), run_time=2.5)
            self.play(FadeIn(authors, shift=UP * 0.2), FadeIn(journal, shift=UP * 0.2), tau.look_at(title))
            self.wait(tr.until(0.8))
            self.play(tau.blink())
        self.play(FadeOut(card), FadeOut(tau))
