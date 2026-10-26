from common import *


class Conclusion(VScene):
    def construct(self):
        title = Text("So what?", font="CMU Serif", font_size=42).to_edge(UP, buff=0.25)
        hl = M(r"\DDt{\rho} = C\,", r"\bb^{\,\alpha}", r"(\rho_i-\rho)", size=44).move_to([-3.3, 1.6, 0])
        hl_l = Text("Herron & Langway (1980)", font="CMU Serif", font_size=22, color=GREY_B).next_to(hl, DOWN, buff=0.2)
        art = Text("derived from a viscous model by\nArthern et al. (2010), assuming:", font="CMU Serif", font_size=24, line_spacing=0.9).move_to([2.9, 2.0, 0])
        assume = VGroup(
            Text("• a steady state", font="CMU Serif", font_size=26),
            Text("• negligibly small surface grains", font="CMU Serif", font_size=26, color=C_R),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).next_to(art, DOWN, buff=0.3).align_to(art, LEFT)
        e = M(r"\rs = 0 \;\Rightarrow\; \dd{\ph}{z} = -\frac{z^{n-1}\ph^m(1-\ph)}{\al}", size=40).move_to([0, -1.2, 0])
        e_l = Text("no dependence on accumulation", font="CMU Serif", font_size=26, color=YELLOW).next_to(e, DOWN, buff=0.3)
        with self.vo("So, let's step back. Remember the Herron and Langway law, with the compaction rate proportional to the accumulation rate? "
                     "Arthern and colleagues showed that a law like that can be derived from a viscous model like this one, by assuming a steady state, "
                     "and negligibly small grains at the surface.") as tr:
            self.play(Write(title), run_time=0.8)
            self.play(Write(hl), FadeIn(hl_l), run_time=1.5)
            self.wait(tr.until(0.45))
            self.play(FadeIn(art), run_time=1)
            self.play(LaggedStart(*[FadeIn(a, shift=RIGHT * 0.2) for a in assume], lag_ratio=0.5), run_time=1.5)
        with self.vo("And that is precisely the case in which beta cancels. The accumulation hidden in the compaction law, and the accumulation driving advection, "
                     "balance each other exactly, and, in that low-density regime, accumulation has no effect on the steady-state density.") as tr:
            self.play(Indicate(assume[1], color=C_R), run_time=1)
            self.play(Write(e), run_time=1.8)
            self.play(FadeIn(e_l), run_time=0.8)
        tau = TauCreature(height=1.4).move_to([6.0, -2.8, 0])
        with self.vo("In other words, models built on that kind of law are quietly making an assumption about grain size, "
                     "and it's one that hasn't been examined in detail.") as tr:
            self.play(Circumscribe(hl[1], color=C_B), FadeIn(tau, shift=UP * 0.3), run_time=1.2)
            self.play(tau.change("surprised", hl), run_time=0.6)
            self.play(tau.blink())
        self.play(*[FadeOut(m) for m in [hl, hl_l, art, assume, e, e_l, title]], tau.change("plain"), run_time=0.8)

        c1 = VGroup(Text("1", font="CMU Serif", font_size=60, color=C_Z),
                    Text("The link between firn thickness and accumulation\ngets stronger as the surface grain size increases.",
                         font="CMU Serif", font_size=30, line_spacing=1.0)).arrange(RIGHT, buff=0.5)
        c2 = VGroup(Text("2", font="CMU Serif", font_size=60, color=C_Z),
                    Text("Assumptions about grain size are hidden inside\nwidely used firn compaction models.",
                         font="CMU Serif", font_size=30, line_spacing=1.0)).arrange(RIGHT, buff=0.5)
        concl = VGroup(c1, c2).arrange(DOWN, aligned_edge=LEFT, buff=0.8).move_to([-0.4, 0.6, 0])
        hd = Text("Conclusions", font="CMU Serif", font_size=42).to_edge(UP, buff=0.3)
        with self.vo("So, the two main conclusions. One: in this model, the strength of the link between firn thickness and accumulation rate "
                     "increases with the surface grain size. Two: assumptions about grain size are hiding inside widely used firn compaction models.") as tr:
            self.play(Write(hd), run_time=0.8)
            self.play(FadeIn(c1, shift=UP * 0.3), run_time=1.2)
            self.wait(tr.until(0.55))
            self.play(FadeIn(c2, shift=UP * 0.3), run_time=1.2)
        self.play(FadeOut(concl), FadeOut(hd), run_time=0.6)

        why = Text("Why it matters", font="CMU Serif", font_size=40).to_edge(UP, buff=0.3)
        rng = np.random.default_rng(9)
        sites = VGroup()
        for k, (x, r) in enumerate([(-4.5, 0.06), (-1.5, 0.12), (1.5, 0.09), (4.5, 0.16)]):
            disc = Circle(radius=0.9).move_to([x, 0.9, 0]).set_fill("#20384F", 1).set_stroke(GREY_B, 2)
            grains = grain_cluster(radius=0.85, rmin=r * 0.7, rmax=r * 1.3, seed=k + 20, porosity=0.55).move_to(disc)
            sites.add(VGroup(disc, grains))
        sens = Text("different surface grains  →  different sensitivity to snowfall", font="CMU Serif", font_size=28, color=C_R).move_to([0, -0.8, 0])
        uses = VGroup(Text("ice-core gas ages & past climate", font="CMU Serif", font_size=28),
                      Text("ice-sheet mass change from satellite altimetry", font="CMU Serif", font_size=28)).arrange(DOWN, buff=0.25).move_to([0, -2.2, 0])
        with self.vo("Why does this matter? Surface grain size isn't constant. It changes with weather and climate, from place to place and from year to year. "
                     "If this mechanism operates in real firn, then the sensitivity of firn thickness to snowfall should vary too. "
                     "That matters for reconstructing past climate from ice cores, and for turning satellite measurements of surface height "
                     "into ice-sheet mass change.") as tr:
            self.play(Write(why), run_time=0.8)
            self.play(LaggedStart(*[FadeIn(s) for s in sites], lag_ratio=0.3), run_time=2)
            self.wait(tr.until(0.42))
            self.play(FadeIn(sens), run_time=1)
            self.wait(tr.until(0.7))
            self.play(LaggedStart(*[FadeIn(u, shift=UP * 0.2) for u in uses], lag_ratio=0.4), run_time=1.5)
        self.play(FadeOut(VGroup(why, sites, sens, uses)), run_time=0.6)

        lim_t = Text("Deliberately simple", font="CMU Serif", font_size=34, color=GREY_A).move_to([-3.6, 2.7, 0])
        lims = VGroup(*[Text("• " + s, font="CMU Serif", font_size=25) for s in [
            "uniform temperature", "dry snow, no meltwater", "purely viscous compaction",
            "one simple grain-size law (1/r²)", "normal grain growth", "zero compaction gradient at the surface"]]
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22).next_to(lim_t, DOWN, buff=0.35).align_to(lim_t, LEFT).shift(LEFT * 0.6)
        fut_t = Text("Next steps", font="CMU Serif", font_size=34, color=GREY_A).move_to([2.6, 2.7, 0])
        futs = VGroup(*[Text("• " + s, font="CMU Serif", font_size=24) for s in [
            "radar measurements of compaction", "grain sizes from ice cores", "lab experiments on artificial firn",
            "firn-thickness compilations"]]
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22).next_to(fut_t, DOWN, buff=0.35).align_to(fut_t, LEFT).shift(LEFT * 0.6)
        with self.vo("Of course, the model is deliberately simple. The temperature is uniform. The snow is dry, with no meltwater. Compaction is purely viscous, "
                     "with one simple dependence on grain size; grains follow normal grain growth; and the compaction gradient vanishes at the surface. "
                     "The authors point towards testing these ideas with radar measurements of firn compaction, grain sizes measured in ice cores, "
                     "and laboratory experiments on artificial firn with controlled grain sizes.") as tr:
            self.play(FadeIn(lim_t), LaggedStart(*[FadeIn(l, shift=RIGHT * 0.2) for l in lims], lag_ratio=0.35), run_time=4)
            self.wait(tr.until(0.6))
            self.play(FadeIn(fut_t), LaggedStart(*[FadeIn(f, shift=RIGHT * 0.2) for f in futs], lag_ratio=0.35), run_time=3)
        self.play(FadeOut(VGroup(lim_t, lims, fut_t, futs)), run_time=0.6)

        self.play(tau.animate.move_to([-2.0, -1.4, 0]).scale(1.5), run_time=1)
        b, anim = tau.says(r"So more snow \emph{does}\\mean thicker firn...", mood="happy")
        with self.vo("So, more snow does mean thicker firn...") as tr:
            self.play(anim, run_time=1.2)
        with self.vo("...but how much thicker depends on the size of the grains. Sometimes, the smallest things decide the biggest questions.") as tr:
            self.play(FadeOut(b), tau.change("pondering"), run_time=0.8)
            g = grain_cluster(radius=1.1, porosity=0.5, seed=4).move_to([2.6, -0.4, 0])
            self.play(LaggedStart(*[GrowFromCenter(x) for x in g], lag_ratio=0.02), tau.look_at(g), run_time=2)
            self.play(tau.blink())
            self.wait(tr.until(0.9))
            self.play(tau.change("happy"), run_time=0.6)
        self.play(FadeOut(g), FadeOut(tau), run_time=0.8)

        ref = VGroup(
            Text("Kingslake, Skarbek, Case & McCarthy (2022)", font="CMU Serif", font_size=32),
            Text("Grain-size evolution controls the accumulation dependence of modelled firn thickness",
                 font="CMU Serif", font_size=24, slant=ITALIC),
            Text("The Cryosphere, 16, 3413–3430   ·   doi:10.5194/tc-16-3413-2022", font="CMU Serif", font_size=24, color=GREY_A),
            Text("Model code: github.com/ldeo-glaciology/firn_accumulation_model   ·   doi:10.5281/zenodo.5213880",
                 font="CMU Serif", font_size=20, color=GREY_B),
        ).arrange(DOWN, buff=0.3).move_to([0, 0.5, 0])
        credits = Text("Animated with Manim  ·  simulations re-run in Python from the paper's model  ·  voice: Kokoro TTS",
                       font="CMU Serif", font_size=18, color=GREY_C).to_edge(DOWN, buff=0.4)
        with self.vo("The paper, and all the code needed to run the model and reproduce its figures, are openly available. Thanks for watching.") as tr:
            self.play(LaggedStart(*[FadeIn(r, shift=UP * 0.2) for r in ref], lag_ratio=0.3), FadeIn(credits), run_time=2.5)
        self.wait(2.0)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)
