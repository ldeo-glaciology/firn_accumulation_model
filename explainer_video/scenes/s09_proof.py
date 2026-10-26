from common import *


class FullModelProof(VScene):
    def construct(self):
        title = Text("Is the cancellation real?", font="CMU Serif", font_size=40).to_edge(UP, buff=0.25)
        sub = Text("Appendix C: the full steady-state model", font="CMU Serif", font_size=26, color=GREY_A).next_to(title, DOWN, buff=0.15)
        tau = TauCreature(height=1.4).move_to([-5.9, -2.8, 0])
        with self.vo("You might worry that this perfect cancellation is just an artefact of all our simplifications. It isn't. "
                     "The paper's appendix shows that it holds in the full steady-state model, with no assumption about the velocity, "
                     "and for any stress exponent.") as tr:
            self.play(Write(title), FadeIn(tau, shift=UP * 0.3), run_time=1.2)
            b, anim = tau.thinks(r"Artefact?", mood="confused")
            self.play(anim, run_time=1)
            self.wait(tr.until(0.6))
            self.play(FadeIn(sub), FadeOut(b), tau.change("pondering"), run_time=0.8)
        self.play(FadeOut(tau), run_time=0.5)
        tau.move_to([6.2, 2.6, 0])

        c1 = M(r"\ww\,\dd{\ph}{z} = -\frac{|\sg|^n\ph^m(1-\ph)}{\al\,\rsq}", size=40).move_to([-2.8, 1.6, 0])
        t1 = eq_tag("C1").next_to(c1, RIGHT, buff=0.3)
        with self.vo("Start from the steady porosity equation.") as tr:
            self.play(Write(c1), FadeIn(t1), run_time=1.5)

        c2 = M(r"\pd{\ph}{t} = \pd{}{z}\big[(1-\ph)\,\ww\big]", size=34)
        c3 = M(r"(1-\ph)\,\ww = C_1 = \bt", size=34)
        c4 = M(r"\ww = \frac{\bt}{1-\ph}", size=40)
        col = VGroup(c2, c3, c4).arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([3.6, 1.2, 0])
        with self.vo("Combining the porosity and velocity equations, in a steady state the ice flux, one minus phi times w, has to be the same at every depth. "
                     "The surface boundary condition fixes that constant: it's beta. So w equals beta over one minus phi.") as tr:
            self.play(Write(c2), run_time=1.3)
            self.wait(tr.until(0.35))
            self.play(Write(c3), run_time=1.3)
            self.wait(tr.until(0.75))
            self.play(Write(c4), run_time=1.2)
        c5 = M(r"|\sg| = \int_0^z (1-\ph)\,\mathrm{d}z'", size=40).move_to([3.6, -1.0, 0])
        with self.vo("The stress is just the integral of one minus phi.") as tr:
            self.play(Write(c5), run_time=1.5)
        c6 = M(r"\dd{\rsq}{z} = \frac{1}{\ww} = \frac{1-\ph}{\bt}", size=36)
        c8 = M(r"\rsq = \rs + \frac{1}{\bt}\int_0^z (1-\ph)\,\mathrm{d}z'", size=38)
        g = VGroup(c6, c8).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to([-2.8, -0.6, 0])
        with self.vo("And for the grains, with delta zero, d r squared by d zed is one over w, which is one minus phi over beta. Integrate, "
                     "and r squared is r s squared, plus one over beta times the very same integral of one minus phi.") as tr:
            self.play(Write(c6), run_time=1.5)
            self.wait(tr.until(0.5))
            self.play(Write(c8), run_time=1.8)
            self.wait(tr.until(0.85))
            self.play(Indicate(c8[0][-12:], color=YELLOW), Indicate(c5[0][-12:], color=YELLOW), run_time=1.2)

        c9 = M(r"\bt", r"\,\frac{1}{1-\ph}\,\dd{\ph}{z} \;=\; -", r"\bt", r"\left[\int_0^z(1-\ph)\,\mathrm{d}z'\right]^{n-1}\frac{\ph^m(1-\ph)}{\al}",
               size=42).move_to([0.3, -2.6, 0])
        t9 = eq_tag("C9").next_to(c9, RIGHT, buff=0.3)
        zero = M(r"\rs = 0", size=36, color=YELLOW).next_to(c8, RIGHT, buff=0.5)
        with self.vo("Now substitute all three into the porosity equation, with r s squared equal to zero. Beta appears on the left, from the advection "
                     "of porosity. And beta appears on the right, from the advection of grain size.") as tr:
            self.play(Write(zero), run_time=0.8)
            self.play(TransformMatchingShapes(VGroup(c1.copy(), c4.copy(), c5.copy(), c8.copy()), c9), FadeIn(t9), run_time=2.5)
            self.wait(tr.until(0.55))
            self.play(c9[0].animate.set_color(C_W).scale(1.3), run_time=0.8)
            self.wait(tr.until(0.8))
            self.play(c9[2].animate.set_color(C_R).scale(1.3), run_time=0.8)
        with self.vo("They cancel.") as tr:
            x1, x2 = Cross(c9[0], stroke_width=5), Cross(c9[2], stroke_width=5)
            self.play(Create(x1), Create(x2), run_time=0.8)
        box = label_box(VGroup(c9[1], c9[3]), color=YELLOW)
        with self.vo("What's left is a differential equation for phi that doesn't contain beta at all, whatever the value of n. "
                     "So, with zero surface grain size, steady-state firn thickness is completely independent of accumulation, "
                     "in the full model too.") as tr:
            self.play(Create(box), run_time=1)
            self.wait(tr.until(0.5))
            self.play(FadeIn(tau, shift=DOWN * 0.3), run_time=0.6)
            self.play(tau.change("happy", c9), run_time=0.8)
            self.play(tau.blink())
        self.wait(0.5)
        self.clear_all()
