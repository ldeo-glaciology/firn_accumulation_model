from common import *


class SteadyODE(VScene):
    def construct(self):
        title = Text("A steady-state shortcut", font="CMU Serif", font_size=40).to_edge(UP, buff=0.3)
        pdes = VGroup(
            M(r"\pd{\ph}{t}", r"= -\frac{1}{\al}\frac{|\sg|^n\ph^m(1-\ph)}{\rsq} - \ww\pd{\ph}{z}", size=36),
            M(r"\pd{\rsq}{t}", r"= (1-\dl\rsq) - \ww\pd{\rsq}{z}", size=36),
            M(r"\pd{\Ag}{t}", r"= 1 - \ww\pd{\Ag}{z}", size=36),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to([0, 0.6, 0])
        zeros = VGroup(*[M(r"0", size=40, color=YELLOW).move_to(p[0]) for p in pdes])
        with self.vo("If we only care about the steady state, there's a shortcut. Set every time derivative to zero.") as tr:
            self.play(Write(title), FadeIn(pdes), run_time=1.5)
            self.play(*[Transform(p[0], z) for p, z in zip(pdes, zeros)], run_time=1.2)
        ode = VGroup(
            M(r"\dd{\ph}{z} = -\frac{|\sg|^n\ph^m(1-\ph)}{\al\,\ww\,\rsq}", size=34),
            M(r"\dd{\sg}{z} = -(1-\ph)", size=34),
            M(r"\dd{\ww}{z} = -\frac{|\sg|^n\ph^m}{\al\,\rsq}", size=34),
            M(r"\dd{\rsq}{z} = \frac{1-\dl\rsq}{\ww}", size=34),
            M(r"\dd{\Ag}{z} = \frac{1}{\ww}", size=34),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22).move_to([-3.6, -0.4, 0])
        tag = eq_tag(22).next_to(ode, RIGHT, buff=0.3).align_to(ode, UP)
        bcs = M(r"\ph(0)=\ph_s,\ \ \sg(0)=0,\ \ \ww(0)=\tfrac{\bt}{1-\ph_s},\ \ \rsq(0)=\rs,\ \ \Ag(0)=0", size=26).next_to(ode, DOWN, buff=0.35).align_to(ode, LEFT)
        with self.vo("The partial differential equations collapse into a set of ordinary differential equations in depth alone, "
                     "which we can integrate straight down from the surface, starting from the surface values.") as tr:
            self.play(FadeOut(pdes), FadeOut(zeros), run_time=0.5)
            self.play(LaggedStart(*[FadeIn(o, shift=RIGHT * 0.2) for o in ode], lag_ratio=0.15), FadeIn(tag), run_time=2.5)
            self.play(FadeIn(bcs), run_time=1)

        d = load("fig2")
        z, zo = d["z"][-1], d["ode_z"]
        plot = DepthPlot([-1, 2.0], [0, 1.0], width=4.4, height=5.0, x_ticks=[-1, 0, 1, 2], z_ticks=[0, 0.5, 1.0], x_dec=0, z_dec=1,
                         x_label=r"\text{nondimensional variable}", z_label=r"z", num_size=20, label_size=28)
        plot.shift(np.array([1.6, 2.2, 0]) - plot.top.get_start())
        items = [("phi", C_PHI, r"\ph"), ("sigma", C_SIG, r"\sg"), ("w", C_W, r"\ww"), ("r2", C_R, r"\rsq"), ("A", C_A, r"\Ag")]
        full = VGroup(*[plot.curve(d[k][-1], z, color=c, width=6) for k, c, _ in items])
        odec = VGroup(*[DashedVMobject(plot.curve(d["ode_" + k], zo, color=BLACK, width=2.5), num_dashes=40) for k, c, _ in items])
        legend = VGroup(*[M(s, size=28) for _, _, s in items]).arrange(RIGHT, buff=0.35).next_to(plot, DOWN, buff=0.25)
        key = VGroup(Line(ORIGIN, RIGHT * 0.5, stroke_width=6, color=GREY_A), Text("full model", font="CMU Serif", font_size=22),
                     DashedLine(ORIGIN, RIGHT * 0.5, stroke_width=2.5, color=BLACK).set_stroke(GREY_A), Text("ODE model", font="CMU Serif", font_size=22))
        key.arrange(RIGHT, buff=0.15).next_to(legend, DOWN, buff=0.2)
        with self.vo("Here's that ODE solution, the thin dashed lines, laid over the steady state of the full time-dependent simulation, "
                     "the thick coloured lines. They're indistinguishable.") as tr:
            self.play(FadeIn(plot), run_time=0.8)
            self.play(LaggedStart(*[Create(c) for c in full], lag_ratio=0.15), FadeIn(legend), run_time=2)
            self.play(LaggedStart(*[Create(c) for c in odec], lag_ratio=0.15), FadeIn(key), run_time=2)
        mis = np.mean([np.mean(np.abs(np.interp(z, zo, d["ode_" + k]) - d[k][-1])) for k, _, _ in items])
        mx = np.max([np.max(np.abs(np.interp(z, zo, d["ode_" + k]) - d[k][-1])) for k, _, _ in items])
        note = M(r"\text{mean mismatch} \approx %.1f\times10^{-3}" % (mis * 1e3), r"\\ \text{max} \approx %.1f\times 10^{-3}" % (mx * 1e3),
                 size=30).move_to([-3.2, -3.3, 0])
        tau = TauCreature(height=1.3).move_to([-6.2, -2.9, 0])
        with self.vo("The mismatch is around a thousandth, and it shrinks as the grid gets finer. So now we have a lightweight model that we trust, "
                     "and that we can start to pull apart.") as tr:
            self.play(FadeOut(bcs), Write(note), run_time=1.5)
            self.wait(tr.until(0.55))
            self.play(FadeIn(tau, shift=UP * 0.3), run_time=0.6)
            self.play(tau.change("happy", plot), run_time=0.6)
            self.play(tau.blink())
        self.wait(0.4)
        self.clear_all()
