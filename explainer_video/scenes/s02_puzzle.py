from common import *


def moving_column(x, rate, n=12, top=2.2, width=1.5, height=4.4):
    """Layers that advect downward at a given rate (updater driven)."""
    layers = VGroup()
    for k in range(n):
        r = Rectangle(width=width, height=height / n, stroke_width=1, stroke_color=config.background_color)
        r.set_fill(layer_color(k, scale=4.0), 1).move_to([x, top - (k + 0.5) * height / n, 0])
        layers.add(r)
    frame = Rectangle(width=width, height=height, stroke_color=WHITE, stroke_width=2).move_to([x, top - height / 2, 0])
    return layers, frame


class Puzzle(VScene):
    def construct(self):
        title = Text("What does extra snowfall do?", font="CMU Serif", font_size=44).to_edge(UP)
        with self.vo("So, what does extra snowfall actually do to the firn? There are two competing effects."):
            self.play(Write(title), run_time=1.5)

        # ---------------- effect 1: advection
        h1 = Text("1. Advection", font="CMU Serif", font_size=34, color=C_W).move_to([-3.6, 2.5, 0])
        colA, frA = moving_column(-5.0, 0.4, top=1.9)
        colB, frB = moving_column(-2.2, 1.0, top=1.9)
        labA = M(r"\text{low } \bb", size=30).next_to(frA, DOWN)
        labB = M(r"\text{high } \bb", size=30).next_to(frB, DOWN)
        arrA = Arrow(UP * 0.35, DOWN * 0.35, color=C_W, buff=0).next_to(frA, RIGHT, buff=0.1)
        arrB = Arrow(UP * 0.7, DOWN * 0.7, color=C_W, buff=0).next_to(frB, RIGHT, buff=0.1)
        zA = DashedLine(frA.get_left() + DOWN * 0.3, frA.get_right() + DOWN * 0.3, color=C_Z, stroke_width=4)
        zB = DashedLine(frB.get_left() + DOWN * 1.2, frB.get_right() + DOWN * 1.2, color=C_Z, stroke_width=4)
        zlabA = M(r"\zc", size=26).next_to(zA, LEFT, buff=0.08)
        zlabB = M(r"\zc", size=26).next_to(zB, LEFT, buff=0.08)

        def advect(speed, frame):
            n = 12
            hgt = frame.height / n

            def upd(group, dt):
                top = frame.get_top()[1]
                for r in group:
                    r.shift(DOWN * speed * dt)
                    if r.get_y() < top - frame.height:
                        r.shift(UP * frame.height)
                # recolour by depth so the column keeps its density structure
                for r in group:
                    s = (top - r.get_y()) / frame.height
                    r.set_fill(layer_color(s * 12, scale=4.0 + 4 * (speed > 0.6)), 1)
                    r.set_opacity(1 if r.get_top()[1] <= top + 1e-6 else 0)
            return upd

        with self.vo("The first is advection. If snow piles up faster, every layer is buried faster. "
                     "Low-density firn from near the surface gets carried deeper before it has time to compact. "
                     "On its own, that makes the firn thicker.") as tr:
            self.play(FadeIn(h1), FadeIn(colA), FadeIn(colB), Create(frA), Create(frB), FadeIn(labA), FadeIn(labB), run_time=1.2)
            colA.add_updater(advect(0.25, frA))
            colB.add_updater(advect(0.7, frB))
            self.play(GrowArrow(arrA), GrowArrow(arrB), run_time=1)
            self.wait(tr.until(0.6))
            self.play(Create(zA), Create(zB), FadeIn(zlabA), FadeIn(zlabB), run_time=1.2)
            self.fill(tr)

        # ---------------- effect 2: compaction / Herron-Langway
        h2 = Text("2. Compaction", font="CMU Serif", font_size=34, color=C_SIG).move_to([3.4, 2.5, 0])
        hl = M(r"\DDt{\rho}", "=", r"C\,", r"\bb^{\,\alpha}", r"(\rho_i - \rho)", size=46).move_to([3.4, 1.0, 0])
        hl_ref = Text("Herron & Langway (1980)", font="CMU Serif", font_size=24, color=GREY_B).next_to(hl, DOWN, buff=0.3)
        hl_tag = eq_tag(1).next_to(hl, RIGHT, buff=0.3)
        weight = VGroup(*[Arrow(UP * 0.5, DOWN * 0.1, color=C_SIG, buff=0).shift(RIGHT * x) for x in [-0.4, 0, 0.4]]).move_to([3.4, -1.0, 0])
        parcel = Square(0.6).set_fill(SNOW, 0.8).set_stroke(WHITE, 2).next_to(weight, DOWN, buff=0.05)
        with self.vo("The second is compaction. More snow means more weight piling onto each parcel of firn, more quickly. "
                     "Many firn models, going back to Herron and Langway in 1980, put the accumulation rate directly into the compaction law. "
                     "Their equation says the rate at which a parcel's density increases is proportional to the accumulation rate, b, "
                     "raised to a power close to one, times how far the density still is from that of ice.") as tr:
            self.play(FadeIn(h2), run_time=0.6)
            self.play(FadeIn(weight, shift=DOWN * 0.3), FadeIn(parcel), run_time=1)
            self.play(parcel.animate.stretch(0.7, 1, about_edge=DOWN), run_time=1)
            self.wait(tr.until(0.38))
            self.play(Write(hl), FadeIn(hl_tag), FadeIn(hl_ref), run_time=2)
            self.wait(tr.until(0.7))
            self.play(Circumscribe(hl[3], color=C_B), run_time=1.5)

        tau = TauCreature(height=1.7).to_corner(DR).shift(UP * 0.0 + LEFT * 0.2)
        with self.vo("But hold on. If firn behaves like a viscous material, compaction should be driven by the stress on a parcel, "
                     "the weight of everything above it, not by how quickly that weight is changing. "
                     "So why should the accumulation rate appear in the compaction law at all?") as tr:
            self.play(FadeIn(tau, shift=UP * 0.4), FadeOut(weight), FadeOut(parcel), run_time=1)
            b, anim = tau.thinks(r"Stress, not the\\rate it builds up?", mood="confused")
            self.wait(tr.until(0.45))
            self.play(anim, run_time=1.2)
            self.wait(tr.until(0.85))
            self.play(tau.blink())
        self.play(FadeOut(b), run_time=0.5)

        with self.vo("There's a practical problem too. With b inside the compaction law, a change in snowfall today would instantly "
                     "change the compaction rate everywhere in the column, even in firn that was buried centuries ago.") as tr:
            colB.clear_updaters()
            flash = Rectangle(width=frB.width, height=frB.height).move_to(frB).set_fill(C_SIG, 0.55).set_stroke(C_SIG, 4)
            now = Text("all at once!", font="CMU Serif", font_size=26, color=C_SIG).next_to(frB, RIGHT, buff=0.5)
            self.wait(tr.until(0.35))
            self.play(FadeIn(flash), FadeIn(now), run_time=0.6)
            self.play(FadeOut(flash), run_time=0.6)
            self.play(FadeIn(flash), run_time=0.6)
            self.play(FadeOut(flash), FadeOut(now), run_time=0.6)
            self.play(tau.change("sad", frB))

        self.play(*[FadeOut(m) for m in [h1, colA, colB, frA, frB, labA, labB, arrA, arrB, zA, zB, zlabA, zlabB, hl_ref]],
                  FadeOut(h2), VGroup(hl, hl_tag).animate.scale(0.7).to_corner(UL).shift(DOWN * 0.9), run_time=1)
        colA.clear_updaters()

        # ---------------- this paper: b only as a boundary condition
        col = firn_column(width=1.7, height=4.4, top=np.array([0.0, 1.8, 0]))
        flakes = VGroup(*[Snowflake(0.09).move_to([x, 3.1 + 0.2 * np.sin(5 * x), 0]) for x in np.linspace(-0.7, 0.7, 6)])
        bc_arrow = Arrow([0, 2.9, 0], [0, 1.85, 0], color=C_B, buff=0, stroke_width=6)
        bc_lab = M(r"\bb", r"\text{ enters only here}", size=34).next_to(bc_arrow, RIGHT, buff=1.0)
        bc_lab[0].scale(1.3)
        with self.vo("This paper takes a cleaner approach. The accumulation rate enters only where it physically belongs: "
                     "as a boundary condition at the surface, setting how fast new snow arrives. Any effect it has on compaction "
                     "has to emerge from the physics itself. And, as we'll see, a lot of that physics runs through grain size.") as tr:
            self.play(FadeOut(title), FadeIn(col), run_time=1)
            self.play(GrowArrow(bc_arrow), FadeIn(flakes, shift=DOWN * 0.3), Write(bc_lab), tau.change("happy", col), run_time=1.5)
            self.wait(tr.until(0.85))
            self.play(tau.blink())
        self.wait(0.4)
        self.clear_all()
