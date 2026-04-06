import matplotlib.pyplot as plt
from typing import Optional
import numpy as np

from Arrhenius_decay import Arrhenius_decay
temps = [60, 65, 70]

class PlotClass:

    model: Arrhenius_decay


    def plot_fig5E(self, out_png: Optional[str] = None) -> None:

        weeks = [0, 1, 2, 3]

        plt.figure(figsize=(8, 6))

        for temp in temps:
            if temp >= 70:
                weeks_D = [0, 1, 2]
            else:
                weeks_D = weeks

            sim_D = [self.model.remaining_dna_frac(temp, False, w) for w in weeks_D]
            plt.plot(weeks_D, sim_D, linestyle=":", marker=None, label=f"D-DNA sim {temp}°C")

            sim_E = [self.model.remaining_dna_frac(temp, True, w) for w in weeks]
            plt.plot(weeks, sim_E, linestyle="-", marker=None, label=f"E-DNA sim {temp}°C")

        plt.yscale("linear")
        plt.xlabel("Time (weeks)")
        plt.ylabel("Remaining fraction (C/C0)")
        plt.title("Fig 5E)")
        plt.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.3, 1))
        plt.tight_layout()

        if out_png:
            plt.savefig(out_png, dpi=200)
        plt.show()

        """
        Fig 5G: half-life vs temperature from Arrhenius equation k(T).
        """
    def plot_fig5G(self, out_png: Optional[str] = None, temp_min: float = -20.0, temp_max: float = 80.0):
        temp_space = np.linspace(temp_min,temp_max,200)

        half_D = [self.model.half_life(float(t),False) for t in temp_space]
        half_E = [self.model.half_life(float(t), True) for t in temp_space]

        marker_temps = [60,65,70]
        markers_dDNA = [self.model.half_life(float(t),False) for t in marker_temps]
        markers_eDNA = [self.model.half_life(float(t), True) for t in marker_temps]


        plt.figure(figsize= (8,6))
        line_dDNA, = plt.plot(temp_space,half_D,label = "D-DNA", linestyle = ":", color = "blue")
        plt.plot(marker_temps,markers_dDNA,"o", color = line_dDNA.get_color(), markersize = 8)

        line_eDNA, = plt.plot(temp_space,half_E, label = "E-DNA", linestyle = "-", color = "red")
        plt.plot(marker_temps, markers_eDNA, "s", color = line_eDNA.get_color(), markersize = 8)

        plt.yscale("log")
        plt.xlabel("Temperature (Celsius)")
        plt.ylabel("Half-life (years)")
        plt.title("Fig 5G")
        plt.legend()
        if out_png:
            plt.savefig(out_png,dpi = 200)
        plt.show()




















#<------------------------------EXTRA's--------------------------------->#

"""    def plot_fig5E(self, xlsx_path: str, out_png: Optional[str] = None) -> None:
        temps = [60, 65, 70]
        weeks = [0, 1, 2, 3]
        tgrid = np.linspace(0.0, 3.0, 200)

        emp_D = self.model.empirical_fractions(xlsx_path, encapsulated=False)
        emp_E = self.model.empirical_fractions(xlsx_path, encapsulated=True)

        plt.figure(figsize=(8, 6))
        ax = plt.gca()

        eps = 1e-12  # for log safety

        for temp in temps:
            # ---- D-DNA simulation line (skip after week 2 if temp >= 70) ----
            if temp >= 70:
                tgrid_D = tgrid[tgrid <= 2.0]  # stop at 2 weeks
            else:
                tgrid_D = tgrid

            sim_D = [self.model.remaining_dna_frac(temp, False, w) for w in tgrid_D]
            line_d, = ax.plot(tgrid_D, sim_D, linestyle=":", label=f"D-DNA sim {temp}°C")
            color_d = line_d.get_color()

            # ---- E-DNA simulation line (0..3 weeks) ----
            sim_E = [self.model.remaining_dna_frac(temp, True, w) for w in tgrid]
            line_e, = ax.plot(tgrid, sim_E, linestyle="-", label=f"E-DNA sim {temp}°C")
            color_e = line_e.get_color()

            # ---- Empirical points (mean ± SEM) ----
            def add_points(emp, isE: bool, plot_color):
                xs, ys = [], []
                yerr_low, yerr_high = [], []

                for w in weeks:
                    vals = emp.get((temp, w), [])
                    if not vals:
                        continue

                    # log safety
                    vals = [max(v, eps) for v in vals]

                    y = float(np.mean(vals))
                    sem = float(np.std(vals, ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else 0.0

                    # asymmetric error bars so lower doesn't go <= 0 on log axis
                    low = min(sem, max(y - eps, 0.0))
                    high = sem

                    xs.append(w)
                    ys.append(max(y, eps))
                    yerr_low.append(low)
                    yerr_high.append(high)

                marker = "o" if isE else "s"
                ax.errorbar(
                    xs, ys,
                    yerr=[yerr_low, yerr_high],
                    fmt=marker, capsize=3, linestyle="none",
                    color=plot_color, ecolor=plot_color,
                    label=f"{'E' if isE else 'D'}-DNA data {temp}°C"
                )

            add_points(emp_D, isE=False, plot_color=color_d)
            add_points(emp_E, isE=True, plot_color=color_e)

        # ---- Axis styling to match the paper ----
        ax.set_yscale("log")
        ax.set_xlabel("Time (weeks)")
        ax.set_ylabel("Remaining fraction (C/C0)")
        ax.set_title("Fig 5E")

        ax.set_ylim(1e-6, 1.2)
        ax.set_yticks([1, 1e-2, 1e-4, 1e-6])
        ax.set_yticklabels(["1", "0.01", "0.0001", "0.000001"])

        ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.3, 1))
        plt.tight_layout()

        if out_png:
            plt.savefig(out_png, dpi=200)
        plt.show()"""





