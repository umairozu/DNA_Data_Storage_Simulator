
from __future__ import annotations

from cProfile import label

import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from typing import Optional, List

from numpy.ma.core import append

from Arrhenius_decay import Arrhenius_decay

def isDna(seq: str)-> bool:
    if not seq:
        return False
    for char in seq.upper():
        if char not in "ACGT":
            return False
    return True

def fig_5E(
        model: Arrhenius_decay, temp_C: list[float], encapsulated: bool,
        weeks: float, graph_points: int, out_png: str
            ):

    mode = "E-DNA" if encapsulated else "D-DNA"
    t_grid = np.linspace(0.0,weeks,graph_points)
    plt.figure(figsize=(7, 5))

    for T in temp_C:
        fractions = np.array([model.remaining_dna_frac(T,encapsulated,w) for w in t_grid], dtype= float)
        k_values = model.k(T,encapsulated)
        half_values = model.half_life(T, encapsulated)
        plt.plot(
            t_grid, fractions, linestyle = "--", label=f"{T:g}°C"
        )


    #plt.plot(t_grid,fractions)
    plt.ylim(1e-6, 1.2)
    plt.yticks([1, 1e-2, 1e-4, 1e-6], ["1", "0.01", "0.0001", "0.000001"])
    plt.yscale("linear")
    plt.xlabel("Weeks")
    plt.ylabel("DNA content C / C0")
    plt.legend(loc = "best", fontsize = 9)
    #plt.title(f"{mode} at {temp_C:g}°C | k={k_values:.3e} s⁻¹ | t½={half_values:.5g} years")
    plt.title(f"{mode}: DNA remaining vs time")
    plt.tight_layout()
    plt.savefig(out_png, dpi = 200)



def fig5G(
        model: Arrhenius_decay, temp_min: float, temp_max: float,
        step_c: float, mark_temp: Optional[list[float]],out_png: str
            ):

    temps = np.arange(float(temp_min), float(temp_max), float(step_c))

    half_D = np.array([model.half_life(float(T), encapsulated=False) for T in temps], dtype=float)
    half_E = np.array([model.half_life(float(T), encapsulated=True)  for T in temps], dtype=float)

    plt.figure(figsize=(7, 5))
    plt.plot(temps, half_D, linestyle="--", label="D-DNA")
    plt.plot(temps, half_E, linestyle="-",  label="E-DNA")

    if mark_temp is not None:
        for T in mark_temp:
            half_D_pt = model.half_life(T, encapsulated=False)
            half_E_pt = model.half_life(T, encapsulated=True)
            plt.scatter([T], [half_D_pt], marker="s", label=f"D @ {T:g}°C")
            plt.scatter([T], [half_E_pt], marker="o", label=f"E @ {T:g}°C")

    plt.yscale("log")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Half-life (years)")
    plt.title(" Half-life vs Temperature (°C)")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

def main():
    p = argparse.ArgumentParser(description="DNA cassette tape simulator CLI")

    #p.add_argument("--xlsx", required= True, help="Path to RawData.xlsx")
    p.add_argument("--seq", required=True, help="DNA sequence (A/C/G/T)")
    p.add_argument("--temp", required= True, action= "append", type= float, help="temperatures in celsius for Fig5E and 5G")
    p.add_argument("--encapsulated", action= "store_true", help= "if encapsulated use E-DNA, else D-DNA")

    p.add_argument("--weeks", type= float, default= 3.0, help= "Simulating upto mentioned number of weeks (default = 3)")

    p.add_argument("--temp_min", type= float, default= -20.0, help="Min temp for fig5G (default = -20.0)")
    p.add_argument("--temp_max", type= float, default= 80.0, help= "Max temp for fig5G (default = 80.0)")

    p.add_argument("--save", default= "output", help= "prefix for output file (e.g run1_fig5G.png)")

    args = p.parse_args()

    if not isDna(args.seq):
        sys.exit("Input sequence must only contain A/C/G/T characters")

    model = Arrhenius_decay.from_xlsx("RawData.xlsx")

    out_fig5E = f"{args.save}_fig5E.png"
    out_fig5G = f"{args.save}_fig5G.png"

    temp_list = args.temp

    fig_5E(
        model = model,
        temp_C = temp_list,
        encapsulated= args.encapsulated,
        weeks = args.weeks if args.weeks > 0.0 else sys.exit("Invalid num weeks"),
        graph_points= 200,
        out_png= out_fig5E
    )

    fig5G(
        model = model,
        temp_min = args.temp_min,
        temp_max = args.temp_max,
        step_c = 1.0,
        mark_temp = temp_list,
        out_png= out_fig5G
    )

    """Terminal outputs"""
    mode = "E-DNA (encapsulated DNA)" if args.encapsulated else "D-DNA (decapsulated DNA)"
    print("=== DcTS --DNA cassette tape simulator ===")
    print(f"Sequence length: {len(args.seq)} nt")
    print(f"Fig 5E mode: {mode}")
    print(f"Temps used: {temp_list}")

    #print(f"Temp for fig 5E: {args.temp} °C" )

    print(f"Saved fig 5E as: {out_fig5E}")
    print(f"Saved fig 5G as: {out_fig5G}")

if __name__ == "__main__":
        main()





