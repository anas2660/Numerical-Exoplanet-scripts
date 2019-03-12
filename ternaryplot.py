#!/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import ternary
import sys

colors = ["r", "y", "g", "b", "g", "y", "r"]


def plotsetup(figure):
    # Draw Boundary and Gridlines
    fig[1].boundary(linewidth=1.0)
    fig[1].gridlines(multiple=0.1, color="black")
    # Axis and tick appear
    fig[1].ticks(axis="lbr", multiple=0.1, linewidth=1.0,
                 tick_formats="%.1f", offset=0.02)
    fig[1].ax.set_aspect("equal", adjustable='box')
    fig[1].ax.axis("off")
    # Set Axis labels and Title
    fig[1].left_axis_label("$\\beta$ (silicate)", offset=0.16)
    fig[1].right_axis_label("$\\gamma$ (water)", offset=0.16)
    fig[1].bottom_axis_label("$\\alpha$ (iron)", offset=0.06)
    fig[1].set_title(arg, pad=16)


def plot(data, figure=None, **kwargs):
    if figure is None:
        fig, tax = ternary.figure()
    else:
        fig, tax = figure

    gammaData = np.ones(len(data[:, 0])) - data[:, 0] - data[:, 1]

    # Plot curve
    tax.plot(np.stack((data[:, 0], gammaData, data[:, 1]), axis=1), **kwargs)


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Too few arguments were passed.", file=sys.stderr)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            fig = ternary.figure()
            plotsetup(fig)
            for i in range(0, 7):
                try:
                    data = np.genfromtxt("out/"+arg+str(i)+".csv",
                                         delimiter=",")
                    plot(data, color=colors[i], figure=fig)
                except Exception as e:
                        print(e, file=sys.stderr)
            fig[1].savefig("out/"+arg+".pdf")
    plt.show()
else:
    print("Program not called directly.", file=sys.stderr)
