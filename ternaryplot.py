#!/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import ternary
import sys

colors = ["r", "y", "g", "b", "g", "y", "r"]

def plot(data, save="", show=False, figure=None, **kwargs):
    
    if figure == None:
        fig, tax = ternary.figure()
    else:
        fig, tax = figure
    
    gammaData = np.ones(len(data[:,0])) - data[:,0] - data[:,1]
    
    # Draw Boundary and Gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(multiple=0.1, color="black")
    
    tax.ticks(axis="lbr", multiple=0.1, linewidth=1.0, tick_formats = "%.1f", offset = 0.02)
    tax.ax.set_aspect("equal",adjustable='box')
    tax.ax.axis("off")
    
    # Set Axis labels and Title
    tax.left_axis_label("$\\beta$ (silicate)", offset= 0.16)
    tax.right_axis_label("$\\gamma$ (water)", offset= 0.16)
    tax.bottom_axis_label("$\\alpha$ (iron)", offset= 0.06)

    tax.plot(np.stack((data[:,0], gammaData, data[:,1]),axis=1), **kwargs)
    
    if save != "":
        fig.savefig(save, bbox_inches="tight")
    

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Too few arguments were passed.", file=sys.stderr)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            fig = ternary.figure()
            fig[1].set_title(arg, pad= 16)
            for i in range(0,7):
                try:
                    data = np.genfromtxt("out/"+arg+str(i)+".csv", delimiter=",")
                    plot(data, show=True, color=colors[i], figure=fig)
                except Exception as e:
                        print(e, file=sys.stderr)
    plt.show()

else:
    print("Program not called directly.", file=sys.stderr)
