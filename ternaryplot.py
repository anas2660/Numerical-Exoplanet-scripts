#!/bin/env python3
import numpy
import ternary
import sys

def plot(data, save="", show=False, figure=None, **kwargs):
    
    if figure == None:
        fig, tax = ternary.figure()
    else:
        fig, tax = figure
    
    # Draw Boundary and Gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(multiple=0.1, color="black")
    
    tax.ticks(axis="lbr", multiple=0.1, linewidth=1.0, tick_formats = "%.1f", offset = 0.02)
    tax.ax.set_aspect("equal",adjustable='box')
    #tax.clear_matplotlib_ticks()
    tax.ax.axis("off")
    
    # Set Axis labels and Title
    tax.left_axis_label("$\\beta$ (silicate)", offset= 0.16)
    tax.right_axis_label("$\\gamma$ (water)", offset= 0.16)
    tax.bottom_axis_label("$\\alpha$ (iron)", offset= 0.06)

    tax.plot(data, **kwargs)
    
    if save != "":
        fig.savefig(save, bbox_inches="tight")
    
    if show:
        tax.show()
    
if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Too few arguments were passed.", file=sys.stderr)
        sys.exit(1)
    else:
        fig = ternary.figure()
        for arg in sys.argv[1:]:
            try:
                data = numpy.genfromtxt(arg, delimiter=",")
                plot(data, show=True, color="red", figure=fig)
            except Exception as e:
                print(e, file=sys.stderr)
else:
    print("Program not called directly.", file=sys.stderr)