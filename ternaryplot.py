import ternary

def plot(data, save="", show=False, **kwargs):

    figure, tax = ternary.figure()

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(multiple=0.1, color="black")

    tax.ticks(axis='lbr', multiple=0.1, linewidth=1.0, tick_formats = "%.1f", offset = 0.02)
    tax.ax.set_aspect("equal",adjustable='box')
    tax.clear_matplotlib_ticks()
    tax.ax.axis("off")

    # Set Axis labels and Title
    tax.left_axis_label("$\\beta$ (silicate)", offset= 0.16)
    tax.right_axis_label("$\\gamma$ (water)", offset= 0.16)
    tax.bottom_axis_label("$\\alpha$ (iron)", offset= 0.06)

    tax.plot(data, **kwargs)
    
    if save != "":
        figure.savefig(save, bbox_inches="tight")
    
    if show:
        tax.show()
        
if __name__ == "__main__":
    plot([[0.8,0.15,0.05],[0.2,0.3,0.5],[0.4,0.4,0.2]], show=True, color="red")
    plot([[0.8,0.15,0.05],[0.2,0.3,0.5],[0.4,0.3,0.3]], show=True)
