import ternary

figure, tax = ternary.figure(scale=1.0)

# Draw Boundary and Gridlines
tax.boundary(linewidth=1.0)
tax.gridlines(multiple=0.1,color="black")

tax.ticks(axis='lbr', multiple=0.125, linewidth=1, tick_formats = "%.1f", offset = 0.02)
tax.clear_matplotlib_ticks()

# Set Axis labels and Title
fontsize = 10
tax.set_title("Various Lines", fontsize=20)

tax.left_axis_label("$\\beta$ (silicate)", offset= 0.12)
tax.right_axis_label("$\\gamma$ (water)", offset= 0.12)
tax.bottom_axis_label("$\\alpha$ (iron)", offset= 0.12)

# Draw lines parallel to the axes
tax.horizontal_line(0.6)
tax.left_parallel_line(0.1, linewidth=2., color='red', linestyle="--")
tax.right_parallel_line(0.2, linewidth=3., color='blue')

tax.show()
