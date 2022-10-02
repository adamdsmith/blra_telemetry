source("R/vis_reception.R")

# Hexagonal centers based on consistent detection range
node_range <- 50
x_range <- 200
y_range <- 200

# Default spaces nodes based on expected consistent, minimum detection distance
vis_reception(x_range, y_range, node_range, res = 1)

# Space about 20% under expected detection range produces decent results
# except at some edges
vis_reception(x_range, y_range, node_range, spacing = 40)
vis_reception(x_range, y_range, node_range, spacing = 30)

node_range <- 75
vis_reception(x_range, y_range, node_range, spacing = 50)

