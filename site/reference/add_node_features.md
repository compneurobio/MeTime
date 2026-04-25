# Function to add features to visnetwork plot from another plotter object

Function to add node features to see the nodes in the network that
affected differently

## Usage

``` r
add_node_features(object, results_indices, which_calculation = 1)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- results_indices:

  indices as a list to define which results to use. Eg:
  list(network=1/"name of the results", guide=2)

- which_calculation:

  index for plot_data to be used. Set to 1 by default.

## Value

network plotter object with new node colors/features
