# Setting a plotting method for the metime_analyser class

Function to plot results of a certain calculation

## Usage

``` r
# S4 method for class 'metime_analyser,ANY'
plot(x, results_index, interactive, plot_type, ...)
```

## Arguments

- x:

  An S4 object of class metime_analyser

- results_index:

  Index/name of the results to be plotted

- interactive:

  logical. Set TRUE for interactive plot

- plot_type:

  to define the type of plot. Accepted inputs are "dot", "tile", "box",
  "forrest", "manhattan"

- ...:

  other parameters to pass color, fill, strat, viz(character vector with
  colnames for interactive)

## Value

plots for a certain set of results

## See also

[mod_generate_plots](https://compneurobio.github.io/MeTime/reference/mod_generate_plots.md)
