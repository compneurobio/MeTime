# Get summary on class edges

Function to check how the different edges in a GGM are associated to
their respective classes(it could be super-pathway or sub-pathway)

## Usage

``` r
get_class_info_from_edges(calc_networks, metadata, phenotypes)
```

## Arguments

- calc_networks:

  a list of calculated networks

- metadata:

  a dataframe containing the metadata of the edges present.

- phenotypes:

  a character vector to define phenotypes that were used for correcting
  the data

## Value

A dataframe with information on different type of edges present
