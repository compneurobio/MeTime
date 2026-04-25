# Calculation function to analyze the metabolite conservation index

Calculation function to analyze the metabolite conservation index

## Usage

``` r
calc_conservation_metabolite(
  object,
  which_data,
  verbose = F,
  cols_for_meta = NULL,
  stratifications = NULL,
  name = "conservation_index_metabolite_1"
)
```

## Arguments

- object:

  a S4 object of class metime_analyser.

- which_data:

  a character specifiying the name of the dataset to be used.

- verbose:

  a logical whether the status of processing should be printed. Default
  set to FALSE, no status information.

- cols_for_meta:

  a list of a character vectors that define column names to be used for
  plotting purposes. The characters should be named the same way as eg:
  list(lipid_data=c(id="id", sub_pathway="sub_pathway"),
  nmr_data=c(id="id", sub_pathway="Group")). The default is set to NULL,
  thereby not adding columns as metadata.

- stratifications:

  a list of parameters by which the data should be stratified. I.e.
  stratifications = list(timepoint=c(0,1,2)). In this case the dataset
  will be filtered to timepoints 0, 1, and 2. Default set to NULL, with
  no stratification applied.

- name:

  a character vector to define the name of the results generated. The
  length should be equal to which_data. The default is set to
  conservation_index_metabolite_1.

## Value

a S4 object with the analysis output saved in the results section
defined by name.

## Details

The metabolite conservation index compares the metabolite profiles
across all subjects at two time points. The metabolite conservation
index quantifies the conservation (stability) of a single metabolite in
comparison to all other measured metabolites. The processed data of the
metime_analyzer are used to calculate the metabolite conservation index
CI(x) for metabolite x as \$\$CI(x)=1-(rank(x)-1)/(N-1)\$\$, where
\$\$rank(x)=N-z and z=\|{y∈⟦1,N⟧≠x∨ρ_p (M_x^bl,M_x^fu )≥ρ_p
(M_x^bl,M_y^fu )}\|\$\$, with N being the number of subjects, M_i^bl
being the aggregated metabolite profile i across all subjects at
baseline, M_i^fu being the aggregated metabolite profile i of all
subjects at follow up, and ρ_p \$\$M_i^bl,M_j^fu\$\$) being the Pearson
correlation coefficient of metabolite i at baseline and metabolite j at
follow up for \$\$x,i,j∈⟦1,N⟧\$\$. The resulting index ranges between 0
and 1, with 0 = no conservation and 1 = maximal conservation. This
approach has been first described by [Yousri et al. 2014, Long-term
conservation of human metabolic phenotypes and link to heritability.
*Metabolomics*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4145193/)

## See also

[calc_conservation_metabotype](https://compneurobio.github.io/MeTime/reference/calc_conservation_metabotype.md),
[get_metadata_for_columns](https://compneurobio.github.io/MeTime/reference/get_metadata_for_columns.md),
[get_stratified_data](https://compneurobio.github.io/MeTime/reference/get_stratified_data.md)
