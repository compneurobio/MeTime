# Choosing Methods

## Choosing Methods

> One-line summary: select the right analysis method for your biological
> question.

### Method map

| Question | Method family | Typical function(s) | Output |
|----|----|----|----|
| Distribution changes? | Distributions | `calc_distribution_*` | summary tables/plots |
| Data representation? | Dimensionality reduction | `calc_dimensionality_reduction_*` | results table/plots |
| Which features matter most? | Feature selection | `calc_featureselection_*` | ranked features/plots |
| Associations and trends? | Regression | `calc_lm`, `calc_lmm`, `calc_gamm` | model coefficients/plots |
| Data driven network? | Networks | `calc_ggm_*`, `calc_temporal_network`, `calc_ggm_multibipartite_lasso` | edge/node tables/plots |
| Metabolite clusters? | Clusters | `calc_clusters_wgcna` | clusters info/plots |
| Correlation, colinearity, and sample similarity? | Distances | `calc_distance_samples`, `calc_correlation_features`, `calc_colinearity_features` | results table/heatmaps |
|  |  |  |  |

### Decision checklist

- Is your target continuous/categorical?
- Do you need repeated-measures handling?
- Is missingness high enough to require imputation?
- Do you need interpretability vs predictive power?
