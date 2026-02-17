# Choosing Methods

> One-line summary: select the right analysis method for your biological question.

## Method map

| Question | Method family | Typical function(s) | Output |
|---|---|---|---|
| Distribution changes? | Distributions | `calc_distribution_*` | summary tables/plots |
| Which features matter most? | Feature selection | `calc_featureselection_*` | ranked features |
| Time association? | Regression | `calc_lm`, `calc_lmm`, `calc_gamm` | model coefficients |
| Network rewiring? | Networks | `calc_ggm_*`, `calc_temporal_network` | edge/node tables |

## Decision checklist

- Is your target continuous/categorical?
- Do you need repeated-measures handling?
- Is missingness high enough to require imputation?
- Do you need interpretability vs predictive power?
