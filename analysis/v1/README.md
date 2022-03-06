# ulcerative colitis and systemic corticosteroid therapy

Project with Andreas Blesl (Medical University Graz) and Konrad Aden

*January 2022*

----

### Preparing data

Cleaning/filtering input data and save it in formats for easy use in later steps.

```shell
Rscript analysis/v1/scripts/00_clean_data.R
```

Map ASV-16S sequences to HRGM 16S rRNA genes

```shell
bash analysis/v1/scripts/01_map-asv_to_HRGM.sh
```



### Run metabolism simulations

From within R:

```r
source("analysis/v1/scripts/02_communityFBA.")
```