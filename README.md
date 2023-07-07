# CREMA (Control of Regulation Extracted from Multiomic Assays)

CREMA is a computational framework for recovering regulatory circuits from single cell multiome data. The framework doesn't rely on chromatin peak calling, which greatly increases its recovery of true regulatory regions missed by peak calling algorithms. You can read the [manuscript](https://doi.org/10.1101/2023.06.23.544355) for more details.

![](fig_crema_method.png?raw=true)

## Install

```
library(devtools)
devtools::install_github("zidongzh/CREMA")
```

## Vignettes

Please see `vignettes/` for detailed usage of CREMA.


## Citation

**Peak-agnostic high-resolution cis-regulatory circuitry mapping using single cell multiome data**

Zidong Zhang, Frederique Ruf-Zamojski, Michel Zamojski, Daniel J. Bernard, Xi Chen, Olga G. Troyanskaya, Stuart C. Sealfon

bioRxiv 2023.06.23.544355; doi: https://doi.org/10.1101/2023.06.23.544355

