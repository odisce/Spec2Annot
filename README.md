
# Spec2Annot <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test coverage](https://codecov.io/gh/odisce/Spec2Annot/branch/main/graph/badge.svg)](https://app.codecov.io/gh/odisce/Spec2Annot?branch=main)
[![R-CMD-check](https://github.com/odisce/Spec2Annot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/odisce/Spec2Annot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Spec2Annot** is an [R](https://www.r-project.org/) package which implement different ways to annotate 
and query mass spectra. An implementation on [Galaxy](https://workflow4metabolomics.usegalaxy.fr/) of those 
functionalities is under progress.

## Installation

You can install the development version of Spec2Annot from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("odisce/Spec2Annot")
```

## Exemples
### Basic functions

**Spec2Annot** has built-in tables containings useful informations on:  
 - common adducts: `Spec2Annot::Adduct_db`  
 - common losses: `Spec2Annot::Losses_db`  
 - common isotopes: `Spec2Annot::Isotopes_db`  
 - common charges: `Spec2Annot::db_monocharge`  
 - a pre-computed list combining the above: `Spec2Annot::ions_database`  
 - atomic elements: `Spec2Annot::Element`  

Those tables are used to search as chemical-knowledge to annotate ions, but the package 
contains a set of functions to calculate mass from strings, like:  
  - Chemical formula: `Spec2Annot::mz_from_string("C6H12O3")`  
  - Chemical loss or adducts: `Spec2Annot::mz_from_string("C6H12O3-H2O+H")`  
  - Chemical notation as described by Damont *et al.* (2019): `Spec2Annot::mz_from_string("[C6H2O3-H2O+H]+")`  
  - And isotopes: `Spec2Annot::mz_from_string("C6H12O3_13C2")`  

Calculate mass:  
  - range using ppm: `Spec2Annot::mz_range(120.2532, ppm = 10)`  
  - ppm from range: `Spec2Annot::mz_ppm(120.2532, 120.2540)`  
  - ion from mass: `Spec2Annot::mz_calc_ion(120.2532, form = "-H")`  


### Elemental composition finder

**Spec2Annot** contains also more advanced function like a combinatory-search algorithm in `C++` 
to find any elemental composition from an exact mass. This algorithm can find common element but 
also any isotops or isotopologs. By default it will use the "seven golden rules" from Kind & Fiehn (2007) 
and the most common atomic element found in organic compounds, but any composition can be used.  

```{r}
Spec2Annot::find_compo_from_mass(
  mass_target = 120.25325,
  ppm = 5,
  use_golden_ratio = TRUE,
  elements_vc = NULL
)
```

For more advanced usage you can access the underlying function here: `Spec2Annot::brute_force_const()`.


### Spectra annotation

**Spec2Annot** contains a wrapper to calculate the elemental composition of a 
mass spectra by adding restriction of a known formula (useful with MS2 spectra of
known precursor) or not. It will also add some chemical scores/metrics like the Double
Bond Equivalent (DBE), Nitrogen rule (nrule) or Senior score (adapted from Morikawa
and Newbold, 2003).  

```{r}
Spec2Annot::annotate_mz(
  input_spectrum = spectra_ms2,
  ppm = 3,
  polarity = 1,
  compo = "C10H12N5O6P1"
)
```
