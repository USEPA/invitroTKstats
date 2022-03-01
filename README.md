# R Package "invitroTKstats"

Standardizes the documentation and statistical analysis for in vitro
Assays that allow prediction of toxicokinetics (that is, the absorption, 
distribution, metabolism, and elimination of chemicals by the body). The assays 
covered include intrinsic clearance hepatocyte incubations; three variants of plasma 
protein binding experiments, CACO-2 membrane permeability, and blood to plasma concentration ratio. 
Analysis methods include frequentist point estimates and, in some cases, Bayesian methods for identifying 
distributions of likely parameter values. Analysis is based on mass spectrometry ratios of analyte to 
internal standard peak areas. Data are formatted for loading into databases.


## Background


## Getting Started

### Dependencies

* Users will need the freely available R statistical computing language: <https://www.r-project.org/>
* Users will likely want a development environment like RStudio: <https://www.rstudio.com/products/rstudio/download/>

### Installing

* Getting Started with R Package bayesmarker from the R command line
```
library(devtools)
install_github("HumanExposure/bayesmarker")
```
* RStudio provides a menu ‘Install Packages’ under ‘Tools’ tab
* Load the bayesmarker data and functions
```
library(bayesmarker)
```
* Check what version you are using 
```
packageVersion(bayesmarker)
```

## Authors

lead package developer John Wambaugh
[@stanfield.zach@epa.gov]

Anna Kreutz
[@kreutz.anna@epa.gov]

Marci Smeltz
[@smeltz.marci@epa.gov]

Barbara Wetmore
[@wetmore.barbara@epa.gov]



## License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>