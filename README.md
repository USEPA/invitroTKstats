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
* Users will need the freely available "Just Another Gibbs Sampler" (JAGS):
<https://mcmc-jags.sourceforge.io/>
* Users will need to have the following package installed in addition to
`invitroTKstats`:
  * `ggplot2`
  * `parallel`
  * `runjags`
  * `stats4`
* Users will likely want a development environment like RStudio: <https://www.rstudio.com/products/rstudio/download/>

### Installing

Getting Started with R Package `invitroTKstats`.

* Installing directly from the GitHub repo from the R console
```
devtools::install_git(
  "https://github.com/jfwambaugh/invitroTKstats.git",
  subdir = "invitroTKstats",
  ref = "main"
)
```

* Installing a local clone of the GitHub repo.
  
  1. Go to the GitHub repo for
  [`invitroTKstats`](https://github.com/jfwambaugh/invitroTKstats)
  2. Choose the "Code" button and copy the repo URL.
  3. In your local command line terminal and navigate to the directory location
  you wish to store your local copy of the repo.
  4. In your terminal type
  `git clone <https://github.com/jfwambaugh/invitroTKstats>`
  4. After cloning completes then open an R session.
  5. In the R console use the following commands.
  
```
devtools::install_local(
  "<file_path_to_invitroTKstats_repo>/invitroTKstats"
)
```

* Installation may also be done via the RStudio provided ‘Install Packages’ menu
under the ‘Tools’ tab


## Loading the Package

To load the `invitroTKstats` data and functions into your local R session in
the R console 

```
library(invitroTKstats)
```

Check the package version installed and in use 

```
packageVersion(invitroTKstats)
```

## Authors

lead package developer John Wambaugh
[@wambaugh.john@epa.gov]

Anna Kreutz
[@kreutz.anna@epa.gov]

Marci Smeltz
[@smeltz.marci@epa.gov]

Barbara Wetmore
[@wetmore.barbara@epa.gov]

Caroline Ring
[@ring.caroline@epa.gov]

Sarah Davidson-Fritz
[@davidsonfritz.sarah@epa.gov]

## License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>
