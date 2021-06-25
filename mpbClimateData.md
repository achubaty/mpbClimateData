---
title: "mpbClimate"
author:
  - "Alex M. Chubaty"
  - "Barry J. Cooke"
  - "Eliot J. B. McIntire"
date: "June 3, 2021"
output:
  html_document:
    keep_md: yes
    toc: true
    toc_depth: 3
    toc_float: true
bibliography:
  - bibliography.bib
editor_options:
  chunk_output_type: console
---



# Overview

Over the last 20 years, and using four indices of MPB climatic suitability, there appears to be a long-term trend in enhanced MPB survival throughout BC and Alberta [@Cooke:2017fem], fuelled by accelerating warmth through the last 20 years.
The effects of temperature on seasonal biology in the US Rocky Mountains have been modelled by Logan et al. (1992).
The Logan suitability index (L) is based on summer temperatures [@Logan:2003fr].
The effects of temperature on winter mortality in the US Rocky Mountains have been modelled by @Regniere:2007ip.
In the Canadian context, the many effects of temperature on MPB recruitment in British Columbia have been modelled by @Safranyik:1975bk.
The Safranyik index (S) was subsequently extended to Alberta by @Carroll:2004rm.
Once the MPB had come into Alberta in large numbers in 2006 the Régnière index (R) was field validated in Alberta [@Cooke:2009ow].
The North America-wide outputs of these models under standard climate change scenarios were first presented in @Nealis:2008gc, and later received peer-review validation through @Bentz:2010bs and @Safranyik:2010ce.
A composite SLR index (G, the default) takes the geometric mean of these three models. These are described in further detail in @Nealis:2008gc and @Nealis:2014re.

![](figures/fig4_Cooke&Carroll_climate_indices.png)

Our simulation models explore 7 MPB climate scenarios (using each of the four indices: S, L, R, G) used as model drivers.
We use the `BioSIM` R package [@Fortin:2021] to get MPB SLR data from BioSim software [@Regniere:1995BioSim] to generate the climate suitability maps [see @Bentz:2010bs; @Logan:2003fr; @Safranyik:2010ce).
The values of each of these indices are bound between 0 and 1, and this value is used to scale the vertical shift of the red-top recruitment curve.
All climate maps are projected using a Lambert Conformal Conic projection and cover all of Canada.
Where possible, all data downloads and preprocessing were scripted for reproducibility from raw, original sources.

# Usage


```r
library(SpaDES.core)

setPaths(modulePath = file.path(".."))
getPaths() # shows where the 4 relevant paths are

times <- list(start = 0, end = 10)

parameters <- list(
  #.progress = list(type = "text", interval = 1), # for a progress bar
  ## If there are further modules, each can have its own set of parameters:
  #module1 = list(param1 = value1, param2 = value2),
  #module2 = list(param1 = value1, param2 = value2)
)
modules <- list("mpbClimateData")
objects <- list()
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
```

# Parameters

Provide a summary of user-visible parameters.


```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|paramName        |paramClass |default    |min  |max  |paramDesc                                                                                                                                            |
|:----------------|:----------|:----------|:----|:----|:----------------------------------------------------------------------------------------------------------------------------------------------------|
|climateScenario  |character  |RCP45      |NA   |NA   |The climate scenario to use. One of RCP45 or RCP85.                                                                                                  |
|suitabilityIndex |character  |R          |NA   |NA   |The MPB climatic suitabilty index to use. One of 'S', 'L', 'R', or 'G'.                                                                              |
|windMonths       |integer    |7          |1    |12   |A vector of length 1 or more, indicating the month(s) of the year from which to extract wind data. This should correspond to the months of dispersal |
|.maxMemory       |numeric    |1e+09      |NA   |NA   |Used to set the 'maxmemory' raster option. See '?rasterOptions'.                                                                                     |
|.plots           |character  |screen     |NA   |NA   |This describes the simulation time at which the first plot event should occur                                                                        |
|.plotInitialTime |numeric    |start(sim) |1981 |2100 |This describes the simulation time at which the first plot event should occur                                                                        |
|.plotInterval    |numeric    |NA         |NA   |NA   |This describes the simulation time interval between plot events                                                                                      |
|.saveInitialTime |numeric    |NA         |NA   |NA   |This describes the simulation time at which the first save event should occur                                                                        |
|.saveInterval    |numeric    |NA         |NA   |NA   |This describes the simulation time interval between save events                                                                                      |
|.tempdir         |character  |           |NA   |NA   |Temporary (scratch) directory to use for transient files (e.g., GIS intermediates).                                                                  |
|.useCache        |logical    |FALSE      |NA   |NA   |Should this entire module be run with caching activated?                                                                                             |

# Events

Describe what happens for each event type.

## Plotting

Write what is plotted.

## Saving

Write what is saved.

# Data dependencies

## Input data

How to obtain input data, and a description of the data required by the module.
If `sourceURL` is specified, `downloadData("mpbClimateData", "..")` may be sufficient.


```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|objectName          |objectClass     |desc                                                                                                                                                                                                                                                                                                                                                                           |sourceURL                                                                                                                                                                                   |
|:-------------------|:---------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|climateMapRandomize |list            |List with 2 elements: srcYears, rmYears. Each of these should be a numeric vector of years (sim$climateSuitabilityMaps must have layer names that contain the full 4-digit numeric year, e.g., X2010), with 'scrYears' being the pool of years to take data from to use in 'rmYears'. This is intended for cases where future projected climate data is unavailable or corrupt |https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing                                                                                                          |
|climateMapFiles     |character       |Vector of filenames correspoding to climate suitablity map layers                                                                                                                                                                                                                                                                                                              |https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing                                                                                                          |
|rasterToMatch       |RasterLayer     |if not supplied, will default to standAgeMap                                                                                                                                                                                                                                                                                                                                   |NA                                                                                                                                                                                          |
|standAgeMap         |RasterLayer     |stand age map in study area, default is Canada national stand age map                                                                                                                                                                                                                                                                                                          |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif |
|studyArea           |SpatialPolygons |The study area to which all maps will be cropped and reprojected.                                                                                                                                                                                                                                                                                                              |NA                                                                                                                                                                                          |

## Output data

Description of the module outputs.


```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|objectName             |objectClass |desc                                                                                                |
|:----------------------|:-----------|:---------------------------------------------------------------------------------------------------|
|climateSuitabilityMaps |RasterStack |A time series of climatic suitablity RasterLayers, each with previx 'X' and the year, e.g., 'X2010' |
|windDirStack           |RasterStack |RasterStack of dominant wind direction maps for every location and year in the study area           |
|windSpeedStack         |RasterStack |RasterStack of wind speed maps (km/h) for every location and year in the study area                 |

# Links to other modules

Mountain Pine Beetle Red Top Growth Model: Short-run Potential for Establishment, Eruption, and Spread.

- `mpbMassAttacksData`
- `mpbPine`
- `mpbRedTopSpread`

# References
