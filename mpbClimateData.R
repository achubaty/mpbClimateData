## package BioSIM not available on CRAN nor GitHub; needs to be installed as follows:
## install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type="source")
## install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest", repos = NULL,  type="source")

## remotes::install_github("CWFC-CCFB/J4R@v1.1.8") ## v1.1.9 broken on ubuntu
## remotes::install_github("RNCan/BioSimClient_R")

defineModule(sim, list(
  name = "mpbClimateData",
  description = "Mountain pine beetle climate sensitivity, and wind layers",
  keywords = c("BioSIM", "ClimaticWind_Monthly", "MPB_SLR"),
  authors = c(
    person(c("Alex", "M"), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre")),
    person(c("Eliot", "JB"), "McIntire", email = "eliot.mcintire@canada.ca", role = "aut")
  ),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "mpbClimateData.Rmd"),
  reqdPkgs = list("achubaty/amc@development",
                  "BioSIM", ## RNCan/BioSimClient_R; needs J4R v1.1.8 (v1.1.9 is broken on ubuntu)
                  "grid", "CircStats",
                  "PredictiveEcology/LandR@development (>= 1.0.4)",
                  "magrittr", # "maptools",
                  "PredictiveEcology/mpbutils (>= 0.1.2)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/SpaDES.core@development (>= 1.0.8.9002)",
                  "quickPlot", "raster", "terra",
                  "PredictiveEcology/reproducible@development (>= 1.2.7.9002)",
                  "sp", "sf", "data.table", "spatialEco"),
  parameters = rbind(
    defineParameter("climateScenario", "character", "RCP45", NA_character_, NA_character_,
                    "The climate scenario to use. One of RCP45 or RCP85."),
    defineParameter("cloudCacheFolderID", "character",
                    "1225VXvtkHdkhLgQbl_CUd-6aVYOy9OD7", ## effective 2023-04-01
                    #"175NUHoqppuXc2gIHZh5kznFi6tsigcOX", ## Eliot's old cache dir
                    NA_character_, NA_character_,
                    "Google Drive folder ID where cached BioSim windmaps can be retrieved."),
    defineParameter("cloudCacheFileIDs", "character",
                    # c("28c428e741e18c6f", "7c98cfb38d1b0783", "75b95f32359cbe0c", "ed8bff9a39b67a94"), ## ABSK (old)
                    c("08cbf5d4fa5d122a", "8c6f0c4965370cf3", "b74daa01256a65a2", "5c958db389b8ff6f"), ## ABSK_9.2
                    ## TODO: need better mechanism to deal with studyArea changes !!
                    NA_character_, NA_character_,
                    "BioSim windmap cache ids to retrieve from `cloudCacheFolderID`."),
    defineParameter("suitabilityIndex", "character", "R", NA_character_, NA_character_,
                    "The MPB climatic suitabilty index to use. One of 'S', 'L', 'R', or 'G'."),
    defineParameter("usePrerun", "logical", TRUE, NA, NA,
                    "Should cached or pre-run calls to BioSIM functions be used?"),
    defineParameter("windMonths", "integer", 7L, 1L, 12L,
                    paste("A vector of length 1 or more, indicating the month(s) of the year from which to extract wind data.",
                          "This should correspond to the months of dispersal")),
    defineParameter(".maxMemory", "numeric", 6e+10, NA, NA,
                    "Used to set the 'maxmemory' raster option. See '?rasterOptions'."),
    defineParameter(".plots", "character", "screen", NA_character_, NA_character_,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInitialTime", "numeric", start(sim), 1981, 2100,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".tempdir", "character", NULL, NA, NA,
                    "Temporary (scratch) directory to use for transient files (e.g., GIS intermediates)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should this entire module be run with caching activated?")
  ),
  inputObjects = bindrows(
    expectsInput("climateMapRandomize", c("list", "logical"),
                 "List with 2 elements: srcYears, rmYears. Each of these should be a numeric vector of years
                 (sim$climateSuitabilityMaps must have layer names that contain the full 4-digit numeric year, e.g., X2010),
                 with 'scrYears' being the pool of years to take data from to use in 'rmYears'. This is
                 intended for cases where future projected climate data is unavailable or corrupt. If TRUE, then it will
                 currently select 2021 as the final year of 'data' and every year beyond that will be randomized",
                 sourceURL = "https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing"),
    expectsInput("climateMapFiles", "character",
                 desc = "Vector of filenames correspoding to climate suitablity map layers",
                 sourceURL = "https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing"),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "if not supplied, will default to standAgeMap",
                 sourceURL = NA),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("studyArea", "SpatialPolygons",
                 desc = "The study area to which all maps will be cropped and reprojected.",
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    # createsOutput("climateMaps", "RasterStack", "Stack of climatic suitablity maps."),
    createsOutput("climateSuitabilityMaps", "RasterStack", "A time series of climatic suitablity RasterLayers, each with previx 'X' and the year, e.g., 'X2010'"),
    createsOutput("windDirStack", "RasterStack",
                  desc = "RasterStack of dominant wind direction maps for every location and year in the study area"),
    createsOutput("windSpeedStack", "RasterStack",
                  desc = "RasterStack of wind speed maps (km/h) for every location and year in the study area"),
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.mpbClimateData <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         "init" = {
           ### check sim init params etc.
           stopifnot(start(sim) > 1981, end(sim) < 2100)

           # do stuff for this event
           sim <- GenerateClimate(sim)
           # sim <- switchLayer(sim)

           # schedule future event(s)
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbClimateData", "save", .last())
         },
         "switchLayer" = {
           sim <- switchLayer(sim)

           sim <- scheduleEvent(sim, time(sim), "mpbClimateData", "switchLayer", .first())
         },
         warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                       "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

GenerateClimate <- function(sim) {
  # CHECK MEMORY -- THIS WORKS BETTER WITH A LOT
  minMemNeeded <- 6e10
  for (i in 1:2) {
    coToIgnore <- capture.output({
      ro <- rasterOptions()
    })
    if (Par$.maxMemory >= minMemNeeded) {
      if (ro$maxmemory < Par$.maxMemory) {
        message("Increasing rasterOptions maxmemory to P(sim)$.maxMemory, which is ", Par$.maxMemory)
        rasterOptions(maxmemory = Par$.maxMemory)
      }
    } else {
      if (ro$maxmemory < minMemNeeded && pemisc::availableMemory() > 2e11) {
        message("Your machine has plenty of memory yet rasterOptions is set to ", ro$maxmemory)
        if (Par$.maxMemory > ro$maxmemory) {
          message(". Increasing this to ", Par$.maxMemory)
          rasterOptions(maxmemory = Par$.maxMemory)
        } else {
          message("Consider increasing this with e.g., set .maxMemory parameter in this module to ", minMemNeeded,
                  " or rasterOptions(maxmemory = ", minMemNeeded,")")
        }
      } else {
        break
      }
    }
  }

  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  ## load the maps from files, crop, reproject, etc.
  # files <- sim$climateMapFiles
  #
  # layerNames <- c("X1981.2010", "X2011.2040", "X2041.2070", "X2071.2100")
  #
  # rawClimateMaps <- raster::stack(files)
  # sim$climateMaps <- Cache(postProcess,
  #                          rawClimateMaps,
  #                          studyArea = sim$studyAreaLarge,
  #                          filename2 = NULL,
  #                          rasterToMatch = sim$rasterToMatch,
  #                          overwrite = TRUE,
  #                          useCache = TRUE) %>%
  #   raster::stack() %>%
  #   setNames(layerNames) %>%
  #   pemisc::normalizeStack()

  # sim$climateSuitabilityMaps <- prepInputs(
  #   url = "https://drive.google.com/file/d/1NIqSv0O5DQbm0nf2YZhfOpKqqws8I-vV/",
  #   destinationPath = inputPath(sim),
  #   fun = "qs::qread"
  # ) |> rast() ## TODO: do MPB_SLR call direcly:
  ## LandR::BioSIM_getMPBSLR(dem, years = years, SLR = "R", climModel = "GCM4", rcp = "RCP45")

  ## Make coarser template
  resClimate <- rep(1e4, 2)
  message("Using Climate resolution (wind, climate suitability) of ", paste(resClimate, collapse = " x "), " m")

  RTMaggregated <- terra::rast(sim$rasterToMatch)
  fact <- sqrt(prod(resClimate)/prod(res(sim$rasterToMatch)))

  fact <- fact * 3
  RTMaggregated <- terra::aggregate(RTMaggregated, fact = fact)

  windModel <- if (!grepl("spades", Sys.info()["nodename"])) {
    for (iii in 1:2) {
      gml <- try(getModelList())
      if (is(gml, "try-error"))
        BioSIM::shutdownClient()
      else
        break
    }

    if (is(gml, "try-error"))
      "ClimaticWind_Monthly"
    else {
      whModel <- grep("ClimaticWind_Monthly", gml)
      gml[whModel] ## "ClimaticWind_Monthly"
    }
  } else {
    message("Not attempting BioSIM because using BorealCloud where it doesn't work")
    try(stop(), silent = TRUE)
  }

  RTMaggregated <- reproducible:::suppressWarningsSpecific(falseWarnings = "raster has no values",
                                                           Cache(aggregateRasByDT,
                                                                 sim$rasterToMatch, RTMaggregated, fn = mean)
  )

  # Make Vector dataset
  cellsWData <- which(!is.na(RTMaggregated[]))
  cells <- terra::xyFromCell(RTMaggregated, cell = cellsWData)
  sps <- sf::st_as_sf(terra::vect(cells, crs = crs(RTMaggregated)))
  sps <- sf::st_transform(sps, crs = 4326)
  locations <- data.table(Name = paste0("ID", 1:NROW(sps)), st_coordinates(sps),
                          cellsWData = cellsWData)

  # Do call to BioSIM
  # Until this gets fixed in J4R; this is the fix
  # stWind <- system.time(
  #   wind <- Cache(getModelOutput, 2010, 2021, locations$Name,
  #                          locations$Y, locations$X, rep(1000, NROW(sps)),
  #                          modelName = windModel,
  #                          rcp = "RCP85", climModel = "GCM4"))
  #
  # windModel <- Cache(getModelList)[17]
  DEM <- Cache(LandR::prepInputsCanDEM,
               rasterToMatch = sim$rasterToMatch,
               studyArea = sim$studyArea,
               destinationPath = dPath)
  aggDEM <- Cache(aggregateRasByDT, DEM, RTMaggregated, mean)

  locations$elevation <- aggDEM[][locations$cellsWData]
  # aggDEM[is.na(aggDEM[])] <- -9999

  splitInd <- ceiling(1:NROW(locations) / 1000)
  locationsList <- split(locations, splitInd)
  windCacheFolderID <-  P(sim)$cloudCacheFolderID
  windCacheIds <- if (isTRUE(P(sim)$usePrerun)) {
    P(sim)$cloudCacheFileIDs
  } else {
    lapply(locationsList, function(x) NULL)
  }

  ## TODO: use CMIP6 climate models that match LandR-fS sims

  modelsThatNeedGT1Yr <- c("BudBurst", "Climate_Mosture_Index_Annual", "EmeraldAshBorerColdHardiness_Annual",
                           "Gypsy_Moth_Seasonality", "HemlockWoollyAdelgid_Annual", "HemlockWoollyAdelgid_Daily",
                           "MPB_Cold_Tolerance_Annual", "MPB_Cold_Tolerance_Daily", "MPB_SLR",
                           "PlantHardinessCanada", "PlantHardinessUSA", "Spruce_Budworm_Biology_Annual",
                           "SpruceBeetle", "Standardised_Precipitation_Evapotranspiration_Index"
  )

  climateSuitabilityModel <- "MPB_SLR" # "TminTairTmax_Daily"
  modelNames <- c(windModel, climateSuitabilityModel)


  if (FALSE) { # used to simplify for dput to send to Mathieu
    for(col in c("X", "Y", "elevation")) {
      set(location, NULL, col, round(location[[col]], 3))
    }
  }

  ## BioSIM call
  ## 43 minutes with 3492 locations
  BioSim <- Map(modelName = c(modelNames),#, windModel),
                function(modelName) {
                  outer <- Map(location = locationsList, cacheId = windCacheIds,
                               function(location, cacheId) {
                                 ## uses cloud cache to retrieve wind maps
                                 strtYrs <- start(sim)
                                 endYrs <- end(sim)
                                 if (modelName %in% modelsThatNeedGT1Yr) {
                                   nyrs <- 1
                                   if (!isTRUE(all.equal(strtYrs, endYrs))) {
                                     strtYrs <- strtYrs - nyrs
                                   }
                                   strtYrs <- seq(strtYrs, endYrs - nyrs)
                                   endYrs <- seq(strtYrs[1] + nyrs, endYrs)
                                 }
                                 inner <- Map(strtYr = strtYrs, endYr = endYrs, function(strtYr, endYr) {
                                   message(modelName, ": yrs: ", strtYr, "-", endYr)
                                   intermediateInner <<-
                                     Cache(generateWeather, fromYr = strtYr, toYr = endYr,
                                           location$Name,
                                           latDeg = location$Y, longDeg = location$X,
                                           elevM = location$elevation,
                                           modelNames = modelName,
                                           rcp = "RCP85", climModel = "GCM4",
                                           #useCloud = TRUE,
                                           #cloudFolderID = windCacheFolderID,
                                           cacheId = cacheId)
                                 })

                                 intermediateInner2 <<- inner2 <-
                                   rbindlist(lapply(inner, function(x) x[[1]])) # only 1 modelName used here

                               })
                  out <- rbindlist(outer)
                })

  nams <- outputObjects(sim)$objectName
  # put correct order
  nams <- nams[pmatch(c("climate", "windDir", "windSpeed"), nams)]

  ## Rasterize BioSIM outputs
  # Climate Suitability
  sim[[nams[1]]] <- rasterizeClimateSuitability(BioSim[[climateSuitabilityModel]], RTMaggregated, cellsWData)
  # Wind Direction
  sim[[nams[2]]] <- rasterizeWindDir(BioSim[[windModel]], RTMaggregated, cellsWData)
  # WindSpeed
  sim[[nams[3]]] <- rasterizeWindSpeed(BioSim[[windModel]], RTMaggregated, cellsWData, whMonths = Par$windMonths)

  ## pre-digest the layers (for multiple Cache calls below)
  digs <- list()
  for (nam in nams) {
    digs[[nam]] <- .robustDigest(list(sim[[nam]], fact))
  }

  ## Visualize the small ones, including incorrect for posterity sake if they are in the future
  if (anyPlotting(Par$.plots)) {
    stNoColons <- gsub(":", "-", format(Sys.time()))

    for (nam in nams) {
      titl <- paste0("Small ", nam, ", pre-randomization")
      fn <- paste0(titl, ", ", start(sim), " to ", end(sim), "_", stNoColons)
      Plots(sim[[nam]], title = titl, filename = fn, maxnl = terra::nlyr(sim[[nam]]),
            plg = list(cex = 2), deviceArgs = list(width = 1500, height = 1200)) |>
        Cache(omitArgs = c("filename", "data"), .cacheExtra = digs[[nam]])
    }
  }

  ## Make resolution match sim$rasterToMatch
  for (nam in nams) {
    sim[[nam]] <-
      disaggregateToStack(sim[[nam]], sim$rasterToMatch, fact = fact) |>
      Cache(.cacheExtra = digs[[nam]], omitArgs = c("x", "y", "fact"))
  }

  ## GIS check
  if (!do.call(compareGeom, unname(mget(nams, envir = envir(sim))))) {
    warning(paste("Wind or climate suitability rasters is not same resolution as rasterToMatch.\n",
                  "This may be due to retrieving cached wind maps for another studyArea.\n",
                  "Please check that you have passed the correct cloudCacheFileIDs parameter values."))
  }

  ## Randomize the "future" years with past years
  if (end(sim) > nowYear()) {
    if (!is.null(sim$climateMapRandomize)) {
      if (isTRUE(sim$climateMapRandomize)) {
        sim$climateMapRandomize <- list(srcYears = start(sim):nowYear(), rmYears = (nowYear() + 1):(end(sim) + 1))

      }
      for (nam in nams) {
        sim[[nam]] <-
          randomizeSomeStackLayers(sim[[nam]], sim$climateMapRandomize,
                                   endTime = end(sim) + 1, inputPath = inputPath(sim)) |>
          Cache(omitArgs = c("stk", "swaps", "endTime", "inputPath"), .cacheExtra = digs[[nam]])
        digs[[nam]] <- .robustDigest(list(sim[[nam]], fact))
      }
    }
  }
  message(crayon::green(mean(sim$climateSuitabilityMaps[[nlyr(sim$climateSuitabilityMaps)]][], na.rm = TRUE)))

  ## Visualize
  if (anyPlotting(Par$.plots)) {
    # Visualize the small ones, including incorrect for posterity sake if they are in the future
    stNoColons <- gsub(":", "-", format(Sys.time()))

    for (nam in nams) {
      titl <- paste0("Full ", nam, ", randomized")
      fn <- paste0(titl, ", ", start(sim), " to ", end(sim), "_", stNoColons)
      Plots(sim[[nam]], title = titl, new = TRUE, filename = fn, maxnl = terra::nlyr(sim[[nam]]),
            plg = list(cex = 2), deviceArgs = list(width = 1500, height = 1200)) |>
        Cache(
          omitArgs = c("filename", "data", "title"),
          .cacheExtra = digs[[nam]])
    }
  }

  ## write to disk; they are enormous
  for (nam in nams) {
    out <-
      writeTo(sim[[nam]],
              writeTo = paste0(figurePath(sim),
                               filenameGenerator(nam, start(sim), end(sim), stNoColons = stNoColons),
                               ".tif")) |>
      Cache(omitArgs = c("from", "writeTo"), .cacheExtra = digs[[nam]])
    sim[[nam]] <- out

  }

  return(sim)
}

switchLayer <- function(sim) {
  stop("switchLayer is not longer functional because bypassing 4-time-sequence climateMap --> direct sim$climateSuitabilityMaps")
  sim$climateSuitabilityMap <- if (time(sim) <= 2010) {  ## TODO: update this to get annual rasters
    sim$climateMaps[[1]]
  } else if (time(sim) <= 2040) {
    sim$climateMaps[[2]]
  } else if (time(sim) <= 2070) {
    sim$climateMaps[[3]]
  } else if (time(sim) <= 2100) {
    sim$climateMaps[[4]]
  } else {
    stop("No climate suitabliity projections available beyond year 2100.")
  }

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  cPath <- cachePath(sim)
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  #mod$prj <- paste("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113",
  #                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  # mod$prj <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
  #                  "+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  prj <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

  ## load study area
  if (!suppliedElsewhere("studyArea")) {
    sim$studyArea <- mpbStudyArea(ecoregions = c(112, 120, 122, 124, 126), prj,
                                  cPath, dPath)
  }

  ## stand age map
  if (!suppliedElsewhere("standAgeMap", sim)) {
    sim$standAgeMap <- LandR::prepInputsStandAgeMap(
      startTime = 2010,
      ageUrl = na.omit(extractURL("standAgeMap")),
      destinationPath = dPath,
      studyArea = sim$studyArea,
      userTags = c("stable", currentModule(sim))
    )
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  ## raster to match
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- sim$standAgeMap
  }

  ## download the climate map files
  if (!suppliedElsewhere("climateMapFiles", sim)) {
    # only download the files; we will load them later during init event
    fileInfo <- Cache(preProcess,
                      targetFile = NULL,
                      archive = asPath(paste0(currentModule(sim), ".zip")),
                      destinationPath = dPath,
                      url = extractURL("climateMapFiles"),
                      fun = "raster::raster",
                      studyArea = sim$studyAreaLarge,
                      rasterToMatch = sim$rasterToMatch,
                      method = "bilinear",
                      datatype = "FLT4S",
                      filename2 = NULL,
                      overwrite = TRUE,
                      userTags = c("stable", currentModule(sim)))

    suffix <- switch(P(sim)$suitabilityIndex,
                     "S" = "_SafP[.]tif",
                     "L" = "_LoganP[.]tif",
                     "R" = "_RegP[.]tif",
                     "G" = "_GeoP[.]tif",
                     stop("suitability index must be one of S, L, R, or G."))

    files <- dir(path = dPath, pattern = suffix, full.names = TRUE)
    files <- c(files[1], grep(P(sim)$climateScenario, files, value = TRUE))

    if (length(files) == 0) {
      stop("mpbClimateData: missing data files")
    }

    sim$climateMapFiles <- files
  }

  if (!suppliedElsewhere("climateMapRandomize", sim)) {

    if (end(sim) > nowYear())
      sim$climateMapRandomize <- TRUE

  }


  # library(ggspatial)
  # library(ggplot2)
  # ggplot() +
  #   geom_spatial(windStk)

  return(invisible(sim))
}

### helper functions
filenameGenerator <- function(titl, start, end, stNoColons = "", sep = "_") {
  paste0(paste(titl, start, end, sep = sep), stNoColons)
}

digNameGenerator <- function(nam) {
  paste0("dig", nam)
}

aggregateRasByDT <- function(ras, newRas, fn = sum) {
  whNonNA <- which(!is.na(ras[]))
  rc2 <- rowColFromCell(ras, whNonNA)
  if (!all(((res(newRas)/res(ras)) %% 1) == 0))
    stop("The resolutions of the original raster and new raster are not integer multiples")
  disaggregateFactor <- unique(res(newRas)/res(ras))
  dt <- data.table(vals = ras[][whNonNA], ceiling(rc2 / disaggregateFactor))
  setnames(dt, old = c("V1", "V2"), new = c("row", "col"))
  dt2 <- dt[, list(vals = fn(vals)), by = c("row", "col")]
  pixes <- cellFromRowCol(newRas, row = dt2$row, col = dt2$col)
  newRasOut <- terra::rast(newRas)
  newRasOut[pixes] <- dt2$vals
  names(newRasOut) <- names(ras)
  newRasOut
}

sumAngles <- function(angles, magnitude) {
  dimsM <- dim(magnitude)
  dimsA <- dim(angles)
  if (!identical(dimsM, dimsA)) {
    if (!length(angles) %in% dimsM) {
      stop("angles must have same number of rows as magnitude matrix")
    }
    a <- rep(angles, NROW(magnitude))
    angles <- matrix(a, nrow = NROW(magnitude), ncol = length(angles), byrow = TRUE)
  }

  x <- round(magnitude * sin(rad(angles)), 5)
  y <- round(magnitude * cos(rad(angles)), 5)
  x <- apply(x, 1, mean)
  y <- apply(y, 1, mean)
  deg(atan2(x, y)) %% 360
}

randomizeSomeStackLayers <- function(stk, swaps, endTime, inputPath) {
  rmYears <- paste0("X", swaps$rmYears)
  if (endTime >= min(swaps$rmYears)) {
    layerNames <- as.numeric(gsub("X", "", layerNames(stk)))
    keepLayers <- endTime >= layerNames
    layerNames <- layerNames[keepLayers]
    srcYears <- paste0("X", swaps$srcYears)
    lenSrcYears <- length(srcYears)
    numNeedYears <- grep(yrNamesPlus1(endTime - 1), rmYears)
    if (length(numNeedYears) == 0) stop("Need end(sim) to be later")
    if (numNeedYears > 0) {
      rmYears <- rmYears[seq(numNeedYears)]
    }

    mapName <- deparse(substitute(stk))
    message("randomizing ", mapName, "; removing: ", paste(rmYears, collapse = ", "))
    newYears <- sample(srcYears, size = numNeedYears, replace = TRUE)
    aa <- stk[[srcYears]]
    bb <- stk[[newYears]]
    names(bb) <- rmYears
    stk <- c(aa, bb) # in terra, these are 2 different objects, but they have the same .tif
                     #   when we put them together, the filename doesn't contain the information
                     #   needed to keep the aa and bb separate, so when wrapSpatRaster is run,
    #   it essentially does c(stk, stk); so, need to rebuild with 1 file backend
    if (any(nzchar(Filenames(stk)))) {
      tf <- file.path(terraOptions(print = FALSE)[["tempdir"]], basename(tempfile(fileext = ".tif")))
      stk1 <- writeRaster(stk, filename = tf)
      fn <- unique(Filenames(stk))
      rm(stk)
      unlink(fn)
      out <- file.copy(Filenames(stk1), fn)
      stk <- terra::rast(fn)
    }
  }
  return(stk)
}

yrNamesPlus1 <- function(yrNames) {
  yrNames <- gsub("X([[:digit:]]{4,4})", "\\1", yrNames)
  paste0("X", as.integer(yrNames) + 1)
}

disaggregateToStack <- function(x, y, fact) {
  x <- terra::disagg(x, fact)
  rr <- terra::crop(x, y)
  # raster::stack(rr)
  rr
}

nowYear <- function() {
  as.numeric(format(Sys.Date(), "%Y")) - 1
}

rasterizeClimateSuitability <- function(BioSimOutput, RTMaggregated, cellsWData) {
  yrsChar <- paste0("X", unique(BioSimOutput$Year))
  climateSuitabilityMaps <- terra::rast(lapply(yrsChar, function(x) RTMaggregated))
  names(climateSuitabilityMaps) <- yrsChar
  # climateSuitabilityMaps[] <-
  m <- matrix(BioSimOutput$CT_Survival,
              ncol = length(names(climateSuitabilityMaps)))
  climateSuitabilityMaps[cellsWData] <- m
  climateSuitabilityMaps
}

rasterizeWindSpeed <- function(wind, RTMaggregated, cellsWData, whMonths = 7) {
  stk <- rastTemplateGenerate(RTMaggregated, Years = wind$Year)

  mnths <- months(as.POSIXct("2021-01-15") + dmonth(1) * 0:11)
  # whMonths <- windMonths
  message("Using only ", crayon::red(paste(mnths[whMonths], collapse = ", ")), " wind directions")

  # Take only the month specified in whMonths
  if (any(colnames(wind) == "Month")) {
    wind <- wind[Month %in% whMonths]
  }

  # Wind Speed
  windSpeed <- wind[, list(WindSpeed = WindSpeed[1]), by = c("KeyID", "Year")]
  windSpeedWide <- dcast(windSpeed, formula = KeyID ~ Year, value.var = "WindSpeed")
  windSpeedWide <- windSpeedWide[match(unique(wind$KeyID), KeyID )]
  colnms <- grep("[[:digit:]]{4,4}", colnames(windSpeedWide), value = TRUE)
  stk[cellsWData] <- as.matrix(windSpeedWide[, ..colnms])
  # stk <- raster::stack(stk) # convert from Brick
  stk
}

rasterizeWindDir <- function(wind, RTMaggregated, cellsWData) {
  stk <- rastTemplateGenerate(RTMaggregated, Years = wind$Year)
  yrsChar <- yearsAsChar(wind$Year)
  windCols <- grep("^W[[:digit:]]", colnames(wind), value = TRUE)

  # There are cases where one or more columns is character
  areAllNumerics <- which(!sapply(wind[, ..windCols], function(x) all(is(x, "numeric"))))
  if (length(areAllNumerics) > 0) {
    chColNames <- names(areAllNumerics)
    wind[, (chColNames) := lapply(.SD, function(x) as.numeric(x)), .SDcols = chColNames]
  }

  cols <- c("KeyID", "Month", windCols)

  # Convert to single main direction -- sum of all vectors (i.e., magnitude and direction)
  # Means converting to x and y dimensions, then summing each of those, then reconverting back
  #  to angles
  angs0To360 <- (seq_along(colnames(wind[, ..windCols])) - 1) * 10

  ang <- sumAngles(angs0To360, wind[, ..windCols])
  set(wind, NULL, "angleMean", ang)

  # For case with multiple months -- this will collapse to 1 datapoint per pixel per year
  wind <- wind[, list(Month = Month[1], Longitude = Longitude[1],
                      Latitude = Latitude[1], angleMean = mean(angleMean)),
               by = c("Year", "KeyID")]

  # Wind direction
  windDirWide <- dcast(wind, formula = KeyID ~ Year, value.var = "angleMean")
  windDirWide <- windDirWide[match(unique(wind$KeyID), KeyID )]
  colnms <- grep("[[:digit:]]{4,4}", colnames(windDirWide), value = TRUE)
  stk[cellsWData] <- as.matrix(windDirWide[, ..colnms]) ## WARNING : In x@data@values[i] <- value : number of items to replace is not a multiple of replacement length

  # windStk <- raster::stack(windStk) # convert from Brick
  # names(windStk) <- yrsChar
  stk
}

yearsAsChar <- function(Year) {
  paste0("X", unique(Year))
}

rastTemplateGenerate <- function(RTMaggregated, Years) {
  yrsChar <- yearsAsChar(Years)
  yrsChar <- paste0("X", unique(Years))
  stk <- terra::rast(lapply(yrsChar, function(x) RTMaggregated))
  names(stk) <- yrsChar
  stk[] <- NA
  stk
}
