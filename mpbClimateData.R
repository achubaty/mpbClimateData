## package BioSIM not available on CRAN nor GitHub; needs to be installed as follows:
## install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type="source")
## install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest", repos = NULL,  type="source")

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
  reqdPkgs = list("achubaty/amc@development", "BioSIM",
                  "grid",
                  "PredictiveEcology/LandR@LCC2010 (>= 1.0.4)",
                  "magrittr", "maptools",
                  "PredictiveEcology/mpbutils (>= 0.1.2)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/SpaDES.core@modifyList2 (>= 1.0.8.9001)",
                  "quickPlot", "raster", "PredictiveEcology/reproducible@development (>= 1.2.7.9002)", "sp", "spatialEco"),
  parameters = rbind(
    defineParameter("climateScenario", "character", "RCP45", NA_character_, NA_character_,
                    "The climate scenario to use. One of RCP45 or RCP85."),
    defineParameter("suitabilityIndex", "character", "R", NA_character_, NA_character_,
                    "The MPB climatic suitabilty index to use. One of 'S', 'L', 'R', or 'G'."),
    defineParameter("windMonths", "integer", 7, 1, 12,
                    paste("A vector of length 1 or more, indicating the month(s) of the year from which to extract wind data.",
                    "This should correspond to the months of dispersal")),
    defineParameter(".maxMemory", "numeric", 1e+9, NA, NA,
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
    expectsInput("climateMapFiles", "character",
                 desc = "Vector of filenames correspoding to climate suitablity map layers",
                 sourceURL = "https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing"),
    expectsInput("windDirStack", "RasterStack",
                 desc = "RasterStack of dominant wind direction maps for every location and year in the study area",
                 sourceURL = ""),
    expectsInput("windSpeedStack", "RasterStack",
                 desc = "RasterStack of wind speed maps (km/h) for every location and year in the study area",
                 sourceURL = ""),
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
    createsOutput("climateSuitabilityMaps", "RasterStack", "A time series of climatic suitablity RasterLayers, each with previx 'X' and the year, e.g., 'X2010'")
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
      sim <- importMaps(sim)
      # sim <- switchLayer(sim)

      # schedule future event(s)
      # sim <- scheduleEvent(sim, start(sim) + 1, "mpbClimateData", "switchLayer", .first())
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbClimateData", "plot", .last() - 1)
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbClimateData", "save", .last())
    },
    "plot" = {
      # do stuff for this event
      # names(sim$climateSuitabilityMap) <- "layer"
      Plot(sim$climateSuitabilityMaps, title = "Climate Suitability Maps", new = TRUE)
      Plot(sim$studyArea, addTo = "sim$climateSuitabilityMaps", gp = gpar(col = "black", fill = 0),
           title = "")

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "mpbClimateData", "plot")
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
  mod$prj <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

  ## load study area
  if (!suppliedElsewhere("studyArea")) {
    sim$studyArea <- mpbStudyArea(ecoregions = c(112, 120, 122, 124, 126), mod$prj,
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

  if (!suppliedElsewhere("windDirStack")) {
    # Make coarser
    aggRTM <- raster::raster(sim$rasterToMatch)
    aggRTM <- raster::aggregate(aggRTM, fact = 40)

    windModel <- try(Cache(getModelList)[17]) ## "ClimaticWind_Monthly"
    workingWindCacheId <- "5f588195a51652d2"
    if (is(windModel, "try-error")) {
      #library(googledrive);
      #driveDL <- Cache(googledrive::drive_download, as_id("16xEX2HVDTT2voLC5doRDEZ-WzNq_69EP"), overwrite = TRUE)
      #wind <- qs::qread(driveDL$local_path)
      windCacheId <- workingWindCacheId
    } else {
      windCacheId <- NULL
    }
    aggRTM <- aggregateRasByDT(sim$rasterToMatch, aggRTM, fn = mean)

    # Make Vector dataset
    cellsWData <- which(!is.na(aggRTM[]))
    cells <- xyFromCell(aggRTM, cell = cellsWData)
    sps <- sf::st_as_sf(SpatialPoints(cells, proj4string = crs(aggRTM)))
    sps <- sf::st_transform(sps, crs = 4326)#"+init=epsg:3857")#  "4326")
    locations <- data.table(Name = paste0("ID", 1:NROW(sps)), st_coordinates(sps))

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
    aggDEM <- aggregateRasByDT(DEM, aggRTM, mean)
    stWind <- system.time({
      ## 43 minutes with 3492 locations
      mess <- capture.output(type = "message", {
        wind <- Cache(getModelOutput, 2010, 2021, locations$Name,
                      locations$Y, locations$X, aggDEM[][cellsWData],
                      modelName = windModel,
                      rcp = "RCP85", climModel = "GCM4", useCloud = TRUE,
                      cloudFolderID = "175NUHoqppuXc2gIHZh5kznFi6tsigcOX", # Eliot's Gdrive: Hosted/BioSIM/ folder
                      cacheId = windCacheId)
        })
    })
    ignore <- lapply(mess, function(m) message(crayon::blue(gsub("^.+: mpbClm", "", m))))
    if (is.null(windCacheId)) {
      windAttr <- attr(wind, "tags")
      curCacheId <- if (!is.null(windAttr)) {
        gsub("cacheId:", "", windAttr)
      } else {
        gsub("^.+\\((.+)\\..+\\)", "\\1", grep("Object to retrieve", mess, value = TRUE))
      }
      if (curCacheId != workingWindCacheId) {
        message(crayon::red("You need to update the windCacheId; it is now ", curCacheId)    )
      }
    }

    # Make RasterStack
    setDT(wind)
    yrsChar <- paste0("X", unique(wind$Year))
    windStk <- stack(lapply(yrsChar, function(x) aggRTM))
    windStk[] <- NA
    windSpeedStk <- windStk

    mnths <- months(as.POSIXct("2021-01-15") + dmonth(1) * 0:11)
    whMonths <- P(sim)$windMonths
    message("Using only ", crayon::red(paste(mnths[whMonths], collapse = ", ")),
            " wind directions")

    if (any(colnames(wind) == "Month")) {
      wind <- wind[Month %in% whMonths]
    }

    # Wind Speed
    windSpeed <- wind[, list(WindSpeed = WindSpeed[1]), by = c("KeyID", "Year")]
    windSpeedWide <- dcast(windSpeed, formula = KeyID ~ Year)
    windSpeedWide <- windSpeedWide[match(unique(wind$KeyID), KeyID )]
    colnms <- grep("[[:digit:]]{4,4}", colnames(windSpeedWide), value = TRUE)
    windSpeedStk[cellsWData] <- as.matrix(windSpeedWide[, ..colnms])
    windSpeedStk <- raster::stack(windSpeedStk) # convert from Brick
    names(windSpeedStk) <- yrsChar

    # Wind Direction
    windCols <- grep("^W[[:digit:]]", colnames(wind), value = TRUE)
    cols <- c("KeyID", "Month", windCols)

    # Convert to single main direction -- sum of all vectors (i.e., magnitude and direction)
    # Means converting to x and y dimensions, then summing each of those, then reconverting back
    #  to angles
    angs0To360 <- (seq_along(colnames(wind[, ..windCols]))-1) * 10

    ang <- sumAngles(angs0To360, wind[, ..windCols])
    set(wind, NULL, "angleMean", ang)

    # For case with multiple months -- this will collapse to 1 datapoint per pixel per year
    wind <- wind[, list(Month = Month[1], Longitude = Longitude[1],
                            Latitude = Latitude[1], angleMean = mean(angleMean)),
                 by = c("Year", "KeyID")]

    # Wind direction
    windDirWide <- dcast(wind, formula = KeyID ~ Year)
    windDirWide <- windDirWide[match(unique(wind$KeyID), KeyID )]
    colnms <- grep("[[:digit:]]{4,4}", colnames(windDirWide), value = TRUE)
    windStk[cellsWData] <- as.matrix(windDirWide[, ..colnms])
    windStk <- raster::stack(windStk) # convert from Brick
    names(windStk) <- yrsChar

    # Visualize
    Plots(windStk, filename = "windDirStack")
    Plots(windSpeedStk, filename = "windSpeedStack")

    windDirStack <- Cache(disaggregate, windStk, fact = 40)
    sim$windDirStack <- raster::stack(crop(windDirStack, sim$rasterToMatch))

    windSpeedStack <- Cache(disaggregate, windSpeedStk, fact = 40)
    sim$windSpeedStack <- raster::stack(crop(windSpeedStack, sim$rasterToMatch))

    if (!compareRaster(sim$windDirStack, sim$windSpeedStack, sim$rasterToMatch, stopiffalse = FALSE)) {
      warning("wind raster is not same resolution as sim$rasterToMatch; please debug")
      browser() ## TODO: remove
    }
  }

  # library(ggspatial)
  # library(ggplot2)
  # ggplot() +
  #   geom_spatial(windStk)

  return(invisible(sim))
}

### helper functions
importMaps <- function(sim) {
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
  sim$climateSuitabilityMaps <- prepInputs(fun = "qs::qread", destinationPath = inputPath(sim),
                   url = "https://drive.google.com/file/d/1NIqSv0O5DQbm0nf2YZhfOpKqqws8I-vV/view?usp=sharing")

  return(sim)
}

aggregateRasByDT <- function(ras, newRas, fn = sum) {
  whNonNA <- which(!is.na(ras[]))
  rc2 <- rowColFromCell(ras, whNonNA)
  if (!all(((res(newRas)/res(ras)) %% 1) == 0))
    stop("The resolutions of the original raster and new raster are not integer multiples")
  disaggregateFactor <- unique(res(newRas)/res(ras))
  dt <- data.table(vals = ras[][whNonNA], ceiling(rc2 / disaggregateFactor))
  dt2 <- dt[, list(vals = fn(vals)), by = c("row", "col")]
  pixes <- cellFromRowCol(newRas, row = dt2$row, col = dt2$col)
  newRasOut <- raster(newRas)
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
