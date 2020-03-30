defineModule(sim, list(
  name = "mpbClimateData",
  description = "Mountain pine beetle climate sensitivity layers",
  keywords = c("insert key words here"),
  authors = c(person(c("Alex", "M"), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "mpbClimateData.Rmd"),
  reqdPkgs = list("amc", "grid", "magrittr", "quickPlot", "raster", "reproducible", "sp"),
  parameters = rbind(
    defineParameter("climateScenario", "character", "RCP45", NA_character_, NA_character_,
                    "The climate scenario to use. One of RCP45 or RCP85."),
    defineParameter("suitabilityIndex", "character", "R", NA_character_, NA_character_,
                    "The MPB climatic suitabilty index to use. One of S, L, R, or G."),
    defineParameter(".maxMemory", "numeric", 1e+9, NA, NA,
                    "Used to set the 'maxmemory' raster option. See '?rasterOptions'."),
    defineParameter(".plotInitialTime", "numeric", start(sim), 1981, 2100,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".tempdir", "character", tempdir(), NA, NA,
                    "Temporary (scratch) directory to use for transient files (e.g., GIS intermediates)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should this entire module be run with caching activated?")
  ),
  inputObjects = bind_rows(
    expectsInput("climateMapFiles", "character",
                 desc = "Vector of filenames correspoding to climate suitablity map layers",
                 sourceURL = "https://drive.google.com/file/d/1u4TpfkVonGk9FEw5ygShY1xiuk3FqKo3/view?usp=sharing"),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "if not supplied, will default to standAgeMap", # TODO: description needed
                 sourceURL = NA),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar"),
    expectsInput("studyArea", "SpatialPolygons",
                 desc = "The study area to which all maps will be cropped and reprojected.",
                 sourceURL = NA),
    expectsInput("studyAreaLarge", "SpatialPolygons",
                 desc = "The larger study area to use for spread parameter estimation.", ## TODO: better desc needed
                 sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput("climateMaps", "RasterStack", "Stack of climatic suitablity maps."),
    createsOutput("climateSuitabilityMap", "RasterLayer", "A climatic suitablity map for the current year.")
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

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "mpbClimateData", "switchLayer", .first())
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbClimateData", "plot", .last())
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbClimateData", "save", .last() + 1)
    },
    "plot" = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      names(sim$climateSuitabilityMap) <- "layer"
      Plot(sim$climateSuitabilityMap, title = "Climate Suitability Map", new = TRUE)
      Plot(sim$studyArea, addTo = "sim$climateSuitabilityMap", gp = gpar(col = "black", fill = 0))

      # schedule future event(s)

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "mpbClimateData", "plot")
      # ! ----- STOP EDITING ----- ! #
    },
    "switchLayer" = {
      sim <- switchLayer(sim)

      sim <- scheduleEvent(sim, time(sim) + 30, "mpbClimateData", "switchLayer") ## TODO: make this work with start times != 2011
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

switchLayer <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  sim$climateSuitabilityMap <- if (time(sim) <= 2010) {
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

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  ## load study area
  if (!suppliedElsewhere("studyArea")) {
    prj <- paste("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113",
                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
    sim$studyArea <- amc::loadStudyArea(dataPath(sim), "studyArea.kml", prj)
  }

  ## stand age map
  if (!suppliedElsewhere("standAgeMap", sim)) {
    standAgeMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.tif")
    sim$standAgeMap <- Cache(prepInputs,
                             targetFile = basename(standAgeMapFilename),
                             archive = asPath(c("kNN-StructureStandVolume.tar",
                                                "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.zip")),
                             destinationPath = dPath,
                             url = na.omit(extractURL("standAgeMap")),
                             fun = "raster::raster",
                             studyArea = sim$studyArea,
                             #rasterToMatch = sim$rasterToMatch,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = paste0(tools::file_path_sans_ext(basename(standAgeMapFilename)), "_cropped"),
                             overwrite = TRUE,
                             userTags = c("stable", currentModule(sim)))
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

  return(invisible(sim))
}

### helper functions
importMaps <- function(sim) {
  ## load the maps from files, crop, reproject, etc.
  files <- sim$climateMapFiles

  layerNames <- c("X1981.2010", "X2011.2040", "X2041.2070", "X2071.2100")

  out <- stack(files)
  #amc::cropReproj(out, studyArea, layerNames = layerNames, filename = amc::tf(".tif"))
  out <- Cache(postProcess,
               out,
               studyArea = sim$studyAreaLarge,
               filename2 = NULL,
               rasterToMatch = sim$rasterToMatch,
               overwrite = TRUE,
               useCache = TRUE) %>%
    stack()

  # ensure all cell values between 0 and 1
  out[out[] < 0.0] <- 0
  out[out[] > 1.0] <- 1

  sim$climateMaps <- setMinMax(out) %>% stack(.) %>% setNames(., layerNames)

  return(sim)
}
