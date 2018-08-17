defineModule(sim, list(
  name = "mpbClimateData",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = c(person(c("Alex", "M"), "Chubaty", email = "alexander.chubaty@canada.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "mpbClimateData.Rmd"),
  reqdPkgs = list("magrittr", "quickPlot", "raster", "reproducible", "sp"),
  parameters = rbind(
    defineParameter("climateScenario", "character", "RCP45", NA_character_, NA_character_, "The climate scenario to use. One of RCP45 or RCP85."),
    defineParameter("suitabilityIndex", "character", "G", NA_character_, NA_character_, "The MPB climatic suitabilty index to use. One of S, L, R, or G."),
    defineParameter(".plotInitialTime", "numeric", start(sim), 1981, 2100, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "numeric", FALSE, NA, NA, "Should this entire module be run with caching activated?")
  ),
  inputObjects = bind_rows(
    expectsInput("studyArea", "SpatialPolygons", "The study area to which all maps will be cropped and reprojected.", sourceURL = NA)
  ),
  outputObjects = bind_rows(
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
      sim <- sim$mpbClimateDataImportMaps(sim)

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
      Plot(sim$studyArea, addTo = "sim$climateSuitabilityMap")

      # schedule future event(s)

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "mpbClimateData", "plot")
      # ! ----- STOP EDITING ----- ! #
    },
    "switchLayer" = {
      sim <- mpbClimateDataSwitchLayer(sim)

      sim <- scheduleEvent(sim, time(sim) + 40, "mpbClimateData", "switchLayer")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template for your event2
mpbClimateDataSwitchLayer <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  sim$climateSuitabilityMap <- if (time(sim) <= 2010) {
    sim$mpbClimateDataMaps[[1]]
  } else if (time(sim) <= 2040) {
    sim$mpbClimateDataMaps[[2]]
  } else if (time(sim) <= 2070) {
    sim$mpbClimateDataMaps[[3]]
  } else if (time(sim) <= 2100) {
    sim$mpbClimateDataMaps[[4]]
  } else {
    stop("No climate suitabliity projections available beyond year 2100.")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  if (!('studyArea' %in% sim$.userSuppliedObjNames)) {
    f <- file.path(modulePath(sim), "mpbClimateData", "data", "studyArea.kml")
    prj <- paste("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113",
                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
    sim$studyArea <- readOGR(f, "studyArea.kml") %>%
      sp::spTransform(., prj)
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### helper functions
mpbClimateDataImportMaps <- function(sim) {
  suffix <- switch(P(sim)$suitabilityIndex,
                   "S" = "_SafP[.]tif",
                   "L" = "_LoganP[.]tif",
                   "R" = "_RegP[.]tif",
                   "G" = "_GeoP[.]tif",
                   stop("suitability index must be one of S, L, R, or G."))
  files <- dir(path = file.path(modulePath(sim), "mpbClimateData", "data"),
               pattern = suffix, full.names = TRUE)
  files <- c(files[1], grep(P(sim)$climateScenario, files, value = TRUE))

  if (length(files) == 0) {
    stop("mpbClimateData: missing data files")
  }

  fn1 <- function(files, studyArea) {
    layerNames <- c("X1981.2010", "X2011.2040", "X2041.2070", "X2071.2100")
    out <- stack(files) %>%
      amc::cropReproj(., studyArea, layerNames = layerNames, filename = amc::tf(".tif"))

    # ensure all cell values between 0 and 1
    out[out[] < 0.0] <- 0
    out[out[] > 1.0] <- 1
    out <- setMinMax(out) %>% stack(.) %>% set_names(., layerNames)
    out
  }
  sim$mpbClimateDataMaps <- Cache(fn1, files, sim$studyArea)

  return(sim)
}
