library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

source("optimSDP/optimSDP.R", echo=FALSE)
sourceCpp('optimSDP/OptimSPD.cpp', showOutput = FALSE)

library(sf); sf::sf_use_s2(FALSE)
library(geosphere)
library(tidyverse)

load("output/Map/eaafMap.rda")
load("output/breedTab_revision.rda")
load("output/mudflatTab.rda")
load("output/tempTab_revision.rda")
load("output/empTrackList.rda")
empTrackList <- empTrackList[c(1,2,3,4,6)]

## species
spParms <- setNames(lapply(c(250, 144, 100, 55, 25), sizeParams), 
                    c("Godwit", "GreatKnot", "RedKnot", "CurlewSandpiper", "RedNeckedStint"))
breedTab <- breedTab %>% filter(species%in%names(spParms)) %>% st_transform(4326)
spCols   <- c("darkblue", "cadetblue", "chartreuse4", "brown3", "darkgoldenrod2", "yellow2")

allSpSim <- lapply(names(spParms), function(sp) {
  
  subBreedTab <- breedTab %>% filter(species == sp)
  
  indSim <- parallel::mclapply(1:nrow(subBreedTab), function(ind) {
    
    #######################
    ### site parameters ###
    #######################
    StartEnd_cell <- tibble(lon = c(subBreedTab[ind,] %>% pull("lon_start"), st_coordinates(subBreedTab)[ind,1]),
                            lat = c(subBreedTab[ind,] %>% pull("lat_start"), st_coordinates(subBreedTab)[ind,2])) %>%
      st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% rowwise() %>%
      mutate(index = which.min(st_distance(geometry, mudflatTab %>% st_transform(4326))))
    
    size     <- (mudflatTab %>% st_centroid() %>% st_transform(4326) %>%
                   mutate(areaHist = ifelse(is.na(histArea_inner), currArea_outer, histArea_outer),
                          areaCurr = currArea_outer,
                          areaMang = mangArea_outer,
                          lake     = ifelse(is.na(lake_area), FALSE, TRUE),
                          pmax     = pmax(areaCurr, areaHist)) %>% 
                   dplyr::select(pmax, areaHist, areaCurr, areaMang, lake) %>%
                   rownames_to_column(var = "index") %>% mutate(index = as.integer(index))) %>%
      filter((areaHist>0 | areaCurr>0 | areaMang>0 | lake)) %>%
      filter(st_coordinates(.)[,1] > 104 | st_coordinates(.)[,1] < -150) %>% 
      arrange(st_distance(geometry, StartEnd_cell[1,])) %>%
      bind_rows(StartEnd_cell[2,] %>% mutate(areaHist = NA, areaCurr = NA, areaMang = NA, lake = FALSE) %>% dplyr::select(names(.))) %>%
      relocate(geometry, .after = last_col()) %>%
      suppressWarnings()
    
    distM  <- st_distance(size, by_element = F)/1000
    bearM  <- abs(distm(st_coordinates(size), fun = bearing))
    
    ### Revision
    reward <- c(subBreedTab[ind,] %>% pull(start_hist), subBreedTab[ind,] %>% pull(start_curr),  subBreedTab[ind,] %>% pull(start_future))
    
    #######################
    ### sdp Objects #######
    #######################
    sdpObjects <- makeSDPobjects(
      list(
        species = subBreedTab$species[ind],
        spParms = spParms[[which(names(spParms)==sp)]],
        minT    = as.numeric(format(as.POSIXct("2012-01-01"), "%j")),
        maxT    = subBreedTab[ind,] %>% pull(start_hist) + 15,
        MaxX    = 100,
        B0      = 3,
        w       = 0.028,
        xc      = 10,
        WindAssist = 0,
        WindProb   = 1,
        ZStdNorm   = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
        PStdNorm   = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092),
        reward     = reward,
        decError   = 250000,
        pred  = c(1e-3, 1e-4, 1e-6, 2, 2),
        sites = size,
        latF  = 1,
        dist  = distM,
        bear  = bearM,
        tTab  = tempTab[size$index,,] - 2.5 ## Wind compensation
      ))
    
    parallel::mclapply(1:3, function(x) {
      model   <- bwdIteration(sdpObjects[[x]])
      simu    <- tryCatch(fwdSimulation(model, 100, start_t = 1, start_site = 1, start_x = c(30,50)), error = function(e) NULL)
      simNetwork(simu, model, crds_ind = mudflatTab %>% st_centroid() %>% st_coordinates() %>% suppressWarnings(), plot = F)
    }, mc.cores = 3)
    
  }, mc.cores = 3)
  
  list(
    Reduce("+", lapply(indSim, function(ind) ind[[1]][[1]])[sapply(indSim, function(ind) class(ind[[1]][[1]])[1])=="matrix"]),
    Reduce("+", lapply(indSim, function(ind) ind[[2]][[1]])[sapply(indSim, function(ind) class(ind[[2]][[1]])[1])=="matrix"]),
    Reduce("+", lapply(indSim, function(ind) ind[[3]][[1]])[sapply(indSim, function(ind) class(ind[[3]][[1]])[1])=="matrix"]),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[1]][[2]], error = function(e) NULL))),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[2]][[2]], error = function(e) NULL))),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[3]][[2]], error = function(e) NULL))),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[1]][[3]], error = function(e) NULL))),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[2]][[3]], error = function(e) NULL))),
    do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[3]][[3]], error = function(e) NULL))),
    do.call("rbind", lapply(1:length(indSim), function(ind) tryCatch(as_tibble(indSim[[ind]][[1]][[4]]) %>% setNames(c("id", "time", "site", "x")) %>% mutate(id = glue::glue("{ind}_{id}")), error = function(e) NULL))),
    do.call("rbind", lapply(1:length(indSim), function(ind) tryCatch(as_tibble(indSim[[ind]][[2]][[4]]) %>% setNames(c("id", "time", "site", "x")) %>% mutate(id = glue::glue("{ind}_{id}")), error = function(e) NULL))),
    do.call("rbind", lapply(1:length(indSim), function(ind) tryCatch(as_tibble(indSim[[ind]][[3]][[4]]) %>% setNames(c("id", "time", "site", "x")) %>% mutate(id = glue::glue("{ind}_{id}")), error = function(e) NULL))),
    do.call("rbind", lapply(1:length(indSim), function(ind) rbind(tryCatch(tibble(t = 1, f = indSim[[ind]][[1]][[5]]), error = function(e) NULL), 
                                                                  tryCatch(tibble(t = 2, f = indSim[[ind]][[2]][[5]]), error = function(e) NULL), 
                                                                  tryCatch(tibble(t = 3, f = indSim[[ind]][[3]][[5]]), error = function(e) NULL))))
  )
})

save(allSpSim, file = "output/allSpSim_R2.rda")
