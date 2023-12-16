library(sf); sf::sf_use_s2(FALSE)
library(raster)
library(tidyverse)
library(stars)

## Basic flyway map
proj <- "+proj=laea +lat_0=15 +lon_0=162 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

xlim <- c(60, 210)
ylim <- c(-55, 82)

poly <- st_as_sf(as(extent(c(xlim, ylim)), "SpatialPolygons")) %>% st_set_crs(4326) %>% st_transform(CRS(proj)) %>% st_bbox() %>%
  st_as_sfc() %>% st_set_crs(proj)

map <- read_sf("resources/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>%
  st_difference(st_as_sf(as(extent(c(-10,10,-90,90)), "SpatialPolygons")) %>% st_set_crs(4326)) %>%
  st_shift_longitude() %>% st_geometry() %>% st_transform(crs = proj) %>% st_buffer(0) %>% st_intersection(poly)

grid <- read_sf("resources/NaturalEarth/10m_physical/ne_10m_graticules_20/ne_10m_graticules_20.shp", ) %>%
  st_shift_longitude() %>% st_geometry() %>% st_transform(crs = CRS(proj)) %>% st_intersection(poly)

eaafMap <- list(map = map, grid = grid, bbox = poly)
# save(eaafMap, file = "output/Map/eaafMap.rda")
# load("output/Map/eaafMap.rda")

hex <- st_make_grid(eaafMap$bbox, cellsize = 400000, square = FALSE) %>% st_set_crs(proj)

ggplot() +
  geom_sf(data = eaafMap$grid, linetype = 3, linewidth = 0.4) +
  geom_sf(data = eaafMap$map, color = "grey30", fill = "grey80") +
  geom_sf(data = eaafMap$bbox, color = "grey20", fill = NA) +
  geom_sf(data = hex, fill = NA, color = "cornflowerblue", alpha = 0.4) +
  theme_void()



###############################################
## Intertidal mudflats                       ##
### Historical mudflat extent (Yellow Sea)   ##
###############################################

hist <- read_sf("Murray_etal_2014/East_Asia_Tidal_Flat_1952_1964.shp") %>%     
  st_transform(st_crs(eaafMap$map)) %>% st_union()

areaHist <- hex %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% mutate() %>% st_intersection(hist %>% st_buffer(0)) %>%
  mutate(area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>% full_join(
    st_centroid(hex) %>% st_buffer(400000) %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% st_intersection(hist %>% st_buffer(0)) %>%
      mutate(area_bfr = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble()
  ) %>% full_join(hex %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% rename(geometry = x)) %>% arrange(id) %>% st_as_sf()


###############################################
### Current mudflat extent                   ##
###############################################

fls  <- list.files("Murray_etal_2018/Extracted polygons", pattern = "*extracted_union.shp", full.names = T)

mudTabs <- do.call("cbind", parallel::mclapply(fls, function(f) {
  mud <- read_sf(f) %>% st_union() %>% st_transform(st_crs(hex))
  areaHist %>% select(id) %>% st_intersection(mud) %>%
    mutate(innter_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>% full_join(
      st_centroid(areaHist) %>% select(id) %>% st_buffer(400000) %>% st_intersection(mud) %>%
        mutate(full_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble()
    ) %>% full_join(tibble(id = 1:length(hex))) %>% arrange(id) %>% select(innter_area, full_area)
}, mc.cores = 3))

areaCurr <- mudTabs %>% setNames(paste0(rep(c("a", "b"), ncol(.)/2), seq_along(.))) %>%
  mutate(currArea_inner = rowSums(select(., starts_with("a")), na.rm = T),
         currArea_outer = rowSums(select(., starts_with("b")), na.rm = T))


###############################################
### All mudflat extent                       ##
###############################################

mudflatTab <- areaHist %>% rename(histArea_inner = area, histArea_outer = area_bfr) %>%
  bind_cols(areaCurr %>% select(c(currArea_inner, currArea_outer))) %>% select(histArea_inner, histArea_outer, currArea_inner, currArea_outer)

# ggplot(eaafMap$map) + geom_sf() +
# geom_sf(mudflatTab, mapping = aes(geometry = geometry, fill = currArea_outer), col = NA)


###############################################
### Mangroves                                ##
###############################################

mang <- read_sf("/Volumes/slisovski/RemoteSensedData/MangroveDistribution/01_Data/14_001_WCMC010_MangroveUSGS2011_v1_3.shp") %>% st_transform(proj) %>%
  st_buffer(0) %>% st_intersection(st_bbox(mudflatTab) %>% st_as_sfc())

mangArea <- mudflatTab %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>%
  st_intersection(mang %>% select(AREA_KM2)) %>% st_drop_geometry() %>% group_by(id) %>% summarise(mangArea_inner = sum(AREA_KM2)) %>%
  full_join(mudflatTab %>% st_centroid() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>% st_buffer(400000) %>%
              st_intersection(mang %>% select(AREA_KM2)) %>% st_drop_geometry() %>% group_by(id) %>% summarise(mangArea_outer = sum(AREA_KM2))) %>%
  full_join(tibble(id = 1:nrow(mudflatTab))) %>% arrange(id)

### join with mudflat table
mudflatTab <- mudflatTab %>% mutate(mangArea_inner = ifelse(is.na(mangArea$mangArea_inner), 0, mangArea$mangArea_inner),
                                    mangArea_outer = ifelse(is.na(mangArea$mangArea_outer), 0, mangArea$mangArea_outer)) %>%
  relocate(geometry, .after = last_col())

###############################################
### Lakes                                    ##
###############################################

mudflatTab <- mudflatTab %>% mutate(lake_area = mudflatTab %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>% st_intersection(
      read_sf("Data/Map/ne_10m_lakes/ne_10m_lakes.shp" )%>% st_transform(proj) %>%
      st_buffer(0) %>% st_intersection(st_bbox(mudflatTab) %>% st_as_sfc()) %>%
      filter(st_coordinates(st_centroid(.) %>% st_transform(4326))[,2]<65 &
               (st_coordinates(st_centroid(.) %>% st_transform(4326))[,1]>100 |
                  st_coordinates(st_centroid(.) %>% st_transform(4326))[,1]<0)) %>%
      mutate(area = as.numeric(st_area(.))) %>% filter(area >= quantile(.$area, probs = 0.9)) %>% st_combine()) %>%
  mutate(lake_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>%
  full_join(tibble(id = 1:nrow(mudflatTab)))  %>% arrange(id) %>% pull(lake_area)) %>% relocate(geometry, .after = last_col())

# save(mudflatTab, file = "output/mudflatTab.rda")
# load("output/mudflatTab.rda")

# ggplot(eaafMap$map) + geom_sf() +
#   geom_sf(mudflatTab, mapping = aes(geometry = geometry, fill = lake_area), col = NA)



###############################################
### Temperature (daily)                      ##
###############################################

path <- "NCEP/Surface_temp/"

flsTemp <- tibble(path = list.files(path)) %>% mutate(year = as.numeric(unlist(lapply(strsplit(path, ".nc"), function(x) substring(x, 12, 15))))) %>%
  filter(year%in%c(1950:2070) | year%in%c(2000:2020))

   ## template
   tmp_r <- brick(paste0(path, flsTemp$path[1]))
   tmp_i <- rasterize(mudflatTab$geometry %>% st_transform(4326) %>% st_shift_longitude() %>% as("Spatial"), tmp_r[[1]])

tempArray <- abind::abind(parallel::mclapply(flsTemp$path, function(x) {

  (brick(paste0(path, x)) %>% st_as_stars() %>% st_set_crs(4326) %>% st_as_sf() %>% st_centroid() %>%
    st_transform(st_crs(hex)) %>%
    mutate(cell = apply(st_intersects(., hex, sparse = FALSE), 1, function(t) ifelse(any(t), which(t), NA))) %>%
    filter(!is.na(cell)) %>% st_drop_geometry() %>% group_split(cell) %>%
    lapply(., function(z) apply(z[,-ncol(z)], 2, function(y) quantile(as.numeric(y), probs = 0.4, na.rm = T))) %>% Reduce("rbind", .))[,1:364]

    tmp_r <- brick(paste0(path, flsTemp$path[1]))
    
    proj4string(tmp_r) <- "+proj=longlat +ellps=WGS84 +pm=-360 +datum=WGS84 +no_defs"
    tt <- data.frame(lon = coordinates(tmp_r)[,1], lat = coordinates(tmp_r)[,2], temp = tmp_r[[1]][])
    tt$lon[tt$lon>180] <-  - (360 - tt$lon[tt$lon>180])
    rast <- raster::rasterize(tt[,1:2], raster(xmn = -179.9, xmx = 179.9, ymn = -90, ymx = 90, res = 2.5), field = tt$temp)
    tempGr <- st_make_grid(eaafMap$bbox, cellsize = 150000) %>% st_transform("+proj=longlat")
    extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
    tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)
    
    tmp_t <- tibble(index = tmp_i[]) %>% bind_cols(tmp_r[]) %>% filter(!is.na(index)) %>% pivot_longer(cols = starts_with("X")) %>%
      group_by(index, name) %>% summarise(temp = quantile(value, probs = 0.4, na.rm = T)) %>%
      pivot_wider(names_from = name, values_from = temp) %>% dplyr::select(-index) %>% as.matrix()
    tmp_t[,1:364]

}, mc.cores = 7), along = 3)

## daily values (two periods)
tempTab <- abind::abind(parallel::mclapply(1:dim(tempArray)[1], function(x) {

  tmp01 <- tibble(y = rep(flsTemp$year, each = dim(tempArray)[2]), doi = rep(1:364, dim(tempArray)[3]), temp = c(tempArray[x,,]) - 273.15)
  tmp02 <- tmp01 %>% bind_rows(tmp01 %>% mutate(doi = doi-364)) %>% bind_rows(tmp01 %>% mutate(doi = doi+364))

  mod0 <- predict(loess(temp~doi, data = tmp02 %>% filter(y > 1960 & y < 1970), span = 0.2), newdata = data.frame(doi = 1:365))
  mod1 <- predict(loess(temp~doi, data = tmp02 %>% filter(y > 2010 & y < 2020), span = 0.2), newdata = data.frame(doi = 1:365))

  array(c(as.numeric(mod0), as.numeric(mod1)), dim = c(1, 365, 2))

}, mc.cores = 8), along = 1)

tempTab <- abind::abind(tempTab, apply(tempTab, 1:2, function(x)
  predict(lm(temp~years, data = tibble(years = c(1960, 2010), temp = c(x))), newdata = data.frame(years = 2060))),
  along = 3)


### CIMP4 AWI SPP2 Scenario (2060s)
cimp  <- (read_stars("Data/CIMP6_SPP2/tas_day_AWI-CM-1-1-MR_ssp245_r1i1p1f1_gn_20600101-20651231_v20190529.nc") %>%
            st_set_crs(4326)) %>% suppressMessages() %>% suppressWarnings()
dates   <- st_get_dimension_values(cimp, 3)
breedID <- unique(apply(st_intersects(breedTab %>% st_transform(st_crs(mudflatTab)), mudflatTab, sparse = F), 1, which))

tempCimp <- matrix(nrow = nrow(mudflatTab), ncol = 365)

for(x in 1:nrow(mudflatTab)) {
  
  cat(sprintf('\n pixel %d if %d', x, nrow(mudflatTab))) 
  
  if((!is.na(sum(mudflatTab[x,] %>% st_drop_geometry(), na.rm = T)) & 
      sum(mudflatTab[x,] %>% st_drop_geometry(), na.rm = T) > 0 &
      all(is.na(tempCimp[x,]))) || x %in% breedID) {
    
    poly <- mudflatTab[x,] %>% st_geometry() %>% st_transform(st_crs(cimp)) %>% st_shift_longitude()
    rst  <- st_extract(cimp, poly) %>% suppressMessages() %>% suppressWarnings()
    arr  <- apply(rst[[1]], 1:2, function(y) c(y))
    
    # tab <- tibble(date = rep(dates, dim(arr)[2]*dim(arr)[3]), temp = c(unlist(arr)) - 273.15) %>% filter(!is.na(temp)) %>%
    #   mutate(doy = as.numeric(format(date, "%j"))) %>% dplyr::select(doy, temp) %>%
    #   group_by(doy) %>% summarise(temp = median(temp))
    
    tab <- tibble(date = dates, temp = c(unlist(arr)) - 273.15) %>% filter(!is.na(temp)) %>%
      mutate(doy = as.numeric(format(date, "%j"))) %>% dplyr::select(doy, temp) %>%
      group_by(doy) %>% summarise(temp = median(temp))
    tabMod    <- rbind(tab %>% mutate(doy = doy - 365), tab, tab %>% mutate(doy = doy + 365))
    loessTemp <- predict(loess(temp~doy, data = tabMod, span = 0.2), newdata = data.frame(doy = 1:365))
    tempCimp[x,] <- as.numeric(loessTemp)
    invisible(gc())
  }
}

tempTab[,,3] <- tempCimp
# save(tempTab, file = "resources/tempTab_revision.rda")
# load("resources/tempTab_revision.rda")


###############################################
### Snow melt timing                         ##
###############################################
library(bbmle)
gaussMLE <- function(day, size, prob, thresh) {

  tab <- data.frame(day = day, size = size, p = prob)

  gauss.curve <- function(parms, intv = 1) {
    t <- seq(1, 366, intv)
    parms <- as.list(parms)
    fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
    fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
    c(fit1, fit2[-1])
  }

  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), 1)
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\r")
    -sum(dbinom(x = round(tab[,3]*tab[,2],0), size = rep(100, length(fit)), prob = fit[day], log=TRUE), na.rm=T)
  }

  mle <- suppressWarnings(mle2(gauss.loglik, method="L-BFGS-B",
                               start=list(a1 = 225, a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
                               lower=list(a1 = 50,  a2 = 5,   a3 = 0.5,a4=5,  a5 = 0.5),
                               upper=list(a1 = 225, a2 =  Inf,  a3 =  Inf, a4 =  Inf, a5 =  Inf),
  ))

  t <- seq(1, 366, 1)
  fit <- gauss.curve(coef(mle), intv = 1)

  start <- t[min(which(fit<thresh))]

  list(result = data.frame(week = t, fit = fit), start = start, end = end)
}
# end maximum likelihood function

geoTab <- readxl::read_excel("output/geolocTab.xlsx") %>% mutate(migDur = as.numeric(difftime(arr, dep, units = "days")))

library(ncdf4)

projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
ncDat <- nc_open("~/Google Drive/My Drive/GeoDat/nhsce_v01r01_19661004_20220103.nc")
t <- ncvar_get(ncDat, "time")
z <- ncvar_get(ncDat, "snow_cover_extent")
ylat <- ncvar_get(ncDat, "longitude")
xlon <- ncvar_get(ncDat, "latitude")
nc_close(ncDat)

## date
start <- as.POSIXct("1966-10-03", "GMT")
date  <- start + ((t-6)*24*60*60)
doy   <- as.numeric(format(date, "%j"))


dat <- matrix(c(z), ncol = dim(z)[3], byrow = F)

pl <- (st_as_sf(data.frame(lat = c(xlon), lon = c(ylat), index = 1:(dim(z)[2]*dim(z)[1])), coords = c("lon", "lat")) %>%
                  st_set_crs(4326) %>% st_transform(projSnow)) %>% mutate(snow = dat[,122]) %>% dplyr::select(snow) %>%
  st_transform("+proj=longlat") %>% st_coordinates()


rast <- raster::rasterize(pl[,1:2], raster(xmn = -179.9, xmx = 179.9, ymn = -90, ymx = 90, res = 5), field = dat[,122])
rast[is.na(rast[])] <- 0
tempGr <- st_make_grid(eaafMap$bbox, cellsize = 150000) %>% st_transform("+proj=longlat")
extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)

extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)

rastInd <- (st_as_sf(data.frame(lat = c(xlon), lon = c(ylat), index = 1:(dim(z)[2]*dim(z)[1])), coords = c("lon", "lat")) %>%
              st_set_crs(4326) %>% st_transform(projSnow)) %>%
  filter(apply(dat, 1, function(x) sum(x, na.rm = T)>0))

crds_breed <- st_as_sf(geoTab, coords = c("lon_end", "lat_end"), crs = 4326) %>%
  st_transform(projSnow) %>% rowwise() %>% mutate(NHWCindex = rastInd$index[which.min(st_distance(geometry, rastInd))])

sm <-  do.call("rbind", lapply(crds_breed$NHWCindex, function(x) {

  tmp <- tibble(year = as.numeric(format(date, "%Y")), doy = doy, y = dat[x,])
  mod1 <- suppressWarnings(with(tmp %>% filter(year%in%c(1975:1985)), gaussMLE(day = doy, size = rep(100, length(doy)), prob = y, thresh = 0.75)))
  mod2 <- suppressWarnings(with(tmp %>% filter(year%in%c(2000:2020)), gaussMLE(day = doy, size = rep(100, length(doy)), prob = y, thresh = 0.75)))

  tibble(start_hist = mod1$start, start_curr = mod2$start)
  }))

breedTab <- crds_breed %>% bind_cols(sm)

### Start 2060s
snowTempTab <- tibble(ID = apply(st_intersects(breedTab %>% st_transform(st_crs(mudflatTab)), mudflatTab, sparse = F), 1, which),
                      start_hist = breedTab$start_hist, start_curr = breedTab$start_curr) %>%
               mutate(temp_hist    = apply(., 1, function(x) tempTab[x[1], x[2], 1]), temp_curr = apply(., 1, function(x) tempTab[x[1], x[3], 2])) %>%
               mutate(start_future = apply(., 1, function(x) min(which(tempTab[as.numeric(x[1]), , 3] >= median(as.numeric(x[4:5])))))) %>%
               mutate(temp_future  = apply(., 1, function(x) tempTab[x[1], x[6], 3])) %>% dplyr::select("ID", "start_hist", "start_curr", "start_future",
                                                                                                        "temp_hist", "temp_curr", "temp_future")

breedTab <- breedTab %>% mutate(start_future = snowTempTab$start_future)

# save(breedTab, file = "ouutput/breedTab_revision.rda")
# load("output/breedTab_revision.rda")



