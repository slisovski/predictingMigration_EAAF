## Supporting Information for
## Predicting resilience of migratory birds to environmental change
## Lisovski et al.

library(tidyverse)
library(sf); sf_use_s2(FALSE)
library(stars)

load("output/eaafMap.rda")
load("output/mudflatTab.rda")
load("output/tempTab_revision.rda")
load("output/breedTab_revision.rda")

#####################################
###### Figures ######################
#####################################

### S1: Model description

## Figure S1
{
w  <- 0.028
xc <- 10 
x  <- seq(0, 100, length = 100)

eqA1 <- 0.5*(((exp(w*(x-xc)) - exp(-w*(x-xc)))/(exp(w*(x-xc)) + exp(-w*(x-xc))))+1)

ggplot(tibble(Rx = eqA1, x = x), aes(x = x, y = Rx)) +
  geom_line(color = "darkred", linewidth = 2)
}

### S2: Environmental parameters

## Figure S2
{
cells <- c(355, 303, 574)

ggplot() +
  geom_sf(data = eaafMap$map, fill = "grey90", color = "grey30") +
  geom_sf(data = eaafMap$bbox, fill = NA, color = "grey10") +
  geom_sf(data = mudflatTab$geometry, fill = NA) +
  geom_sf(data = mudflatTab$geometry[cells], fill = c("orange", "darkred", "darkblue"), alpha = 0.75) +
  theme_void()
}

## Figure S3
{
temp_path <- "NCEP/Surface_temp/" ## path to NCEP dataset

flsTemp <- tibble(path = list.files(temp_path)) %>% 
  mutate(year = as.numeric(unlist(lapply(strsplit(path, ".nc"), function(x) substring(x, 12, 15))))) %>%
  filter(year%in%c(1960:2070) | year%in%c(2010:2020))

   ## template
   tmp_r <- brick(paste0(path, flsTemp$path[1]))
   tmp_i <- rasterize(mudflatTab$geometry %>% st_transform(4326) %>% st_shift_longitude() %>% as("Spatial"), tmp_r[[1]])
   
   #### To be done with access to NCEP files
   
   
## CMIP4 for future scenario  
# Copernicus Toolbox Request
#' {
#'      import cdstoolbox as ct
#'      
#'      @ct.application(title='Download data')
#'      @ct.output.download()
#'      def download_application():
#'        data = ct.catalogue.retrieve(
#'          'projections-cmip6',
#'          {
#'            'temporal_resolution': 'daily',
#'            'experiment': 'ssp1_2_6',
#'            'variable': 'near_surface_air_temperature',
#'            'model': 'awi_cm_1_1_mr',
#'            'year': [
#'              '2060', '2061', '2062',
#'              '2063', '2064', '2065',
#'              '2066', '2067', '2068',
#'              '2069', '2070',
#'            ],
#'            'month': [
#'              '01', '02', '03',
#'              '04', '05', '06',
#'              '07', '08', '09',
#'              '10', '11', '12',
#'            ],
#'            'day': [
#'              '01', '02', '03',
#'              '04', '05', '06',
#'              '07', '08', '09',
#'              '10', '11', '12',
#'              '13', '14', '15',
#'              '16', '17', '18',
#'              '19', '20', '21',
#'              '22', '23', '24',
#'              '25', '26', '27',
#'              '28', '29', '30',
#'              '31',
#'            ],
#'          }
#'        )
#'      return data
#'    }  

cimp_file <- "tas_day_AWI-CM-1-1-MR_ssp245_r1i1p1f1_gn_20600101-20651231_v20190529.nc"

cimp  <- (read_stars(cimp_file) %>%
          st_set_crs(4326)) %>% suppressMessages() %>% suppressWarnings()
dates <- st_get_dimension_values(cimp, 3)

cimpDat <- lapply(cells, function(c) {
  poly <- mudflatTab[c,] %>% st_geometry() %>% st_transform(st_crs(cimp)) %>% st_shift_longitude()
  rst  <- st_extract(cimp, poly) %>% suppressMessages() %>% suppressWarnings()
  arr  <- apply(rst[[1]], 1:2, function(y) c(y))
  
  tab <- tibble(date = dates, temp = c(unlist(arr)) - 273.15) %>% filter(!is.na(temp)) %>%
    mutate(doy = as.numeric(format(date, "%j"))) %>% dplyr::select(doy, temp) %>%
    mutate(cell = c, .before = doy)
}) %>% Reduce("rbind",.)

loessDat <- cimpDat %>% filter(cell == 355) %>% 
  bind_rows(cimpDat %>% filter(cell == 355) %>% mutate(doy = doy - 365)) %>% 
  bind_rows(cimpDat %>% filter(cell == 355) %>% mutate(doy = doy + 365))
  
loessLine <- predict(loess(temp~doy, data = loessDat, span = 0.2), newdata = data.frame(doy = 1:365))

ggplot(cimpDat %>% filter(cell == 355), aes(x = doy, y = temp)) +
  geom_point() +
  xlab("doy") + ylab("temperature") +
  geom_line(data = tibble(loess = loessLine, doy = 1:365), 
            mapping = aes(x = doy, y = loess), color = "orange", linewidth = 2) +
  theme_light()
   
}

## Figure S4
{
  tempScen <- tibble(scenarios = rep(c("1960s", "2010s", "2060s"), each = 365),
                     doy      = rep(1:365, 3),
                     temp     = unlist(lapply(1:3, function(x) tempTab[303,,x])))
  
  ggplot(tempScen, aes(x = doy, y = temp, color = scenarios)) +
    geom_line(linewidth = 1) +
    ylab("temperature") +
    theme_light()

}
   
## Figure S5
{
  ## Tidal mudflat with buffer
  cell <- 303
  poly <- mudflatTab$geometry[303]
  
  mudflat_folder <- "Murray_etal_2018/Extracted polygons"
  
  fls  <- list.files(mudflat_folder, pattern = "*extracted_union.shp", full.names = T)
  bf   <- st_centroid(poly) %>% st_buffer(400000) %>% st_as_sf()
  
  pols <- lapply(fls, function(f) {
      mud <- read_sf(f) %>% st_union() %>% st_transform(st_crs(bf)) %>% st_intersection(bf)
      if(length(mud)>0) mud else NULL
  })

    ggplot() +
      geom_sf(data = bf) +
      geom_sf(data = pols %>% Reduce("rbind",.) %>% st_as_sfc() %>% st_set_crs(st_crs(bf)), fill = "darkred", color = "transparent") +
      geom_sf(data = poly, fill = "transparent") +
      theme_light()
}
   
## Figure S6
{
  ## Snowmelt
  cell <- 574
  {
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
    ## end maximum likelihood function
    } ## likelihood function
  
  library(ncdf4)

  ncFile <- "nhsce_v01r01_19661004_20220103.nc"
  
  projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
  ncDat <- nc_open(ncFile)
  t <- ncvar_get(ncDat, "time")
  z <- ncvar_get(ncDat, "snow_cover_extent")
  ylat <- ncvar_get(ncDat, "longitude")
  xlon <- ncvar_get(ncDat, "latitude")
  nc_close(ncDat)
  
  ## date
  start <- as.POSIXct("1966-10-03", "GMT")
  date  <- start + ((t-6)*24*60*60)
  doy   <- as.numeric(format(date, "%j"))

  pl <- (st_as_sf(data.frame(lat = c(xlon), lon = c(ylat), index = 1:(dim(z)[2]*dim(z)[1])), coords = c("lon", "lat")) %>%
                    st_set_crs(4326) %>% st_transform(projSnow))
  
  ind <- which((pl %>% st_intersects(mudflatTab$geometry[cell] %>% st_transform(projSnow), sparse = F))[,1])
  snowDat <- lapply(1:dim(z)[3], function(x) tibble(date = date[x], snow = c(z[,,x])[ind])) %>% Reduce("rbind",.) %>%
    mutate(doy = as.numeric(format(date, "%j"))) %>% filter(as.numeric(format(date, "%Y")) %in% c(1975:1985, 2010:2020)) %>%
    mutate(scenario = if_else(as.numeric(format(date, "%Y"))<2000, "1960s", "2010s"))
  
  sm <-  lapply(unique(snowDat$scenario), function(x) {
    suppressWarnings(with(snowDat %>% filter(scenario==x), gaussMLE(day = doy, size = rep(100, length(doy)), prob = snow, thresh = 0.75)))
  })
  
  ggplot(snowDat, aes(x = doy, y = jitter(snow, 0.1), color = scenario)) +
    geom_point(shape = 16, alpha = 0.7) +
    scale_color_manual(values = c("orange4", "darkgreen"), name = "Scenarios") +
    geom_line(data = sm[[1]]$result, mapping = aes(x = week, y = fit), color = "orange4", linewidth = 1.2) +
    geom_point(data = tibble(doy = sm[[1]]$start, sm = 0.75), mapping = aes(x = doy, y = sm), color = "orange4", shape = 18, size = 6) +
    geom_line(data = sm[[2]]$result, mapping = aes(x = week, y = fit), color = "darkgreen", linewidth = 1.2) +
    geom_point(data = tibble(doy = sm[[2]]$start, sm = 0.75), mapping = aes(x = doy, y = sm), color = "darkgreen", shape = 18, size = 6) +
    ylab('Snow cover (1 = snow, 0 = snow free")') +
    theme_light()
  
  
}

## Figure S7
{
  snowTempTab <- tibble(ID = apply(st_intersects(breedTab %>% st_transform(st_crs(mudflatTab)), mudflatTab, sparse = F), 1, which),
                        start_hist = breedTab$start_hist, start_curr = breedTab$start_curr) %>%
                 mutate(temp_hist    = apply(., 1, function(x) tempTab[x[1], x[2], 1]), 
                        temp_curr = apply(., 1, function(x) tempTab[x[1], x[3], 2])) %>%
                 mutate(start_future = apply(., 1, function(x) min(which(tempTab[as.numeric(x[1]), , 3] >= median(as.numeric(x[4:5])))))) %>%
                 mutate(temp_future  = apply(., 1, function(x) tempTab[x[1], x[6], 3])) %>% dplyr::select("ID", "start_hist", "start_curr", "start_future",
                                                                                                          "temp_hist", "temp_curr", "temp_future")
  
  ggplot(snowTempTab[,1:4] %>% setNames(c("ID", "1960s", "2010s", "2060s")) %>% pivot_longer(-ID), aes(x = value, y = ..density.., fill = as.factor(name))) +
    geom_histogram(show.legend = FALSE) +
    geom_density(bw = 6, fill = NA, linewidth = 0.4) +
    scale_fill_manual(values = c("grey90", "grey50", "grey20")) +
    ylab("Frequency") + xlab("Snowmelt [doy]") +
    facet_wrap(~name, nrow = 3) +
    theme_light()

}

### S3: Species parameters

spParms <- tibble(species = c("Godwit", "GreatKnot", "RedKnot", "CurlewSandpiper", "RedNeckedStint"),
                  LBM = c(250, 144, 100, 55, 25))

## Figure S9
{
  lbm <- seq(0.01, 0.3, length = 100)
  fdr <- 2.37 * lbm^-0.27
  
  opar <- par(mar = c(4,6,1,1))
  plot(lbm, fdr, type = 'l', ylim = c(0, 9), xlab = "LBM [kg]",
       ylab = "Fuel deposition rate \n(% of lean body mass /day)")
  abline(v = spParms$LBM/1000, lty = 3)
  par(opar)
}



## Figure S10 BMR
empTab <- tibble(species = c("Oystercatcher", "Littel ringed plover", 
                             "Grey plover", "Sanderling", "Redshank", "Turnstone"),
                 BM = c(554, 36, 226, 50, 149, 114),
                 BMR_W = c(2.91, 0.41, 1.78, 0.56, 1.56, 0.99)) %>% mutate(id = 1:nrow(.))
   
plot(empTab$BM, empTab$BMR_W, xlim = c(30,600), ylim = c(0.35, 3.5), las = 1,
     xlab = "Lean Body Mass (g)", ylab = "Basal Metabolic Rate (Watts)",
     pch = 21, cex = 4, bg = adjustcolor("cornflowerblue", alpha.f = 0.4))
text(empTab$BM, empTab$BMR_W, empTab$id)
lines(seq(30, 600, length = 100), 5.06 * (seq(30, 600, length = 100)/1000)^(0.729), lty = 2)
points(spParms$LBM, 5.06 * (spParms$LBM/1000)^(0.729), pch = 16, cex = 2, col = adjustcolor("firebrick", alpha.f = 0.4))
text(spParms$LBM+50, 5.06 * (spParms$LBM/1000)^(0.729), spParms$species, pos = 4)
