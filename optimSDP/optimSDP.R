setClass(
  "SDPMig",
  slots = c(Init    = "list",
            Species = "list",
            Sites   = "list",
            Results = "list")
)

sizeParams     <- function(lbm) {
  
  out <- list(lbm   = lbm,
              mbm   = approx(c(15, 300), c(23.05, 643.4), xout = lbm)$y)
  
  ### Fuel deposition rate
  out$X1gX  <- (out$mbm-out$lbm)/100
  out$X1xkJ <- out$X1gX*32.77
  out$FDR   <- 2.37*(out$lbm/1000)^-0.27
  out$FDRg  <- (out$FDR*(out$lbm))/100
  out$FDRx  <- out$FDRg/out$X1gX
  out$daysToRefuel <- (out$mbm-out$lbm)/out$FDRg
  
  ### Metabolic rate
  out$DER    <- 912*(lbm/1000)^0.704
  out$BMR    <- 5.06*(lbm/1000)^0.729
  out$cond   <- 0.0067*(lbm)^0.4384
  out$Kesm   <- 41-(out$BMR/(1.2*out$cond))
  out$W      <- 57.3*(lbm/1000)^0.81
  
  ### Flight parameters
  x1  <- c(250, 144, 110, 20)
  dst <- c(14000, 7500, 6500, 3000) 
  fit1 <- lm(log(dst)~log(x1))
  pr   <- exp(predict(fit1))
  out$FlightCap <- exp(predict(fit1, newdata = data.frame(x1 = lbm))) + 1000
  
  # Energy exp ~ temp
  t  <- seq(-35, 45, length = 100)
  te <- 41 - (out$BMR/(1.2*out$cond))
  b = out$BMR - (-out$cond)*te
  
  tm <- data.frame(tm = seq(-35, te, length = 100), W = -out$cond*seq(-25, te, length = 100) + b)
  tm <- rbind(tm, data.frame(tm = seq(te, 45, length = 100), W = rep(out$BMR, 100)))
  
  out$EEFnc  <- suppressWarnings(approxfun(x = tm[,1], y = 2*tm[,2]*86.4, rule = 3))
  
  out$speed  <- 1440
  out$c      <- out$FlightCap / (1-(1+99/100)^-0.5)
  
  return(out)
}

Mat3Array      <- function(x) {
  
  rows = length(x)
  cols = length(x[[1]])
  d3 = length(x[[1]][[1]])
  
  out <- array(dim = c(rows, cols, d3))
  
  for(r in 1:rows) {
    for(c in 1:cols) {
      for(d.1 in 1:d3) {
        out[r, c, d.1] <- x[[r]][[c]][[d.1]]
      }
    }
  }
  out 
}

CalcTR         <- function(x, t) {
  TR <- approx(parms$xFTReward, parms$yFTReward, t, rule = 3)$y
  s  <- parms$w * (x - parms$xc)
  
  if(x == 0) { SR <- 0 } else {
    SR <- (((exp(s) - exp(-s)) / (exp(s) + exp(-s))) + 1)/2
  }
  TR * SR + parms$B0
}

makeSDPobjects <- function(parms) {
  
  
  pmax <- apply(parms$sites %>% st_drop_geometry(), 1, function(x) pmax(x[3]+((x[5]*(1000000))/1000)*0.3, x[4]+((x[5]*(1000000))/1000)*0.3))
  hist <- apply(parms$sites %>% st_drop_geometry(), 1, function(x) x[3]+((x[5]*(1000000))/1000)*0.3)
  curr <- apply(parms$sites %>% st_drop_geometry(), 1, function(x) x[4]+((x[5]*(1000000))/1000)*0.3)
  ## intake rates
  intake <- tibble(index = parms$site$index,
                   hist = approx(range(pmax, na.rm = T), 
                                 c(0.5, 1.5), hist)$y,
                   curr = approx(range(pmax, na.rm = T), 
                                 c(0.5, 1.5), curr)$y) %>%
    mutate(hist = ifelse(parms$species=="RedNeckedStint" & hist<1.25 & parms$sites$lake, 1.25, hist), 
           curr = ifelse(parms$species=="RedNeckedStint" & curr<1.25 & parms$sites$lake, 1.25, curr),
           intHist = hist * parms$spParms$FDRx + parms$spParms$EEFnc(parms$spParms$Kesm)/parms$spParms$X1xkJ,
           intCurr = curr * parms$spParms$FDRx + parms$spParms$EEFnc(parms$spParms$Kesm)/parms$spParms$X1xkJ)
  
  # future intake
  intake$intFut <- apply(intake %>% dplyr::select(intHist, intCurr) %>% as.matrix(), 1, function(x) {
    if(any(is.na(x))) max(x, na.rm = T) else predict(lm(intake~years, data = tibble(years = c(1960, 2010), intake = c(x))), newdata = data.frame(years = 2060))
  })
  # 
  ### Latitude EIR relationship
  intake <- intake %>% mutate(intHist = intHist * approx(seq(0,90, length = 100), as.numeric(parms$latF)^c(seq(0,1, length = 100)), abs((parms$site %>% st_coordinates())[,2]))$y,
                              intCurr = intCurr * approx(seq(0,90, length = 100), as.numeric(parms$latF)^c(seq(0,1, length = 100)), abs((parms$site %>% st_coordinates())[,2]))$y,
                              intFut  = intFut  * approx(seq(0,90, length = 100), as.numeric(parms$latF)^c(seq(0,1, length = 100)), abs((parms$site %>% st_coordinates())[,2]))$y)
  
  # opar <- par(mfrow = c(1,3))
  # plot(eaafMap$map, col = "grey90", border = "grey30")
  # plot(parms$site[parms$site$index%in%intake$index,] %>% st_transform(st_crs(eaafMap$map)) %>% st_geometry(),
  #      col = rev(viridis::inferno(100, alpha = 0.6))[cut(intake %>% pull("intHist"), seq(min(intake[-nrow(intake),4:6], na.rm = T), max(intake[-nrow(intake),4:6], na.rm = T), length = 99), labels = F)],
  #      pch = 18, cex = 3, add = T)
  # plot(eaafMap$map, col = "grey90", border = "grey30")
  # plot(parms$site[parms$site$index%in%intake$index,] %>% st_transform(st_crs(eaafMap$map)) %>% st_geometry(),
  #      col = rev(viridis::inferno(100, alpha = 0.6))[cut(intake %>% pull("intCurr"), seq(min(intake[-nrow(intake),4:6], na.rm = T), max(intake[-nrow(intake),4:6], na.rm = T), length = 99), labels = F)],
  #      pch = 18, cex = 3, add = T)
  # plot(eaafMap$map, col = "grey90", border = "grey30")
  # plot(parms$site[parms$site$index%in%intake$index,] %>% st_transform(st_crs(eaafMap$map)) %>% st_geometry(),
  #      col = rev(viridis::inferno(100, alpha = 0.6))[cut(intake %>% pull("intFut"), seq(min(intake[-nrow(intake),4:6], na.rm = T), max(intake[-nrow(intake),4:6], na.rm = T), length = 99), labels = F)],
  #      pch = 18, cex = 3, add = T)
  # par(opar)
  
  ## energy expenditure
  expend <- apply(parms$tTab[,parms$minT:parms$maxT,], 1:3,function(x) parms$spParms$EEFnc(x)/parms$spParms$X1xkJ)
  
  ## predation
  pred <- tibble(b0 = rep(parms$pred[1], nrow(intake)), b1 = parms$pred[2], b2 = parms$pred[3])
  
  ### start site
  intake[1,4:5]  <- (parms$spParms$FDRx + parms$spParms$EEFnc(parms$spParms$Kesm)/parms$spParms$X1xkJ)
  pred[1,1]      <- pred[1,1]*0.1
  pred[1,2]      <- pred[1,2]*0.1
  pred[1,3]      <- pred[1,3]*0.1
  
  lapply(1:3, function(x) {
    new(
      "SDPMig",
      Init  = list(
        MinT   = parms$minT,
        MaxT   = parms$maxT,
        NSites = nrow(intake)-1,
        MaxX   = parms$MaxX
      ),
      Species = list(
        B0    = parms$B0,
        w     = parms$w,
        xc    = parms$xc,
        c     = as.numeric(parms$spParm$c),
        speed = parms$spParm$speed,
        WindAssist = parms$WindAssist,
        WindProb   = parms$WindProb,
        ZStdNorm = parms$ZStdNorm,
        PStdNorm = parms$PStdNorm,
        xFTReward  = as.numeric(c(0, c(unlist(c(parms$reward[x]-1, parms$reward[x], parms$reward[x]+10, parms$reward[x]+11)), parms$maxT) - parms$minT)),
        yFTReward  = c(0,0,2,2,0,0),
        decError   = parms$decError
      ),
      Sites = list(
        crds  = parms$site %>% st_coordinates(),
        index = parms$sites$index,
        dist  = parms$dist,
        bear  = parms$bear,
        b0    = pred %>% pull(b0),
        b1    = pred %>% pull(b1),
        b2    = pred %>% pull(b2),
        pred_a1 =  parms$pred[4],
        pred_a2 =  parms$pred[5],
        expend  =  as.matrix(expend[,,x]),
        gain    =  as.numeric(c(unlist(intake[,x+3])))
      ),
      Results = list(
        FitnessMatrix     = NA,
        DecisionMatrix    = NA,
        ProbMatrix        = NA
      )
    )
  })
  
}

##########################
### Backward iteration ###
##########################

bwdIteration <- function(obj, pbar = FALSE) {
  
  Init(obj@Init$MinT, 
       obj@Init$MaxT, 
       obj@Init$NSites, 
       obj@Init$MaxX,
       obj@Species$w,
       obj@Species$xc,
       obj@Species$B0,
       obj@Sites$b0,
       obj@Sites$b1,
       obj@Sites$b2,
       obj@Sites$pred_a1,
       obj@Sites$pred_a2,
       obj@Species$c,
       obj@Species$speed,
       obj@Species$WindAssist,
       obj@Species$WindProb,
       obj@Species$ZStdNorm,
       obj@Species$PStdNorm,
       as.numeric(obj@Species$xFTReward),
       obj@Species$yFTReward,
       obj@Species$decError,
       as.matrix(obj@Sites$dist),
       as.matrix(obj@Sites$bear),
       obj@Sites$gain,
       as.matrix(obj@Sites$expend))
  
  out <- BackwardIteration(pbar = pbar)
  
  obj@Results$FitnessMatrix <- out[[1]]
  DM <- array(dim = c(dim(out[[2]]),2))
  DM[,,,1] <- out[[2]]
  DM[,,,2] <- out[[3]]
  obj@Results$DecisionMatrix <- DM
  PM <- array(dim = c(dim(out[[4]]),2))
  PM[,,,1] <- out[[4]]
  PM[,,,2] <- out[[5]]
  obj@Results$ProbMatrix <- PM
  
  obj
}

############################
#### Forward Simulation#####
############################

fwdSimulation <- function(model, NrInd, start_t, start_site, start_x) {
  
  InitSim(model@Init$MinT, 
          model@Init$MaxT, 
          model@Init$NSites, 
          model@Init$MaxX,
          model@Species$w,
          model@Species$xc,
          model@Species$B0,
          model@Sites$b0,
          model@Sites$b1,
          model@Sites$b2,
          model@Sites$pred_a1,
          model@Sites$pred_a2,
          model@Species$c,
          model@Species$speed,
          model@Species$WindAssist,
          model@Species$WindProb,
          model@Species$ZStdNorm,
          model@Species$PStdNorm,
          model@Species$xFTReward,
          model@Species$yFTReward,
          model@Species$decError,
          model@Sites$dist,
          model@Sites$bear,
          model@Sites$gain,
          model@Sites$expend)
  
  
  x <- round(runif(NrInd, start_x[1], start_x[2]),0)
  
  if(length(start_site)>1 & length(start_site)<start_x[1]) {
    stop("start_site must have same length as numbers of individuals or a single site.")
  }
  if(length(start_site)==1) start_site <- rep(start_site, NrInd)
  
  SimOut = array(dim = c(length(x), 6, dim(model@Results$FitnessMatrix)[1]))
  
  ### First entry
  for(i in 1:dim(SimOut)[1]) {
    SimOut[i, ,start_t] <- c(start_t, start_site[i], x[i], 0, 0, 0)
  }
  
  
  ## SimOut: 1 = time, 2 = site, 3 = x, 4 = decision, 5 = flying, 6 = dead {
  for(time in 1:(dim(SimOut)[3]-1)) {
    
    for(ind in 1:dim(SimOut)[1]) {
      
      ## Not dead, not arrived, not flying
      if(!SimOut[ind, 6, time] & 
         sum(SimOut[ind, 2, ] >= nrow(model@Sites$crds), na.rm = T)<1 & !SimOut[ind, 5, time]) {
        
        ## Decision
        if(runif(1) <  model@Results$ProbMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]) {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]
        } else {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 2]
        }
        
        ## Action
        if(decision>=0) { ## Flying
          
          fl_help = simFlying(decision, time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          nextt = fl_help[1] + 1
          if(nextt<=time) nextt <- time+1
          if(nextt>dim(SimOut)[3]) time <- dim(SimOut)[3]
          
          nextx = fl_help[2] + 1
          if(nextx <  0) {
            nextx = 0
            dead  = 1 } else  dead = 0
          if(nextx > model@Init$MaxX) nextx = model@Init$MaxX
          
          SimOut[ind,,nextt] = c(nextt, decision+1, nextx, NA, 0, dead)
          if(nextt>(time+1)) SimOut[ind,5:6,(time+1):(nextt-1)] = cbind(1,0)
          if(SimOut[ind, 6, nextt])  SimOut[ind, 6, nextt:dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
          if(SimOut[ind, 2, nextt]==nrow(model@Sites$crds)) {
            SimOut[ind,2:6, nextt:dim(SimOut)[3]] <- SimOut[ind, 2:6, nextt]
            SimOut[ind,1,   nextt:dim(SimOut)[3]] <- seq(nextt, dim(SimOut)[3])
          }
          
        } else { ## Feeding
          
          fo_help = simForaging(abs(decision+1.0), time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          newx = fo_help[1]+1
          dead = fo_help[2]
          if(newx<=0) dead = 1
          
          if(newx > model@Init$MaxX) newx = model@Init$MaxX
          
          SimOut[ind,,time+1] = c(time+1, SimOut[ind, 2, time], newx, abs(decision+1.0), 0, dead)
          
          if(SimOut[ind, 6, time+1])  SimOut[ind, 6, (time+1):dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
        }
        
      }
      
    } ## Ind loop
  } ## time loop
  
  SimOut       <- SimOut[,-5,] 
  SimOut[,2,] <- SimOut[,2,]-1
  
  SimOut
}

pltNetwork <- function(simu, model, map = eaafMap) {
  
  sitesCrds <- model@Sites$crds
  dead      <- apply(simu[,5,], 1, sum, na.rm = T)
  trkS      <- simu[dead<1,2,]+1
  
  diagT <-   apply(trkS, 1, function(x) {
    tt <- as.data.frame(table(x[!is.na(x)]))
    m  <- merge(data.frame(Var1 = 1:nrow(sitesCrds)), tt, all.x = T)
    m[nrow(m),2] <- NA
    ifelse(is.na(m[,2]), 0, m[,2])
  })
  tmp.out   <- matrix(0, ncol = nrow(sitesCrds), nrow = nrow(sitesCrds)) 
  diag(tmp.out) <- apply(diagT, 1, sum)
  
  transT <- as.data.frame(table(apply(do.call("rbind", apply(trkS, 1, function(x) {
    tmp01 <- cbind(x[!is.na(x)][-sum(!is.na(x))], x[!is.na(x)][-1])
    as.data.frame(tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),])
  })), 1, function(y) paste(y, collapse = "_"))))
  trans <- cbind(t(apply(transT, 1, function(z) as.numeric(strsplit(z[1], "_")[[1]]))), transT[,2])
  
  
  lwd.f <- approxfun(c(0, max(trans[,3])), c(0.1, 5), rule = 3)
  
  opar <- par(mfrow = c(1,1))
  plot(map$map, col = "grey90", border = "grey50")
  pr <- st_as_sf(model@Sites$crds %>% as_data_frame(), coords = c("X", "Y")) %>% st_set_crs(4326) %>% st_transform(st_crs(map$map)$input) %>% st_coordinates()
  segments(pr[trans[,1],1], pr[trans[,1],2], pr[trans[,2],1], pr[trans[,2],2], lwd = lwd.f(trans[,3]), col = "grey20")
  
  site01 <- cbind(sitesCrds, d = diag(tmp.out))
  cex.f <- approxfun(x = c(1, max(site01[-1,3])), y = c(0.1, 15), rule = 3)
  
  cex <- cex.f(site01[,3])
  cex[model@Sites$start] <- 2
  
  pch = rep(21, nrow(site01))
  pch[model@Sites$start] <- 23
  
  pts <- st_as_sf(site01[,1:2] %>% as_data_frame(), coords = c("X", "Y")) %>% st_set_crs(4326) %>% st_transform(st_crs(map$map)$input) %>% st_coordinates()
  
  points(pts, cex = cex, pch = pch, col  = "grey10",
         bg = c("white", rep(adjustcolor("orange", alpha.f = 0.55), nrow(site01)-1)), lwd = 1.4)
  points(pr[1,1], pr[1,2], pch = 21, bg = adjustcolor("white", alpha.f = 0.3))
  par(opar)
  
  cat(sum(dead<1)/length(dead))
  
}

simNetwork <- function(simu, model, crds_ind, map = eaafMap, plot = FALSE) {
  
  dead_orNA <- apply(cbind(apply(simu[,5,], 1, sum, na.rm = T)==0, 
                           apply(simu[,2,], 1, function(x) any(x[!is.na(x)]==max(simu[,2,], na.rm = T)))), 1, 
                     function(x) all(x))
  trkS      <- simu[dead_orNA,2,]+1
  
  sites <- tibble(index = model@Sites$index, ts = apply(trkS, 1, function(y) {
    seqTmp <- y[!is.na(y)][max(which(y[!is.na(y)]==1)):min(which(y[!is.na(y)]==max(y[!is.na(y)])))]
    tibble(sites = seqTmp) %>% group_by(sites) %>% summarise(count = n()) %>%
      full_join(tibble(sites = 1:length(model@Sites$index)), by = "sites") %>% arrange(sites) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>% pull(count)
  }) %>% rowSums()) %>% group_by(index) %>% summarise(ts = sum(ts)) %>%
    full_join(tibble(index = 1:nrow(crds_ind)), by = "index") %>%
    arrange(index) %>% mutate(ts = ifelse(is.na(ts), 0, ts)) %>% pull(ts)
  
  if(plot) {
    opar <- par(mar = c(0,0,0,0))
    plot(eaafMap$grid, col = "grey80", lty = 3)
    plot(eaafMap$map, col = "grey90", border = "grey60", add = T)
    plot(eaafMap$bbox, add = T, border = "grey60")
  }
  
  ### Relative Site Use
  relSite <- do.call("rbind", lapply(1:nrow(trkS), function(y) {
    tmp <- trkS[y,]
    tibble(ts = model@Sites$index[tmp[!is.na(tmp)][(max(which(tmp[!is.na(tmp)]==1))+1):(min(which(tmp[!is.na(tmp)]==max(tmp[!is.na(tmp)])))-1)]]) %>%
      group_by(ts) %>% summarise(days = n()) %>% arrange(desc(days)) %>% mutate(site = 1:nrow(.)) %>%
      dplyr::select(site, ts, days)
  }))
  migrPhen <- do.call("rbind", lapply(1:nrow(trkS), function(y) {
    tmp <- trkS[y,]
    tibble(dep = (max(which(!is.na(tmp) & tmp==min(tmp, na.rm=T))+1))+model@Init$MinT,
           arr = (min(which(!is.na(tmp) & tmp==max(tmp, na.rm=T))))+model@Init$MinT)
  })) %>% mutate(dur = arr-dep)
  
  
  
  transTab <- apply(trkS, 1, function(y) {
    tmp00 <- model@Sites$index[y[!is.na(y)][max(which(y[!is.na(y)]==1)):min(which(y[!is.na(y)]==max(y[!is.na(y)])))]]
    tmp01 <- cbind(tmp00[!is.na(tmp00)][-sum(!is.na(tmp00))], tmp00[!is.na(tmp00)][-1])
    as.data.frame(tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),])
  }) %>% do.call("rbind",.) %>% setNames(c("a", "b")) %>% group_by(a, b) %>% summarise(count = n())
  
  trans <- matrix(0, ncol = length(sites), nrow = length(sites))
  trans[transTab %>% dplyr::select(a, b) %>% as.matrix()] <- transTab %>% pull(count)
  
  
  if(plot) {
    transT   <- tibble(a = rep(1:nrow(trans), nrow(trans)), b = rep(1:nrow(trans), each = nrow(trans)), count = c(trans)) %>% filter(count > 0)
    trans_sf <- lapply(1:nrow(transT), function(t) st_linestring(crds_ind[unlist(transT[t,1:2]),])) %>% st_sfc() %>% st_set_crs(st_crs(eaafMap$map)) %>%
      st_sf() %>% mutate(trans = transT[,3]) %>% dplyr::select(trans, geometry)
    plot(trans_sf %>% st_geometry(), add = T, lwd = approx(c(1, 100), c(0, 5), transT$count)$y, col = adjustcolor("grey30", alpha.f = 0.7))
    points(crds_ind, cex = approx(range(sites), c(0, 8), sites)$y,
           pch = 21, lwd = 1.5, col = "grey40", bg = adjustcolor("orange", alpha.f = 0.6))
  }
  
  diag(trans) <- sites
  
  ### fuelling
  trkS  <- simu[,2,]+1
  trkX  <- simu[,3,]
  trkD  <- simu[,5,]
  
  fueling <- lapply(1:dim(simu)[1], function(i) {
    if(all(trkD[i,]==0)){
      cbind(i, model@Init$MinT:model@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
    } else NULL
  }) %>% Reduce("rbind", .)
  
  
  fitness <- lapply(1:dim(simu)[1], function(x) {
    if(all(simu[x,5,]==0)) {
      t <- min(which(simu[x,2,]==max(simu[x,2,], na.rm = T)))
      matrix(c(t, simu[x,3,t]), ncol = 2)
    } else NULL
  }) %>% do.call("rbind",.) %>% apply(., 1, function(f) model@Results$FitnessMatrix[f[1], dim(model@Results$FitnessMatrix)[2],f[2]])
  
  
  list(trans, relSite, migrPhen, fueling, fitness)
  
} 

condProfile <- function(simu, model) {
  trkS  <- simu[,2,]+1
  trkX  <- simu[,3,]
  trkD  <- simu[,5,]
  
  cls <- rainbow(max(trkS, na.rm = T))[sample(1:max(trkS, na.rm = T))]
  
  opar <- par(mfrow = c(1,1))
  matplot(model@Init$MinT:model@Init$MaxT, t(trkX), type = "n", 
          xaxt = "n", xlab  = "", ylab = "Body condition", las = 1)
  for(i in 1:dim(simu)[1]) {
    if(all(trkD[i,]==0)){
      tmp <- cbind(model@Init$MinT:model@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
      points(tmp[,1], tmp[,3], type = "l", col = "grey90")
      points(tmp[,1], tmp[,3], pch = 16, cex = 0.5, col = cls[tmp[,2]])
    }
  }
  axis(1, at = seq(model@Init$MinT,model@Init$MaxT, by = 15), labels = format(as.POSIXct("2012-01-01")+seq(model@Init$MinT,model@Init$MaxT, by = 15)*24*60*60, "%b-%d"))
  par(opar)
}

