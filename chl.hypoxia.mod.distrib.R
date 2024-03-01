## clean and commented
## incorporate new model for spring turnover (6.3.2022)
## edit scripts to include only code useful for distribution
dochlmod <- function(dat.merge, dat.merge2, modat0, fitout = NULL, runmod = T, varout = NULL) {

    require(mgcv)

#    dat.merge$site.id.c <- factor(dat.merge$site.id.c)

    ## parameters for criteria calculation
    lat.site <- 45
    lon.site <- -91
    elev.site <- 1775
    crittemp <- 17
    doc.site <- 5.2
    depth.site <- 16

    credint <- 75
    refuge.depth <- 0.3
    DOtarget <- 5
    area.site <- 850
    tdepth.site <- 1

    ## function to compute saturation DO concentration
    DO.sat <- function (temp.C, elevation.m = NULL) {
        tk <- 273.15 + temp.C
        bar.press.atm <- exp(-9.80665*0.0289644*elevation.m/(8.31447*288.15))
        A1 <- -139.34411
        A2 <- 157570.1
        A3 <- 66423080
        A4 <- 1.2438e+10
        A5 <- 862194900000
        DO <- exp(A1 + (A2/tk) - (A3/(tk^2)) + (A4/(tk^3)) -
            (A5/(tk^4)))

        theta <- 0.000975 - temp.C * 1.426e-05 + (temp.C^2) * 6.436e-08
        u <- exp(11.8571 - (3840.7/tk) - (216961/(tk^2)))
        Fp <- ((bar.press.atm - u) * (1-(theta * bar.press.atm)))/
            ((1 - u) * (1 - theta))
        Cp <- DO *  Fp
        return(Cp)
    }

    ## Prepare raw data
    # drop littoral duplicates
    dat.merge <- subset(dat.merge, sample.type == "MICX")
    dat.merge.sav <- dat.merge
    # select variables used in the model
    dat.merge <- dat.merge[, c("site.id.c", "doc.result", "therm.sav", "chl",
                               "domean",
                               "do.hyp.fail", "do.hyp.pass",
                               "yday.x",  "tmax",
                               "index.lat.dd", "index.lon.dd",
                               "index.site.depth",
                               "elevation", "area.ha",  "year")]
    dat.merge$nhyp <- dat.merge$do.hyp.fail + dat.merge$do.hyp.pass
    dat.merge <- dat.merge[, c("site.id.c", "doc.result", "therm.sav", "chl",
                               "domean", "nhyp",
                               "yday.x",  "tmax",
                               "index.lat.dd", "index.lon.dd",
                               "index.site.depth",
                               "elevation", "area.ha", "year")]


    dat.merge2 <- dat.merge2[, c("unique.id", "doc", "therm.sav", "chla",
                                 "domean","nhyp",
                                 "yday", "tmax",
                                 "lat.dd83", "lon.dd83",
                                 "index.site.depth",
                                 "elevation", "area.ha")]
    dat.merge2$year <- "17"
    names(dat.merge2) <- names(dat.merge)
    dat.merge <- rbind(dat.merge, dat.merge2)

    sitemean <- function(x, site, varname,dolog = F) {
        if (! dolog) {
            xm <- tapply(x, site, mean, na.rm = T)
        }
        else {
            xm <- tapply(x, site, function(x) mean(log(x), na.rm = T))
        }
        return(xm)
    }
    ## compute mean values of lake characteristics and merge
    ##"tmean.pt", "tmin.pt",
    varlist <- c("doc.result", "index.lat.dd",
                 "index.lon.dd", "elevation", "area.ha", "index.site.depth")
    matout <- matrix(NA, ncol = length(varlist),
                     nrow = length(levels(dat.merge$site.id.c)))
    dolog <- c(T, rep(F, times = 7))
    for (i in 1:length(varlist)) {
        matout[,i] <- sitemean(dat.merge[, varlist[i]],
                               dat.merge$site.id.c, dolog = dolog[i])
    }
    dat.mean <- data.frame(levels(dat.merge$site.id.c), matout)
    names(dat.mean) <- c("site.id.c", "docmean", "Lat",
                         "Lon", "Elev", "Area", "Depth")

    dat.mean$LakeRatio <- (dat.mean$Area*10000)^0.25/dat.mean$Depth
    ## convert area to km2 to match dimensions of the priors
    ## for turnover date
    dat.mean$Area <- dat.mean$Area*0.01

    interpairtemp <- function(xlat, xlon, pdat.all) {
        lat <- pdat.all$lat
        lon <- pdat.all$lon
        prismdat <- pdat.all$prismdat
        airtemp <- rep(NA, times = length(xlat))
        for (i in 1:length(xlat)) {
            ilat <- which.min(abs(xlat[i]-lat))
            ilon <- which.min(abs(xlon[i]-lon))
            airtemp[i] <- prismdat[ilat, ilon]
        }
        return(airtemp)
    }

    ## merge in mean data
    dat.merge<-merge(dat.merge, dat.mean,by = "site.id.c")

    ## compute mean depth below thermocline
    dat.merge$hyp.depth <-dat.merge$Depth - dat.merge$therm.sav
    incvec <- dat.merge$hyp.depth < 0
    incvec[is.na(incvec)] <- F
    dat.merge$hyp.depth[incvec] <- NA # drop sites with negative hyp depth
    hypdepth.mn <- tapply(dat.merge$hyp.depth, dat.merge$site.id.c, mean,
                          na.rm = T)
    dftemp <- data.frame(site.id.c = names(hypdepth.mn),
                         hypdepth.mn = as.vector(hypdepth.mn))
    dat.merge <- merge(dat.merge, dftemp, by = "site.id.c")

    incvec <- ! is.na(dat.merge$domean) & ! is.na(dat.merge$hypdepth) &
        ! is.na(dat.merge$chl) & ! is.na(dat.merge$docmean)
    dat.merge <- dat.merge[incvec,]

    ## get beginning DO concentration at turnover as DO with T = 4
    dat.merge$maxDO <- DO.sat(4, elevation.m = dat.merge$Elev)

    ## select lakes that seasonally stratify based on lake ratio
    dat.merge <- subset(dat.merge, LakeRatio < 3)

    ## drop 7 partial profiles depth > 25 m and number of measurements < 16
    ## another temporary drop of hypolimnion samples
#    incvec <- (dat.merge$do.hyp.pass + dat.merge$do.hyp.fail < 16) &
#        dat.merge$hyp.depth > 25
#    dat.merge <- dat.merge[!incvec,]

    ## compute adjusted latitude for finding dimictic lakes
    lat0 <- seq(0, 90, by = 10)
    adj <- c(0.27, 0.31, 0.34, 0.39, 0.46, 0.54, 0.68, 0.89, 1.3, 2.4)
    dat.merge$adjloc <- approx(lat0, adj, dat.merge$Lat)$y
    dat.merge$adjlat <- dat.merge$Lat + dat.merge$adjloc*dat.merge$Elev/100

    ## using air4 now to define dimictic gives about 80 more lakes
##    dat.merge <- subset(dat.merge, adjlat > 40) # drop monomictic lakes


    ## get mean air temperature
    load("pdat.all.rda")
    dat.merge$Temp <-interpairtemp(dat.merge$Lat, dat.merge$Lon,
                                   pdat.all)

    ## standardize variables for model
    dat.merge$chl <- log(dat.merge$chl)
    varlist <- c("hypdepth.mn", "docmean", "chl", "yday.x", "Temp")
    mnsav <- apply(dat.merge[, varlist], 2, mean)
    sdsav <- apply(dat.merge[, varlist], 2, sd)

    for (i in varlist) {
        dat.merge[,paste(i, "sc", sep = ".")] <- (dat.merge[,i] - mnsav[i])/sdsav[i]
    }

    ## get day4 (day of the year in which mean air temp is 4
    dat.merge$air4 <- interpairtemp(dat.merge$Lat, dat.merge$Lon, airtemp.4)
    dat.merge$air4 <- (dat.merge$air4 - mnsav["yday.x"])/sdsav["yday.x"]

    ## drop sites with missing air4
    incvec <- is.na(dat.merge$air4)
    dat.merge <- dat.merge[!incvec,]

    ## select NLA sites with at least 2 profiles
    numvis <- table(dat.merge$site.id.c)
    idsav <- names(numvis)[numvis >= 2]
    incvec <- rep(F, times = nrow(dat.merge))
    for (i in idsav) incvec <- dat.merge$site.id.c == i | incvec
    dat.merge.sing <- dat.merge[!incvec,]
#    dat.merge <- dat.merge[incvec,]

    dat.merge$site.id.c <- factor(dat.merge$site.id.c)
    dat.merge$sitenum <- as.numeric(dat.merge$site.id.c)

    ## sitedata
    dfsite <- unique.data.frame(dat.merge[, c("sitenum", "hypdepth.mn.sc",
                                           "docmean.sc","maxDO",
                                              "Lat",  "air4",
                                              "Area", "Depth", "Temp.sc")])
    dfsite <- dfsite[order(dfsite$sitenum),]
    print(nrow(dfsite))

    ## split data into missing (DO < 2) or observed
    incvec <- dat.merge$domean < 2
    dftemp1 <- dat.merge[!incvec,]
    dftemp2 <- dat.merge[incvec,]

    df1 <- modat0

    df1$Elev <- df1$elev*0.3048
    df1$air4 <- interpairtemp(df1$Lat, df1$Long, airtemp.4)
    df1$air4 <- (df1$air4 - mnsav["yday.x"])/sdsav["yday.x"]
    df1$Temp <- interpairtemp(df1$Lat, df1$Long, pdat.all)
    df1$Temp.sc <- (df1$Temp - mnsav["Temp"])/sdsav["Temp"]

    ## convert area from acres to km2
    df1$resarea <- df1$resarea*0.004047

    ## standardize MO data by nla vals
    df1$chlmean.sc <- (df1$chlmean - mnsav["chl"])/sdsav["chl"]
    df1$docmean.sc <- (df1$docmean - mnsav["docmean"])/sdsav["docmean"]
    df1$yday.sc <- (df1$yday - mnsav["yday.x"])/sdsav["yday.x"]
    df1$depth.hm.sc<- (df1$depth.hm - mnsav["hypdepth.mn"])/sdsav["hypdepth.mn"]

    df1$mo.maxDO2 <- DO.sat(4, elevation.m = 250) #use mean elev of 250 for now

    ## drop two MO sites with non-monotonic changes
    incvec <- df1$id == "152"# | df1$id == "74"
    df1 <- df1[!incvec,]

    ## set up site ids
    df1$id <- factor(as.character(df1$id))
    df1$sitenum <- as.numeric(df1$id)
    df1$idyear <- factor(df1$idyear)
    df1$siteyrnum <- as.numeric(df1$idyear)

    lookup0 <- unique.data.frame(df1[, c("siteyrnum", "air4", "Temp.sc",
                                         "dmax.all", "resarea")])
    lookup0 <- lookup0[order(lookup0$siteyrnum),]

    dfsite0 <- unique.data.frame(df1[, c("sitenum", "depth.hm.sc", "chlmean.sc",
                                         "docmean.sc",  "mo.maxDO2",
                                         "air4", "dmax.all", "resarea")])
    dfsite0 <- dfsite0[order(dfsite0$sitenum),]
    print(dfsite0)

    datstan <- list(n = c(nrow(dftemp1), nrow(dftemp2), nrow(dat.merge)),
                    nsite = nrow(dfsite),
                    chl = dat.merge$chl.sc,
                    air4 = dfsite$air4,
                    area = log(dfsite$Area),
                    temp = dfsite$Temp.sc,
                    depthall = log(dfsite$Depth),
                    depth = dfsite$hypdepth.mn.sc,
                    doc = dfsite$docmean.sc,
                    dmax = dfsite$maxDO,
                    t1 = dftemp1$yday.x.sc,
                    t2 = dftemp2$yday.x.sc,
                    domean = dftemp1$domean,
                    sitenum = dat.merge$sitenum,
                    sitenum1 = dftemp1$sitenum,
                    sitenum2 = dftemp2$sitenum,
                    n0 = nrow(df1),
                    t0 = df1$yday.sc,
                    domean0 = df1$domean.all,
                    dmax0 = dfsite0$mo.maxDO2,
                    nsite0 = max(dfsite0$sitenum),
                    chl0 = dfsite0$chlmean.sc,
                    depth0 = dfsite0$depth.hm.sc,
                    air40 = lookup0$air4,
                    temp0 = lookup0$Temp.sc,
                    depthall0 = log(lookup0$dmax.all),
                    area0 = log(lookup0$resarea),
                    doc0 = dfsite0$docmean.sc,
                    sitenum0 = df1$sitenum,
                    nsiteyr0 = max(df1$siteyrnum),
                    siteyrnum0 = df1$siteyrnum)


    print(str(datstan))

    ## bayesian model code
    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 6
        options(mc.cores = nchains)

        fit <- stan(file = "hypmodel.new.stan",  data = datstan,
                    iter = 900, chains = nchains, warmup = 300,
                    control = list(max_treedepth = 15,
                                   adapt_delta = 0.9), thin = 2)

        return(fit)
    }
    else {

        outputforshiny <- TRUE
        if (outputforshiny) {
            ## prepare data file for shiny that includes chlsite and tstart
            tstart <- apply(varout$tstart, 2, mean)
            chlsite <- apply(varout$chlsite, 2, mean)
            dftemp <- data.frame(sitenum = 1:length(tstart), tstart = tstart,
                                 chlsite = chlsite)
            ## make sure everything in sent to shiny app is unscaled
            dftemp$tstart <- dftemp$tstart*sdsav["yday.x"] + mnsav["yday.x"]
            dftemp$chlsite <- dftemp$chlsite*sdsav["chl"] + mnsav["chl"]
            print(nrow(dat.merge))
            dat.merge <- merge(dat.merge, dftemp, by = "sitenum")
            print(nrow(dat.merge))
            print(names(dat.merge))
            natdodat <- dat.merge[, c("site.id.c", "yday.x", "year",
                                      "Elev", "Area", "Lat", "Depth",
                                      "Lon", "LakeRatio", "domean",
                                      "tstart", "chlsite", "docmean",
                                      "hypdepth.mn","air4", "maxDO")]
            natdodat$air4 <- natdodat$air4*sdsav["yday.x"] +mnsav["yday.x"]

            print(summary(natdodat))
            save(natdodat, mnsav, sdsav, file = "natdodat.2017.rda")

            ## drop big variables from varout
            ip <- which("tstart" == names(varout))
            varout <- varout[-ip]
            ip <- which("chlsite" == names(varout))
            varout <- varout[-ip]
            save(varout, file = "varout.2017.rda")
#            return()
        }

        ## Computations to derive a Chl criteria.

        ## find temperature release day, subset temperature
        ## prediction to after the hottest day
        ## (wider range of days modeled initially to provide
        ## a meaningful plot
        ## load water temperature mode
        load("mod.watertemp.2017.rda")
        new.data <- data.frame(yday.x = 122:290)
        new.data$Lat = lat.site
        new.data$Lon = lon.site
        new.data$elevation = elev.site
        new.data$year = "12"
        predtemp <- predict(mod.watertemp, new.data,se.fit = T)
        ydays <- seq(122, 290, by = 1)
        n <- length(predtemp$fit)
        ip <- 83:n
        daycrit <- approx(predtemp$fit[ip],ydays[ip], crittemp)$y
        daycrit.up <- approx(predtemp$fit[ip] + predtemp$se.fit[ip],
                             ydays[ip], crittemp)$y
        daycrit.dn <- approx(predtemp$fit[ip] - predtemp$se.fit[ip],
                             ydays[ip],  crittemp)$y
        endday <- c(daycrit, 0.5*(daycrit.up - daycrit.dn))
        cat("Release day:, ",endday[1], "\n")

        scalevar <- function(x, mnsav, sdsav, forward = T) {
            ifelse (forward, y <- (x-mnsav)/sdsav, y <- x*sdsav + mnsav)
            return(y)
        }

        ## get starting day (when temperature reaches 4)
        ## uses model parameters
        load("airtemp.4.lores.rda")
        ilat <- which.min(abs(lat.site-airtemp.4$lat))
        ilon <- which.min(abs(lon.site-airtemp.4$lon))
        air4pt <- airtemp.4$prismdat[ilat, ilon]
        air4pt.sc <- scalevar(air4pt, mnsav["yday.x"],
                              sdsav["yday.x"])
        ydaystart <- scalevar(varout$b[,1] +
                                  varout$b[,4]*air4pt.sc +
                                      varout$b[,2]*log(area.site) +
                                          varout$b[,3]*log(depth.site),
                              mnsav["yday.x"],
                              sdsav["yday.x"], FALSE)
        cat("Initial day", mean(ydaystart), "\n")

        airtemp.raw <- air4pt
        ## get probability distribution for the length of time
        ## for DO depletion (endday - ydaystart)
        drange <- rnorm(length(ydaystart), mean = endday[1],
                               sd = endday[2]) - ydaystart
        drange <- as.vector(drange)

        ## predict mean DO as a function of chl and site parameters
        ## that are provided
        x <- seq(min(natdodat$chlsite),
                 max(natdodat$chlsite), length = 40)
        x.sc <- scale(x, mnsav["chl"], sdsav["chl"])

        predout <- matrix(NA, ncol = 3, nrow = length(x))
        maxDO <- DO.sat(4, elevation.m = elev.site)
        # convert input to scaled values
        doc.pick <- scalevar(log(doc.site),
                             mnsav["docmean"],
                             sdsav["docmean"])
        depth.pick <- scalevar(depth.site - tdepth.site,
                               mnsav["hypdepth.mn"],
                               sdsav["hypdepth.mn"])
        maxDO <- DO.sat(4, elevation.m=elev.site)

        for (k in 1:length(x)) {
            y <- maxDO + (-exp(varout$mud[,1]) +
                              varout$mud[,2]*x.sc[k]+
                                  varout$mud[,3]*depth.pick +
                                      varout$mud[,4]*doc.pick)*drange/
                                          sdsav["yday.x"]
            predout[k,] <- quantile(y, prob = c(0.5*(1-credint/100), 0.5,
                                           1-0.5*(1-credint/100)))
        }

        ## find mean DO at the moment of temperature release
        ## profile based on surface temperature right before
        ## release from temperature constraint
        DOtarg <- 0.5*maxDO^2/(maxDO-DOtarget)*refuge.depth/
            (depth.site - tdepth.site)
        if (maxDO < DOtarg) DOtarg <- maxDO

        critlo <- approx(predout[,1],x,DOtarg)$y
        cat("Chl criterion:", exp(critlo), "\n")
        return()

    }
}

## runmod variable set to T to run simulation and set to F to
##  run post processing.
load("pdat.all.rda")
modat0 <- read.csv("modat0.csv")
load("airtemp.4.lores.rda")
load("dat.merge.17.rda")
load("dat.merge.cross.rda")

fitout <- dochlmod(dat.merge.cross, dat.merge.17, modat0, runmod = T)
varout.temp <- extract(fitout, pars = c("mud", "b", "tstart", "chlsite"))

## post process
dochlmod(dat.merge.cross, dat.merge.17, modat0, fitout, runmod = F,
         varout = varout.temp)
