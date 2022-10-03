## clean and commented
## incorporate new model for spring turnover (6.3.2022)
dochlmod <- function(dat.merge, dat.merge2, modat0, fitout = NULL, runmod = T, varout = NULL) {

    require(mgcv)

#    dat.merge$site.id.c <- factor(dat.merge$site.id.c)

    ## parameters for criteria calculation
    lat.site <- 41.41
    lon.site <- -84.32
    elev.site <- 1350
    crittemp <- 19
    doc.site <- 7.2
    depth.site <- 13
    credint <- 0.8
    refuge.depth <- 0.3
    DOtarget <- 4

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
            stop()
        }

        ## so, in my first model, I could get a high r2 because
        ## lots of the variability was absorbed by the site specific
        ## estimate of tstart
        ## compute predictions for single samples to test model

        b <- apply(varout$b, 2, mean)
        tstart <- apply(varout$tstart, 2, mean)
        if (dim(varout$b)[2] == 4 ) {
            tstart.mod <- b[4]*dfsite$air4 + b[1] + b[2]*log(dfsite$Area) +
                b[3]*log(dfsite$Depth)
        }
        else {
            tstart.mod <- b[1] + b[2]*dfsite$Temp.sc
        }
        dftemp <- data.frame(sitenum = 1:length(tstart.mod),
                             tstart.mod = as.vector(tstart.mod),
			     tstart = tstart)

        dat.merge <- merge(dat.merge, dftemp, by = "sitenum")

        mud <- apply(varout$mud, 2, mean)
        d2 <- -exp(mud[1]) + mud[2]*dat.merge$chl.sc +
            mud[3]*dat.merge$hypdepth.mn.sc +
            mud[4]*dat.merge$docmean.sc
        print(names(dat.merge))
        pred <- dat.merge$maxDO + d2*(dat.merge$yday.x.sc - dat.merge$tstart)
        incvec <- dat.merge$domean >2

        plot(pred[incvec], dat.merge$domean[incvec])
        abline(0,1)
        rms <- function(x,y) sqrt(sum((x-y)^2)/length(x))
        print(rms(pred[incvec], dat.merge$domean[incvec]))


        b <- apply(varout$b, 2, mean)
	tstart <- apply(varout$tstart, 2, mean)
        tstart.raw <- tstart*sdsav["yday.x"] + mnsav["yday.x"]
        air4.raw <- dfsite$air4*sdsav["yday.x"] + mnsav["yday.x"]

        png(width = 3, height = 2.5, pointsize = 7, units = "in",
            res = 600, file = "airwater.png")
        par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
        plot(air4.raw, tstart.raw, xlab = "AIR4 (day)", ylab = "WATER4 (day)",
             pch = 21, col = "grey39", bg = "white")
        points(laketemp.sat$airtemp, laketemp.sat$day4, pch = 21, bg = "grey39", col = "black")
        dev.off()
        stop()
        area.mn <- mean(log(dfsite$Area))
	depth.mn <- mean(log(dfsite$Depth))


        laketemp.sat$day4.sc <- (laketemp.sat$day4 - mnsav["yday.x"])/sdsav["yday.x"]
        laketemp.sat$airtemp.sc <- (laketemp.sat$airtemp - mnsav["yday.x"])/sdsav["yday.x"]



        mod <- lm(day4.sc ~ airtemp.sc + log(area) + log(depth), data = laketemp.sat)
        print(summary(mod))

        stop()
        dfnew <- dftemp
        dfnew$air4 <- mean(dfnew$air4)
#        dfnew$Area <- exp(mean(log(dfnew$Area)))
        dfnew$Depth <- exp(mean(log(dfnew$Depth)))
        predout <- predict(mod, dfnew)
        plot(log(dftemp$Area), dftemp$tstart)
        iord <- order(dftemp$Area)
        lines(log(dftemp$Area)[iord], predout[iord])
        points(log(dftemp$Area)[dftemp$sat], dftemp$tstart[dftemp$sat],
               pch = 16, col = "red")

        stop()
        mod <- lm(day4.sc ~ airtemp.sc + log(area) + log(depth), data = laketemp.sat)
        print(summary(mod))
        print(b)

        dfsite$air4.raw <- dfsite$air4*sdsav["yday.x"] + mnsav["yday.x"]
        plot(dfsite$air4.raw, tstart.raw)
        points(laketemp.sat$airtemp, laketemp.sat$day4, pch = 16,
               col = "red", cex = 0.6)

        stop()
        par(mar = c(4,4,1,1), mfrow = c(1,2))

        plot(log(dfsite$Depth), tstart.raw, xlim = range(c(log(dfsite$Depth),
                                                          log(laketemp.sat$depth))))
        points(log(laketemp.sat$depth), laketemp.sat$day4,
               pch = 16, col = "red", cex = 0.6)
        plot(log(dfsite$Area), tstart.raw, xlim = range(c(log(dfsite$Area),
                                                          log(laketemp.sat$area))))
        points(log(laketemp.sat$area), laketemp.sat$day4,
               pch = 16, col = "red", cex = 0.6)
        stop()
        plot(dfsite$air4*sdsav["yday.x"] + mnsav["yday.x"],
             tstart*sdsav["yday.x"] + mnsav["yday.x"])
        stop()
#        sigt <- mean(varout$sigt)

        ## predict with uncertainty
        pred <- FALSE
        if (pred) {
            xnew <- seq(min(dfsite$air4), max(dfsite$air4), length = 40)
            predout <- matrix(NA, ncol =3, nrow = length(xnew))
            nsamp <- length(varout$sigt)
            for (i in 1:length(xnew)) {
                y <- varout$b[,1] + xnew[i] + varout$b[,2]*area.mn +
                    varout$b[,3]*depth.mn
                noise <- rnorm(nsamp, mean = 0, sd = sigt)
                predout[i,] <- quantile(y + noise, prob = c(0.05, 0.5, 0.95))
            }
            xnew.raw<- xnew*sdsav["yday.x"] + mnsav["yday.x"]
            predout.raw <- predout*sdsav["yday.x"] + mnsav["yday.x"]
        }



        ## look at temporal history of sites with similar VOD
        print(names(dat.merge))

        stop()

        predout <- dat.merge$maxDO + d2*(dat.merge$yday.x.sc -
                                             dat.merge$tstart)
        incvec <- dat.merge$domean > 2#& dat.merge$year != "17"

        deltDO <- (dat.merge$domean - dat.merge$maxDO)/(dat.merge$yday.x.sc -
                                                            dat.merge$tstart)


        docpick <- 1
        depthpick <- mean(dat.merge$hypdepth.mn.sc)

        d2 <- -exp(mud[1]) + mud[2]*dat.merge$chl.sc +
            mud[3]*depthpick +
                mud[4]*docpick

        selvec <- sqrt((docpick - dat.merge$docmean.sc)^2 +
                           (depthpick - dat.merge$hypdepth.mn.sc)^2) < 0.5
        par(mar = c(4,4,1,1), mfrow = c(1,2))
        plot(dat.merge$docmean.sc, dat.merge$hypdepth.mn.sc, col = "grey")
        points(docpick, depthpick, pch = 16, col = "black")
        points(dat.merge$docmean.sc[selvec],
               dat.merge$hypdepth.mn.sc[selvec], col = "red")

        plot(dat.merge$chl.sc[incvec], deltDO[incvec], col = "grey")
        points(dat.merge$chl.sc[incvec &selvec],
               deltDO[incvec & selvec], col = "red")
        iord <- order(dat.merge$chl.sc)
        lines(dat.merge$chl.sc[iord], d2[iord])
#        abline(0,1)
        stop()

        dev.new()
        par(mfrow = c(1,2), mgp = c(2.3,1,0), mar = c(4,4,1,1))
        plot(predout[incvec], dat.merge$domean[incvec])
        abline(0,1)
        stop()
        err <- dat.merge$domean[incvec] - predout[incvec]
        err[err > 12] <- NA
        dat.merge$err <- rep(NA, times = nrow(dat.merge))
        dat.merge$err[incvec] <- err
        plot(dat.merge$Lat, dat.merge$err)
        plot(dat.merge$Lon, dat.merge$err)
        plot(log(dat.merge$Area), dat.merge$err)
        mod <- lm(err ~ log(Area), data = dat.merge)
        print(summary(mod))
        plot(log(dat.merge$Depth), dat.merge$err)
        plot(dat.merge$air4, dat.merge$err)
        plot(log(dat.merge$LakeRatio), dat.merge$err)

        mod <- lm(err ~ log(Depth), data = dat.merge)
        print(summary(mod))
        abline(mod)

        print(sqrt(sum((predout[incvec] - dat.merge$domean[incvec])^2)/
                       sum(incvec)))
        stop()

        tstart <- apply(varout$tstart, 2, mean)
        d2n <- apply(varout$d2n, 2, mean)
        dev.new()
        par(mar = c(4,4,1,1), mfrow = c(3,3))

        nvis <- table(dat.merge$sitenum)
        ip <- as.numeric(names(nvis))[nvis >= 3]

        for (i in ip[10:18]) {
            incvec <- dat.merge$sitenum == i
            plot(dat.merge$yday.x.sc[incvec], dat.merge$domean[incvec])
            int0 <- dfsite$maxDO[i] - d2n[i]*tstart[i]
            abline(int0, d2n[i])
        }
        stop()

        tstart0 <- apply(varout$tstart0, 2, mean)
        dftemp <- data.frame(siteyrnum = 1:length(tstart0),
                             tstart0 = tstart0)

        df1 <- merge(df1, dftemp, by = "siteyrnum")

        d2 <- apply(varout$d2, 2, mean)
        dftemp <- data.frame(sitenum = 1:length(d2),
                             d2 = d2)

        df1 <- merge(df1, dftemp, by = "sitenum")

        df1$predout <- df1$mo.maxDO2  + df1$d2*(df1$yday.sc - df1$tstart0)

        dev.new()
        plot(df1$predout, df1$domean.all - df1$predout)
        abline(h=0)

        incvec <- df1$predout < 1
        print(df1[incvec,])

        rms <- function(x) sqrt(sum(x^2)/length(x))
        errsite <- tapply(df1$domean.all - df1$predout, df1$sitenum, rms)
        print(errsite)
        print(mean(errsite))

        dev.new()
        par(mar = c(4,4,1,1), mfrow = c(3,3), mgp = c(2.3,1,0))
        for (i in 1:9) {
            incvec <- df1$sitenum == i
            dftemp <- df1[incvec,]
            plot(dftemp$yday.sc, dftemp$domean.all, type = "n",
                 xlim = range(c(dftemp$yday.sc, dftemp$tstart0)),
                 ylim = range(c(dftemp$domean.all, dftemp$mo.maxDO2)))
            xnew <- seq(min(dftemp$yday.sc), max(dftemp$yday.sc),
                        length = 40)
            j0 <- unique(dftemp$siteyrnum)
            for (j in j0) {
                incvec <- dftemp$siteyrnum == j
                dftemp2 <- dftemp[incvec,]
                text(dftemp2$yday.sc, dftemp2$domean.all, lab = paste(j))
                dftemp3 <- unique.data.frame(dftemp2[, c("tstart0", "d2",
                                                        "mo.maxDO2")])
                predout <- dftemp3$mo.maxDO2 + dftemp3$d2*(xnew - dftemp3$tstart0)

                lines(xnew, predout)
                text(dftemp3$tstart0, dftemp3$mo.maxDO2, lab = paste(j),
                     col = "red")
            }
        }
        stop()

        ## plot median posterior prediction vs. observation
#        varout.p <- extract(fitout, pars = c("do0p", "dop"))
#        varout <- extract(fitout, pars = c("b", "mud"))
        dopmed <- apply(varout.p$dop, 2, median)
        do0pmed <- apply(varout.p$do0p, 2, median)
        dev.new()
        par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0),
            bty = "l")
        plot(dopmed, dftemp1$domean, xlab = "Predicted DO",
             ylab = "Observed DO", pch = 21, col = "grey39",
             bg = "white", main = "NLA")
        abline(0,1)
        plot(do0pmed, df1$domean.all, xlab = "Predicted DO",
             ylab = "Observed DO", pch = 21, col= "grey39",
             bg = "white", main = "MO")
        abline(0,1)

        ## Computations to derive a Chl criteria.
        ## Requires lake location to estimate stratification start
        ## date, release of temperature constraint in surface waters
        ##
        ## compute spring stratification data from air temp based on
        ## lake location and model parameters
        new.data <- data.frame(Lat = lat.site, Lon = lon.site,
                               Elev = elev.site)
        airtemp.raw<- predict(modairtemp, new.data)
        airtemp <- (airtemp.raw  - mnsav["Temp"])/sdsav["Temp"]
        ## compute estimate of spring strat day
        stratstart.sc <- varout$b[,1] +varout$b[,2]*as.vector(airtemp)

        ## get predicted surface water temperature by day for selected location
        ## only computing temps when temps start declining
        yday.all <- seq(204, 290, by = 1)
        new.data <- data.frame(yday.x = yday.all)
        new.data$Lat <- lat.site
        new.data$Lon <- lon.site
        new.data$Elev <- elev.site
        predout.all <- predict(mod.watertemp, new.data, se.fit = T)
        up <- predout.all$fit + predout.all$se.fit
        dn <- predout.all$fit - predout.all$se.fit
        endday <- approx(predout.all$fit,yday.all, crittemp)$y
        endday.up <- approx(up, yday.all, crittemp)$y
        endday.dn <- approx(dn, yday.all, crittemp)$y
        ## estimate of SD from half of +/- se.fit
        val <- c(endday, 0.5*(endday.up - endday.dn))

        endday.sc <- (endday -  mnsav["yday.x"])/sdsav["yday.x"]
        endday.se <- val[2]/sdsav["yday.x"]
        ## scaled number of days between start of stratification
        ## and temperature release in surface layer
        ## Distribution of values available from posterior coefs for
        ## stratstart.  Create a normal distribution for
        ## predicted temperature release based on mean = endday.sc
        ## and sd = endday.se
        nsamp <- length(stratstart.sc)
        drange <- rnorm(nsamp, mean = endday.sc, sd = endday.se) - stratstart.sc

        ## compute initial DO from location
        new.data <- data.frame(Lat = lat.site, Lon = lon.site,
                               Elev = elev.site)
        tempmin <- predict(modmintemp, new.data)
        maxDO <- DO.sat(as.vector(tempmin), elevation.m = elev.site)

        # convert other model inputs to scaled values
        doc.pick <- (log(doc.site) - mnsav["docmean"])/sdsav["docmean"]
        depth.pick <- (depth.site-mnsav["hypdepth.mn"])/
            sdsav["hypdepth.mn"]

        x.sc <- seq(min(dat.merge$chl.sc), max(dat.merge$chl.sc), length = 40)
        predout <- matrix(NA, ncol = 3, nrow = length(x.sc))

        for (k in 1:length(x.sc)) {
            y <- maxDO + (-exp(varout$mud[,1]) + varout$mud[,2]*x.sc[k]+
                               varout$mud[,3]*depth.pick +
                                   varout$mud[,4]*doc.pick)*drange

            predout[k,] <- quantile(y, prob = c(0.5*(1-credint), 0.5,
                                           1-0.5*(1-credint)))
        }

        ## compute depth-averaged DO target from DO threshold
        ## and refugia depth
        ## DO in surface layer based on critical temperature
        maxDO <- DO.sat(crittemp, elev.site)
        DOtarg.da <- 0.5*maxDO^2/(maxDO-DOtarget)*refuge.depth/depth.site
        DOtarg.da

        dev.new()
        par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
        xraw <- x.sc*sdsav["chl"] + mnsav["chl"]
        plot(xraw, predout[,2], type = "n", ylim = c(0, max(predout)),
             axes = F, xlab = expression(Chl~italic(a)~(mu*g/L)),
             ylab = expression(DO[m]~(mg/L)))
        logtick.exp(0.001, 10, c(1), c(F,F))
        axis(2)
        polygon(c(xraw, rev(xraw)),
                c(predout[,1], rev(predout[,3])), col = "grey", border = NA)
        lines(xraw, predout[,2])
        abline(h = DOtarg.da, lty = "dashed")
        chlcrit <- approx(predout[,1], xraw, DOtarg.da)$y
        cat("Chl criterion:", round(exp(chlcrit), digits = 1), "\n")


    }
}

## runmod variable set to T to run simulation and set to F to
##  run post processing.
#fitout <- dochlmod(dat.merge.cross, dat.merge.17, modat0, runmod = T)
#varout.temp <- extract(fitout, pars = c("mud", "b", "tstart"))

## post process
dochlmod(dat.merge.cross, dat.merge.17, modat0, fitout, runmod = F,
         varout = varout)
