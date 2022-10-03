## 12.9.2019
## Cleaned and commented
## Inputs:
## dat.merge.all: master data file
## 5.5.2021: Add in 2017 data: df2: 2017 zoop and phyt data
## 11.21.2021: running again with revised zoop data
## Factor to account for different years of the survey included
## More than 4 depth categories gives non-monotonic change in criteria across
## the groups, which indicates uncertainty due to lack of data
chlzoopmodel <- function(df1, df2, dftempmax, varout = NULL, runmod = T) {

    require(rstan)
    require(mapproj)
    require(maps)
    require(mgcv)
    nchains <- 3          # number of chains for mcmc simulation

    ## specify parameters here to calculate criteria
    credint <- 0.95       # credible interval for criteria derivation
    targ <- 0.0              # targeted threshold
    j0 <- 1:4            # select depth class for output (1-3)

    ## combine 2017 data with 2012 (2007 data will drop out
    ## because no phytoplankton data are available
    df1 <- df1[, c("site.id.c", "visit.no", "year",
                   "totl.bio", "chl",  "taxonomist", "biov.tot",
                   "index.site.depth", "index.lon.dd", "index.lat.dd")]
    df1$year <- as.character(df1$year)
    df2$year <- "17"
    df2 <- df2[, c("unique.id", "visit.no", "year",
                            "totl.bio", "chla", "analyst", "biov",
                            "index.site.depth", "lon.dd83", "lat.dd83")]
    names(df2)  <- names(df1)
    df1 <- rbind(df1, df2)

    ## merge in seasonal max lake temperature
    df1 <- merge(df1, dftempmax, by = "site.id.c")

    ## drop one sample from analyst XX in 2017
    incvec <- df1$taxonomist == "XX"
    incvec[is.na(incvec)] <- F
    print(sum(incvec))
    df1 <- df1[!incvec,]
    df1$taxonomist <- factor(df1$taxonomist)

    ## drop missing chl and zero chl
    incvec <- df1$chl > 0
    incvec[is.na(incvec)] <- F
    df1 <- df1[incvec,]
    df1$logchl <- log(df1$chl)

    ## calculate mean depth by site and then drop lakes with depth < 3
    depth.mn <- tapply(df1$index.site.depth, df1$site.id.c, mean, na.rm = T)
    lat.mn <- tapply(df1$index.lat.dd, df1$site.id.c, mean, na.rm = T)
    lon.mn <- tapply(df1$index.lon.dd, df1$site.id.c, mean, na.rm = T)
    dftemp <- data.frame(site.id.c = names(depth.mn),
                         depth.mn = as.vector(depth.mn),
                         lat.mn = as.vector(lat.mn),
                         lon.mn = as.vector(lon.mn))
    df1 <- merge(df1, dftemp, by = "site.id.c")
    incvec <- df1$depth.mn > 3
    incvec[is.na(incvec)] <- F
    df1 <- df1[incvec,]

    df1 <- na.omit(df1[, c("site.id.c", "visit.no", "totl.bio",
                           "depth.mn", "logchl", "biov.tot",
                           "taxonomist", "year", "tempmax", "lat.mn",
                           "lon.mn")])

    print(summary(df1))

    ## log transform and center phytoplankton variables
    df1$biov.log <- log(df1$biov.tot)
    mn.biov <- mean(df1$biov.log)
    df1$biov.log <- df1$biov.log - mn.biov

    mn.chl <- mean(df1$logchl)
    df1$logchl <- df1$logchl - mn.chl
    mnval <- c(mn.biov, mn.chl)
    names(mnval) <- c("biov", "chl")
    print(mnval)
    mnval <- c(0, mnval) # zbiomass is not centered in this model
    ## so adding a 0 mnval so this still works in the shiny app
    names(mnval) <- c("zbiomass", "biov", "chl")
    print(mnval)

    ## set sitenums
    df1$site.id.c <- factor(df1$site.id.c)
    df1$sitenum <- as.numeric(df1$site.id.c)
    dflookup <- unique.data.frame(df1[, c("site.id.c", "sitenum")])

    ## set up taxonomist factor
    df1$taxonomist <- factor(df1$taxonomist)
    df1$taxonnum <- as.numeric(df1$taxonomist)

    ## compute the mean zbiomass and zdensity at each
    ## site.id.c and visit no. Cleans out duplicate entries of zoop
    ## at littoral and index samples while
    ## retaining samples in which only littoral zone phytoplankton is
    ## available with an index zoop sample
    df1$sitevis <- paste(df1$site.id.c, df1$visit.no, df1$year, sep = "--")
    zbiomass.mn <- tapply(df1$totl.bio, df1$sitevis, mean)
    df2 <- data.frame(sitevis = names(zbiomass.mn),
                      zbiomass = as.vector(zbiomass.mn))

    ## recover site.id.c and visit.no in df2
    w <- regexpr("--", as.character(df2$sitevis))
    df2$site.id.c <- substring(as.character(df2$sitevis), 1, w-1)
    nch0 <- nchar(as.character(df2$sitevis))
    df2$year <- substring(as.character(df2$sitevis),
                                     nch0-1, nch0)
    df2$visit.no <- substring(as.character(df2$sitevis),
                              nch0-4, nch0-4)
    df2 <- df2[, c("site.id.c", "visit.no", "year", "zbiomass")]

    ## set up site-yr as identifier
    df2$siteyr <- paste(df2$site.id.c, df2$year, sep = "--")
    df2$siteyrnum <- as.numeric(factor(df2$siteyr))

    ## merge sitenums back into zoop data
    df2 <- merge(df2, dflookup[, c("site.id.c", "sitenum")], by = "site.id.c")

    ## assign index to year
    df2$year <- factor(df2$year)
    print(table(df2$year))

    ## define 4 temperature classes and merge into df2
    dfdeep <- unique.data.frame(df1[, c("sitenum", "depth.mn", "tempmax",
                                        "lon.mn", "lat.mn")])
    cutp <- quantile(dfdeep$tempmax, prob = seq(0,1, length = 5))
    cutf2 <- cut(dfdeep$tempmax, cutp, include.lowest = T)

    print("Classes and N within each")
    dfdeep$grp <- cutf2
    print(table(dfdeep$grp))

    dfdeep$depthnum <- as.numeric(dfdeep$grp)
    dfdeep <- dfdeep[order(dfdeep$sitenum),]

    df2 <- merge(df2, dfdeep[, c("sitenum", "grp", "depthnum",
                                 "lon.mn", "lat.mn")], by = "sitenum")
    df2$yeargrp <- interaction(df2$grp, df2$year)
    df2$yeargrpnum <- as.numeric(df2$yeargrp)
    print(table(df2$yeargrp))

    ## output data for mapping temperature zones
    dftempout <- unique.data.frame(df2[, c("sitenum", "grp", "lon.mn", "lat.mn")])
    save(dftempout, file = "dftempout.rda")

    ## map classes
    require(maps)
    require(mapproj)
    map("state", proj = "albers",par = c(30,40))
    pout <- mapproject(df2$lon.mn, df2$lat.mn, proj = "")
    col0 <- brewer.pal(4, "PuBuGn")
#    col0 <- c("blue", "green4", "orange", "red")
    for (i in 1:4) {
        incvec <- df2$grp == levels(df2$grp)[i]
        points(pout$x[incvec], pout$y[incvec], pch = 21,
               col = "grey", bg = col0[i])
    }
    stop()

    dflookup2 <- unique.data.frame(df2[, c("siteyrnum", "sitenum",
                                           "depthnum", "yeargrpnum")])
    dflookup2 <- dflookup2[order(dflookup2$siteyrnum),]

    zoopdat <- df2
    zoopdat$zbiomass <- log(zoopdat$zbiomass)
    ## save files for shiny
    save(zoopdat, mnval, cutp, file = "shinydat.2017.rda")

    xnew <- seq(from = min(df1$logchl), to = max(df1$logchl), length = 60)

    datstan <- list(nsamp = nrow(df1),
                    nsite = max(df1$sitenum),
                    sitenum = df1$sitenum,
                    biov = df1$biov.log,
                    chl = df1$logchl,
                    nz = nrow(df2),
                    zoopb = log(df2$zbiomass),
                    sitenumz = df2$siteyrnum,
                    ngrp = max(dfdeep$depthnum),
                    depthf = dflookup2$depthnum,
                    taxonnum = df1$taxonnum,
                    ntaxon = max(df1$taxonnum),
                    nsiteyrnum = max(df2$siteyrnum),
                    ix = dflookup2$sitenum,
                    yeargnum = dflookup2$yeargrpnum)

    print(str(datstan))

    if (runmod) {

        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)

        fit <- stan(file = "zoopmod.stan", data = datstan,
                    chains = nchains, iter = 1400, warmup = 400,
                    control = list(adapt_delta = 0.95,
                        max_treedepth = 12), thin = 2)
        return(fit)
    }

    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    ## compute mean estimated site biov
    chl.site.obs <- tapply(df1$logchl, df1$sitenum, mean)
    biov.site.obs <- tapply(df1$biov.log, df1$sitenum, mean)

    dobiom <- F  # set to true to plot Chl and biov comparisons (Fig 5)
    if (dobiom) {
#        png(width = 4.5, height = 2.25, pointsize = 7, units = "in", res = 600,
#            file = "chl.p.png")
        dev.new()
        biov.site <- apply(varout$biov_site, 2, mean)
        par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
        plot(biov.site + mnval["biov"] + log(1.e-6),
             biov.site.obs + mnval["biov"] + log(1.e-6),
             pch = 21, col = "grey", bg = "white", axes = F,
         xlab = expression(Mean~biovolume~(mm^3/L)),
             ylab = expression(Measured~biovolume~(mm^3/L)))
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        abline(0,1)
        plot(biov.site + mnval["biov"] + log(1.e-6),
             chl.site.obs + mnval["chl"],
             pch = 21, col = "grey", bg = "white", axes = F,
             xlab = expression(Mean~biovolume~(mm^3/L)),
             ylab = expression(Chl~italic(a)~(mu*g/L)))
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        abline(-mnval["biov"] + mnval["chl"] - log(1.e-6), 1)
        stop()
        dev.off()
    }

    b <- apply(varout$b, 2, median)
    cp <- apply(varout$cp,2, median)

    biov.site <- apply(varout$biov_site, 2, mean)
    dftemp <- data.frame(sitenum = 1:length(biov.site),
                         biov.site = biov.site)
    dftemp <- merge(dftemp, df2, by = "sitenum")

    nit <- nrow(varout$cp)

    # compute binned data for each group
    ngrp <- max(df2$depthnum)
    bout1 <- as.list(rep(NA, times = ngrp))
    bout2 <- as.list(rep(NA, times = ngrp))
    print(mnval)
    for (j in 1:ngrp) {
        incvec <- dftemp$depthnum == j & dftemp$year == "12"
        bout1[[j]] <- binv(dftemp$biov.site[incvec] + mnval["chl"],
                           log(dftemp$zbiomass)[incvec], 10) # average 5 measurements for each bin
        incvec <- dftemp$depthnum == j & dftemp$year == "17"
        bout2[[j]] <- binv(dftemp$biov.site[incvec] + mnval["chl"],
                           log(dftemp$zbiomass)[incvec], 10) # average 5 measurements for each bin

    }

    xlim1 <- range(sapply(bout1, function(x) range(x$xb)))
    ylim1 <- range(sapply(bout1, function(x) range(x$yb)))


#    png(width = 6.5*2/3, height = 2, pointsize = 8, units = "in", res = 600,
#        file = "zoop.new.all.png")
#    dev.new()  # Figure 6 in document
#    pdf(file = "plotout2.pdf", width = 9, height = 6, pointsize = 10)

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(4,2), mgp = c(2.3,1,0))

    ## average f0 for two years into one
    ## reformat varout format to match format in 2012 only model
    fnew <- rep(NA, times = 3*ngrp*nrow(varout$f0))
    dim(fnew) <- c(nrow(varout$f0),ngrp,3)
    for (j in 1:ngrp) {
        fnew[,j,1] <- 0.5*(varout$f0[,j] + varout$f0[,j+ngrp])
        fnew[,j,2] <- varout$f[,j,1]
        fnew[,j,3] <- varout$f[,j,2]
    }
    print(str(fnew))

    ## output varout.zoop after replacing f with
    ## averaged f0 and removing biov.site
    varout.zoop <- varout
    print(names(varout.zoop))
    ip <- which("biov_site" == names(varout.zoop))
    varout.zoop <- varout.zoop[-ip]
    print(names(varout.zoop))
    varout.zoop$f <- fnew
    ip <- which("f0" == names(varout.zoop))
    varout.zoop <- varout.zoop[-ip]
    print(names(varout.zoop))
#    save(varout.zoop, biov.site, file = "varout.zoop.2017.rda")
#    stop()

    biomasspred <- matrix(NA, ncol = nit, nrow = length(xnew))
    biomasspred2 <- matrix(NA, ncol = nit, nrow = length(xnew))
    for (j in 1:ngrp) {
        # compute predictions
        for (i in 1:nrow(biomasspred)) {
            biomasspred[i,] <- fnew[,j,1] +
                fnew[,j,3]*xnew[i] +
                (fnew[,j,3] - fnew[,j,2])*exp(varout$b[,j])*
                    log(1+exp(-(xnew[i] - varout$cp[,j])/exp(varout$b[,j])))
        }


        zoopb.q <- apply(biomasspred, 1, quantile, prob = c(0.5*(1-credint),
                                               0.5, 1- 0.5*(1-credint)))
        ## estimated relationship for second year
        for (i in 1:nrow(biomasspred)) {
            ## f0 for second year is j+ngrp
            biomasspred2[i,] <- varout$f0[,j+ngrp] +
                varout$f[,j,2]*xnew[i] +
                (varout$f[,j,2] - varout$f[,j,1])*exp(varout$b[,j])*
                    log(1+exp(-(xnew[i] - varout$cp[,j])/exp(varout$b[,j])))
        }
        zoopb.q2 <- apply(biomasspred2, 1, quantile, prob = c(0.5*(1-credint),
                                          0.5, 1- 0.5*(1-credint)))

        ylim2 <- range(c(ylim1, zoopb.q, zoopb.q2))
        ynew.raw <- xnew +  mnval["chl"]
        xnew.raw <- xnew + mnval["chl"]

        xlim2 <- range(c(xlim1, ynew.raw))

        ## plot zoop biomass vs. Chl
        plot(bout1[[j]]$xb, bout1[[j]]$yb,xlim = xlim2,ylim = ylim2,
             xlab = expression(Chl~italic(a)~(mu*g/L)),
             ylab = expression(Zooplankton~biomass~(mu*g/L)), axes = F,
             pch = 21, col = "grey", bg = "white", type = "p")
        points(bout2[[j]]$xb, bout2[[j]]$yb, pch = 16, cex=0.5)
        logtick.exp(0.00001, 10, c(1,2), c(F,F))
        polygon(c(ynew.raw,rev(ynew.raw)), c(zoopb.q[1,],rev(zoopb.q[3,])),
                col = grey.t, border = NA)

        lines(ynew.raw,zoopb.q[2,])
        lines(ynew.raw,zoopb.q2[2,], lty = "dashed")
        incvec <- dftemp$depthnum == j
        dftemp2 <- dftemp[incvec,]
        mod <- gam(log(zbiomass) ~ s(biov.site, k = 3),
                   data = dftemp2)
        predout <- predict(mod)
        iord <- order(dftemp2$biov.site)
        lines(dftemp2$biov.site[iord] + mnval["chl"], predout[iord],
              col = "red")
#        print(summary(mod))

        ## plot slope vs.Chl
        ## slopes do not vary among years because year effect
        ## is strictly additive
        delt <- xnew.raw[2]-xnew.raw[1]
        zoopgrad <- matrix(NA, ncol = nit, nrow = length(xnew)-1)
        for (i in 1:nrow(zoopgrad)) {
            zoopgrad[i,] <- (biomasspred2[i+1,] - biomasspred2[i,])/delt
        }

        xmean <- 0.5*(xnew.raw[-1] + xnew.raw[-length(xnew.raw)])
        ymean <- 0.5*(ynew.raw[-1] + ynew.raw[-length(ynew.raw)])
        zoopgq <- apply(zoopgrad, 1, quantile,prob = c(0.5*(1-credint),
                                                  0.5, 1-0.5*(1-credint)))

        plot(ymean, zoopgq[2,],  ylim = range(zoopgq),
             type = "n", axes = F,
             ylab = expression(Slope~(Delta~log(zoop)/Delta~log(phyt))),
             xlab = expression(Chl~italic(a)~(mu*g/L)))
        logtick.exp(0.0001, 10, c(1), c(F,F))
        axis(2)
        lines(ymean, zoopgq[1,], lty = "dashed")
        lines(ymean, zoopgq[3,], lty = "dashed")
        lines(ymean,zoopgq[2,])
        abline(h = targ, lty = "dotted")


        critout <- exp(approx(zoopgq[1,], ymean, targ)$y)

        cat("Criterion value:", critout, "\n")

    }
    stop()
    dev.off()
    return()

}


## uncomment next statement to run MCMC sim
fitout <- chlzoopmodel(dat.merge.cross, dat.merge.17, dftempmax, NULL, runmod = T)

## uncomment next two statements to post process
varout <- extract(fitout, pars = c( "f",  "cp", "biov_site", "b", "f0"))
chlzoopmodel(dat.merge.cross, dat.merge.17, dftempmax, varout, runmod = F)

