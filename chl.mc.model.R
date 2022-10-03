## 12.10.2019: Cleanded and commented
## 8.02.2018: simplified MC model to speed up
## 9.18.2020: Eqn numbers refer to numbers in criterion document
## 8.23.2021: add in 2017 data and rerun
##          : add in DOC and depth model
## 6.17.2022: Look at ecoregions in cyano-mc relationship
chlrat <- function(df1, df2, qa, lookup, varout = NULL, runmod = T,
                   orig.mod = T, title0 = "") {

    require(rstan)
    require(mgcv)
#    source("logtick.exp.R")
#    source("binv.R")

    nchains <- 4  # number of concurrent chains in MCMC simulation

    ## set parameters here for calculating criteria
    p.exceed <- 0.01  # set allowable exceedance rate
    thold <- 8         # set MC threshold
    credint <- 0.6    # set credible interval for candidate criterion

    ## combine 12 and 17 data
    df1 <- df1[, c("site.id.c", "visit.no", "year", "sample.type", "biov.cyano",
                   "prop.cyano", "biov.tot",
                   "chl", "result",
                   "st.nla2012", "us.l3code", "doc.result","cond.result",
                   "index.site.depth", "taxonomist")]
    names(df1)[9] <- "mc"
    df1$mc[is.na(df1$mc)] <- 0  # set mc = NA in 2007 - 2012 data to zero

    df2$sample.type <- "MICX"
    df2$year <- "17"
    df2 <- df2[, c("unique.id", "visit.no", "year", "sample.type",
                   "biov.cyano",
                   "prop.cyano", "biov", "chla",
                   "micx", "pstl.code", "us.l3code", "doc","cond",
                   "index.site.depth", "analyst")]
    names(df2) <- names(df1)

    ## rbind 2007/2012 with 2017
    df1 <- rbind(df1, df2)

    ## make mc into integer times 10 to work with negbin distribution
    df1$mc <- round(df1$mc*10)

    ## set below detection limit chl to NA
    incvec <- df1$chl == 0
    incvec[is.na(incvec)] <- F
    df1$chl[incvec] <- NA

    ## log transforms
    df1$biov.cyano <- log(df1$biov.cyano*1e-6)
    df1$biov.tot <- log(df1$biov.tot*1e-6)
    df1$chl <- log(df1$chl)
    df1$doc.result <- log(df1$doc.result)

    ## validation data are the 2007 data because no phytoplankton
    ## count data in 2007
    ## put these aside before centering and classifying
    dfvalid <- df1[df1$year == "07",]
    df1 <- df1[df1$year != "07",]

    ## samples with zero cyano biov set to missing
    df1$biov.cyano[is.infinite(df1$biov.cyano)] <- NA

    ## drop samples missing cyano or chl
    incvec <- !is.na(df1$biov.cyano) & ! is.na(df1$chl) &
        ! is.na(df1$index.site.depth) &  ! is.na(df1$mc)
    df1 <- df1[incvec,]

    ## outliers: two cyano biov measurements that are very low omitted
    incvec <- df1$biov.cyano < -3 & df1$mc > 50
    df1 <- df1[!incvec,]

    ## center chl and biov
    mnval.chl <- mean(df1$chl)
    df1$chl <- df1$chl - mnval.chl
    mnval.biov <- mean(df1$biov.tot)
    df1$biov.tot <- df1$biov.tot - mnval.biov
    df1$biov.cyano <- df1$biov.cyano - mnval.biov
    mnval <- c(mnval.chl, mnval.biov)
    names(mnval) <- c("chl","biov")

    ## prepare QA data
    ## convert QA sample ids to same format as dat.merge.cross
    qa$idnew <- factor(qa$idnew)
    str <- as.character(qa$idnew)
    w <- regexpr("--", str)
    qa$site.id <- substring(str, 1, w-1)
    qa <- merge(qa, lookup, by = "site.id")
    qa$year <- "12"
    qa$sample.type <- paste("MIC", substring(str, nchar(str), nchar(str)),
                            sep = "")
    qa$visit.no <- as.numeric(substring(str, nchar(str)-3, nchar(str)-3))
    qa$biov.cyano1 <- log(qa$biov1.cyano*1.e-6)
    qa$biov.cyano2 <- log(qa$biov2.cyano*1.e-6)
    qa$prop.cyano2 <- qa$biov2.cyano/qa$biov2
    qa$biov.tot2 <- log(qa$biov2*1.e-6) - mnval.biov

    ## only MSU resample is a duplicate so only use that one
    selvec <- ! is.na(qa$biov.cyano2)
    qa <- qa[selvec,]
    ## merge qa data with main data set to find qa samples that
    ## match recorded samples (122 out of 129)
    df2 <- merge(df1, qa[, c("site.id.c", "visit.no", "year","sample.type",
                             "biov.cyano2", "prop.cyano2", "biov.tot2")],
                 by = c("site.id.c", "year", "visit.no", "sample.type"))

    ## set up matched qa data to combine with full dataset
    df1b <- df2[, c("site.id.c", "visit.no", "year", "sample.type",
                    "biov.cyano2", "prop.cyano2", "biov.tot2")]
    names(df1b) <- c("site.id.c", "visit.no", "year", "sample.type",
                     "biov.cyano", "prop.cyano", "biov.tot")

    ## add in blank columns so that I can rbind QA data with main data
    df1b$chl <- NA
    df1b$mc <- NA
    df1b$st.nla2012 <- NA
    df1b$us.l3code <- NA
    df1b$doc.result <- NA
    df1b$cond.result <- NA
    df1b$index.site.depth <- NA
    df1b$taxonomist <- NA

    df1 <- rbind(df1, df1b)

    ## drop one sample identified by taxonomist "XX"
    incvec <- df1$taxonomist == "XX"
    incvec[is.na(incvec)] <- FALSE
    df1 <- df1[!incvec,]

    ## reset site id factors
    df1$idnew <- factor(paste(df1$site.id.c, df1$year, df1$visit.no,
                              df1$sample.type, sep = "--"))

    df1$sampnum <- as.numeric(df1$idnew)

    df1$site.id.c <- factor(df1$site.id.c)
    df1$st.nla2012 <- factor(df1$st.nla2012)
    df1$istate <- as.numeric(df1$st.nla2012)
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)
    df1$taxonomist <- factor(df1$taxonomist)
    df1$itax <- as.numeric(df1$taxonomist)

    ## specify dfsamp as samples with mc (not replicate cyano)
    incvec <- ! is.na(df1$mc)
    dfsamp <- df1[incvec,]
    dfsamp <- dfsamp[order(dfsamp$sampnum),]

    ## set up four depth groups
    cutp1 <- quantile(dfsamp$index.site.depth, prob = seq(0, 1, length = 5))
    depthfac <- cut(dfsamp$index.site.depth, cutp1, include.lowest = T)
    fac <- depthfac
    dfsamp$idepth <- as.numeric(fac)

    dfsamp$chl.loge <- dfsamp$chl # add chl.loge to match old field name for shiny
    ## write out dfsamp for shiny app
    save(dfsamp, mnval, cutp1, file = "shinydat.2017.rda")

    ## lookup table matching sample to site
    lookup2 <- unique.data.frame(dfsamp[, c("sampnum", "site.id.c")])
    lookup2$sitenum <- as.numeric(lookup2$site.id.c)
    lookup2 <- lookup2[order(lookup2$sampnum),]

    datstan <- list(n = nrow(df1),
                    nsamp = nrow(dfsamp),
                    chl = dfsamp$chl,
                    chl2 = dfsamp$chl^2,
                    nsite = max(lookup2$sitenum),
                    sitenum = lookup2$sitenum,
                    propcyano = round(df1$prop.cyano*100),
                    biov = df1$biov.tot,
                    id = df1$sampnum,
                    mc = dfsamp$mc,
                    nstate = length(levels(dfsamp$st.nla2012)),
                    istate = as.numeric(dfsamp$st.nla2012),
                    ndepth = max(dfsamp$idepth),
                    idepth = dfsamp$idepth,
                    neco = max(dfsamp$econum),
                    econum = dfsamp$econum,
                    ntax = max(dfsamp$itax),
                    itax = dfsamp$itax)

    print(str(datstan))

    if (runmod) {

        rstan_options(auto_write= TRUE)
        options(mc.cores = nchains)
        fitout <- stan(file = "micro.ecodepth.stan", data = datstan,
                       chains = nchains, iter = 1200, warmup = 400, thin =1,
                       control = list(max_treedepth = 15, adapt_delta = 0.9))
        return(fitout)
    }
    else {

        grey.t <- adjustcolor("grey39", alpha.f = 0.5)

        plot1 <- F
        if (plot1) {
            ## ECOREGIONAL EFFECTS
            ## compare coef for cyano-mc to alloc DOC levels
            ## from explore.doc model (u2)
            ## this is evidence that some of the differences among
            ## ecoregion are due to the average amount of alloc DOC
            ## loading
            d <- apply(varout$d, c(2,3), mean)
            dftemp <- data.frame(econum = 1:nrow(d), d)
            names(dftemp) <- c("econum", "d1", "d2")
            dftemp <- merge(dftemp,
                            unique.data.frame(df1[, c("econum", "us.l3code")]),
                            by = "econum")

            dftemp <- merge(dftemp, dfu2, by = "us.l3code")

            dev.new()
            par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0), bty = "l")
            plot(dftemp$u2, dftemp$d1)
            require(mgcv)
            mod <- gam(d1 ~ s(u2, k =3), data = dftemp)
            iord <- order(dftemp$u2)
            lines(dftemp$u2[iord], predict(mod)[iord])
            print(summary(mod))
            plot(dftemp$u2, dftemp$d2)
            mod <- lm(d2 ~ u2, data = dftemp)
            print(summary(mod))
            abline(mod)
        }

        plot2 <- FALSE
        if (plot2) {
            ## CHL VS. CYANO RELATIVE ABUNDANCE
            ## plot chl vs. cyano relative abundance (prop.cyano)
            ## prepare prop.cyano for plotting by converting 1s and 0s to
            ## values close to 1 and 0, so they can be logit transformed
            incvec <- dfsamp$prop.cyano == 1
            dfsamp$prop.cyano[incvec] <- NA
            maxval <- max(dfsamp$prop.cyano, na.rm = T)
            dfsamp$prop.cyano[incvec] <- 1 - (1-maxval)*0.5
            incvec <- dfsamp$prop.cyano == 0
            dfsamp$prop.cyano[incvec] <- NA
            minval <- min(dfsamp$prop.cyano, na.rm = T)
            dfsamp$prop.cyano[incvec] <- 0.5*minval

            sig <- apply(varout$sig, 2, mean)
            sig.pcyano <- sqrt(sig[2]^2 + sig[3]^2)
            print(sig.pcyano)
            png(width = 5, height = 4, pointsize = 10, units = "in",
                res = 600, file = "cyanoabund.png")
            par(mar = c(4,4,1,1), mfrow = c(2,2), mgp = c(2.3,1,0))
            for (i in 1:4) {
                incvec <- dfsamp$idepth == i
                dftemp <- dfsamp[incvec,]
                ## average prop.cyano in bins of about 20 samples
                bout <- binv(dftemp$chl,
                             qlogis(dftemp$prop.cyano), 10)

                ## compute posterior predictions of mean relationship
                np <- 80   # posterior predictions at 80 points
                ## set up storage for median and upper and lower credible intervals
                predout <- rep(NA, times = np*3)
                dim(predout) <- c(np, 3)
                xrange <- range(bout$xb)
                x <- seq(xrange[1], xrange[2], length = np)
                nsamp <- length(varout$muf[,i,1])
                for (j in 1:np) {
                    y <- rep(varout$muf[,i,1] + varout$muf[,i,2]*x[j],
                             times = 5) +
                        rnorm(nsamp*10, mean = 0, sd = sig.pcyano)
                    ## compute median prediction and 90% credible intervals
                    predout[j,] <- plogis(quantile(y,prob = c(0.4,0.5,0.6)))
                }

                plot(bout$xb+mnval["chl"], plogis(bout$yb), axes = F,
                     pch = 21,  col = "grey39",
                     bg = "white", ylab = "Cyanobacteria rel. abund.",
                     xlab = expression(Chl~italic(a)~(mu*g/L)),
                     ylim = c(0, 1),
                     xlim = quantile(dfsamp$chl+mnval["chl"], prob = c(0.01, 0.99)))
                logtick.exp(0.001, 10, c(1), c(F,T))
                axis(2)
                polygon(c(x, rev(x)) + mnval["chl"],
                        c(predout[,1], rev(predout[,3])),
                        border= NA, col = grey.t)
                lines(x + mnval["chl"], predout[,2])

            }
            dev.off()
        }

        ## PREDICTED VS OBSERVED MC PLOT
        ## drop misssing data from validation data
        incvec <- ! is.na(dfvalid$chl) & ! is.na(dfvalid$index.site.depth)
        dfvalid <- dfvalid[incvec,]

        ## set up depth classes in validation data and merge in ecoregion
        dfvalid$depthfac <- cut(dfvalid$index.site.depth, cutp1,
                                include.lowest = T)
        dfvalid$idepth <- as.numeric(dfvalid$depthfac)
        ## drop validation sites with depths outside the range of
        ## calibration data
        dfvalid <- dfvalid[! is.na(dfvalid$idepth),]
        dfvalid$idnew <- paste(dfvalid$site.id.c, dfvalid$year,
                               dfvalid$visit.no, sep = "--")
        dfvalid <- merge(dfvalid,
                         unique.data.frame(df1[,c("us.l3code", "econum")]),
                         by = "us.l3code")

        ## center chl
        dfvalid$chl <- dfvalid$chl - mnval.chl
#        dfvalid <- dfsamp

        ## calculate mean predicted MC (original model and revised model)
        phi <- mean(varout$phi)
        if (! orig.mod) {
            muf <- apply(varout$muf, c(2,3), mean)
            d <- apply(varout$d, c(2,3),mean)
            alpha <- muf[dfvalid$idepth, 1] + muf[dfvalid$idepth,2]*dfvalid$chl +
                muf[dfvalid$idepth,3]*dfvalid$chl^2
            cyano <- dfvalid$chl + log(plogis(alpha))
            mcmean <- d[dfvalid$econum,1] + d[dfvalid$econum,2]*cyano
        }
        else {
            ## original model
            muf <- apply(varout$muf, 2, mean)
            cp <- mean(varout$cp)
            d <- apply(varout$d, 2, mean)
            alpha = muf[1] + muf[2] * dfvalid$chl +
                muf[3]*(dfvalid$chl^2)
            cyano_mean = dfvalid$chl + log(plogis(alpha));
            mcmean = d[1] + d[2]*cyano_mean +
                (d[3]-d[2])*(cyano_mean - cp) * plogis(10*(cyano_mean-cp));
        }

        dfout1 <- data.frame(idnew = dfvalid$idnew, chl = dfvalid$chl,
                             mcmean = mcmean,  mc = dfvalid$mc)

        nb = 40 # number of bins for estimating mean MC
        cutp <- quantile(dfout1$mcmean, prob = seq(0, 1,
                              length = round(length(dfout1$mcmean)/nb)))
        cutf <- cut(dfout1$mcmean, cutp, include.lowest = T)
        cutm <- 0.5*(cutp[-1] + cutp[-length(cutp)])
        val <- rep(NA, times = length(levels(cutf)))
        val.sd <- rep(NA, times = length(levels(cutf)))
        for (i in 1:length(levels(cutf))) {
            incvec <- cutf == levels(cutf)[i]
            if (! all(dfout1$mc[incvec] == 0)) {
                mod <- gam(dfout1$mc[incvec] ~ 1, family = negbin(theta = phi))
                val[i] <- summary(mod)$p.coeff
                val.sd[i] <- summary(mod)$se
            }
        }
        cutm <- cutm - log(10)
        cutp <- cutp - log(10)
        val <- val - log(10)

        up <- val + 1.96*val.sd # 95 confidence limits on estimated mean
        dn <- val - 1.96*val.sd # for each bin

        plot(cutm, val, ylim = range(c(up,dn), na.rm = T),
             axes = F, xlab = expression(Predicted~MC~(mu*g/L)),
             ylab = expression(Observed~MC~(mu*g/L)))
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        mtext(title0, side = 3, line = -1, cex = 1.5)
        segments(cutm, up, cutm,dn)
        segments(cutp[-length(cutp)], val,
                 cutp[-1], val)
        points(cutm, val, pch = 21, bg = "white")

        rms <- function(x,y) {
            incvec <- ! is.na(x)
            sqrt(sum((x[incvec]-y[incvec])^2/sum(incvec)))
        }

        bias <- mean(cutm - val, na.rm = T)
        cat("Mean bias:", bias,"\n")
        ## calculate prediction error after correcting for mean bias
        abline(0,1, lty = "dashed")
        rms0 <- rms(val, cutm - bias)
        cat("RMS error:", rms0, "\n")

        return()

        ## simulate mc given chl using model output

        # set up chl values where predicted output will be computed
        # 50 values ranging from minimum to maximum chl
        chlnew <- seq(min(dfsamp$chl.loge),
                      max(dfsamp$chl.loge),length = 50)
        nit <- length(varout$phi)  # number of samples in MCMC sim

        ## store sampled posterior distribution of model parameters
        ## repeat them fac times to increase sample size
        ## for plotting
        fac <- 2
        f1 <- rep(varout$muf[,1], times = fac)
        f2 <- rep(varout$muf[,2], times = fac)
        f3 <- rep(varout$muf[,3], times = fac)
        sig1 <- rep(varout$sig[,1], times = fac)
        sig2 <- rep(varout$sig[,2], times = fac)
        cp <- rep(varout$cp, times = fac)
        d1 <- rep(varout$d[,1], times = fac)
        d2 <- rep(varout$d[,2], times = fac)
        d3 <- rep(varout$d[,3], times = fac)
        sigchl <- rep(varout$sigchl[,2], times = fac)

        ## set up storage location for predicted mc
        mcpred <- matrix(NA, nrow=length(chlnew), ncol = nit*fac)
        ## compute predicted MC
        for (i in 1:length(chlnew)) {
            ## get error associated with estimating seasonal mean chl
            chlsamp <- rnorm(fac*nit, mean = chlnew[i], sd = sigchl)
            ## predict proportion cyano for each chl
            alpha <- rnorm(fac*nit, mean = f1 + f2*chlsamp + f3*chlsamp^2,
                           sd = sig1)
            ## estimate of cyano abundance
            mn0 <- rnorm(fac*nit, mean = log(plogis(alpha)) + chlsamp,
                         sd = sig2)
            ## predict mean MC given cyano
            y.sc <- d1 + d2*mn0 + (d3 - d2)*mn0*plogis(10*(mn0 - cp))
            # sample from negative binomial to get distribution of
            # observed MC given mean mc at the selected percentile
            mcpred[i,] <- qnbinom(1-p.exceed, mu = exp(y.sc),
                                  size = varout$phi)
        }

        # compute mean, 90th and 99th percentile of predicted distribution
        # of microcysin
        mcmean0 <- apply(mcpred, 1, median)
        ## compute selected credible intervals
        mcup <- apply(mcpred, 1, quantile, prob = 1-(1-credint)*0.5)
        mclo <- apply(mcpred, 1, quantile, prob = (1-credint)*0.5)

        # compute candidate chl criterion
        crit <- approx(mcup, chlnew, thold*10)$y

        dev.new()
        layout(matrix(c(1,2), nrow = 2, ncol = 1),
               heights = c(1, 0.5))
        par(mar = c(0,5,1,1), mgp = c(2.3,1,0))

        ## identify observed microcystin samples that are zero
        ## these needed to be plotted differently in log-transformed coord
        iszero <- hbdat$dfsamp$mc == 0
        hbdat$dfsamp$mc[iszero] <- NA
        xlim <- range(c(hbdat$dfsamp$chl.loge + mnval["chl"]))
        ylim <- c(log(0.1), max(log(mcup*0.1)))

        ymin <- ylim[1] - 0.04*diff(ylim)
        xmin <- xlim[1] - 0.04*diff(xlim)

        ## set zero values of lower CL to ymin so that error bound is drawn
        ## correctly in log-transformed coord
        mclo[mclo ==0] <- exp(ymin)
        ## plotting observed data for grab chl vs mc, but
        ## modeled relationship is based on
        ## seasonal mean chl, because the difference isn't that great.
        ## Effect of seasonal mean chl is manifested in the credible
        ## intervals
        plot(hbdat$dfsamp$chl.loge + mnval["chl"],
             log(hbdat$dfsamp$mc*0.1),
             xlab = expression(Chl~italic(a)~(mu*g/L)),
             ylab = expression(MC~(mu*g/L)), axes = F,
             xlim  = xlim, ylim = ylim,  type = "n")
        points(hbdat$dfsamp$chl.loge + mnval["chl"],
           log(hbdat$dfsamp$mc*0.1), pch = 21, col = "grey39",
           bg = "white", cex = 0.8)
#            points(bout$xb, log(bout$yb), pch = 21, col = "grey39",
#                   bg = "white")
        logtick.exp(0.0001, 10, c(2), c(F,F))
        tickloc <-  c(seq(0.01, 0.09, 0.01), seq(0.1, 0.9, 0.1),
                      seq(1,9, 1), seq(10,90,10), seq(100,900, 100))
        axis(1, at = log(tickloc), lab = NA, tcl = -0.2)
        axis(1, at = log(c(0.01, 0.1, 1, 10, 100, 1000)), lab = NA)

        # add predicted percentile lines
        polygon(c(chlnew + mnval["chl"], rev(chlnew + mnval["chl"])),
                c(log(mcup*0.1), rev(log(mclo*0.1))), border = NA,
                col = grey.t)
        lines(chlnew + mnval["chl"], log(mcmean0*0.1))

        # draw crit line
        if (!is.na(crit)) {
            segments(xmin, log(thold), crit + mnval["chl"], log(thold),
                     col = "red")
            segments(crit+mnval["chl"], log(thold),
                     crit+mnval["chl"], ymin, col = "red")
        }
        else cat("Criteria could not be calculated.\n")

        ## put plot of non-detects at the bottom
        par(mar = c(4,5,1,1))
        cutp <- quantile(dfsamp$chl.loge, prob = seq(0,1, length = 21))
        cutm <- (cutp[-1] + cutp[-length(cutp)])*0.5
        cutf <- cut(dfsamp$chl.loge, cutp, include.lowest = T)
        pnond <- tapply(as.numeric(iszero), cutf, mean, include.lowest = T)

        plot(cutm + mnval["chl"], pnond, xlim = xlim, axes = F,
             xlab = expression(Chl~italic(a)~(mu*g/L)),
             ylab = "Proportion\n non-detect")
        axis(2)
        logtick.exp(0.001, 10, c(1), c(F,F))

    }

    return()

}

#fitout <- chlrat(dat.merge.cross, dat.merge.17, qa.biov, nars.cross, runmod = T)

#png(width = 6, height = 3, pointsize = 6, units = "in", res = 600,
#    file = "pred.v.obs.png")
#par(mar = c(4,4,1,1), mgp = c(2.3,1,0), mfrow = c(1,2))
#chlrat(dat.merge.cross, dat.merge.17, qa.biov, nars.cross, varout = varout.orig,
#       runmod = F, orig.mod = T, title0 = "Original")
#chlrat(dat.merge.cross, dat.merge.17, qa.biov, nars.cross, varout = varout,
#       runmod = F, orig.mod = F, title = "Revised")
#dev.off()

