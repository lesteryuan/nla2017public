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
        ## set parameters here for calculating a criterion
        p.exceed <- 0.05  # set allowable exceedance rate
        thold <- 8         # set MC threshold
        credint <- 0.75    # set credible interval for candidate criterion
        depth <- 10
        ecoreg <- 51

        ## find bin for depth
        idepth0 <- which(depth < cutp1[-1] & depth > cutp1[-length(cutp1)])
        econum0 <- which(levels(dfsamp$us.l3code) == ecoreg)

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
        f1 <- rep(varout$muf[,idepth0, 1], times = fac)
        f2 <- rep(varout$muf[,idepth0, 2], times = fac)
        f3 <- rep(varout$muf[,idepth0, 3], times = fac)
        sig1 <- rep(varout$sig[,1], times = fac)
        sig2 <- rep(varout$sig[,2], times = fac)
        d1 <- rep(varout$d[,econum0,1], times = fac)
        d2 <- rep(varout$d[,econum0,2], times = fac)
        sigchl <- rep(varout$sigchl[,2], times = fac)

        ## set up storage location for predicted mc
        mcpred <- matrix(NA, nrow=length(chlnew), ncol = nit*fac)
        ## compute predicted MC
        for (i in 1:length(chlnew)) {
            ## get error associated with estimating seasonal mean chl
            chlsamp <- rnorm(fac*nit, mean = chlnew[i], sd = sigchl)
            ## predict proportion cyano for each chl
            alpha <- rnorm(fac*nit, mean = f1 + f2*chlsamp + f3*chlsamp^2,
                           sd = sig2)
            ## estimate of cyano abundance
            mn0 <- rnorm(fac*nit, mean = log(plogis(alpha)) + chlsamp,
                         sd = sig1)
            ## predict mean MC given cyano
            y.sc <- d1 + d2*mn0 
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
        crit <- approx(mcup, chlnew, thold*10)$y + mnval.chl
        cat("Criterion value:", round(exp(crit), digits = 1), "\n")

        dev.new()
        layout(matrix(c(1,2), nrow = 2, ncol = 1),
               heights = c(1, 0.5))
        par(mar = c(0,5,1,1), mgp = c(2.3,1,0))

        ## identify observed microcystin samples that are zero
        ## these needed to be plotted differently in log-transformed coord
        iszero <- dfsamp$mc == 0
        dfsamp$mc[iszero] <- NA
        xlim <- range(c(dfsamp$chl.loge + mnval["chl"]))
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
        plot(dfsamp$chl.loge + mnval["chl"],
             log(dfsamp$mc*0.1),
             xlab = expression(Chl~italic(a)~(mu*g/L)),
             ylab = expression(MC~(mu*g/L)), axes = F,
             xlim  = xlim, ylim = ylim,  type = "n")
        points(dfsamp$chl.loge + mnval["chl"],
           log(dfsamp$mc*0.1), pch = 21, col = "grey39",
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
            segments(xmin, log(thold), crit, log(thold),
                     col = "red")
            segments(crit, log(thold),
                     crit, ymin, col = "red")
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
#varout <- extract(fitout, pars = c("d", "sig", "phi", "muf", "sigchl"))
## run next statement to calculate criteion
#chlrat(dat.merge.cross, dat.merge.17, qa.biov, nars.cross, varout = varout,
#       runmod = F)


