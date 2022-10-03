## 11.6.2019: Production version of TP-chl model
## 12.17.2019Cleaned and commented

## 4.29.2022: Added error model for Chl and really slows down
## model fitting and only yielded a minor change in the fitted limit
## relationship. Don't think it's worth it.

## 5.2.2022. linear model for u-depth was not converging easily. Maybe not
## linear enough. Went back to 10 depth categories.

ntumodel <- function(df1, df2, varout = NULL, runmod = T) {
    require(rstan)

    nchains <- 3    # select number of chains

    depthsel <- 3.2              # lake depth
    ecosel <- 65                # Level III ecoregion code
    credint <- 0.80             # Credible interval for criterion calc
    chltarg <- 10               # Chl target

    varlist1 <- c("ptl.result", "chl",  "us.l3code", "turb.result",
                  "index.site.depth")
    varlist2 <- c("ptl", "chla", "us.l3code", "turb", "index.site.depth")

    df2 <- df2[, c("unique.id", "visit.no", varlist2)]

    df2$year <- "17"
    df2 <- df2[, c("unique.id", "year", "visit.no", varlist2)]
    names(df2) <- c("site.id.c", "year", "visit.no", varlist1)

    ## select index sites from 12-17 data
    df1 <- subset(df1, sample.type == "MICX")
    df1$site.id.c <- factor(df1$site.id.c)
    df1 <- rbind(df1[, c("site.id.c", "year", "visit.no", varlist1)],
                 df2)
    ## run with only 07-12 data
    #df1 <- df1[, c("site.id.c", "year", "visit.no", varlist1)]

    ## omit records that are missing data
    incvec <- ! is.na(df1$ptl.result) &
        ! is.na(df1$chl) & ! is.na(df1$turb.result) &
            ! is.na(df1$index.site.depth) &
                ! is.na(df1$us.l3code)
    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]

    df1$site.id.c <- factor(df1$site.id.c)
    print(nrow(df1))
    print(table(table(df1$site.id.c)))

    ## drop below detection limit turb and tp
    norig <- nrow(df1)
    incvec <- df1$turb.result > 0.01
    df1 <- df1[incvec,]
    incvec <- df1$ptl.result >= 3
    df1 <- df1[incvec,]
    incvec <- df1$chl >0
    df1 <- df1[incvec,]
    cat("N dropped for detection limit:", norig - nrow(df1), "\n")
    cat("Final N:", nrow(df1), "\n")

    ## scale chl and TP
    chlmn <- mean(log(df1$chl))
    chlsc <- exp(chlmn)
    print(chlsc)

    df1$chl.sc <- df1$chl/chlsc
    tpmn <- mean(log(df1$ptl.result))
    tpsc <- exp(tpmn)
    df1$tp.sc <- df1$ptl.result/tpsc
    print(tpsc)

    print(summary(df1$index.site.depth))

    ## define 10 depth classes based on quantiles
    cutp.depth <- quantile(log(df1$index.site.depth),
                           prob = seq(0, 1,length = 11))
    cutm <- 0.5*(cutp.depth[-1] + cutp.depth[-length(cutp.depth)])
    df1$dclass <- cut(log(df1$index.site.depth), cutp.depth, include.lowest = T)
    names(cutm) <- levels(df1$dclass)
    df1$dclassnum <- as.numeric(df1$dclass)

    ##drop HI, only dealing with conterminous US
    ## HI is dropped because of NA ecoregion

    ## reset L3 ecoregion codes
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    ## save data out to disk
    tpchldat <- df1
    save(tpchldat, tpsc, chlsc, cutp.depth, file = "tpchldat.rda")

    ## split tp data into two error regimes
    incvec <- df1$ptl.result < 20
    n1 <- sum(!incvec)
    n2 <- sum(incvec)
    print(n1)
    print(n2)
    ip1 <- (1:nrow(df1))[!incvec]
    ip2 <- (1:nrow(df1))[incvec]

    incvec.c <- df1$chl < 15
    n1c <- sum(!incvec.c)
    n2c <- sum(incvec.c)
    print(n1c)
    print(n2c)
    ip1c <- (1:nrow(df1))[!incvec.c]
    ip2c <- (1:nrow(df1))[incvec.c]

    datstan <- list(n = nrow(df1), n1 = n1, n2 = n2,
                    ndepth = max(df1$dclassnum),
                    depthnum = df1$dclassnum,
                    depth = log(df1$index.site.depth),
                    neco = max(df1$econum), econum = df1$econum,
                    ntu = log(df1$turb.result),
                    tp1 = log(df1$tp.sc)[!incvec],
                    tp2 = df1$tp.sc[incvec],
                    chl = log(df1$chl.sc),
                    ip1 = ip1,
                    ip2 = ip2)


    print(str(datstan))
    print(median(table(df1$econum)))
    stop()


    modstan <- '
        data {
            int n;                 // number of samples
            int n1;
            int n2;
            int ndepth;
            int depthnum[n];
//            vector[n] depth;       // measured depth
            int neco;              // number of L3 ecoregions
            int econum[n];         // ecoregion assigment
            vector[n] ntu;         // NTU
            vector[n1] tp1;          // TP
            vector[n2] tp2;
            vector[n] chl;         // Scaled Chl
            int ip1[n1];
            int ip2[n2];

        }
        parameters {
//            real a[2];      // coefficients for sediment-depth relationship

            real muu;
            vector[ndepth] eta_u1;
            real<lower= 0> sigu[2];  // standard deviation of ntu_np
            vector[n] eta_u2;

            real<lower = 0> signtu;  // measurement error of ntu

            real mub;                // mean coef of chl-ntu relationship
            real<lower = 0> sigb;    // standard deviation of b among ecoregions
            vector[neco] etab;
            real<lower = 0> muk[3];  // exponents on chl and u in models
            vector[2] mud;           // mean coef for tp model

            real<lower = 0> sigd[2];

            vector[neco] etad1;
            vector[n] etad2;
            real<lower = 0> sigtp[2];    // measurement error of tp


        }
        transformed parameters {
            vector[neco] b;
            vector[ndepth] u1;
            vector[n] u;
            vector[neco]  d1;
            vector[n] d1a;

            b = mub + etab*sigb;
            u1 = muu + sigu[1]*eta_u1;
            u = u1[depthnum] + sigu[2]*eta_u2;
//            u = a[1] + a[2]*depth + eta_u2*sigu;

            d1 = mud[1] + sigd[1]*etad1;
            d1a = d1[econum] + sigd[2]*etad2;
        }
        model {
            vector[n] ntu_mn;
            vector[n] tp_mn;

//            a ~ normal(0,5);
            muu ~ normal(0,3);
            sigu ~ cauchy(0,3);

            eta_u1 ~ normal(0,1);
            eta_u2 ~ normal(0,1);
            signtu ~normal(0.1,0.002);
            mub ~ normal(0,4);
            sigb ~ cauchy(0,3);
            etab ~ normal(0,1);

            muk ~ normal(1,0.3);    // muk should be somewhere around 1
            mud ~ normal(0,4);
            sigd ~ cauchy(0,3);

            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);


//            sigtp ~ cauchy(0,4);  // loose prior for sigtp
            sigtp[1] ~ normal(0.1, 0.002); //prior for log-transformed error
            sigtp[2] ~ normal(0.05, 0.003); // scaled prior for tp error
                                           // 2 ug/L divided overall mean of 40


            for (i in 1:n)
                tp_mn[i] = log_sum_exp(d1a[i] + muk[2]*u[i],
                                          mud[2] + muk[3]*chl[i]);
            // Eqn 25 from document
           for (i in 1:n)
               ntu_mn[i] = log_sum_exp(b[econum[i]] + muk[1]*chl[i], u[i]);

           ntu ~ student_t(4,ntu_mn, signtu);

           // Eqn 30 from document
            tp1 ~ student_t(4,tp_mn[ip1], sigtp[1]);
            tp2 ~ student_t(4, exp(tp_mn[ip2]), sigtp[2]);
        }
    '


    if (runmod) {

        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1100, chains = nchains,
                    warmup = 300, thin = 1, control = list(max_treedepth = 13))
        return(fit)
    }


    ## post processing
    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    muk <- apply(varout$muk, 2, mean)
#    d1 <- apply(varout$d1, 2, mean)
#    d2 <- apply(varout$d2, 2, mean)

#    muumn <- mean(varout$muu_mn)
#    d2raw <- d2 - muk[3]*log(chlsc)

    mud <- apply(varout$mud, 2, mean)

    mub <- mean(varout$mub)
    mubraw <- mub - muk[1]*log(chlsc)

    ## PLOT: relationship between Chl and turbidity (Fig 24)
#    png(width = 3, height = 2.5, pointsize = 8, units = "in", res = 600,
#        file = "chlturb.png")
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,1), mgp = c(2.3,1,0))
    plot(log(df1$chl), log(df1$turb.result), axes = F,
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = "Turbidity (NTU)", pch = 21, col = "grey",
         bg="white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    abline(mubraw, muk[1])
#    dev.off()

    dev.new()
    plot(log(df1$chl.sc), log(df1$ptl.result),
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), pch= 21, col = "grey",
         bg = "white", axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    abline(mud, muk[3])
    stop()

    b <- apply(varout$b, 2, mean)
    braw <- b - muk[1]*log(chlsc)

    ## PLOT: predicted values for Pdiss and muu (Fig 25)
    png(width = 6, height = 2.5, pointsize = 8, units = "in", res = 600,
        file = "Pdiss.turb.depth.png")
    par(mar = c(4,4,1,1), mgp = c(2.3, 1, 0), mfrow= c(1,2))
    plot(cutm, exp(muu), xlab = "Depth (m)", axes = F,
         ylab = expression(Turb[np]~(NTU)),
         pch = 21, col = "grey39", bg="white")
    axis(2)
    logtick.exp(0.001, 10, c(1), c(F,F))

    plot(cutm, exp(d1), xlab = "Depth (m)", axes = F,
         ylab = expression(P[diss]~(mu*g/L)),
         pch = 21, col = "grey39", bg="white")
    axis(2)
    logtick.exp(0.001, 10, c(1), c(F,F))
    dev.off()

    ## load in ntu_np mean values into main data
    df1$umean <- umean

    ## merge ecoregion and depth specific coefficients into main data
    dfd <- data.frame(num = 1:length(muu),  muu)
    names(dfd) <- c("dclassnum",  "muu")
    print(nrow(df1))
    df1 <- merge(df1, dfd, by = "dclassnum")
    print(nrow(df1))
    dfd3 <- data.frame(num = 1:length(d3), d3, d3raw, d2, d2raw)
    print(dim(dfd3))
    names(dfd3) <- c("econum", "d3", "d3raw", "d2", "d2raw")
    df1 <- merge(df1, dfd3, by = "econum")

    ## save ecoregion coefficients to file for mapping
    dfd3 <- merge(dfd3, unique.data.frame(df1[, c("econum", "us.l3code")]),
                  by = "econum")
    save(dfd3, file= "dfd3.rda")  # output to different script for mapping

        ## compute mean predicted TP
    df1$predout <- exp(mud[1]) + exp(df1$d2)*exp(df1$umean - muumn)^muk[2] +
        exp(df1$d3)*df1$chl.sc^muk[3]


    png(width = 4, height = 4, pointsize = 10, units = "in", res = 600,
        file = "pred.v.obs.all.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot(log(df1$predout), log(df1$ptl.result), xlab = "Predicted TP",
         ylab = "Observed TP",  pch = 21, col = "grey39",
         bg = "white", axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    abline(0,1)
    dev.off()

    rmsfnc <- function(x,y) {
        return(sqrt(sum((x-y)^2)/length(x)))
    }
    cat("Whole data RMS:",rmsfnc(log(df1$predout), log(df1$ptl.result)), "\n")

    ## print numerical summaries of parameters
    print("mub")
    print(exp(quantile(varout$mub - varout$muk[,1]*log(chlsc),
                   prob = c(0.05, 0.5, 0.95))))
    print("k distribution")
    print(apply(varout$muk, 2, quantile, prob = c(0.05, 0.5, 0.95)))
    mud <- apply(varout$mud, 2, mean)

    print("mud[2]")
    print(exp(quantile(varout$mud[,2] - varout$muk[,2]*varout$muu_mn,
                       prob = c(0.05, 0.50, 0.95))))
    print("mud[3]")
    print(exp(quantile(varout$mud[,3] - varout$muk[,3]*log(chlsc),
                       prob = c(0.05, 0.5, 0.95))))

    # PLOT: NTU_np vs TP, Chl vs TP (fig 27)
#    png(width = 6, height = 2.5, pointsize = 8, units = "in", res = 600,
#        file = "tpchl.png")

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(df1$umean, log(df1$ptl.result), xlab = expression(Turb[np]~(NTU)),
         ylab = expression(TP~(mu*g/L)), pch = 21, col = "grey",
         bg="white", axes = F)
    logtick.exp(0.000001, 10, c(1,2), c(F,F))
    x <- seq(min(df1$umean), max(df1$umean), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    for (i in 1:length(x)) {
        y <- varout$mud[,2] + varout$muk[,2]*x[i] - varout$muu_mn*varout$muk[,2]
        predout[i,] <- quantile(y, prob = c(0.05, 0.5, 0.95))
    }
    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,2])

        plot(log(df1$chl), log(df1$ptl.result),
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), pch= 21, col = "grey",
         bg = "white", axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,F))

    x <- seq(min(log(df1$chl)), max(log(df1$chl)), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(x))

    ## calculate estimate based on MO data
    predout.mo <- matrix(NA, ncol = 3, nrow = length(x))
    load("mn.val.mo.rda")

    d3 <- apply(varout$d3, 2, mean)

    print(summary(d3))
#    ip <- which(d3 == sort(d3)[2])
    ip <- 40
    print(ip)
    print(d3[ip])

    nsamp <- nrow(varout$mud)

    for (i in 1:length(x)) {
#        y <- varout$mud[,3] - varout$muk[,3]*log(chlsc) +
#            varout$muk[,3]*x[i]
        y <- rnorm(nsamp, mean = varout$mud[,3], sd = varout$sigd[,3])  -
            varout$muk[,3]*log(chlsc) +
                varout$muk[,3]*x[i]
        mnval <- varout.mo$d2[,3]
        y2 <- log(mn.val["tp"]) +
             mnval -
            varout.mo$k[,3]*log(mn.val["chl"]) +
                varout.mo$k[,3]*x[i]
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        predout.mo[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }
    print(predout)
    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,2])
    lines(x, predout.mo[,2], lty = "dashed")
    lines(x, predout.mo[,1], lty = "dotted")
    lines(x, predout.mo[,3], lty = "dotted")
    stop()

    dev.off()

    ## Plot criterion derivation figure

    ## find ecoregion index number
    ieco <- which(levels(df1$us.l3code) == ecosel)

    ## find correct depth class
    idepth <- 1
    while(exp(cutp.depth[idepth]) < depthsel & (idepth < length(cutp.depth)))
        idepth <- idepth + 1
    print(idepth)

    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    ## calculate mean parameter values
    muk <- apply(varout$muk, 2, mean)
    d3 <- apply(varout$d3, 2, mean)
    d3 <- d3 - muk[3]*log(chlsc)

    incvec <- df1$econum == ieco

    ## define regularly spaced values along chl gradient
    ## for computing predictions
    x <- seq(min(log(df1$chl)), max(log(df1$chl)), length = 50)
    xsc <- x - log(chlsc)

    ## compute predicted ambient and limiting values
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    predout2 <- matrix(NA, ncol = 3, nrow = length(x))
    for (i in 1:length(x)) {
        y <- varout$d3[,ieco] - varout$muk[,3]*log(chlsc) +
            varout$muk[,3]*x[i]
        predout[i,] <- quantile(y, prob = c(0.5*(1-credint),
                                       0.5, 1- 0.5*(1-credint)))
        y2 <- exp(varout$mud[,1]) +
            exp(varout$d2[,ieco])*
                (exp(varout$muu[, idepth] -
                         varout$muu_mn))^varout$muk[,2] +
                             exp(varout$d3[, ieco])*exp(xsc[i])^varout$muk[,3]
        predout2[i,] <- quantile(log(y2), prob = c(0.5*(1-credint),
                                              0.5, 1-0.5*(1-credint)))
    }

    dev.new()
    par(mgp = c(2.3,1,0), bty = "l", mar = c(4,4,1,1))
    plot(log(df1$chl), log(df1$ptl.result), type = "p",
         axes = F, xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), col = "grey80",
         pch = 21, bg = "white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    ## highlight points in the selected ecoregion
    points(log(df1$chl)[incvec], log(df1$ptl.result)[incvec],
           pch = 21, col = "black", bg = "grey")

    cat("Median depth:", median(df1$index.site.depth[incvec]), "\n")

    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col  = grey.t, border = NA)
    lines(x, predout[,2])
    polygon(c(x, rev(x)), c(predout2[,1], rev(predout2[,3])),
            col  = grey.t, border = NA)
    lines(x, predout2[,2], lty = "dashed")

    crit1 <- approx(x, predout[,1], log(chltarg))$y
    crit2 <- approx(x, predout2[,1], log(chltarg))$y

    ylo <- min(log(df1$ptl.result)) - 0.04*diff(range(log(df1$ptl.result)))
    xlo <- min(log(df1$chl)) - 0.04*diff(range(log(df1$chl)))
    segments(xlo, crit1,
             log(chltarg), crit1,  col = "red")
    segments(xlo, crit2,
             log(chltarg), crit2,  col = "red")
    segments(log(chltarg), crit2,
             log(chltarg), ylo, col = "red")
    segments(xlo, crit1,
             xlo, crit2, col = "red", lwd = 3)

    cat("TP criterion:", round(exp(crit2)), "\n")


}

## runmod variable set to T to run simulation and set to F to
##  run post processing.

fitout <- ntumodel(dat.merge.cross, dat.merge.17, runmod = T)
#ntumodel(dat.merge.cross, dat.merge.17, varout = varout, runmod = F)
