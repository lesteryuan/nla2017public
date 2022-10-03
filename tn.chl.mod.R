## TN model using NLA data

tn.model <- function(df1, df2, runmod = F) {
    require(rstan)
    nchains <- 4   # number of chains

    varlist1 <- c("ntl.result", "chl",  "us.l3code", "doc.result",
                  "no3no2.result")
    varlist2 <- c("ntl", "chla", "us.l3code", "doc",
                  "nitrate_nitrite_n")

    df2 <- df2[, c("unique.id", "visit.no", varlist2)]
    df2$year <- "17"
    df2 <- df2[, c("unique.id", "year", "visit.no", varlist2)]
    names(df2) <- c("site.id.c", "year", "visit.no", varlist1)


    ## adjust units on TN and NOx in 2017 data
    df2$ntl.result <- 1000*df2$ntl.result
    df2$no3no2.result <- 1000*df2$no3no2.result

    ## assume missing no3.no2 values are zero
    incvec <- is.na(df2$no3no2.result)
    df2[incvec, "no3no2.result"] <- 0

    ## select index sites from 12-17 data
    df1 <- subset(df1, sample.type == "MICX")
    df1$site.id.c <- factor(df1$site.id.c)
    df1 <- rbind(df1[, c("site.id.c", "year", "visit.no", varlist1)],
                 df2)
    ## start with just 07-12 data to figure out sampling error
    ## for tn
  #  df1 <- df1[, c("site.id.c", "year", "visit.no", varlist1)]

    ## omit records that are missing data
    incvec <- ! is.na(df1$ntl.result) & ! is.na(df1$no3no2.result) &
        ! is.na(df1$chl) & ! is.na(df1$doc.result) & !is.na(df1$us.l3code)

    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]
    print(nrow(df1))

    ## drop chl = 0
    incvec <- df1$chl > 0
    df1 <- df1[incvec,]

    ## compute TN-DIN and drop values that are <= 0
    df1$tkn <- df1$ntl.result - df1$no3no2.result
    incvec <- df1$tkn <= 0
    cat("TKN <= 0:", sum(incvec), "\n")
    df1 <- df1[!incvec,]

    ## reset state and L3 ecoregion factors
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    ## center chl and doc
    ## use same chl scale as TP model
    #load("chlsc.rda")
    #chlsc <- exp(mean(log(df1$chl)))
    load("tpchldat.717.rda")
    df1$chl.sc <- df1$chl/chlsc

    docsc <- exp(mean(log(df1$doc.result)))
    df1$doc.sc <- df1$doc.result/docsc

    ## center tn and nox
    tnsc <- exp(mean(log(df1$ntl.result)))

    ## nox error is 10/tnsc, or 10/615 which is 0.016
    print(tnsc)

    df1$tn.sc <- df1$ntl.result/tnsc
    df1$nox.sc <- df1$no3no2.result/tnsc

    ## split nox into log and normal error
    incvec <- df1$no3no2.result < 100
    ipnox1 <- (1:nrow(df1))[incvec]
    ipnox2 <- (1:nrow(df1))[!incvec]

    ## split tn into log and normal error
    incvec <- df1$ntl.result < 100
    iptn1 <- (1:nrow(df1))[incvec]
    iptn2 <- (1:nrow(df1))[!incvec]

#    save(tnsc, file = "tnsc.rda")

    ## drop 5 outliers in tn doc relationship
#    incvec <- (log(df1$tn.sc - df1$nox.sc) - log(df1$doc.sc)) > 3 |
#        (log(df1$tn.sc - df1$nox.sc) - log(df1$doc.sc)) < -2
#    print(sum(incvec))
#    df1 <- df1[!incvec,]

    tnchldat <- df1[, c("chl", "chl.sc", "tkn", "doc.sc", "doc.result",
                         "ntl.result", "tn.sc",
                        "no3no2.result", "nox.sc",
                        "econum", "us.l3code")]
#    save(tnchldat, docsc,tnsc, file = "tnchldat.717.rda")
#    stop()
    datstan <- list(n = nrow(df1), nnox = c(length(ipnox1), length(ipnox2)),
                    ntn = c(length(iptn1), length(iptn2)),
                    neco = max(df1$econum),econum = df1$econum,
                    tn = log(df1$tn.sc),
#                    nox = log(df1$nox.sc),
#                    tn1 = df1$ntl.result[iptn1],
#                    tn2 = log(df1$ntl.result[iptn2]),
                    nox1 = df1$nox.sc[ipnox1],
                    nox2 = log(df1$nox.sc[ipnox2]),
                    doc = log(df1$doc.sc),
                    chl = log(df1$chl.sc),
                    ipnox1 = ipnox1,
                    ipnox2 = ipnox2,
                    iptn1 = iptn1,
                    iptn2 = iptn2)

    print(str(datstan))
    stop()

    modstan <- '
        data {
            int n;                // number of samples
            int nnox[2];         // split nox into two error
            int neco;            // number of ecoregions
            int econum[n];       // ecoregion of each sample
            vector[n] tn;        // scaled TN
            vector[n] chl;       // scaled chl
//            vector[n] nox;       // scaled NOx
            vector[nnox[1]] nox1;
            vector[nnox[2]] nox2;

            vector[n] doc;       // scaled DOC
            int ipnox1[nnox[1]];
            int ipnox2[nnox[2]];
        }
        parameters {
            real muk;                // mean value of exponent on chl

            real mud[2];              // mean value of model coefficients
            real<lower = 0> sigd[2]; // SD of model coefficients among ecoregions
            vector[neco] etad2;
            vector[n] etad2a;

            real<lower = 0> sigtn;  // measurement error of tn
            vector<lower = 0>[n] noxmean;

        }
        transformed parameters {
            vector[neco] d2;
            vector[n] d2a;

            d2 = mud[2] + sigd[1]*etad2;
            d2a = d2[econum] + sigd[2]*etad2a;

        }
        model {
            matrix[n,3] temp;
            vector[n] tnmean;

            mud ~ normal(0,4);
            muk ~ normal(1,1);

            sigd ~ cauchy(0,4);
            etad2 ~ normal(0,1);
            etad2a ~ normal(0,1);

            sigtn ~ normal(0.1, 0.002);

            noxmean ~ cauchy(0,10);

            nox1 ~ normal(noxmean[ipnox1], 0.016);
            nox2 ~ normal(log(noxmean)[ipnox2], 0.1);

           temp[,1] = mud[1] + muk*chl;
           temp[,2] = d2a + doc;
           temp[,3] = log(noxmean);
           for (i in 1:n) tnmean[i] = log_sum_exp(temp[i,]);

            tn ~ student_t(4, tnmean, sigtn);
        }
    '
    if (runmod) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1100, chains = nchains,
                    warmup = 300, thin = 1)
        return(fit)
    }

    ## from the biogeochm paper. coef is 9.22, slope is 1.01
    b <- log(9.22/tnsc*(chlsc^1.06))
    plot(log(df1$chl.sc), log(df1$tn.sc-df1$nox.sc), col = "grey")
    abline(-1.88, 1.01)
    abline(-1.82, 0.97, col = "red")
    abline(b, 1.06, lty = "dashed")
    abline(v = log(100/chlsc))
    abline(h = log(1000/tnsc))
}

## save extracted variables to varout to post-process
fitout <- tn.model(dat.merge.cross,dat.merge.17, runmod = T)
#tn.model(dat.merge.cross, dat.merge.17, runmod = F)
#varout.n.limnat <- extract(fitout, pars = c("u","muk", "mud",  "d2", "d2a",
#                                       "sigd", "u", "muu"))


