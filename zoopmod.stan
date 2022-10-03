        data {
           int nsite;              // number of distinct sites
           int nsamp;              // number of distinct samples
           vector[nsamp] chl;      // chl
           vector[nsamp] biov;     // biovolume
           int sitenum[nsamp];     // sitenums associated with each sample
           int nz;                 // number of zoop samples
           int ntaxon;             // number of taxonomists
           int ngrp;               // number of depth classes
           int taxonnum[nsamp];    // taxonomist
           int nsiteyrnum;         // number of distinct site-yr ids
           int depthf[nsiteyrnum];      // depth class for each site
           int ix[nsiteyrnum];         // site-yr index
           int sitenumz[nz];       // zoop sitenum
           int yeargnum[nsiteyrnum];        // year factor
           vector[nz] zoopb;       // zoop biomass
        }
        parameters {
           real<lower = 0> sigk;     // sd of chl about biov
           real<lower = 0> sig[2];   // biov sd among site, biov sd within site
           vector[nsite] eta1;

           real<lower = 0> phi; // dispersion parameter for zoop biomass

           vector[ngrp] cp;
           vector[ngrp*2] f0;  // intercept term allows for year differences
           matrix[ngrp,2] f;   // other two coefficients for model
           vector[ngrp] b;

           vector[ntaxon] etataxon;
           real<lower = 0> sigtaxon;
        }
        transformed parameters {
           vector[nsite] biov_site;   // estimate mean phyto biov at site
           vector[nsiteyrnum] zoopbmean;   // estimate mean site z biomass
           vector[ngrp] g;                 // difference of f2 and f1
           vector[ntaxon] btaxon;         // taxonomist-specific adjustment

           btaxon = etataxon*sigtaxon;

           biov_site = sig[1]*eta1;

           g = f[,2] - f[,1];

          // model for biomass for each site and year
          for (i in 1:nsiteyrnum) zoopbmean[i] = f0[yeargnum[i]] +
                      f[depthf[i], 2]*biov_site[ix[i]] +
                      g[depthf[i]]*exp(b[depthf[i]])*
                      log(1+exp(-(biov_site[ix[i]]-
                        cp[depthf[i]])/exp(b[depthf[i]])));
        }
        model {
           b ~ normal(-1,0.5);   // curvature prior is relatively low
           sigk ~ cauchy(0,3);
           eta1 ~ normal(0,1);
           sig ~ normal(0,2);

           etataxon ~ normal(0,1);
           sigtaxon ~ cauchy(0,3);

           cp ~ normal(0,0.5); //priors for breakpoints try to keep cp inside range of data

           f0 ~ normal(0,3);
           f[,1] ~ normal(1, 0.4); // prior for initial slope is based on
                                   // z/p = k in oligo lakes
           f[,2] ~ normal(0,0.2); // prior to push slope to zero at high chl
           phi ~ cauchy(0,3);

           biov ~ student_t(4,biov_site[sitenum] + btaxon[taxonnum], sig[2]);  // Eqn (3)
           chl ~ normal(biov_site[sitenum], sigk);  // Eqn (1)
           zoopb ~ normal(zoopbmean[sitenumz], phi); // Eqn (8)
        }
