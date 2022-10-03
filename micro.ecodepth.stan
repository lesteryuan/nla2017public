// stan model with ecoregion classes for d and depth for f
data {
    int n;                  // total number of samples
    int nsamp;              // number of samples without QA
    vector[nsamp] chl;      // observed log chl
    vector[nsamp] chl2;     // chl squared
    int nsite;              // number of lakes
    int sitenum[nsamp];     // lookup table for sitenum
    int propcyano[n];       // proportion cyanobacteria
    vector[n] biov;         // total phytoplankton biovolume
    int id[n];              // sample number for full dataset
    int mc[nsamp];          // observed MC
    int ndepth;             // number fo depth classes
    int idepth[nsamp];
    int neco;
    int econum[nsamp];
    int ntax;
    int itax[nsamp];
}
parameters {
    // parameters for chl-cyano relationship
    matrix[ndepth,3] muf;    // coef for propcyano

    real<lower = 0> sig[4];  // variance components
    vector[nsamp] eta1;      // overparametrization variables
    vector[nsamp] eta2;
    vector[n] eta3;

    // parameters for cyano-mc relationship
    real mud[2];
    real<lower = 0> sigd[2];
    matrix[neco, 2] etad;

    real<lower = 0> phi; //over dispersion for neg bin

    // parameters for site mean chl
    real<lower = 0> sigchl[2];
    vector[nsite] etachl;

    vector[ntax] etataxon;
    real<lower = 0> sigtaxon;
 }
transformed parameters {
    vector[nsamp] biov_mean;
    vector[nsamp] biov_obs;
    vector[nsamp] alpha;
    vector[nsamp] cyano_mean;
    vector[n] alpha2;
    vector[nsamp] mcmean;
    vector[nsite] chlsite;
    matrix[neco, 2] d;

    for (i in 1:2) d[,i] = mud[i] + etad[,i]*sigd[i];

    chlsite = etachl*sigchl[1];
    // quadratic relationship for chl-propcyano (Eqn 18)
    alpha = muf[idepth, 1] + muf[idepth,2] .* chl  +
            muf[idepth,3] .* chl2 + sig[2]*eta2;
    alpha2 = alpha[id] + sig[3]*eta3;

    // compute mean cyano (Eqn 21)
    biov_mean = chl + eta1*sig[1];
    biov_obs = biov_mean + etataxon[itax]*sigtaxon;
    cyano_mean = biov_mean + log(inv_logit(alpha));

    // piecewise linear relationship for MC (Eqn 22)
    mcmean = d[econum, 1] + d[econum, 2] .* cyano_mean;
}

model {
    // weakly informative priors
    sigchl ~ cauchy(0,3);
    etachl ~ normal(0,1);
    sig ~ cauchy(0,3);
    eta1 ~ normal(0,1);
    eta2~ normal(0,1);
    eta3 ~ normal(0,1);

    for (i in 1:3) muf[,i] ~ normal(0, 4);

    for (i in 1:2) etad[,i] ~normal(0,1);
    mud ~ normal(0,4);
    sigd ~ cauchy(0,3);

    phi ~ normal(0,3);
    etataxon ~ normal(0,1);
    sigtaxon ~ cauchy(0,3);

    // sampling statements
    chl ~ normal(chlsite[sitenum], sigchl[2]);
    biov ~ normal(biov_obs[id], sig[4]);
    propcyano ~ binomial_logit(100, alpha2);
    mc ~ neg_binomial_2_log(mcmean, phi);
}
