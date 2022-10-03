         data {
             int n[3];                    // N above 2, below 2, and all
             int nsite;                   // number of samples
             vector[n[3]] chl;            // sample Chl

             vector[nsite] depth;         // depth below thermocline
             vector[nsite] doc;           // mean doc
             vector[nsite] dmax;          // maximum DO
             vector[nsite] air4;
             vector[nsite] area;
             vector[nsite] depthall;

             vector[n[1]] t1;             // sample day above 2s
             vector[n[2]] t2;             // sample day below 2s
             vector[n[1]] domean;         // mean DO
             int sitenum[n[3]];           // sitenums for different samples
             int sitenum1[n[1]];
             int sitenum2[n[2]];

             int n0;                      // N samples for MO
             vector[n0] t0;               // sampling day for MO
             vector[n0] domean0;          // mean DO for MO
             int nsite0;                  // N sites
             int nsiteyr0;                // N site-year
             vector[nsite0] dmax0;                  // maximum DO for MO
             vector[nsiteyr0] air40;      // Air temperature
             vector[nsiteyr0] depthall0;      // Air temperature
             vector[nsiteyr0] area0;      // Air temperature
             vector[nsite0] chl0;         // Chl
             vector[nsite0] depth0;       // depth below thermocline
             vector[nsite0] doc0;         // mean DOC
             int sitenum0[n0];
             int siteyrnum0[n0];

         }
         parameters {
             real muchl;                       // model for mean site chl
             real<lower =0> sigchl[2];         // variance terms for chl model
             vector[nsite] etachl;

             real mud[4];                      // Coefficients for VOD model
             vector<upper = 2>[n[2]] domean2;  // Estimates for censored DO

             real b[4];             // Coefficients for strat day model
             real<lower = 0> sigt;             // Variance in strat day model
             vector[nsite] etat;
             real<lower = 0> sigt0;	
             vector[nsiteyr0] etat0;

             real<lower = 0> sig;              // variance in predicted DO
         }
         transformed parameters {
             vector[nsite] chlsite;
             vector[nsite] tstart;
             vector[n[1]] dopred1;
             vector[n[2]] dopred2;
             vector[nsiteyr0] tstart0;

             vector[n0] dopred0;
             vector[nsite] d2n;
             vector[nsite0] d2;

             chlsite = muchl + etachl*sigchl[1];

             // strat day models
             tstart = b[4]*air4 + b[1] + b[2]*area + b[3]*depthall + etat*sigt;

             tstart0 = b[4]*air40 + b[1] + b[2]*area0 + b[3]*depthall0 + etat0*sigt0;

             // VOD model in NLA
             d2n  = -exp(mud[1]) + mud[2]*chlsite + mud[3]*depth+ mud[4]*doc;
             // VOD model in MO
             d2 = -exp(mud[1]) + mud[2]*chl0 + mud[3]*depth0+ mud[4]*doc0;

             // Predictions of mean DO in MO
             dopred0 = dmax0[sitenum0] +d2[sitenum0] .* (t0 - tstart0[siteyrnum0]);
             // Eqn 10 and 12, model for VOD and DO
             dopred1 = dmax[sitenum1] + d2n[sitenum1] .*(t1 - tstart[sitenum1]);
             // Same model equation for DO < 2mg/L
             dopred2 = dmax[sitenum2] + d2n[sitenum2] .*(t2 - tstart[sitenum2]);
         }
         model {
             muchl ~ normal(0,3);
             sigchl ~ cauchy(0,3);
             etachl ~ normal(0,1);
             sig ~ cauchy(0,3);

    // turnover time model coefficients:
             b[1] ~ normal(-0.986, 0.245);
	     b[2] ~ normal(0.114, 0.018);
	     b[3] ~ normal(0.151, 0.019);
	     b[4] ~ normal(1.022, 0.047);
             sigt ~ cauchy(0,3);
             etat ~ normal(0,1);
             etat0 ~ normal(0,1);

             mud ~ normal(0,2);

             // sampling statements
             chl ~ normal(chlsite[sitenum], sigchl[2]);
             domean ~ normal(dopred1,  sig);
             domean2 ~ normal(dopred2, sig);
             domean0 ~ normal(dopred0, sig);
         }
