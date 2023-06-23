##################################
#### configurations for plots ####
##################################
twocolors <- c("#fec44f", "#43a2ca")    
palcolors <- colorRampPalette(brewer.pal(11, "Spectral"))(16)
alpha_hgh <- 0.8
alpha_med <- 0.2 
alpha_low <- 0.1
greycol     <- "grey50"
greycol_lgt <- "grey80"
fontfam     <- "Segoe UI Light"
iv_type <- "median_qi"
iv_lims <- c(0.5,0.9)

#################################
#### miscellaneous functions ####
#################################
# factor coding and ordering in one function (https://rdrr.io/github/svmiller/stevemisc/src/R/misc.R)
fct_reorg <- function(fac, ...) {
  dots <- unname(list(...))
  fac <- do.call("fct_relevel", c(list(fac), dots))
  fac <- fct_recode(fac, ...)
  fac
}

#############################
##### modelling function ####
#############################
run_model <- function(
  data,
  outcome_var         = NULL,
  intact_vars         = NULL,
  fin_inst_env        = NULL,
  macroeco_env        = NULL,
  orgnz_charct        = NULL,
  nominal_vars        = NULL,
  random_itcpt        = NULL,
  model_phi_group     = FALSE,
  ranslopes           = FALSE,
  priorset            = 1,
  modelname           = NULL,
  warmup              = NULL,
  iter                = NULL,
  chains              = NULL, 
  withincores         = NULL,
  cores               = NULL,
  seed                = NULL,
  get_prior_sims_too  = FALSE
){
  
  if(isTRUE(get_prior_sims_too) & priorset != 1) stop("priors are only simulated for the main prior set with current code")
  if(isTRUE(model_phi_group) & priorset != 1) stop("phi is only modeled for the main prior set with current code")
  if(is.null(modelname)) stop("specify a model name")
  
  ##############################################
  # order interacted predictors alphabetically #
  ##############################################
  intact_vars <- sort(intact_vars)
  
  ###################
  # create formulas #
  ####################
  if(is.null(random_itcpt)){
    formula_eta <- paste0(outcome_var, " ~ share_grp + cooperation + ", 
                          
                          paste(fin_inst_env, collapse = " + "), " + ",
                          paste(macroeco_env, collapse = " + "), " + ",
                          paste(orgnz_charct, collapse = " + "), " + ",
                          paste(nominal_vars, collapse = " + "), " + ",
                          paste(paste0("share_grp:", intact_vars), collapse = " + "))
    
  } else{
    if(isFALSE(ranslopes)){
      formula_eta <- paste0(outcome_var, " ~ share_grp + cooperation + ", 
                            "(1|", paste0(random_itcpt, collapse = "/"), ") + ",
                            paste(fin_inst_env, collapse = " + "), " + ",
                            paste(macroeco_env, collapse = " + "), " + ",
                            paste(orgnz_charct, collapse = " + "), " + ",
                            paste(nominal_vars, collapse = " + "), " + ",
                            paste(paste0("share_grp:", intact_vars), collapse = " + ")
      )
    } else{
      formula_eta <- paste0(outcome_var, " ~ share_grp + cooperation + ", 
                            "(1|ling_cluster) + (1 + share_grp|ling_cluster:country.name) +",
                            paste(fin_inst_env, collapse = " + "), " + ",
                            paste(macroeco_env, collapse = " + "), " + ",
                            paste(orgnz_charct, collapse = " + "), " + ",
                            paste(nominal_vars, collapse = " + "), " + ",
                            paste(paste0("share_grp:", intact_vars), collapse = " + ")
      )
    }
  }
  if(isTRUE(model_phi_group)) formula_phi <- "phi ~ (1|country.name)"
  
  ###############
  # filter data #
  ###############
  data_model_unscaled <- data %>% 
    drop_na(all_of(c("cooperation",
                     "share_grp",
                     outcome_var,
                     intact_vars,
                     fin_inst_env,
                     macroeco_env,
                     orgnz_charct,
                     nominal_vars,
                     random_itcpt)))
  
  ###################
  # create clusters #
  ###################
  ld <- read_csv("data/LingDistance.csv")
  ld <- ld[,-1]
  country <- countrycode(names(ld),origin = "iso2c", destination = "country.name", custom_match = c("XK" = "Kosovo"))
  names(ld) <- country
  
  ld <- ld[names(ld) %in% data_model_unscaled$country.name, names(ld) %in% data_model_unscaled$country.name]
  
  clusterings <- list(NA)
  meanS       <- as.numeric(NA)
  nr_singles  <- as.numeric(NA)
  m           <- seq(1, nrow(ld)-1, 1)
  
  for(i in m){
    temp <- pam(x = ld,
                k = i,
                diss = TRUE)
    
    clusterings[[i]]  <- temp
    try(meanS[i]      <- mean(silhouette(temp)[,"sil_width"]), silent = TRUE)
    try(nr_singles[i] <- sum(table(silhouette(temp)[,"cluster"])==1), silent = TRUE)
  }
  
  clusts_ld <- clusterings[[which.max(meanS * (nr_singles==0))]]$clustering
  df_ld <- tibble(
    country.name = names(clusts_ld),
    ling_cluster = paste0("cluster_", sprintf("%02d", unname(clusts_ld)))
  )
  data_model_unscaled <- data_model_unscaled %>% 
    select(-ling_cluster) %>% 
    left_join(df_ld, by = "country.name")
  
  #####################################################################################
  # save the data tables which go into the model prior to scaling and transformations #
  #####################################################################################
  write_csv(data_model_unscaled, paste0("modelfiles/data_model_unscaled_", modelname, ".csv"))
  write_csv(df_ld, paste0("modelfiles/df_ld_", modelname, ".csv"))
  write_csv(ld, paste0("modelfiles/ld_", modelname, ".csv"))
  
  ########################
  # transform predictors #
  ########################
  data_model <- data_model_unscaled %>% 
    mutate(across(.cols = all_of(c(random_itcpt,
                                   nominal_vars)),
                  factor),
           across(.cols = all_of(c("cooperation",
                                   intact_vars,
                                   fin_inst_env,
                                   macroeco_env,
                                   orgnz_charct)),
                  scale),
           hdi_sq = hdi^2,
           wdi_gdp_sq = wdi_gdp^2) %>% 
    select(all_of(c("cooperation",
                    "share_grp",
                    outcome_var,
                    intact_vars,
                    fin_inst_env,
                    macroeco_env,
                    orgnz_charct,
                    nominal_vars,
                    random_itcpt))
    )
  
  ###########################
  # set prior distributions #
  ###########################
  if(priorset == 1){
    # coefficients
    scal_coef_prior  <- 1.5
    # variance components
    scal_sd_prior    <- 1 
    df_sd_prior      <- 3
    # global intercept
    p_0_mean  <- round(qlogis(median(data_model %>% pull(!!as.symbol(outcome_var)))), 1)

    # prior for group lending share (because SD != 1)
    p_1_sd    <- round(scal_coef_prior / sd(data_model$share_grp), 1)
    p_1_prior <- c(prior_string(paste0("normal(0, ", p_1_sd, ")"), class = "b", coef = "share_grp"))
    
    # prior for squared terms (because SD != 1)
    if("hdi_sq" %in% macroeco_env){
      p_2_sd       <- round(scal_coef_prior / sd(data_model$hdi_sq), 1)
      p_2_prior    <- c(prior_string(paste0("normal(0, ", p_2_sd, ")"), class = "b", coef = "hdi_sq"))
      
    } else if("wdi_gdp_sq" %in% macroeco_env){
      p_2_sd       <- round(scal_coef_prior / sd(data_model$wdi_gdp_sq), 1)
      p_2_prior    <- c(prior_string(paste0("normal(0, ", p_2_sd, ")"), class = "b", coef = "wdi_gdp_sq"))
    }
    
    # prior for interactions  (because SD != 1)
    intact_terms <- paste0("share_grp:", intact_vars)
    p_3_sd       <- round(scal_coef_prior / apply(data_model$share_grp * data_model[,intact_vars], 2, sd) /2, 1) # "/2" for more stringent restriction on interactions
    p_3_prior    <- c(prior_string(paste0("normal(0, ", p_3_sd, ")"), class = "b", coef = intact_terms))
  } 
  
  cat("observations going into the model:", nrow(data_model), "\n")
  if("mfiname" %in% random_itcpt) cat("mfi intercepts going into the model:", length(unique(data_model$mfiname)), "\n")
  cat("country intercepts going into the model:", length(unique(data_model$country.name)), "\n")
  cat("cultural cluster intercepts going into the model:", length(unique(data_model$ling_cluster)), "\n")
  
  if(priorset == 1){
    if(is.null(random_itcpt)){
      extra_priors <- 
        c(p_1_prior,
          p_2_prior,
          p_3_prior
        )
    } else{
      if ("mfiname" %in% random_itcpt | isTRUE(ranslopes)) {
        sd_prior <- c(prior_string(paste0("student_t(", df_sd_prior, ", 0, ",  scal_sd_prior/2, ")"), class = "sd"))
        if(isTRUE(model_phi_group)){
          sd_prior <- c(sd_prior, prior_string(paste0("student_t(", df_sd_prior, ", 0, ",  scal_sd_prior/2, ")"), class = "sd", dpar = "phi"))
        }
        
      } else {
        sd_prior <- c(prior_string(paste0("student_t(", df_sd_prior, ", 0, ",  scal_sd_prior, ")"), class = "sd"))
        if(isTRUE(model_phi_group)){
          sd_prior <- c(sd_prior, prior_string(paste0("student_t(", df_sd_prior, ", 0, ",  scal_sd_prior, ")"), class = "sd", dpar = "phi"))
        }
      }
      
      extra_priors <- 
        c(sd_prior, 
          p_1_prior,
          p_2_prior,
          p_3_prior
        )
    }
    
    #################
    # fit the model #
    #################
    if(!model_phi_group){
      modelobj  <- ordbetareg(
        # model and data
        bf(formula_eta, center = TRUE),
        data = data_model,
        # technical
        control   = list(adapt_delta = 0.975),
        threads   = threading(withincores),
        warmup    = warmup,
        iter      = iter,
        thin      = 1,
        chains    = chains, 
        cores     = cores, 
        seed      = seed,
        backend   = "cmdstanr",
        # priors
        dirichlet_prior      = c(1,1,1),        
        coef_prior_mean      = 0,               
        coef_prior_SD        = scal_coef_prior, 
        intercept_prior_mean = p_0_mean,        
        intercept_prior_SD   = 5,               
        phi_prior            = 0.1,             
        extra_prior = extra_priors
      )
    } else{
      modelobj  <- ordbetareg(
        # model and data
        bf(formula_eta, formula_phi, center = TRUE),
        data = data_model,
        # technical
        control   = list(adapt_delta = 0.975),
        threads   = threading(withincores),
        warmup    = warmup,
        iter      = iter,
        thin      = 1,
        chains    = chains, 
        cores     = cores, 
        seed      = seed,
        backend   = "cmdstanr",
        phi_reg   = "both",
        # priors
        dirichlet_prior      = c(1,1,1),        
        coef_prior_mean      = 0,               
        coef_prior_SD        = scal_coef_prior, 
        intercept_prior_mean = p_0_mean,        
        intercept_prior_SD   = 5,               
        phi_intercept_prior_mean = 0,           
        phi_intercept_prior_SD = 5,             
        extra_prior = extra_priors
      )
    }
    
  } else if(priorset == 2){
    
    modelobj  <- ordbetareg(
      # model and data
      bf(formula_eta, 
         center = TRUE),
      data = data_model,
      # technical
      control   = list(adapt_delta = 0.975),
      threads   = threading(withincores),
      warmup    = warmup,
      iter      = iter,
      thin      = 1,
      chains    = chains, 
      cores     = cores, 
      seed      = seed,
      backend   = "cmdstanr",
      # priors
      coef_prior_SD = 15,
      extra_prior = 
        c(prior_string(paste0("student_t(", 1, ", 0, ",  5, ")"), class = "sd"))
    )
  }
  
  if(isTRUE(get_prior_sims_too)){
    priorobj  <- ordbetareg(
      # model and data
      bf(formula_eta, 
         center = TRUE),
      data = data_model,
      # technical
      control   = list(adapt_delta = 0.85),
      threads   = threading(withincores),
      warmup    = warmup,
      iter      = iter,
      thin      = 1,
      chains    = chains, 
      cores     = cores, 
      seed      = seed,
      backend   = "cmdstanr",
      sample_prior = "only",
      # priors
      dirichlet_prior      = c(1,1,1),         
      coef_prior_mean      = 0,                
      coef_prior_SD        = scal_coef_prior,  
      intercept_prior_mean = p_0_mean,         
      intercept_prior_SD   = 5,                
      phi_prior            = 0.1,              
      extra_prior = extra_priors
    )
    saveRDS(priorobj, file = paste0("modelfiles/prior_m_", modelname, ".rds"))
  }
  
  ###########################
  #### save model object ####
  ###########################
  saveRDS(modelobj, file = paste0("modelfiles/m_", modelname, ".rds"))
}

##################################################
#### slightly edited function from ordbetareg ####
##################################################
# (to allow for varying intercept-only model for precision parameter)
.load_ordbetareg_custom <- function(beta_prior=NULL,
                                    intercept_prior=NULL,
                                    phi_reg=NULL,
                                    phireg_beta_prior=NULL,
                                    phireg_intercept_prior=NULL,
                                    dirichlet_prior_num=NULL,
                                    phi_prior=NULL,
                                    extra_prior=NULL,
                                    suffix="",
                                    formula=NULL) {
  
  # custom family
  ord_beta_reg <- custom_family("ord_beta_reg",
                                dpars=c("mu","phi","cutzero","cutone"),
                                links=c("identity","log",NA,NA),
                                lb=c(NA,0,NA,NA),
                                type="real")
  
  # stan code for density of the model
  
  stan_funs <- "

    real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu  - thresh[2]);
    } else {
      return log_diff_exp(log_inv_logit(mu   - thresh[1]), log_inv_logit(mu - thresh[2])) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }"
  
  stanvars <- stanvar(scode=stan_funs,block="functions")
  
  # For pulling posterior predictions
  
  posterior_predict_ord_beta_reg <- function(i, draws, ...) {
    
    get_args <- list(...)
    
    if(!is.null(get_args$ntrys) && get_args$ntrys>5) {
      
      type <- "continuous"
      
    } else {
      
      type <- "combined"
      
    }
    
    mu <- brms::get_dpar(draws, "mu", i = i)
    phi <- brms::get_dpar(draws, "phi", i = i)
    cutzero <- brms::get_dpar(draws, "cutzero", i = i)
    cutone <- brms::get_dpar(draws, "cutone", i = i)
    N <- draws$ndraws
    
    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)
    
    pr_y0 <- 1 - plogis(mu - thresh1)
    pr_y1 <- plogis(mu - thresh2)
    pr_cont <- plogis(mu-thresh1) - plogis(mu - thresh2)
    out_beta <- rbeta(n=N,plogis(mu)*phi,(1-plogis(mu))*phi)
    
    # now determine which one we get for each observation
    outcomes <- sapply(1:N, function(i) {
      sample(1:3,size=1,prob=c(pr_y0[i],pr_cont[i],pr_y1[i]))
    })
    
    if(type=="combined") {
      
      final_out <- sapply(1:length(outcomes),function(i) {
        if(outcomes[i]==1) {
          return(0)
        } else if(outcomes[i]==2) {
          return(out_beta[i])
        } else {
          return(1)
        }
      })
      
    } else if (type=="continuous") {
      
      final_out <- out_beta
      
    }
    
    final_out
    
  }
  
  # for calculating marginal effects/conditional expectations
  
  posterior_epred_ord_beta_reg<- function(draws) {
    
    cutzero <- brms::get_dpar(draws, "cutzero")
    cutone <- brms::get_dpar(draws, "cutone")
    
    mu <- brms::get_dpar(draws, "mu")
    
    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)
    
    low <- 1 - plogis(mu - thresh1)
    middle <- plogis(mu-thresh1) - plogis(mu - thresh2)
    high <- plogis(mu - thresh2)
    
    low*0 + middle*plogis(mu) + high
  }
  
  # for calcuating LOO and Bayes Factors
  
  log_lik_ord_beta_reg <- function(i, draws) {
    
    mu <- brms::get_dpar(draws, "mu", i = i)
    phi <- brms::get_dpar(draws, "phi", i = i)
    y <- draws$data$Y[i]
    cutzero <- brms::get_dpar(draws, "cutzero", i = i)
    cutone <- brms::get_dpar(draws, "cutone", i = i)
    
    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)
    
    if(y==0) {
      out <- log(1 - plogis(mu - thresh1))
    } else if(y==1) {
      out <- log(plogis(mu - thresh2))
    } else {
      out <- log(plogis(mu-thresh1) - plogis(mu - thresh2)) + dbeta(y,plogis(mu)*phi,(1-plogis(mu))*phi,log=T)
    }
    
    out
    
  }
  
  ###### Code declaring induced dirichlet prior ####
  # code from Michael Betancourt/Staffan Betner
  # discussion here: https://discourse.mc-stan.org/t/dirichlet-prior-on-ordinal-regression-cutpoints-in-brms/20640
  dirichlet_prior <- "
  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    if(cutnum==1) {

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    // divide in half for the two cutpoints

    // don't forget the ordered transformation

      return   dirichlet_lpdf(p | alpha)
           + log_determinant(J) + cut2;

    } else {

      return(0);

    }


  }

  real induced_dirichlet_rng(vector alpha, real phi, int cutnum, real cut1, real cut2) {

    int K = num_elements(alpha);
    vector[K] p;
    vector[K-1] cutpoints;

    // need to reverse the steps
    // first get the dirichlet probabilities conditional on alpha

    p = dirichlet_rng(alpha);

    // then do the *reverse* transformation to get cutpoints

    for(k in 1:(K-1)) {

       if(k==1) {

          cutpoints[k] = phi - logit(1 - p[k]);

       } else {

          cutpoints[k] = phi - logit(inv_logit(phi - cutpoints[k-1]) - p[k]);

       }

    }

    return  cutpoints[cutnum];
  }
"
  dirichlet_prior_stanvar <- stanvar(scode = dirichlet_prior, block = "functions")
  
  # stanvar(scode = "ordered[2] thresh;
  #             thresh[1] = cutzero;
  #             thresh[2] = cutzero+exp(cutone);",
  #         block = "tparameters") -> # there might be a better way to specify this
  #   dirichlet_prior_ordbeta_stanvar
  
  stanvars <- stanvars + dirichlet_prior_stanvar
  
  # Feel free to add any other priors / change the priors on b,
  # which represent regression coefficients on the logit
  # scale
  
  
  # Set priors --------------------------------------------------------------
  
  if(length(suffix)>1) {
    
    # multiple outcomes
    
    cutzero <- paste0("cutzero",suffix)
    cutone <- paste0("cutone",suffix)
    
    priors <- set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                            collapse=","),"]', 0, 1,", cutzero[1],",", cutone[1],")"),
                        class="cutzero",resp=substring(suffix[1],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero[1],",", cutone[1],")"),
                class="cutone",resp=substring(suffix[1],2)) +
      set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b",resp=substring(suffix[1],2)) +
      set_prior(paste0("exponential(",phi_prior,")"),class="phi",resp=substring(suffix[1],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 1,", cutzero[2],",", cutone[2],")"),
                class="cutzero",resp=substring(suffix[2],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero[2],",", cutone[2],")"),
                class="cutone",resp=substring(suffix[2],2)) +
      set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b",resp=substring(suffix[2],2)) +
      set_prior(paste0("exponential(",phi_prior,")"),class="phi",resp=substring(suffix[2],2))
    
  } else {
    
    
    cutzero <- paste0("cutzero",suffix)
    cutone <- paste0("cutone",suffix)
    
    priors <- set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                            collapse=","),"]', 0, 1,", cutzero,",", cutone,")"),
                        class="cutzero") +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero,",", cutone,")"),
                class="cutone")
    
    # only do phi reg priors for univariate models
    
    if(phi_reg=='both') {
      
      # commented out, sg
      # priors <- priors + set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b") +
      #   set_prior(paste0("normal(",phireg_beta_prior[1],",",phireg_beta_prior[2],")"),class="b",dpar="phi")
      
    } else if(phi_reg=='none') {
      
      priors <- priors + set_prior(paste0("exponential(",phi_prior,")"),class="phi") +
        set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b")
      
    } else if(phi_reg=='only') {
      
      priors <- priors + set_prior(paste0("normal(",phireg_beta_prior[1],",",phireg_beta_prior[2],")"),class="b",dpar="phi")
      
    }
    
  }
  
  
  
  if(!is.null(extra_prior)) {
    
    priors <- priors + extra_prior
    
  }
  
  if(!is.null(intercept_prior)) {
    
    # different priors with/without centering
    
    if(attr(formula$formula, "center")) {
      
      priors <- priors + set_prior(paste0("normal(",intercept_prior[1],",",intercept_prior[2],")"),
                                   class="Intercept")
      
    } else {
      
      priors <- priors + set_prior(paste0("normal(",intercept_prior[1],",",intercept_prior[2],")"),
                                   coef="Intercept",class="b")
      
    }
    
  }
  
  if(!is.null(phireg_intercept_prior) && phi_reg %in% c("only","both")) {
    
    if(attr(formula$formula, "center")) {
      
      priors<- priors + set_prior(paste0("normal(",phireg_intercept_prior[1],",",phireg_intercept_prior[2],")"),
                                  class="Intercept",dpar="phi")
    } else {
      
      priors<- priors + set_prior(paste0("normal(",phireg_intercept_prior[1],",",phireg_intercept_prior[2],")"),
                                  class="Intercept",
                                  dpar="phi")
      
    }
    
  }
  
  return(list(priors=priors,
              stanvars=stanvars,
              log_lik=log_lik_ord_beta_reg,
              posterior_epred=posterior_epred_ord_beta_reg,
              stan_funs=stan_funs,
              family=ord_beta_reg))
  
  
}
environment(.load_ordbetareg_custom) <- asNamespace('ordbetareg')
assignInNamespace(".load_ordbetareg", .load_ordbetareg_custom, ns = "ordbetareg")