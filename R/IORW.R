report_RR_CI <- function(mdl, digits=2){
  cbind(broom::tidy(mdl),
        broom::confint_tidy(mdl)) %>%
    mutate(output = sprintf("%.2f ( %.2f, %.2f)",
                            estimate %>% exp,
                            conf.low %>% exp,
                            conf.high %>% exp
    )
    ) %>%
    filter(term=="raceBlack race") %>%
    pull(output)
}

# Mediation
#' @import survey
IORW <- function(
  otcm,
  xpsr,
  mdtr,
  covar,
  data,
  method = c("IORW", "IOW"),
  btsp_it = 1000
){
   # covar_lst %>%
  # map_dfr(.f = function(covar, N_it = 1000, seed=1) {

    covar <- setdiff(covar, "race")


    cmplt_dat <- dat %>%
      select(outcome, covar, expsr, mdtr, "smp_weight", "inc_htn", "race", "gender") %>%
      filter(complete.cases(.))

    lgt_dgn <- svydesign(~1, weights = ~ smp_weight,
                         strata=~race*gender,
                         data = cmplt_dat %>% mutate(inc_htn = inc_htn %>% as.numeric()-1))


    # Logistic Model to calculate odd
    lgt_fmlr <- paste0(expsr, " ~ ", mdtr,
                       ifelse(!is.null(covar), " + ", ""),
                       paste0(covar, collapse=" + "))



    lgt_mdl <- svyglm(lgt_fmlr,
                      design = lgt_dgn, family = binomial(link="logit"))


    wgt <- 1/(predict(lgt_mdl) %>% as.numeric %>% exp)

    if(nrow(cmplt_dat)!=length(wgt)) stop("Un-equal number of observations")

    cmplt_dat$new_wgt <- cmplt_dat$smp_weight*wgt

    new_dgn <- svydesign(~1, weights = ~ new_wgt,
                         strata=~race*gender,
                         data = cmplt_dat %>% mutate(inc_htn = inc_htn %>% as.numeric()-1))

    mdl_fmlr <- paste0(outcome, " ~ ", expsr,
                       ifelse(!is.null(covar), " + ", ""),
                       paste0(covar, collapse=" + "))


    ttl_mdl <- svyglm(
      formula = mdl_fmlr,
      design = lgt_dgn, family=quasipoisson
    )

    mdtr_mdl <- svyglm(
      formula = mdl_fmlr,
      design = new_dgn, family=quasipoisson
    )

    #N <- ttl_mdl$data %>% nrow()
    # set.seed(seed)
    boot_res <- map_dfr(1:N_it, .f=function(x){
      boot_dat <- cmplt_dat %>% mutate(strata = paste(gender, race, sep="_")) %>% split(.$strata) %>%
        map_dfr(.f = function(dat){
          N_strat <- nrow(dat)
          boot_indices <- sample(1:N_strat, N_strat, replace = T)
          dat[boot_indices,,drop=FALSE]
        })

      boot_lgt_dgn <- svydesign(~1, weights = ~ smp_weight,
                                strata=~race*gender,
                                data = boot_dat %>% mutate(inc_htn = inc_htn %>% as.numeric()-1))

      boot_lgt_mdl <- svyglm(lgt_fmlr,
                             design = boot_lgt_dgn, family = binomial(link="logit"))


      boot_wgt <- 1/(predict(boot_lgt_mdl) %>% as.numeric %>% exp)

      if(nrow(boot_dat)!=length(boot_wgt)) stop("Un-equal number of observations")

      boot_dat$new_wgt <- boot_dat$smp_weight*boot_wgt

      boot_new_dgn <- svydesign(~1, weights = ~ new_wgt,
                                strata=~race*gender,
                                data = boot_dat %>% mutate(inc_htn = inc_htn %>% as.numeric()-1))



      boot_ttl_mdl <- svyglm(
        formula = mdl_fmlr,
        design = boot_lgt_dgn, family=quasipoisson
      )

      boot_mdtr_mdl <- svyglm(
        formula = mdl_fmlr,
        design = boot_new_dgn, family=quasipoisson
      )


      data.frame(
        ttl_beta = tidy(boot_ttl_mdl)%>%
          mutate(beta = estimate) %>%
          filter(term=="raceBlack race") %>%
          pull(beta),
        mdtr_beta = tidy(boot_mdtr_mdl)%>%
          mutate(beta = estimate) %>%
          filter(term=="raceBlack race") %>%
          pull(beta)
      )
    })



    ttl_beta <- tidy(ttl_mdl)%>%
      mutate(beta = estimate) %>%
      filter(term=="raceBlack race") %>%
      pull(beta)

    direct_beta <- tidy(mdtr_mdl)%>%
      mutate(beta = estimate) %>%
      filter(term=="raceBlack race") %>%
      pull(beta)

    indirect_beta <- ttl_beta-direct_beta

    perc_change_CI <- boot_res %>%
      mutate(
        indirect_eff = (ttl_beta-mdtr_beta) ,
        pct_change = (ttl_beta-mdtr_beta)/(ttl_beta)
      ) %>%
      summarize(
        indirect_p = sum(indirect_eff>=indirect_beta)/n(),
        indirect_eff_CI_LL = quantile(indirect_eff, 0.025),
        indirect_eff_CI_UL = quantile(indirect_eff, 0.975),
        CI_LL = quantile(pct_change, 0.025)*100,
        CI_UL = quantile(pct_change, 0.975)*100
      )

    data.frame(
      #N,
      ttl_eff = report_RR_CI(ttl_mdl),
      ttl_p = sprintf("%.3f", broom::tidy(ttl_mdl) %>%
                        filter(term=="raceBlack race") %>%
                        pull(p.value)),

      direct_eff = report_RR_CI(mdtr_mdl),
      direct_p = sprintf("%.3f", broom::tidy(mdtr_mdl) %>%
                           filter(term=="raceBlack race") %>%
                           pull(p.value)),

      indirect_eff = sprintf("%.2f ( %.2f, %.2f)",
                             (ttl_beta-direct_beta) %>% exp,
                             perc_change_CI$indirect_eff_CI_LL %>% exp,
                             perc_change_CI$indirect_eff_CI_UL %>% exp),

      # TODO: Start here, calculate P-value in the bootstrap
      # indirect_p = sprintf("%.3f", perc_change_CI$indirect_p),
      # TODO: Add CI
      Mediation_Percentage = sprintf("%.1f%%", (ttl_beta-direct_beta)/ttl_beta*100#,
                                     # perc_change_CI$CI_LL,
                                     # perc_change_CI$CI_UL
                                     ),
      stringsAsFactors = F
    )}#,
    # .id = "Model") %>%
  # kable(format="pandoc")
}
