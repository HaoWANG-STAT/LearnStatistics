# https://pages.github.roche.com/adaptR/adaptR-tutorials/questionsAnswers.html
# https://cran.r-project.org/web/packages/gscounts/vignettes/Introduction_gscounts.html
library(gscounts)
SS <- gscounts::design_gsnb(rate1 = 1.2, rate2 = 1.5, dispersion = 0.3, ratio_H0 = 1, 
                    power = 0.89, sig_level = 0.025, timing = c(0.5, 1), 
                    esf = obrien, followup_max = 1, random_ratio = 1)
SS$n1


SS2 <- gscounts::design_gsnb(rate1 = 0.6*0.58, rate2 = 0.6, dispersion = 0.0001, ratio_H0 = 1, 
                            power = 0.8, sig_level = 0.025, timing = c(0.0001, 1), 
                            esf = obrien, followup_max = 1, random_ratio = 1)
SS2
