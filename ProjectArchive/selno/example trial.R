n <- 2*c(10)
sigma <- 0.96
rho <- c(0.7)
reduCtrl <- c(0)
reduTest <- seq(0.5, 0.6, 0.05)

scenarios <- expand.grid(n = n, sigma=sigma, rho = rho, reduCtrl=reduCtrl, reduTest=reduTest) |> 
  asplit(MARGIN = 1) 

argsvec <- scenarios[[1]]
 set.seed(123456)
 list_of_defs <- s_define()
 dTrial <- s_generate(list_of_defs, argsvec)
 results <- s_model(dTrial)$results
 #s_single_rep(list_of_defs, argsvec)
 
 
 gm_mean = function(a){prod(a)^(1/length(a))}
 
 dTrial |> mutate(BL=exp(Y0),
                  EP=exp(Y1)) |> 
   ggplot() +
   #geom_point(aes(x=TRT, y=BL)) +
   geom_boxplot(aes(x=TRT, y=BL), alpha=0.8, fill="green") +
   stat_summary(aes(x=TRT, y=BL), fun = mean, geom = "point", col = "darkgreen", size=5) +
   stat_summary(aes(x=TRT, y=BL), fun = gm_mean, geom = "point", col = "darkgreen", size=5, shape =2) +
   #geom_point(aes(x=TRT, y=EP), color="red", shape=2) +
   geom_boxplot(aes(x=TRT, y=EP), alpha=0.2, fill="red")+
   stat_summary(aes(x=TRT, y=EP), fun = mean, geom = "point", col = "darkred", size=5) +
   stat_summary(aes(x=TRT, y=EP), fun = gm_mean, geom = "point", col = "darkred", size=5, shape =2) +
   labs(x = "treatment group", y = "Il1Beta", title = "Data on original scale")
   
 
 
 
   dTrial |> mutate(BL=exp(Y0),
                    EP=exp(Y1)) |> 
     ggplot() +
     #geom_point(aes(x=TRT, y=BL)) +
     geom_boxplot(aes(x=TRT, y=Y0), alpha=0.8, fill="green") +
     stat_summary(aes(x=TRT, y=Y0), fun = mean, geom = "point", col = "darkgreen", size = 5) +
     #geom_point(aes(x=TRT, y=EP), color="red", shape=2) +
     geom_boxplot(aes(x=TRT, y=Y1), alpha=0.2, fill="red")+
     stat_summary(aes(x=TRT, y=Y1), fun = gm_mean, geom = "point", col = "darkred", size = 5) +
   labs(x = "treatment group", y = "log-Il1Beta", title = "Data on log scale")

 