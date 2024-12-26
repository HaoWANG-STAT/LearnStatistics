
x1 <- rnorm(10, 1,1)
x0 <- rnorm(10, 2,1)

# differentiate percentage of mean change from baseline vs. mean of percentage change from baseline
tibble(x1, x0) %>% 
  mutate(cfb=x1-x0,
         pctcfb=cfb/x0) %>% 
  summarise(mbl=mean(x0), mcfb=mean(cfb), 
            mpctcfb=mean(pctcfb)) %>% 
  mutate(pctmcfb=mcfb/mbl)


# calculating change first then averaging is same with calculate mean first than difference
tibble(x1, x0) %>% summarise(m1=mean(x1), mbl=mean(x0)) %>% mutate(mcfb=m1-mbl,
                                                                  pctmcfb=mcfb/mbl)
