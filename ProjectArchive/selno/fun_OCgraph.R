#' Summarise the Operating Characteristics into a 3-Outcome (Go/NoGo/Eval) Figure
#' 
#' @param d the input data frame
#' @param NoGoThreshold (`integer`) The vector of No-Go threshold
#' @param GoThreshold (`integer`) The vector of Go threshold
library(plotly)
createOCSummaryFigure_3ODM <- function(
    d
) {
  
    f <-  d |> 
      select(n, reduCtrl, reduTest, ProbEval, ProbGoAc, effType) |> 
      mutate(NoGo = round(100-ProbEval-ProbGoAc, 1),
             reduTest = 100 * reduTest,
             reduCtrl = 100 * reduCtrl,
             n = as.factor(n))  |> 
      rename(Go = ProbGoAc, 
             Eval = ProbEval,
             nTotal = n) |>
      pivot_longer(cols = c(Go, NoGo, Eval),
                   values_to = "Prob",
                   names_to= "Decision"
      ) %>% 
      #mutate(Decision = factor(Decision, levels = c("Go", "Eval", "NoGo"), ordered = T)) |> 
      ggplot() +
      geom_line(aes(x = reduTest, y = Prob, colour = Decision, linetype = effType)) +
      scale_color_manual(values = c(Go = "darkgreen", Eval = "blue", NoGo = "red")) +
      #facet_wrap(vars(NTotal), labeller = as_labeller(function(z) paste0("N=",z), label_both)) +
      labs(
        x = "True reduction from baseline (%) for Selnoflast",
        colour = "Decision",
        y = "Probability",
        title = paste0("True reduction from baseline (%) for Placebo: ", d |> pull(reduCtrl) |> unique()*100, "%")
      ) +
      theme_minimal()
    
    ggplotly(f, tooltip = c("y", "x")) %>% layout(hovermode = 'x')
  }
  
