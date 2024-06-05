#####################
## Suppl. Figure 1 ##
#####################
library(MOFA2)
library(ggplot)


plot_variance_explained(MOFA2_model.trained, plot_total = TRUE)[[2]] + 
    scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70)) +
    theme_bw() + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title = element_text(size=25),
      axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 25)
    ) 