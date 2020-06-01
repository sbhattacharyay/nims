plot_metric_stripcharts <- function(metric_set,metric_name,file_name,height = 500, width = 650){
  
  rootDir<-generateRootDir('metric_stripcharts')
  
  if ("glmnet" %in% metric_set$classifier){
    metric_set[metric_set == "glmnet"] = "glm"
  }
  
  theme_set(theme_classic())
  
  viz<-ggplot(metric_set, aes(x = classifier, y = values,color = pred_space,shape=pred_space)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = .3,
                                      jitter.height = 0,
                                      dodge.width = .5),
      size = 1.5
    ) +
    labs(x = "Classifier Type", y = metric_name) +
    theme(axis.text.x = element_text(angle = 0, hjust = .5),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill=NA, size=1)) + 
    stat_summary(
      aes(color = pred_space),
      shape = 9,
      fun.data="mean_sdl",  fun.args = list(mult=1), 
      geom = "pointrange",  size = .6,
      position = position_dodge(0.5),
    ) + 
    expand_limits(y = c(0, 1))
  
  file_string <- file.path(rootDir,file_name)
  png(file = file_string,
      width = width,
      height = height)
  print(viz)
  dev.off()
}