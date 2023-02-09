resid_plot_fun <- function(pl) {
  grid.arrange(pl[[1]], pl[[2]], nrow=1) # fitted value & raw residual
  grid.arrange(pl[[1]], pl[[3]], nrow=1) # fitted value & pearson residuals
  grid.arrange(pl[[4]], pl[[5]], nrow=1) # raw resid vs fitted, pears resid vs fitted
  grid.arrange(pl[[6]], pl[[7]], nrow=1) # pearson residuals vs oak_prob, elev
  grid.arrange(pl[[8]], pl[[9]], pl[[10]], pl[[11]], pl[[12]], pl[[13]], nrow=2)
}