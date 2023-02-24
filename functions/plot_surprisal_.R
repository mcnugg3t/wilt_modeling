#'
#'
#'
plot_surprisal_ <- function(plot.dat) {
  var.v <- plot.dat$var |> unique()
  # adjust infinite variables
  plot.dat <- plot.dat |> 
    mutate(surprisal = if_else(is.infinite(surprisal), -1*log(1e-7), surprisal))
  # loop over variables
  for(i in seq_along(var.v)) {
    { # prep for plot
      var.tmp <- var.v[i]
      cat("\nViewing : ")
      cat(crayon::bgBlue(var.tmp))
      p1.dat <- plot.dat |> 
        filter(var == var.tmp)
      
      surp.v <- p1.dat$surprisal
      mean.IC = mean(surp.v)
      med.IC = median(surp.v)
    }
    
    
    { # plot
      plot.tmp <- p1.dat |> 
        ggplot(aes(y= surprisal, group=bw, fill=bw )) + 
        theme(
          plot.title = element_text(size=24),
          axis.text.x = element_text(size=14, angle=0),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
        ) + 
        geom_boxplot() +
        labs(title = paste0(var.tmp, " mean IC = ", mean.IC), x = "bandwidth (m)", y="log(Ix)" )
      plot(plot.tmp)
    }
    
    
    cat("\nPress any key for next...")
    prpt <- readline(prompt=" ")
  }
}
