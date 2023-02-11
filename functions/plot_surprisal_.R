#'
#'
#'
plot_surprisal_ <- function() {
  var.v <- return.dat$var |> unique(); var.v
  
  # infinite values
  return.dat |> 
    group_by(var, bw) |> 
    summarise(count_inf = sum(is.infinite(surprisal))) |> 
    arrange(desc(count_inf))
  
  # max non-infinite value
  max.non.inf <- return.dat |> 
    filter(!is.infinite(surprisal)) |> 
    select(surprisal) |> 
    max(); max.non.inf
  
  replace.inf <- max.non.inf
  
  # adjust infinite values to slightly above max
  plot.dat <- return.dat |> 
    mutate(surprisal = if_else(is.infinite(surprisal), replace.inf + runif(n=1, min=2, max=5), surprisal))
  
  for(i in seq_along(var.v)) {
    
    { # prep for plot
      var.tmp <- var.v[i]
      cat("\nViewing : ")
      cat(crayon::bgBlue(var.tmp))
      p1.dat <- plot.dat |> 
        filter(var == var.tmp) |> 
        mutate(
          plt.x = bw + runif( n=nrow(p1.dat), min=-12.5, max=12.5 ),
          plt.y = surprisal )
      
      surp.v <- p1.dat$surprisal
      mean.surprisal <- mean(surp.v, na.rm=T); mean.surprisal
      median.surprisal = median(log(surp.v), na.rm=T); median.surprisal
      
      central.dat <- p1.dat |> 
        group_by(bw) |> 
        summarise(mean = mean(plt.y),
                  median = median(plt.y)) |>
        mutate(skew = (mean-median)^2) |> 
        ungroup() |> 
        pivot_longer(cols=2:3, names_to="centrality", values_to = "val")
    }
    
    
    { # plot
      clrs <- c("mean" = "violet", "median" = "deeppink")
      plot.tmp <- p1.dat |> 
        ggplot(aes(x=plt.x, y= log(plt.y) )) + 
        theme(
          plot.title = element_text(size=24),
          axis.text.x = element_text(size=14, angle=0),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="black")
        ) + 
        geom_hex(bins=30) +
        geom_point(data=central.dat, aes(x=bw, y=log(val), size=skew, color=centrality)) + 
        labs(title = paste0(var.tmp, " mean IC = ", round(mean.surprisal, 1) ), x = "bandwidth", y="log(Ix)" ) +
        scale_color_manual(values = clrs) +
        scale_fill_viridis_c(option="D", name="frequency")
      plot(plot.tmp)
    }
    
    
    cat("\nPress any key for next...")
    prpt <- readline(prompt=" ")
  }
}
