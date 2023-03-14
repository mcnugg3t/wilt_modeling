require(tidyverse)
require(ggpubr)
#'
#'
#'
var_vs_response_ <- function(dat.df, resp.name, var.name ) {
  dat.subs = dat.df |> 
    select(all_of(c(resp.name, var.name) ) ) |> 
    mutate(tile5 = ntile(!!sym(var.name), 5),
           tile10 = ntile(!!sym(var.name), 10))
  
  left.plot = dat.subs |> 
    group_by(tile5) |> 
    summarise(rate_infect = !!sym(resp.name) / n() ) |> 
    ggplot(aes(x=tile5, y=rate_infect )) + geom_col()
  
  right.plot = dat.subs |> 
    group_by(tile10) |> 
    summarise(rate_infect = !!sym(resp.name) / n() ) |> 
    ggplot(aes(x= tile10, y= rate_infect )) + geom_col()
      
  fig <- ggarrange(left.plot, right.plot, labels=c(paste0(var.name, " 5tile"), paste0(var.name, " 10tile")))
}