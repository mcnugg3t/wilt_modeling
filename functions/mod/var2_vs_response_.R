#'
#'
#'
var2_vs_response_ <- function(dat.df, resp.name, var1, var2 ) {
  dat.subs = dat.df |> 
    select( all_of( c(resp.name, var1, var2) ) )
  
  pts.dat = dat.df |> 
    dplyr::filter( sym(resp.name) > 0 )
  
  ggplot(data=dat.subs, aes({{var1}}, {{var2}})) +
    geom_hex() +
    geom_point(data=pts.dat,
               color="red")
}