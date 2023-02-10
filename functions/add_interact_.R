#'
#'
#'
add_interact_ <- function(dat.df, interact.v, verbose=T, DBG=F) {
  
  for(i in seq_along(interact.v)) {
    term.tmp = interact.v[i]
    term.v <- interact.v[i] |> str_split(pattern=" x ", simplify=T)
    term1 = term.v[1]; term2 = term.v[2]
    t1.v = dat.df[[term1]]; t2.v = dat.df[[term2]]
    t1.v = t1.v - min(t1.v); t2.v = t2.v - min(t2.v) # adjust by minimum
    
    dat.df <- dat.df |> 
      add_column(!!sym(term.tmp) := t1.v * t2.v )
    
  }
  
  return(dat.df)
}