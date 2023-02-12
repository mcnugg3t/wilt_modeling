require(assertthat)
require(crayon)
# {# debug
#   dat.df = extr.join
#   i = 6
# }
#'
#'
#'
add_interact_ <- function(dat.df, interact.v, verbose=T, DBG=F) {
  if(verbose) cat(crayon::bgGreen("\n\tFUNCT : add_interact_"))
  for(i in seq_along(interact.v)) {
    if(DBG) cat(paste0("\n\t\ti = ", i, "\t", interact.v[i] ) )
    term.tmp = interact.v[i] |> str_split(pattern=" x ", simplify=T) # extract term, split on 
    term1 = term.tmp[1]; term2 = term.tmp[2]
    
    t1.v = dat.df[[term1]]; t2.v = dat.df[[term2]] # extract vectors
    t1.v[is.na(t1.v)] <- 0; t2.v[is.na(t2.v)] <- 0
    # check numeric
    assert_that( (class(t1.v) == "numeric") & 
                 (class(t2.v) == "numeric"))
    t1.v = t1.v - min(t1.v, na.rm=T); t2.v = t2.v - min(t2.v, na.rm=T) # adjust each by minimum
    
    dat.df <- dat.df |> 
      add_column(!!sym(interact.v[i]) := t1.v * t2.v )
    
  }
  
  return(dat.df)
}