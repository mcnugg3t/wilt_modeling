# filter to all vars and top 10 interaction terms (using median ic)
filter_topN_interacts_ <- function(dat.in, n.top) {
  var.unique <- dat.in$var |> # unique variables
    unique()
  int.ind <- var.unique |>   # indices of those variables with interact _ x _
    str_detect(" x ")
  interact.vars <- var.unique[int.ind]; interact.vars # subset interact terms
  noninteract.vars <- var.unique[!int.ind]; noninteract.vars # and non-interact terms
  top.n.interact <- dat.in |> # get var names for top 10 interact terms
    filter(var %in% interact.vars) |> 
    group_by(var) |> 
    summarise(med_ic = median(surprisal, na.rm=T)) |> 
    arrange(desc(med_ic)) |> 
    select(var) |> 
    slice(1:n.top) |> 
    pull(var)
  return.dat <- dat.in |> # filter - either a plain covariate or one of the top 10 interaction terms
    filter(var %in% noninteract.vars | var %in% top.n.interact) |> 
    filter(!(var %in% c("SLOPE_CLAS", "SOIL_NAME", "EROSION_DI")))
  return(return.dat)
}

