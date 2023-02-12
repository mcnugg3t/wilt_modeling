require(tidyverse)
require(assertthat)
#'
#'
calc_surprisal_tab_ <- function(tab.in, sample.v.in) {
  assert_that(class(tab.in) == "table")
  tab.tbl <- tab.in |> 
    as_tibble() |> 
    mutate(prob = n/sum(n))
  # if(class(sample.v.in) == "factor") {
  #   tab.tbl <- tab.tbl |> 
  #     mutate(val.v = as.numeric(val.v))
  # }
  surprisal.v <- sample.v.in |> 
    vapply(FUN = function(x) {
      ind.tmp <- which(tab.tbl[[1]] == x)
      prob <- tab.tbl[["prob"]][ind.tmp]
      return(log(prob))
    }, FUN.VALUE = numeric(1))
  return(surprisal.v)
}

# test
# {
# tab.in <- readRDS("tmp/tabl_tst.Rds")
# sample.v.in <- ow.pts.dat |> select(SOIL_NAME) |> as.vector()
# tab.in |> str()
# }
