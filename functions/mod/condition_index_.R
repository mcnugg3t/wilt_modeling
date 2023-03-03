#'
#'
condition_index_ <- function(dat.mat.in) {
  cor.tmp <- cor(dat.mat.in)
  eig.tmp <- cor.tmp |> eigen()
  cond.ind <- sqrt(max(eig.tmp$values)/min(eig.tmp$values))
}
