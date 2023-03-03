vis_1var_ <- function(in.dat, in.var, median.ic.in) {
  return <- in.dat |> 
    filter(var == in.var) |> 
    group_by(as.factor(bw)) |> 
    ggplot(aes(y=surprisal, group=bw, fill=bw)) +
    geom_boxplot() + geom_hline(yintercept=median.ic.in, color="red")
}
