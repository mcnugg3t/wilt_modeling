require(corrplot)
require(RColorBrewer)
pear_corr_ <- function(dat.in) {
  M <- cor(dat.in)
  p1 <- corrplot(M, type="upper", order="hclust",
                 col=brewer.pal(n=8, name="RdYlBu"))
}
