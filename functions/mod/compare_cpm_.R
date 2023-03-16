#' takes 2 CPMs, performs chisq test using delta DOF
#'
#'
compare_cpm <- function(cpm_ha, cpm_h0, DEBUG=F) {
  delta.dof <- as.integer(attr(logLik(cpm_ha), "df")) - as.integer(attr(logLik(cpm_h0), "df"))
  if(DEBUG) {cat(paste0("\n\n\ndelta.dof = ", as.integer(delta.dof)) )}
  t.stat <- 2*(as.numeric(logLik(cpm_ha)) - as.numeric(logLik(cpm_h0)))
  if(DEBUG) {cat(paste0("\nt-stat = ", t.stat))}
  p.val <- pchisq(t.stat, df=delta.dof, lower.tail=F)
  if(DEBUG) cat(paste0("\np = ", p.val))
  return(p.val)
}