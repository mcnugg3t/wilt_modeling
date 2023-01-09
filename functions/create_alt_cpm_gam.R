#' create CPM with no improvement method - allows GAM with method_str = "clik2" and gam_bool=T
#'
#'
create_alt_cpm_gam <- function(var_str, method_str, cluster_type, gam_bool, save_name) {
  mod.tmp <- kppm(formula(var_str), 
                  clusters=cluster_type, 
                  data=im.list, 
                  method=method_str,
                  use.gam=gam_bool)
  var.name <- paste0("cpm_hA_",cluster_type,"_",save_name)
  assign(var.name, mod.tmp)
  do.call("save", list(var.name, file=paste0(mod.res.path, "cpm/", "cpm_hA_", cluster_type, "_", save_name, ".Rds")))
  #return(summary(get(var.name)))
}