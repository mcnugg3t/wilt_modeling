create_train_cpm_gam <- function(var_str, method_str, gam_bool, cluster_type, hyp_str, save_str) {
  mod.tmp <- kppm(formula(var_str), 
                  clusters=cluster_type, 
                  data=im.list, 
                  method=method_str,
                  use.gam=gam_bool)
  var.name <- paste0("cpm_train_", hyp_str, "_", cluster_type, "_", save_str)
  assign(var.name, mod.tmp)
  do.call("save", list(var.name, file=paste0(mod.res.path, "cpm/", var.name, ".Rds")))
}