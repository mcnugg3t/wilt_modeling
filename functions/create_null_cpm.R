create_null_cpm <- function(method_str, improve_str, cluster_type, save_name) {
  mod.tmp <- kppm(sa.quad ~ surv_bool + pine_bool + oak_prob, 
                  clusters=cluster_type, 
                  data=im.list, 
                  method=method_str, 
                  improve.type=improve_str)
  var.name <- paste0("cpm_h0A_",cluster_type,"_",save_name)
  assign(var.name, mod.tmp)
  do.call("save", list(var.name, file=paste0(mod.res.path, "cpm/", "cpm_h0A_", cluster_type, "_", save_name, ".Rds")))
}