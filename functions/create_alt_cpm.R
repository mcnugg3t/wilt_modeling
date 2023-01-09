create_alt_cpm <- function(var_str, method_str, improve_str, cluster_type, save_name) {
  mod.tmp <- kppm(formula(var_str), 
                  clusters=cluster_type, 
                  data=im.list, 
                  method=method_str, 
                  improve.type=improve_str)
  var.name <- paste0("cpm_hA_",cluster_type,"_",save_name)
  assign(var.name, mod.tmp)
  do.call("save", list(var.name, file=paste0(mod.res.path, "cpm/", "cpm_hA_", cluster_type, "_", save_name, ".Rds")))
  #return(summary(get(var.name)))
}