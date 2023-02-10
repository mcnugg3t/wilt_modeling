#'
#'
create_dirs_ <- function(){
  path.0 ="mid_data"
  dir.create(path.0)
  
  mid.dat.0 <- c("10", "30", "wilt")
  mid.dat.1 <- c("grid", "wiscland2", "bedrock", "wilt", "polaris", "manage")
    polaris.dat <- c("alpha", "bd", "clay", 
                     "hb", "ksat", "lambda", 
                     "n", "om", "ph", "sand", 
                     "silt", "theta_r", "theta_s")
  
  
  for(v0 in mid.dat.top) {
    path.1 <- paste0(path.0, "/", v0)
    dir.create(path.1)
    for(v1 in mid.dat.1) {
      if(v1 == "wilt") next
      path.2 <- paste0(path.1, "/", v1)
      dir.create(path.2)
      if(v1 == "polaris") {for(p1 in polaris.dat) path.3 <- paste0(path.2, "/", p1)}
    }
  }
}