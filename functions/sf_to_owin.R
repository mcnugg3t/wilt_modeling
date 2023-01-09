#***
#* function that takes sf polygon as input, returns spatstat owin object
#* 
sf_to_owin <- function(sf_obj, DEBUG=F) {
  cat(paste0("\nCRS of sf_obj : ", st_crs(sf_obj) ) )
  bnd.tbl <- st_geometry(sf_obj)[[1]] %>% 
    as.matrix() %>% 
    as_tibble() %>% 
    setNames(c("x","y")) %>% 
    dplyr::mutate(ind = row_number() )
  # check if clockwise or counterclockwise currently
  # if (DEBUG) {
  #   p1 <- bnd.tbl %>% mutate(ind=as.factor(ind)) %>% ggplot() +
  #     geom_point(aes(x=x,y=y, color=ind)) +
  #     geom_text(aes(x=x,y=y,label=ind)) # clockwise
  #   print(p1)
  # }
  bnd.tbl.flip <- bnd.tbl %>% dplyr::arrange(desc(ind))
  bnd.tbl.flip$ind <- 1:nrow(bnd.tbl)
  # check if now counter-clockwise
  if (DEBUG) {
    p1 <- bnd.tbl.flip %>% 
      mutate(ind=as.factor(ind)) %>% 
      ggplot() +
      geom_point(aes(x=x,y=y, color=ind)) +
      geom_text(aes(x=x,y=y,label=ind)) # counter-clockwise
    print(p1) 
  }
  bnd.tbl.subs <- bnd.tbl.flip %>% 
    dplyr::slice(-nrow(bnd.tbl)) %>% 
    dplyr::select(x,y) 
  # bnd.tbl.subs
  # 
  # define owin
  require(spatstat)
  wilt.owin <- owin(poly=bnd.tbl.subs)
  return(wilt.owin)
}
