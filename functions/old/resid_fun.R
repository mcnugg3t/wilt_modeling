#' model object either mgcv::gam (mb = 0) or spatstat::kppm (mb = 1)
#' data tibble with variables of interest - x, y ...
#' 
#' returns multiplot of residuals mapped, residuals vs fitted, residuals vs important predictors
resid_fun <- function(mod, is_kppm, dat_tbl, pts_shp, p_slice, dbg) {
  require(tidyverse)
  require(gridExtra)
  #
  # spatstat:kppm branch
  #
  if(dbg) cat(paste0("\nModel is kppm: ", (is_kppm == 1), "\n\tcalculating residuals..."))
  if(is_kppm) {
    if(dbg) cat(paste0("\n\t\tkppm branch..."))
    # calc fitted
    fit.tmp <- predict(mod, locations=sa.quad$dummy) # calc fitted values
    fit.tbl <- tibble(x = sa.quad$dummy$x, y= sa.quad$dummy$y, fitted=fit.tmp) # make tibble with fitted values
    # calc residuals
    rd.tmp <- residuals(mod, type="raw", quad=sa.quad) # raw
    rp.tmp <- residuals(mod, type="pearson", quad=sa.quad) # pearson
    resid.tbl <- tibble(x=as.vector(rd.tmp$loc$x), y=rd.tmp$loc$y, rd=rd.tmp$val, rp=rp.tmp$val) # tibble
    
    # join fitted with residuals
    dat.mid <- left_join(fit.tbl, resid.tbl, by=c("x", "y"))
    dat.tmp <- left_join(dat.mid, dat_tbl, by=c("x", "y"))
  #
  # mgcv::gam branch
  #
  } else {
    if(dbg) cat(paste0("\n\t\tgam branch..."))
    dat.tmp <- dat_tbl %>% 
      add_column(rd = residuals(mod, type="deviance")) %>% 
      add_column(rp = residuals(mod, type="pearson")) %>% 
      add_column(fitted = mod$fitted.values) %>% 
      mutate(x = (x*attr(x, "scaled:scale"))+attr(x, "scaled:center"),
             y = (y*attr(y, "scaled:scale")) + attr(y, "scaled:center"))
  }
  # calc 5%
  rd_p05 <- quantile(dat.tmp$rd, probs=c(0.05), na.rm=T)
  rp_p05 <- quantile(dat.tmp$rp, probs=c(0.05), na.rm=T)
  # calc 95%
  rd_p95 <- quantile(dat.tmp$rd, probs=c(0.95), na.rm=T)
  rp_p95 <- quantile(dat.tmp$rp, probs=c(0.95), na.rm=T)
  if(dbg) cat(paste0("\n\tDev Resid 5% = ", rd_p05, "  &  95% = ", rd_p95, "\n\tPears Resid 5% = ", rp_p05, "  &  95% = ", rp_p95, "\n"))
  
  if(rd_p05 <= 0 & rd_p95 <= 0 & (abs(rd_p05) > abs(rd_p95)) ) {
    if(dbg) {cat(paste0("\nDev residuals inverted, adjusting & recalculating..."))}
    dat.tmp$rd <- abs(dat.tmp$rd)
    # recalc quantiles
    rd_p05 <- quantile(dat.tmp$rd, probs=c(0.05), na.rm=T)
    rd_p95 <- quantile(dat.tmp$rd, probs=c(0.95), na.rm=T)
    if(dbg) {cat(paste0("\n\tDev resid 5% = ", rd_p05, "  &  95% = ", rd_p95))}
    }
  if(rp_p05 <= 0 & rp_p95 <= 0 & (abs(rp_p05) > abs(rp_p95)) ) {
    if(dbg) {cat(paste0("\nPears residuals inverted, adjusting & recalculating..."))}
    dat.tmp$rp <- abs(dat.tmp$rp)
    # recalc quantiles
    rp_p05 <- quantile(dat.tmp$rp, probs=c(0.05), na.rm=T)
    rp_p95 <- quantile(dat.tmp$rp, probs=c(0.95), na.rm=T)
    if(dbg) {cat(paste0("\n\tPears resid 5% = ", rp_p05, "  &  95% = ", rp_p95))}
    }
  
  #
  # plots
  #
  if(dbg) cat(paste0("\n\tplotting...\n\t\tspatial..."))
  # spatial
  p1.dat <- dat.tmp %>% filter(rd >= rd_p95)
  p2.dat <- dat.tmp %>% filter(rp >= rp_p95)
  pears.top <- dat.tmp %>% arrange(desc(rp)) %>% slice(1:p_slice)
  # p0 - fitted values with infection sites
  p0 <- ggplot() +
    geom_tile(data=dat.tmp, aes(x=x,y=y,color=fitted)) +
    geom_sf(data=pts_shp, color="red", fill=NA, size=1.5, alpha=0.1) +
    labs(title="fitted values with infection sites") +
    scale_color_viridis_c(option="magma", direction=1)
  # p1 - top 5% dev residuals
  p1 <- ggplot() + 
    geom_point(data=p1.dat, aes(x=x,y=y, color=rd)) + 
    #geom_sf(data=pts_shp, size=0.3, color="red", alpha=0.5) + 
    labs(title="top 5% dev resid") +
    scale_color_viridis_c(option="magma", direction=1)
  # p2 - top 5% pearson residuals
  p2 <- ggplot() + 
    geom_point(data=p2.dat, aes(x=x,y=y, color=rp)) + 
    geom_point(data=pears.top, aes(x=x,y=y), color="red") +
    #geom_sf(data=pts_shp, size=0.3, color="red", alpha=0.5) + 
    labs(title=paste0("top 5% pears resid + top ", p_slice, " (red)") ) +
    scale_color_viridis_c(option="magma", direction=1)

  if(dbg){cat(paste0("\n\t\tresid vs fitted..."))}
  # fitted
  # p3 - deviance residuals vs fitted
  p3 <- dat.tmp %>% ggplot(aes(x=fitted, y=rd)) + 
    geom_point() +
    labs(title="dev resid vs fitted")
  # p4 - pearson residuals vs fitted
  p4 <- dat.tmp %>% ggplot(aes(x=fitted, y=rp)) +
    geom_point() +
    labs(title="pears resid vs fitted")
  
  if(dbg){cat(paste0("\n\t\tresid vs oak_prob, elev..."))}
  # oak_prob and elev
  p5 <- dat.tmp %>% mutate(oak_prob_rank = percent_rank(oak_prob)) %>% 
    ggplot(aes(x=oak_prob_rank, y=rp, color=rp)) + geom_point() +
    labs(title="pears resid vs oak_prob (rank)")
  
  p6 <- dat.tmp %>% mutate(elev_rank = percent_rank(elev)) %>% 
    ggplot(aes(x=elev_rank, y=rp, color=rp)) + geom_point() +
    labs(title="pears resid vs elev (rank)")
  
  if(dbg){cat(paste0("\n\t\tresid vs other predictors..."))}
  # vs predictors
  p7 <- dat.tmp %>% 
    mutate(plan_curv_rank = percent_rank(plan_curv)) %>% 
    ggplot(aes(x=plan_curv_rank,y=rp, color=rp)) + geom_point() +
    labs(title="plan_curv")
  
  p8 <- dat.tmp %>% 
    mutate(ch_net_dist_rank = percent_rank(ch_net_dist)) %>% 
    ggplot(aes(x=ch_net_dist_rank,y=rp, color=rp)) + geom_point() +
    labs(title="ch_net")
  
  p9 <- dat.tmp %>% 
    mutate(bd_rank = percent_rank(bd)) %>% 
    ggplot(aes(x=bd_rank, y=rp, color=rp)) + geom_point() +
    labs(title="bd")
  
  p10 <- dat.tmp %>% 
    mutate(ph_rank = percent_rank(ph)) %>%
    ggplot(aes(x=ph_rank,y=rp, color=rp)) + geom_point() +
    labs(title="ph")
  
  p11 <- dat.tmp %>% 
    mutate(hs_rank = percent_rank(hs)) %>%
    ggplot(aes(x=hs_rank,y=rp, color=rp)) + geom_point() +
    labs(title="hs")
  
  p12 <- dat.tmp %>% 
    mutate(sand_rank = percent_rank(sand)) %>%
    ggplot(aes(x=sand_rank,y=rp, color=rp)) + geom_point() +
    labs(title="sand")
  
  # grid.arrange(p5,p6,p7,p8,p9,p10, nrow=2)
  
  return( list(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) )
}
