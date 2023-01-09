
#' Takes output of map_wilt function and returns tbl that contains info for desired precision/recall plots
#'
#' 1 - compute possible threshold values
#' 2 - loop over threshold values and compute precision + recall stats
#'
prec_rec_fun <- function(mod_pred_tbl, scale_bool, scale_px=101703, scale_wilt=108, DEBUG) {
  #
  # 1 - compute possible threshold values and other stats
  #
  min.pred <- min(mod_pred_tbl$pred); max.pred <- max(mod_pred_tbl$pred) # get min and max pred values
  thresh.v <- seq(from=min.pred, to=max.pred, length.out=300) # construct vector of 300 threshold vals btw min and max
  
  if(DEBUG){cat(paste0("\npopulating threshold values..."))} # debug print
  if(DEBUG){cat(paste0("\n\t\tmin pred : ", min.pred, "\tmax pred : ", max.pred))} # debug print
  
  wilt.ct <- mod_pred_tbl %>% # calc total # wilt in model prediction tibble
    filter(wilt>=1) %>%  
    summarise(wilt.ct = sum(wilt)) %>% 
    as.numeric()
  
  wilt.n.px <- mod_pred_tbl %>% 
    filter(wilt >= 1) %>% 
    summarise(wilt.n.px = n()) %>% 
    as.numeric()
  
  if(scale_bool) { # for calculating % of study area surveyed
    n.total <- scale_px # if scale=T set n.total to scale_px
    wilt.ct <- scale_wilt
  } else {
    n.total <- nrow(mod_pred_tbl) # otherwise n.total <- # pixels in model prediction tibble
  }
  
  if(DEBUG) {cat(paste0("\n\ttotal # wilt cases = ", wilt.ct, " \t in ", wilt.n.px, " pixels."))}
  if(DEBUG) {cat(paste0("\n\tscale_bool is : ", scale_bool, " n.total = ", n.total))}
  #
  # 2 - compute precision/recall stats
  #
  # initialize storage vectors
  recall.v <- c() # empty vector to store calculated recall values
  prec.v <- c() # empty vector to store precision
  area.v <- c() # empty vector to store proportion total area surveyed values (positive recall)
  bacc.v <- c() # empty vector to store balanced accuracy values
  # debug print
  if(DEBUG){cat(paste0("\n\ncomputing precision & recall stats...\n\n"))} 
  # loop over threshold values
  for (t in thresh.v) { 
    # if(DEBUG) cat(paste0("\nthresh = ", t))
    # filter data to locations with model predictions >= threshold
    dat.filt <- mod_pred_tbl %>% filter(pred >= t)
    dat.excl <- mod_pred_tbl %>% filter(pred < t)
    
    n.pred <- dat.filt %>% nrow() # number pixels predicted at this threshold
    n.ppos <- dat.filt %>% filter(wilt >= 1) %>% nrow() # count of pred pixels where wilt is present
    
    n.excl <- dat.excl %>% nrow()
    n.npos <- dat.excl %>% filter(wilt >= 1) %>% nrow() # count of excl pixels where wilt is present
    
    # TP = # wilt correctly identified
    tp <- dat.filt %>% 
      filter(wilt >= 1) %>% 
      summarise(tp = sum(wilt)) %>% 
      as.numeric()
    # FP = count pixels predicted wilt AND is not wilt
    fp <- n.pred - n.ppos
    # TN = count of pixels predicted not wilt AND is not wilt
    tn <- n.excl
    # if(DEBUG & !scale_bool) {assert_that(n.pred == (fp+tn))}
    # FN = # wilt predicted not wilt and is actually wilt
    fn <- dat.excl %>% 
      filter(wilt >= 1) %>% 
      summarise(fn = sum(wilt)) %>% 
      as.numeric()
    # debug print
    if (DEBUG) {cat(paste0("\nt: ", round(t, 4), 
                           "  npred: ", n.pred, 
                           "  out of: ", n.total, 
                           "\tn.excl: ", n.excl, 
                           "  n.ppos: ", n.ppos, 
                           "  n.npos: ", n.npos, 
                           "  tp: ", tp, 
                           "  fp: ", fp, 
                           "  tn: ", tn, 
                           "  fn: ", fn))}
    # recall = TP/(TP+FN)  = proportion of actual cases present in data that the model correctly labels
    rec.tmp <- tp/wilt.ct
    recall.v <- c(recall.v, rec.tmp)
    # precision = proportion of predicted cases that are actual wilt cases = TP/(TP+FP)
    prec.tmp <- tp/n.pred
    prec.v <- c(prec.v, prec.tmp)
    # area surveyed = (TP+FP)/N
    area.tmp <- (n.pred)/n.total
    area.v <- c(area.v, area.tmp)
    # balanced accuracy = 1/2 ( TP/ (TP+FN) + TN/ (TN+FP) )
    bacc.tmp <- 0.5 * ( rec.tmp + prec.tmp )
    bacc.v <- c(bacc.v, bacc.tmp)
  }
  # create tibble to return
  if(DEBUG){cat(paste0("\n\ncreating ROC tibble...\n\n"))}
  roc.tbl <- tibble(thresh = thresh.v, 
                    recall = recall.v, 
                    precision = prec.v, 
                    area_surveyed = area.v, 
                    balanced_accuracy = bacc.v) %>% 
            mutate(rel_eff = recall/area_surveyed)
  # return
  return(roc.tbl)
}
