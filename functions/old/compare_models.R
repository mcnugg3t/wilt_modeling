#' Takes list of model objects and a vector indicating which ones are hp3 models (whose prediction tibbles need help)
#'
#' 1. form prediction tibble, adding and unscaling spatial coords if necessary
#' 2. extract test.pts to prediction tibble
#'
compare_models <- function(model_list, test_pts_tbl, hp3_bool_v, hp3_scale_v, DEBUG_SUB, DEBUG) {
  if(DEBUG_SUB) {d.sub <- DEBUG_SUB} else {d.sub <- F}
  # init tibble to store results
  tbl.tmp <- tibble() 
  # for each model in model_list
  for(i in seq_len(length(model_list))) { 
    # extract model
    mod <- model_list[[i]] 
    # extract model name
    mod.name <- names(model_list)[i] 
    # debug print
    if(DEBUG) {cat(paste0("\nmod : ", mod.name, "\n\tcomputing model predictions..."))}

    # if model is an hp3 model, then need to unscale x & y
    if(hp3_bool_v[i]) { # hp3 models may require addition of spatial coord columns, then un-scaling of spatial coords
      if(DEBUG)cat(paste0("\n\tHP3 branch"))
      pred.tmp <- predict(mod) %>% 
        as_tibble() %>%  
        rename(pred=value)
      
      if(DEBUG){cat(paste0("\n\tunscaling hp3..."))}
      pred.unsc <- pred.tmp %>% 
        add_column(x = oak.dat.scale$x, y = oak.dat.scale$y) %>%  # add x and y from model data
        mutate(x = (x*attr(x, "scaled:scale")) + attr(x, "scaled:center"), # unscale x
               y = (y*attr(y, "scaled:scale")) + attr(y, "scaled:center")) # unscale y
      join.tmp <- left_join(pred.unsc, test_pts_tbl, by=c("x", "y"))
    } else {
      if(DEBUG)cat(paste0("\n\tCPM branch"))
      pred.tmp <- predict(mod, locations=sa.quad$dummy) %>% 
        as_tibble() %>% 
        rename(pred = value) %>% 
        add_column(x = sa.quad$dummy$x) %>% 
        add_column(y = sa.quad$dummy$y) %>% 
        filter(!is.na(pred))
      join.tmp <- map_wilt(pred.tmp, test.pts)
    }
    
    if(DEBUG) cat(paste0("\n\tmapping..."))
    # call map function - identifies which cells have wilt
    
    
    
    if(DEBUG) cat(paste0("\n\tmapped ", sum(join.tmp$wilt, na.rm=T), " wilt points"))
    if(DEBUG_SUB) cat(paste0("\njoin.tmp is now ", nrow(join.tmp), " rows"))
    if(DEBUG) {cat(paste0("\n\tcalculating precision and recall..."))}
    
    # if it's an hp3 model and we want results scaled, call prec_rec with scale bool = T
    if(hp3_scale_v[i]) { 
      prec.rec.tmp <- prec_rec_fun(join.tmp, scale_bool=T, DEBUG=d.sub)
    } else { # otherwise doesn't need to be scaled
      prec.rec.tmp <- prec_rec_fun(join.tmp, scale_bool=F, DEBUG=d.sub)
    }
    # add model name to precision / recall tibble
    prec.rec.tmp <- prec.rec.tmp %>% 
      add_column(mod_name = mod.name) %>% 
      relocate(mod_name, .before=thresh)
    # add this model's precision and recall tibble to the larger one
    tbl.tmp <- rbind(prec.rec.tmp, tbl.tmp)
  }
  
  # pivot longer and scale by model threshold - specific value is arbitrary, want [0,1]
  mod.comp.tbl <- tbl.tmp %>% 
    pivot_longer(cols=recall:rel_eff, names_to="measure", values_to="val") %>% 
    group_by(mod_name) %>% 
    mutate(scale_thresh = (thresh-min(thresh))/ (max(thresh)-min(thresh)) ) %>% # scale threshold value (min-max) for plotting on same axis 
    ungroup() 
  # ggplot(aes(x=thresh, y=val, color=mod_name)) + 
  #   geom_line(size=1.1) +
  #   facet_wrap(~measure) +
  #   labs(title="Model Classification Metrics") +
  #   theme(plot.title = element_text(size=36), strip.text.x = element_text(size=24))
  return(mod.comp.tbl)
}