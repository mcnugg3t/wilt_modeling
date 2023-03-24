require(mgcv)
require(spatstat)
source("functions/mod/map_wilt_.R")
source("functions/mod/prec_rec_fun_.R")
#' Takes a list of model objects; a vector with 1 for GAM, 2 for CPM; a vector with resolution; a list with model input data
#'
#' 1. form prediction tibble
#' 2. extract test.pts to prediction tibble
#'
compare_models_tbl_ <- function(model_list, model_vect, grid_vect, model_in_data, verbose, DBG) {
  # init tibble to store results
  tbl.tmp <- tibble() 
  
  # for each model in model_list...
  for(i in seq_along(model_list) ) { 
    if(DBG) cat(paste0("\nmodel ", i))
    
    # extract model object, info, & dbg print
    mod <- model_list[[i]]
    mod.type = model_vect[i]
    mod.name <- names(model_list)[i]
    mod.grid = grid_vect[i]
    if(verbose) {cat(paste0("\n\tmod name: ", mod.name, "  on grid: ", mod.grid, "-m", "\n\tcomputing model predictions..."))}

    # GAM branch
    if(mod.type == 1) { 
      if(DBG)cat(paste0("\n\tGAM branch"))
      
      # predict -> tibble
      pred.tmp = mgcv::predict.bam(mod) |>  
        as_tibble() |>  
        rename(raw_pred=value) 
      
      # 30-m branch
      if(mod.grid == 30) {
        pred.tmp = pred.tmp |> # add missing spatial info & wilt test (validation points)
          add_column(x = model_in_data[["GAM data 30"]][["x"]],
                     y = model_in_data[["GAM data 30"]][["y"]],
                     wilt = model_in_data[["GAM data 30"]][["wilt_test"]])
        
      # 10-m branch
      } else if(mod.grid == 10) {
        pred.tmp = pred.tmp |> # add missing spatial info & wilt test (validation points)
          add_column(x = model_in_data[["GAM data 10"]][["x"]],
                     y = model_in_data[["GAM data 10"]][["y"]],
                     wilt = model_in_data[["GAM data 10"]][["wilt_test"]])
      }
        
      
      # convert raw to predict, scale by max
      if(str_detect(mod.name, "logit")) {
        pred.tmp = pred.tmp |> 
          mutate(pred = plogis(raw_pred)) |> 
          mutate(pred = pred/max(pred, na.rm=T))
      } else {
        pred.tmp = pred.tmp |> 
          mutate(pred = exp(raw_pred)) |> 
          mutate(pred = pred/max(pred, na.rm=T))
      }
    
    # CPM branch
    } else {
      if(DBG)cat(paste0("\n\tCPM branch"))
      if(mod.grid == 30) {
        quad.tmp = model_in_data[["quads 30"]]
      } else if(mod.grid == 10) {
        quad.tmp = model_in_data[["quads 10"]]
      }
      
      # predict
      if(DBG)cat(paste0("\n\t\tpredicting..."))
      pred.tmp = predict(mod, locations = quad.tmp[["dummy"]]) |>  
        as_tibble() |>  
        rename(raw_pred = value) |>  
        add_column(x = quad.tmp[["dummy"]][["x"]]) |> 
        add_column(y = quad.tmp[["dummy"]][["y"]]) |> 
        mutate(pred = raw_pred/max(raw_pred, na.rm=T)) |> 
        select(x, y, pred)
      #if(i == 3) browser()
      if(mod.grid == 30) {
        pred.tmp = map_wilt_30_(pred.tmp, DBG=T)
      } else if (mod.grid == 10) {
        pred.tmp = map_wilt_10_(pred.tmp, DBG=T)
      }
      
    }
    prec.rec.tmp = prec_rec_fun_(pred.tmp, DEBUG=T)
    
    # add model name to precision / recall tibble
    prec.rec.tmp = prec.rec.tmp |>  
      add_column(mod_name = mod.name) |>  
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

#### test ##### 
# model_list = mod.list; model_vect = mod.v; model_in_data = model.in.data; verbose=T; DBG=T









