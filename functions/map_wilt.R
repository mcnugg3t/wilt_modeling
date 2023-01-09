#' map wilt
#' 
#' takes (1) simple model predictions tibble and (1) sf points object, returns model predictions tibble with added column ( $wilt ) that rasterizes points in test_pts to the model prediction grid. 
#' 
#' Assumes test_pts have already been filtered to the study area. Output can then be fed into associated "prec_rec_fun" function. Details are efficient but not very readable.
map_wilt <- function(mod_pred_tbl, test_pts, DEBUG=T) {
  # convert test_pts to tibble with x, y coords
  if(DEBUG) {cat(paste0("\nconverting test points to tibble..."))}
  test.pts.tbl <- do.call(rbind, st_geometry(test_pts)) %>% as_tibble() %>% setNames(c("x", "y"))
  # identify grid cells with test infection
  if(DEBUG){cat(paste0("\n\nidentifying grid cells with observed infection..."))} # debug print
  res.tbl <- mod_pred_tbl %>% add_column(wilt = 0) # add column of 0's to results tibble
  for(i in seq_len(nrow(test.pts.tbl))) { # for each test point in test.pts.tbl
    row.tmp <- test.pts.tbl[i,] # subset relevant row from tibble
    dist.v <- (row.tmp$x- mod_pred_tbl$x)^2 + (row.tmp$y - mod_pred_tbl$y)^2 # calc euclidian distance from infection point to each grid cell center point
    min.ind <- which(dist.v == min(dist.v)) # get index of grid center with minimum distance
    assert_that(length(min.ind) == 1) # assert only one minimum distance
    wilt.ct <- res.tbl$wilt[min.ind] # get the current count value stored in res.tbl$wilt
    if (wilt.ct >= 1) { # if current count is >= 1
      res.tbl$wilt[min.ind] <- wilt.ct + 1 # bump it up one
      if(DEBUG) {cat(paste0("\n\trepeat wilt assignment found"))} # debug print
    } else {
      res.tbl$wilt[min.ind] <- 1 # otherwise set wilt value to 1
    }
  }
  return(res.tbl) # return results tibble that now has added column rasterizing test_pts
}

#########recycling bin##########

# # calc cell size of mod_pred_tbl - just for reference
# if(DEBUG) {cat(paste0("\n\ncalculating cell size..."))} 
# cell.size <- abs(mod_pred_tbl$y[2]-mod_pred_tbl$y[1]); if(DEBUG) {cat(paste0("\n\tcell size : ", cell.size))}

# if(DEBUG){cat(paste0("\n"))} # debug print template



# a <- tibble(x=c(1, 2, 2, 3, 3, 3, 4, 4, 4), y=c(1, 1, 2, 1, 2, 3, 1, 2, 3))
# b <- tibble(x=c(3.1), y=c(2.1))
# dist.v <- (b$x - a$x)^2 + (b$y - a$y)^2
# min.ind <- which(dist.v == min(dist.v))
# min.dat <- a[min.ind, ]
# 
# ggplot() + 
#   geom_point(data=a, aes(x=x,y=y)) + 
#   geom_point(data=b, aes(x=x,y=y), size=2, color="red") + 
#   geom_point(data=min.dat, aes(x=x,y=y), size=2, color="blue")

# Notes

# for each test infection, we want to know which grid cell it's within (if any) - can we just get the nearest grid point?

# what's the most useful intermediate form of data?
# model prediction tibble + whether or not each cell had a case

# prediction is a tibble

# test dat is a sf points

# many possible ways of squaring the circle:
# 1 - 

# we care about the recall and precision of the model given a certain threshold
# model predictions might be at different resolutions
# test data is in spatial points format

# for a given threshold
# recall = proportion of cases predicted by model
# precision = proportion of predicted cases that are actual cases
# specificity (true negative rate) = proportion of non-cases that are accurately assigned as non-cases

# F1 = harmonic mean of precision and recall = (2 (precision * recall)) / (precision + recall)

# i think what we're interested in, for a given model, 
# is the % of cases you could find relative to the area surveyed
# e.g. 80% of cases by surveying 20% of study area vs 100% of cases by surveying 100%