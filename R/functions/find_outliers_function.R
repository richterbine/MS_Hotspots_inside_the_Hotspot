# adapted from https://rawgit.com/valentinitnelav/plotbiomes/master/html/Whittaker_biomes_examples.html
# function get_outliers

find.outliers <- function (points, poly) {
  (points <- data.table::as.data.table(points))
  (data.table::setnames(x = points, c("x", "y")))
  (points[, `:=`(row_idx, 1:.N)])
  (points_valid <- points[complete.cases(points)])
  (points_valid_sp <- sp::SpatialPoints(coords = points_valid[, -3], proj4string = poly@proj4string))
  (sp_overlay <- sp::over(x = points_valid_sp, y = as(poly, 
                                                 "SpatialPolygons")))
  (outliers <- points_valid[is.na(sp_overlay)])
  data.table::setcolorder(x = outliers, neworder = c("row_idx", 
                                                     "x", "y"))
  return(outliers)
}
