#' Create a single-level voronoi map from summarized data
#' @param values the weights for each area
#' @param groups the grouping variable
#' @param shape the initial shape to start with (ignore for unit circle)
#' @param seed random seed to make images reproducible
#' @param iter maximum iterations to achieve desired accuracy
#' @param accuracy stop after achieving this value of accuracy
#'
#' @return a voronoi tassellation within a given shape
#' @export
#'
#' @examples
#' voronoi_map(c("A", "B", "C"), 1:3)
voronoi_map <- function(values, groups, shape, seed,
                        iter = 200, accuracy = 0.01) {

  if(missing(seed)) seed <- NULL
  if(is.null(seed)) seed <- length(groups)
  if(missing(groups)) groups <- as.character(seq_along(values))
  if(missing(shape)) {
    shape <- terra::buffer(terra::vect(cbind(0, 0),
                                       crs = "+proj=utm +zone=1"), 1, 30)
  }

  set.seed(seed)
  tass <- terra::spatSample(shape, length(groups))
  res <- terra::crop(terra::voronoi(tass, bnd = shape), shape)
  res <- structure(
            list(geom = res,
              sites = tass,
              areas = terra::expanse(res)/sum(terra::expanse(res)),
              groups = groups,
              values = values,
              shape = shape),
            class = "voronoi_map")
  for(i in 1:5) res <- centralize(res)
  for(i in seq(iter)) {
    res <- improve_weights(res)
    if(rmse(res) < accuracy) break;
  }
  res
}

#-------------------------------------------------------------------------------
#' Plotting method for voronoi maps
#'
#' @param x A valid voronoi_map object
#' @param y Ignored
#' @param ... Ignored
#'
#' @return Doesn't return anything; draws plot as side effect
#' @export
#'
#' @examples
#' plot.voronoi_map(obj)
plot.voronoi_map <- function(x, y, ...) {
  terra::plot(x$geom, ...)
  coords <- terra::crds(terra::centroids(x$geom))
  text(x = coords[, 1], y = coords[, 2], label = x$groups)
}

#-------------------------------------------------------------------------------
#' Summary method for Voronoi maps
#'
#' @param x a voronoi map object
#'
#' @return a data frame containing the sites and tile areas
#' @export
#'
#' @examples
#' summary(voronoi_map(1:5))
summary.voronoi_map <- function(x) {
  cbind(data.frame(group = x$groups,
                   value = x$values),
        setNames(as.data.frame(terra::crds(x$sites)), c("x", "y")),
        tile_area = terra::expanse(x$geom) * x$areas,
        relative_area = x$areas)
}

#-------------------------------------------------------------------------------
#' Print method for Voronoi maps
#'
#' @param x a Voronoi map object
#' @param ...
#'
#' @return invisibly returns a summary data frame
#' @export
#'
#' @examples
#' print(voronoi_map(1:5))
print.voronoi_map <- function(x, ...) {
  cat("A voronoi map with", length(x$groups), "sites\n")
  print(summary(x))
}

#-------------------------------------------------------------------------------
# Given a set of points to act as sites, this updates the voronoi map object

update_from_points <- function(tass, pts) {
  tass$sites <- pts
  tass$geom  <- terra::crop(terra::voronoi(pts, bnd = tass$shape), tass$shape)
  sub <- apply(terra::relate(tass$geom, tass$sites, "contains"), 2, which)
  tass$geom <- tass$geom[sub]
  tass$areas <- terra::expanse(tass$geom)/sum(terra::expanse(tass$geom))
  tass
}

#-------------------------------------------------------------------------------
# Lloyd's algorithm to make Voronoi equally spaced initially

centralize <- function(tass, ratio = 1) {
  cent  <- terra::crds(terra::centroids(tass$geom))
  pts <- terra::crds(tass$sites)
  pts <- (1 - ratio) * pts + (ratio * cent)
  cent <- terra::vect(pts, type = 'points', crs = terra::crs(tass$shape))
  update_from_points(tass, cent)
}

#-------------------------------------------------------------------------------
# Generates a matrix of a single index point repeated so that it is the same
# dimensions as an array of all site indices (allows point comparisons)

coords_from_index <- function(tass, index) {
  all_pts <- terra::crds(tass$sites)
  matrix(all_pts[index, ], nrow = nrow(all_pts), ncol = 2, byrow = TRUE)
}

#-------------------------------------------------------------------------------
# Find the vector between an index point and all other points

vector_to_points <- function(tass, index) {
  all_pts <- terra::crds(tass$sites)
  all_pts - coords_from_index(tass, index)
}

closest_edge <- function(tass, index) {
  terra::crds(terra::nearest(tass$sites[index], terra::as.lines(tass$shape)))
}


vector_to_centroid <- function(tass) {
  cent <- terra::crds(terra::centroids(tass$shape))
  all_pts <- terra::crds(tass$sites)
  cent <- matrix(cent, nrow = nrow(all_pts), ncol = 2, byrow = TRUE)
  cent - all_pts
}

centroid_rays <- function(tass) {
  vc <- -100 * vector_to_centroid(tass)
  cent <- terra::crds(terra::centroids(tass$shape))
  cent <- matrix(cent, nrow = nrow(vc), ncol = 2, byrow = TRUE)
  lines <- do.call(rbind, lapply(seq(nrow(cent)), function(x) {
    rbind(c(object = x, x = cent[x, 1], y = cent[x, 2]),
          c(object = x, x = vc[x, 1] + cent[x, 1], y = vc[x, 2] + cent[x, 2]))
  }))
  terra::vect(lines, type = 'lines', crs = terra::crs(tass$shape))
}

#-------------------------------------------------------------------------------
# Find the vectors connecting each point and the boundary, pointing from an
# index point through all other points.

push_vec <- function(tass, index) {
  orig <- coords_from_index(tass, index)
  vec <- 1e6 * vector_to_points(tass, index)
  lines <- do.call(rbind, lapply(seq(nrow(orig)), function(x) {
    rbind(c(object = x, x = orig[x, 1], y = orig[x, 2]),
          c(object = x, x = vec[x, 1] + orig[x, 1], y = vec[x, 2] + orig[x, 2]))
    }))
  lines <- terra::vect(lines, type = 'lines', crs = terra::crs(tass$shape))

  proj <- terra::crds(terra::intersect(lines, terra::as.lines(tass$shape)))
  vec[-index,] <- proj - terra::crds(tass$sites)[-index,]
  vec
}

#-------------------------------------------------------------------------------
# Find the most disproportionate area and push or pull the other points from it

improve_weights <- function(tass) {
  orig <- terra::crds(tass$sites)
  weights <- tass$values/sum(tass$values)
  areas <- tass$areas
  diff <- areas - weights
  ind <- which.max(abs(diff))
  pushvec <- push_vec(tass, ind) * -diff[ind]

  pts <- orig + pushvec
  pts <- terra::vect(pts, 'points', crs = terra::crs(tass$shape))
  update_from_points(tass, pts)
}

#-------------------------------------------------------------------------------
# Finds the root mean square error of relative weights and relative areas

rmse <- function(tass) {
  w <- tass$values / sum(tass$values)
  a <- tass$area
  sqrt(mean((a - w)^2))
}

#' Get polygons from voronoi map
#'
#' @param tass Voronoi Map
#'
#' @return data frame of polygons
#' @export
#'
#' @examples
#' get_polygons(voronoi_map(1:3))
get_polygons <- function(tass) {
  d <- as.data.frame(terra::geom(tass$geom))[-c(2, 5)]
  d$group <- tass$groups[d$geom]
  d$value <- tass$values[d$geom]
  d
}


make_vmt <- function(formula, data, summary_fun = sum, shape, seed) {
  mf <- stats::model.frame(formula = formula, data = data)
  terms <- attr(mf, "terms")
  pred_vars <- attr(terms, "term.labels")
  value_var <- attr(mf, "names")[1]
  mf[pred_vars] <- lapply(mf[pred_vars], function(x) factor(x, unique(x)))
  first_tab <- tapply(mf[[value_var]], mf[[pred_vars[1]]], summary_fun)
  vm <- voronoi_map(as.vector(first_tab), names(first_tab), shape = shape,
                    seed = seed)
  if(length(pred_vars) == 1) return(vm)

  vmt <- list(root = vm,
              children = Map(function(x, y) {
                f <- reformulate(pred_vars[-1], response = value_var)
                make_vmt(f, x, summary_fun = summary_fun, shape = vm$geom[y])
              },
              x = split(mf[-(which(names(mf) == pred_vars[1]))], mf[[pred_vars[1]]]),
              y = seq_along(vm$geom)))

  return(vmt)
}


parse_vmt <- function(vmt, parent = NA, level = 0) {

  if(inherits(vmt, "voronoi_map")) {
    return(cbind(get_polygons(vmt), parent, level))
  }

  nm <- names(vmt)

  do.call(rbind, Map(function(a, b) parse_vmt(a, b, level + 1), vmt, nm))

}


#' Get the polygons for a Voronoi treemap
#'
#' @param formula a formula with the weights on the left hand side and
#' grouping variables in order on the right
#' @param data the data frame where the variables are to be found
#' @param summary_fun in case of multiple values within a single group,
#' a function that will summarize them (defaults to sum)
#' @param shape a SpatVector to act as a boundary. If missing will default
#' to a unit circle
#'
#' @return a data frame of polygons
#' @export
#'
#' @examples
#' voronoi_treemap(Sepal.Length ~ Species, iris)
voronoi_treemap <-  function(formula, data, summary_fun = sum, shape,
                             seed = NULL) {
  if(missing(shape)) {
    shape <- terra::buffer(terra::vect(cbind(0, 0),
                                       crs = "+proj=utm +zone=1"), 1, 30)
  }
  vmt <- make_vmt(formula, data, summary_fun, shape, seed = seed)
  res <- parse_vmt(vmt)
  rownames(res) <- NULL
  res$level <- match(res$level, unique(res$level))
  res[-1]
}
