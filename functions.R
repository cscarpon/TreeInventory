process_inventory <- function(las_str, buffer, infrastructure, epsg_crs, output_dir) {
  
  file_name <- stringr::str_split(basename(las_str), "\\.")[[1]][1]
  
  # Read the LAS file and assign point IDs
  laz <- lidR::readLAS(las_str)
  laz@data$P_ID <- 1:nrow(laz@data)
  
  # Classify and remove noise points
  laz <- lidR::classify_noise(laz, sor(15, 7))
  noise_points <- laz@data$P_ID[laz@data$Classification == 18]  # Assuming '18' is the noise class
  laz <- laz[!laz@data$P_ID %in% noise_points, ]
  
  # Create bounding box for LAS file and intersect it with buffer and infrastructure layers
  bbox <- sf::st_bbox(laz)
  bbox_sf <- sf::st_as_sfc(bbox, crs = epsg_crs)
  
  # Intersect buffer with LAS bounding box
  buffer_intersect <- sf::st_intersection(buffer, bbox_sf)
  
  
  # Clip LAS to ROI if intersects exist
  if (length(buffer_intersect) > 0) {
    laz_roi <- clip_roi(laz, buffer_intersect)
  } else {
    message("Infrastructure intersect not found for file: ", las_str)
    return(NULL)
  }
  
  # Remove ground points
  laz_ground <- filter_ground(laz_roi)
  if (!is.null(laz_ground)) {
    laz_roi <- laz_roi[!laz_roi@data$P_ID %in% laz_ground@data$P_ID, ]
  } else {
    message("Ground points not found for file: ", las_str)
    return(NULL)
  }
  ground_las <- classify_ground(laz_roi, algorithm  = csf())
  ground_las <- filter_ground(ground_las)
  
  laz_roi <- laz_roi[!laz_roi@data$P_ID %in% ground_las@data$P_ID, ]
  
  # Infrastructure filtering
  
  # Intersect infrastructure with LAS bounding box
  infrastructure_intersect <- sf::st_intersection(infrastructure, bbox_sf)
  infrastructure_count <- if (!is.null(infrastructure_intersect)) nrow(as.data.frame(infrastructure_intersect)) else 0
  if (infrastructure_count >= 1) {
    laz_infrastructure <- clip_roi(laz_roi, infrastructure_intersect)
    if (!is.null(laz_infrastructure)) {
      laz_roi <- filter_poi(laz_roi, !(laz_roi@data$P_ID %in% laz_infrastructure@data$P_ID))
    }
  }
  
  laz_roi <- segment_shapes(laz_roi, shp_line(th1 = 10, k = 8), attribute = "Shape")
  laz_roi <- filter_poi(laz_roi, Shape == FALSE)
  laz_roi@data$Shape <- NULL
  
  laz_roi <- filter_poi(laz_roi, Classification != 8)
  
  # Save the processed file
  if (!file.exists(paste0(output_dir, "/laz/laz_", file_name, ".laz"))) {
    writeLAS(laz_roi, paste0(output_dir, "/laz/laz_", file_name, ".laz"))
  } else {
    laz_roi <- readLAS(paste0(output_dir, "/laz/laz_", file_name, ".laz"))
  }
  
  # Create the DTM
  print(paste0("Creating the DTM for tile ", file_name))
  if (file.exists(paste0(output_dir, "/dtm/dtm_", file_name, ".tif"))) {
    print("The DTM already exists, loading that file")
    dtm <- terra::rast(paste0(output_dir, "/dtm/dtm_", file_name, ".tif"))
  } else {
    writeLAS(laz_ground, paste0(output_dir, "/laz/ground_", file_name, ".laz"))
    dtm <- lidR::rasterize_terrain(laz_ground, res = 0.25, algorithm = tin())
    dtm <- terra::mask(dtm, vect(buffer_intersect))
    terra::writeRaster(dtm, paste0(output_dir, "/dtm/dtm_", file_name, ".tif"), gdal = c("COMPRESS=LZW"), overwrite = TRUE)
  }
  
  # Create normalized LAS
  print("Creating the normalized LAS")
  if (file.exists(paste0(output_dir, "/nlas/nlas_", file_name, ".laz"))) {
    print("The normalized LAS already exists, loading that file")
    nlas <- lidR::readLAS(paste0(output_dir, "/nlas/nlas_", file_name, ".laz"))
  } else {
    dtm <- dtm + 4
    nlas <- laz_roi - dtm
    nlas <- lidR::filter_poi(nlas, Z < 30)
    lidR::writeLAS(nlas, paste0(output_dir, "/nlas/nlas_", file_name, ".laz"))
  }
  
  # Create the CHM
  print("Creating the CHM")
  if (file.exists(paste0(output_dir, "/chm/chm_", file_name, ".tif"))) {
    print("The CHM already exists, loading that file")
    chm_clamp <- terra::rast(paste0(output_dir, "/chm/chm_", file_name, ".tif"))
  } else {
    chm <- lidR::rasterize_canopy(nlas, res = 0.5, p2r(0.2, na.fill = tin()))
    chm <- terra::mask(chm, vect(buffer_intersect))
    chm_smooth <- terra::focal(chm, w = matrix(1, 5, 5), fun = median, na.rm = TRUE)
    chm_clamp <- terra::clamp(chm_smooth, 0, 30)
    terra::writeRaster(chm_clamp, paste0(output_dir, "/chm/chm_", file_name, ".tif"), overwrite = TRUE)
  }
  
  # Locate tree tops
  print("Creating the Big Tops")
  if (file.exists(paste0(output_dir, "/tops/topsBig_", file_name, ".shp"))) {
    print("The tops already exists, loading that file")
    ttops <- sf::st_read(paste0(output_dir, "/tops/topsBig_", file_name, ".shp"))
  } else {
    ttops <- lidR::locate_trees(nlas, lmf(ws = 7.5, hmin = 2.0))
    ttops <- sf::st_zm(ttops)  # Remove ZM values if necessary
    sf::st_write(ttops, paste0(output_dir, "/tops/topsBig_", file_name, ".shp"))
  }
}


process_tree_group <- function(tree_group) {
  tryCatch({
    # Find the highest point (tree top)
    tree_top <- tree_group %>%
      dplyr::slice_max(Z) %>%
      select(X, Y, Z)
    
    # Convert the highest point to an sf object
    tree_top_sf <- st_as_sf(tree_top, coords = c("X", "Y"), crs = st_crs(26917))
    
    # Get distinct points for the crown
    distinct_points <- tree_group %>%
      select(X, Y) %>%
      distinct() %>%
      st_as_sf(coords = c("X", "Y"), crs = 26917)
    
    # Count the number of points
    point_count <- nrow(distinct_points)
    
    # Combine points into a MULTIPOINT geometry
    distinct_multipoint <- st_union(distinct_points)
    
    # Create a concave hull
    tree_crowns_sf <- st_concave_hull(distinct_multipoint, ratio = 0.8)
    
    # Ensure the resulting geometry is valid and a polygon
    tree_crowns_sf <- st_make_valid(tree_crowns_sf)
    
    
    # Calculate the area of the tree crown polygon
    crown_area <- as.numeric(st_area(tree_crowns_sf))
    
    # Calculate the density (points per unit area)
    density <- point_count / crown_area
    
    # Add attributes to the crown
    tree_crowns_sf <- st_sf(geometry = tree_crowns_sf, max_height = tree_top$Z, 
                            point_count = point_count, density = density)
    
    # Return the tree crown with height info and the tree top point
    return(list(tree_crowns_sf, tree_top_sf))
    
  }, error = function(e) {
    message("Error processing tree group: ", e)
    return(NULL)
  })
}

count_time <- function(expr) {
  start.time <- Sys.time()  # Start timer
  output <- eval.parent(substitute(expr))  # Evaluate the expression in the parent environment
  end.time <- Sys.time()  # End timer
  time.taken <- end.time - start.time  # Calculate time difference
  
  # Print the time taken with a more readable message
  print(paste0("The task took ", round(time.taken, 2), " seconds to complete."))
  
  return(output)  # Return the output from the evaluated expression
}

compute_pca_metrics <- function(xyz, k = 5) {
  if (!is.matrix(xyz) || nrow(xyz) == 0) {
    stop("Error: XYZ data is not a valid matrix or is empty.")
  }
  
  knn <- get.knn(xyz, k)  # Find k-nearest neighbors
  if (any(is.na(knn$nn.index))) {
    stop("Error: Some nearest neighbor indices are NA, possibly due to empty or malformed data.")
  }
  
  # Initialize metric storage
  linearity <- numeric(nrow(xyz))
  planarity <- numeric(nrow(xyz))
  scatter <- numeric(nrow(xyz))
  omnivariance <- numeric(nrow(xyz))
  anisotropy <- numeric(nrow(xyz))
  sum_eigenvalues <- numeric(nrow(xyz))
  curvature <- numeric(nrow(xyz))
  sphericity <- numeric(nrow(xyz))
  eigenentropy <- numeric(nrow(xyz))
  density <- numeric(nrow(xyz))
  height_variation <- numeric(nrow(xyz))
  
  for (i in 1:nrow(xyz)) {
    # Extract neighbors
    neighbors <- xyz[knn$nn.index[i, ], ]
    
    # Compute covariance matrix
    cov_mat <- cov(neighbors)
    eigen_decomp <- eigen(cov_mat)  # Eigen decomposition
    eigen_vals <- sort(eigen_decomp$values, decreasing = TRUE)  # Sort eigenvalues
    
    # Avoid zero or negative eigenvalues due to numerical issues
    if (any(eigen_vals <= 0)) next  
    
    l1 <- eigen_vals[1]  # Largest eigenvalue
    l2 <- eigen_vals[2]  
    l3 <- eigen_vals[3]  # Smallest eigenvalue
    
    # Compute geometric features
    linearity[i] <- (l1 - l2) / l1
    planarity[i] <- (l2 - l3) / l1
    scatter[i] <- l3 / l1
    omnivariance[i] <- (l1 * l2 * l3)^(1/3)
    anisotropy[i] <- (l1 - l3) / l1
    sum_eigenvalues[i] <- l1 + l2 + l3
    curvature[i] <- l3 / sum_eigenvalues[i]
    sphericity[i] <- l3 / l1
    eigenentropy[i] <- -sum(eigen_vals * log(eigen_vals))
    
    # Compute point density (inverse of mean nearest neighbor distance)
    density[i] <- 1 / mean(knn$nn.dist[i, ])
    
    # Compute local height variation (standard deviation of Z values in neighborhood)
    height_variation[i] <- sd(neighbors[, "Z"])
    
  }
  
  return(data.frame(linearity, planarity, scatter, omnivariance, anisotropy, 
                    sum_eigenvalues, curvature, sphericity, eigenentropy, 
                    density, height_variation))
}

# Function to apply PCA and additional metrics to an LAS object
eigen_features <- function(nlas, k) {
  # Extract XYZ coordinates
  xyz <- as.matrix(nlas@data[, c("X", "Y", "Z")])
  
  if (!is.matrix(xyz) || nrow(xyz) == 0) {
    stop("Error: LAS data is empty or incorrectly formatted.")
  }
  
  start.time <- Sys.time()
  
  # Compute all metrics
  metrics <- compute_pca_metrics(xyz, k = k)
  
  end.time <- Sys.time()
  pca_time_taken <- end.time - start.time
  print(paste0("Time taken to process: ", pca_time_taken))
  
  # Append metrics to LAS data
  nlas@data$linearity <- metrics$linearity
  nlas@data$planarity <- metrics$planarity
  nlas@data$scatter <- metrics$scatter
  nlas@data$omnivariance <- metrics$omnivariance
  nlas@data$anisotropy <- metrics$anisotropy
  nlas@data$sum_eigenvalues <- metrics$sum_eigenvalues
  nlas@data$curvature <- metrics$curvature
  nlas@data$sphericity <- metrics$sphericity
  nlas@data$eigenentropy <- metrics$eigenentropy
  nlas@data$density <- metrics$density
  nlas@data$height_variation <- metrics$height_variation
  
  return(nlas)
}

install.packages(c("lidR", "terra", "sf", "dplyr", "dbscan", "FNN"))



