library(lidR)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(units)



###########################################################################

## Inventories and testing

sample_sf <- st_read("F:/Thesis/Data/Toronto/Sunnybrook/Boundary/Sunnybrook.shp")
sample_sf <- st_transform(sample_sf, 26917)
sample_vect <- vect(sample_sf)

las_trees <- readLAS("F:/Thesis/Data/Toronto/Sunnybrook/Laz/round2/Edits/SB_2019_Veg_clean.laz")

format_number_safe <- function(x) {
  return(gsub("\\.", "_", as.character(x)))  # Replace '.' with '_'
}
# Treetops and Validation
tree_validates <- st_read("F:/Thesis/Data/Toronto/Sunnybrook/Field_Survey/tree_2018.geojson")
tree_validates <- st_transform(tree_validates, 26917)


tree_validates <- st_read("F:/Thesis/Data/Toronto/Sunnybrook/Field_Survey/tallet_survey2018.shp")
# tree_validates <- tree_validates %>%
#   mutate(geometry = st_centroid(geometry))

# Define your functions
lm_small <- function(x) { x * 0.3 + 3}
lm_mid <- function(x) { x * 0.5 + 5 }
lm_large <- function(x) { x * 0.9 + 9 }

# lm_exp <- function(x) {
#   y <- 2.6 * (-(exp(-0.08*(x - 2)) - 1)) + 3
#   y[x < 2] <- 3
#   y[x > 20] <- 5
#   return(y)
# }

tree_func <- list(
  list(func = lm_small, name = "small"),
  list(func = lm_mid, name = "mid"),
  list(func = lm_large, name = "xl"),
  # list(func = lm_exp, name = "exp"),
  list(func = function(x) x + 3, name = "lmf3"),
  list(func = function(x) x + 5, name = "lmf5"),
  list(func = function(x) x + 7, name = "lmf7"),
  list(func = function(x) x + 9, name = "lmf9"),
  list(func = function(x) x + 15, name = "lmf15"),
  list(func = function(x) x + 19, name = "lmf19")
)

tree_func <- (list(func = "survey", name = "survey"))

chm_res <- c(0.5, 0.75, 1, 3, 5)
p2r_res <- c(0.5, 0.75, 1, 3, 5)
lfm_ws <- c(5, 7, 9, 15, 25)
tol_list <- c(1.0, 1.3, 1.5, 2.5, 3.0, 5, 9)
ext_list <- c(3, 5, 7, 9, 11, 15, 19)

# Initialize lists and data frame to store results
ws_df <- data.frame()
results_df <- data.frame()  # Data frame to store results

work_dir <- "F:/Thesis/Data/Toronto/Sunnybrook/ITS/2019/run2/"
results_dir <- paste0(work_dir, "results/")
chm_dir <- paste0(work_dir,"chm/")
ws_dir <- paste0(work_dir,"watershed/")
dalponte_dir <- paste0(work_dir,"dalponte/")
silva_dir <- paste0(work_dir,"silva/")
li_dir <- paste0(work_dir,"li/")

log_file <- paste0(work_dir, "results/process_log.txt")


log_message <- function(message) {
  write(paste(Sys.time(), "-", message), log_file, append = TRUE)
}

log_message("Starting processing...")

# Loop through each combination of CHM resolution and smoothing window size
for (i in 1:length(chm_res)) {
  for (j in 1:length(lfm_ws)) {
    for (q in 1:length(p2r_res)) {
      tryCatch({
        
        file_path <- file.path(paste0(chm_dir, "res_", format_number_safe(chm_res[i]),"_p2r_", format_number_safe(p2r_res[q]), "_ws_", lfm_ws[j], ".tif"))
        
        if (!file.exists(file_path)) {
          start.time <- Sys.time()
          
          # Create the CHM
          chm <- lidR::rasterize_canopy(las_trees, res = chm_res[i], p2r(p2r_res[q]))
          
          # Set NA values to 0 after rasterization
          chm[is.na(chm)] <- 0
          
          # veg_proxy_resample <- terra::resample(veg_proxy, chm, method = "bilinear")
          # chm_proxy <- (chm + 0.1) * veg_proxy_resample
          
          # Mask CHM
          chm_mask <- terra::mask(chm, sample_vect)
          
          # Save the CHM
          file_name <- strsplit(basename(file_path), "\\.")[[1]][1]
          names(chm_mask) <- file_name
          terra::writeRaster(chm_mask, file_path, gdal = c("COMPRESS=LZW"))
          
          end.time <- Sys.time()
          chm_time_taken <- end.time - start.time
          log_message(paste0("Resolution: ", chm_res[i], ", p2r: ", p2r_res[q], ", CHM created in", chm_time_taken, "seconds"))
          
        } else {
          chm_mask <- rast(file_path)
        }
        
        for (n in seq_along(tol_list)) {
          for (m in seq_along(ext_list)) {
            
            file_path_ws_crowns <- file.path(paste0(ws_dir, "crown/res_",
                                                    format_number_safe(chm_res[i]),
                                                    "_ws_", lfm_ws[j],
                                                    "_p2r_", format_number_safe(p2r_res[q]),
                                                    "_tol_", format_number_safe(tol_list[n]),
                                                    "_ext_", ext_list[m], ".shp"))
            
            if (!file.exists(file_path_ws_crowns)) {
              
              ws_start.time <- Sys.time()
              
              # Segment trees using watershed segmentation
              ws_las <- lidR::segment_trees(las_trees,
                                            lidR::watershed(chm_mask, th_tree = 2.0,
                                                            tol = tol_list[n], ext = ext_list[m]))
              
              ws_end.time <- Sys.time()
              ws_time_taken <- ws_end.time - ws_start.time
              log_message(paste0("Resolution: ", chm_res[i],
                                 ", p2r: ", p2r_res[q],
                                 ", Watershed segmentation completed in", ws_time_taken, "seconds"))
              
              # Calculate crown metrics for watershed segmentation
              crowns_ws <- lidR::crown_metrics(ws_las, func = .stdtreemetrics, geom = "convex")
              
              tree_max <- ws_las@data %>%
                group_by(treeID) %>%
                dplyr::slice_max(Z) %>%
                mutate(highest_z = Z) %>%
                select(treeID, highest_z) %>%
                distinct() %>%
                as.data.frame()
              
              crowns_ws_joined <- crowns_ws %>%
                left_join(tree_max, by = "treeID")
              
              # 1. Calculate Area and Perimeter
              crowns_ws_joined$area_m2 <- as.numeric(st_area(crowns_ws_joined))
              crowns_ws_joined$perimeter_m <- as.numeric(st_length(crowns_ws_joined))
              
              # 2. Calculate IPQ
              crowns_ws_joined$IPQ <- (4 * pi * crowns_ws_joined$area_m2) / (crowns_ws_joined$perimeter_m^2)
              
              # 3. Calculate Mean Number of Returns for Each TreeID
              returns_per_tree <- ws_las@data %>%
                group_by(treeID) %>%
                summarize(mean_returns = mean(NumberOfReturns, na.rm = TRUE))
              
              # 4. Normalize Mean Number of Returns (0-1 scale)
              returns_per_tree <- returns_per_tree %>%
                mutate(mean_returns_norm = (mean_returns - min(mean_returns)) / (max(mean_returns) - min(mean_returns)))
              
              # 5. Join back to crown polygons
              crowns_ws_joined <- crowns_ws_joined %>%
                left_join(returns_per_tree, by = "treeID")
              
              # 6. Create Munzinger-style score
              crowns_ws_joined <- crowns_ws_joined %>%
                mutate(tree_score = IPQ + mean_returns_norm)
              
              # 7. Define a threshold value
              threshold_value <- 1.6 # <-- you can adjust this later
              
              # 8. Filter crowns based on threshold
              crowns_ws_joined_filtered <- crowns_ws_joined %>%
                filter(tree_score > threshold_value)
              
              st_write(crowns_ws_joined_filtered, file_path_ws_crowns)
            }
          }
        }
        
        # Tree top identification
        for (p in 1:length(tree_func)) {
          tree_function <- tree_func[[p]][[1]]
          func_name <- tree_func[[p]][[2]]
          
          silva_crowns_file_path <- file.path(paste0(silva_dir, "crown/res_", 
                                                     format_number_safe(chm_res[i]), 
                                                     "_p2r_", format_number_safe(p2r_res[q]), 
                                                     "_ws_", lfm_ws[j], 
                                                     "_func_", func_name, ".shp"))
          
          if (!file.exists(silva_crowns_file_path)) {
            
            silva_start.time <- Sys.time()
            
            # ttops <- lidR::locate_trees(chm_mask, lmf(tree_function))
            tree_vect <- vect(tree_validates)
            tree_vect <- tree_vect[,1]
            tree_vect$treeID <- 1:nrow(tree_vect) 
            
            tree_z <- terra::extract(chm_mask, tree_vect, fun = max, na.rm = TRUE)
            names(tree_z) <- c("treeID", "z")
            tree_vect$z <- tree_z$z
            
            tree_sf <- st_as_sf(tree_vect)
            tree_sf <- tree_sf %>%
              mutate(geometry = st_centroid(geometry))
            
            ttops <-   tree_sf
            
            
            silva_las <- segment_trees(las_trees, 
                                       silva2016(chm_mask,ttops))
            
            silva_end.time <- Sys.time()
            silva_time_taken <- silva_end.time - silva_start.time
            log_message(paste0("Resolution: ", chm_res[i], 
                               ", p2r: ", p2r_res[q], 
                               ", function name: ", func_name, 
                               ", Silva segmentation completed in", silva_time_taken, "seconds"))
            
            # Calculate crown metrics
            # Calculate crown metrics
            crowns_silva <- crown_metrics(silva_las, func = .stdtreemetrics, geom = "convex")
            
            tree_max_silva <- silva_las@data %>%
              group_by(treeID) %>%
              slice_max(Z) %>%
              mutate(highest_z = Z) %>%
              select(treeID, highest_z) %>%
              distinct() %>%
              as.data.frame()
            
            crowns_silva_joined <- crowns_silva %>%
              left_join(tree_max_silva, by = "treeID")
            
            # 1. Area and Perimeter
            crowns_silva_joined$area_m2 <- as.numeric(st_area(crowns_silva_joined))
            crowns_silva_joined$perimeter_m <- as.numeric(st_length(crowns_silva_joined))
            
            # 2. IPQ
            crowns_silva_joined$IPQ <- (4 * pi * crowns_silva_joined$area_m2) / (crowns_silva_joined$perimeter_m^2)
            
            # 3. Mean Number of Returns
            returns_per_tree_silva <- silva_las@data %>%
              group_by(treeID) %>%
              summarize(mean_returns = mean(NumberOfReturns, na.rm = TRUE))
            
            # 4. Normalize
            returns_per_tree_silva <- returns_per_tree_silva %>%
              mutate(mean_returns_norm = (mean_returns - min(mean_returns)) / (max(mean_returns) - min(mean_returns)))
            
            # 5. Join
            crowns_silva_joined <- crowns_silva_joined %>%
              left_join(returns_per_tree_silva, by = "treeID")
            
            # 6. Score
            crowns_silva_joined <- crowns_silva_joined %>%
              mutate(tree_score = IPQ + mean_returns_norm)
            
            # 7. Threshold
            threshold_value <- 1.6
            crowns_silva_joined_filtered <- crowns_silva_joined %>%
              filter(tree_score > threshold_value)
            
            # 8. Save
            st_write(crowns_silva_joined_filtered, silva_crowns_file_path)
          }
          
          ## Dalponte processing
          dalponte_crowns_file_path <- file.path(paste0(dalponte_dir, "crown/res_", 
                                                        format_number_safe(chm_res[i]), 
                                                        "_p2r_", format_number_safe(p2r_res[q]), 
                                                        "_ws_", lfm_ws[j], 
                                                        "_func_", func_name, ".shp"))
          
          if (!file.exists(dalponte_crowns_file_path)) {
            dalponte_start.time <- Sys.time()
            
            dalponte_las <- segment_trees(las_trees, 
                                          dalponte2016(chm_mask, 
                                                       ttops,
                                                       th_tree = 2,
                                                       th_seed = 0.45,
                                                       th_cr = 0.55,
                                                       max_cr = 10))
            
            dalponte_end.time <- Sys.time()
            dalponte_time_taken <- dalponte_end.time - dalponte_start.time
            log_message(paste0("Resolution: ", chm_res[i], 
                               ", p2r: ", p2r_res[q], 
                               ", function name: ", func_name, 
                               ", Silva segmentation completed in", silva_time_taken, "seconds"))
            
            # Calculate crown metrics
            crowns_dalponte <- crown_metrics(dalponte_las, func = .stdtreemetrics, geom = "convex")
            
            tree_max_dalponte <- dalponte_las@data %>%
              group_by(treeID) %>%
              slice_max(Z) %>%
              mutate(highest_z = Z) %>%
              select(treeID, highest_z) %>%
              distinct() %>%
              as.data.frame()
            
            crowns_dalponte_joined <- crowns_dalponte %>%
              left_join(tree_max_dalponte, by = "treeID")
            
            # 1. Area and Perimeter
            crowns_dalponte_joined$area_m2 <- as.numeric(st_area(crowns_dalponte_joined))
            crowns_dalponte_joined$perimeter_m <- as.numeric(st_length(crowns_dalponte_joined))
            
            # 2. IPQ
            crowns_dalponte_joined$IPQ <- (4 * pi * crowns_dalponte_joined$area_m2) / (crowns_dalponte_joined$perimeter_m^2)
            
            # 3. Mean Number of Returns
            returns_per_tree_dalponte <- dalponte_las@data %>%
              group_by(treeID) %>%
              summarize(mean_returns = mean(NumberOfReturns, na.rm = TRUE))
            
            # 4. Normalize
            returns_per_tree_dalponte <- returns_per_tree_dalponte %>%
              mutate(mean_returns_norm = (mean_returns - min(mean_returns)) / (max(mean_returns) - min(mean_returns)))
            
            # 5. Join
            crowns_dalponte_joined <- crowns_dalponte_joined %>%
              left_join(returns_per_tree_dalponte, by = "treeID")
            
            # 6. Score
            crowns_dalponte_joined <- crowns_dalponte_joined %>%
              mutate(tree_score = IPQ + mean_returns_norm)
            
            # 7. Threshold
            threshold_value <- 1.6
            crowns_dalponte_joined_filtered <- crowns_dalponte_joined %>%
              filter(tree_score > threshold_value)
            
            # 8. Save
            st_write(crowns_dalponte_joined_filtered, dalponte_crowns_file_path)
            
          }
          
          
        }
        
        # Save results to results_df
        results_df <- rbind(results_df, data.frame(
          CHM_res = chm_res[i],
          CHM_ws = lfm_ws[j],
          p2r_res = p2r_res[q],
          function_name = func_name,
          Watershed_Crowns_Shapefile = file_path_ws_crowns,
          Dalponte_Crowns_Shapefile = dalponte_crowns_file_path,
          Silva_Crowns_Shapefile = silva_crowns_file_path
        ))
        
      }, error = function(e) {
        log_message(paste("Error occurred at iteration for CHM res", 
                          chm_res[i], "_p2r_", p2r_res[q], 
                          "and WS", lfm_ws[j], ":", e$message))
      })
    }
  }
}

head(ws_df)
head(results_df)

# Save the results to CSV after the loop
write.csv(ws_df, paste0(results_dir, "watershed_results.csv"), row.names = TRUE)
write.csv(results_df, paste0(results_dir, "tree_top_results.csv"), row.names = TRUE)


# Print final results (data frame containing all results)
print(results_df)
print(ws_df)


## li model
li_df <- data.frame()
dt1 <- c(1, 2, 3, 4, 5, 6, 7,8, 9, 11, 15)
radius_li <- c(1, 2, 3, 4, 5, 6, 7,8, 9, 11, 15)


for (i in 1:length(dt1)) {
  for (j in 1:length(radius_li)) {
    
    start.time <- Sys.time()
    
    # File paths for LAS and shapefile outputs
    li_crowns_file_path <- file.path(paste0(li_dir,"crown/li23_dt1_", dt1[i], "_r_", radius_li[j], ".shp"))
    
    # Segment trees using Li 2012 model
    li_las <- segment_trees(las_trees, li2012(
      dt1[i],
      R = radius_li[j],
      Zu = 15,
      hmin = 2,
      speed_up = 20
    ))
    
    end.time <- Sys.time()
    li_time_taken <- end.time - start.time
    
    # Crown metrics
    li_crowns <- crown_metrics(li_las, func = .stdtreemetrics, geom = "convex")
    
    # Fix invalid geometries
    li_crowns <- li_crowns %>%
      st_make_valid() %>%
      filter(!st_is_empty(.))
    
    # Highest Z per tree
    tree_max_li <- li_las@data %>%
      group_by(treeID) %>%
      slice_max(Z) %>%
      mutate(highest_z = Z) %>%
      select(treeID, highest_z) %>%
      distinct() %>%
      as.data.frame()
    
    li_crowns_joined <- li_crowns %>%
      left_join(tree_max_li, by = "treeID")
    
    # 1. Area and Perimeter
    li_crowns_joined$area_m2 <- as.numeric(st_area(li_crowns_joined))
    li_crowns_joined$perimeter_m <- as.numeric(st_length(li_crowns_joined))
    
    # 2. IPQ
    li_crowns_joined$IPQ <- (4 * pi * li_crowns_joined$area_m2) / (li_crowns_joined$perimeter_m^2)
    
    # 3. Mean Number of Returns
    returns_per_tree_li <- li_las@data %>%
      group_by(treeID) %>%
      summarize(mean_returns = mean(NumberOfReturns, na.rm = TRUE))
    
    returns_per_tree_li <- returns_per_tree_li %>%
      mutate(mean_returns_norm = (mean_returns - min(mean_returns)) / (max(mean_returns) - min(mean_returns)))
    
    # 4. Join returns to crowns
    li_crowns_joined <- li_crowns_joined %>%
      left_join(returns_per_tree_li, by = "treeID")
    
    # 5. Munzinger Score
    li_crowns_joined <- li_crowns_joined %>%
      mutate(tree_score = IPQ + mean_returns_norm)
    
    # 6. Save raw crowns
    li_crowns_file_path_raw <- gsub(".shp", "_raw.shp", li_crowns_file_path)
    st_write(li_crowns_joined, li_crowns_file_path_raw, delete_layer = TRUE)
    
    # 7. Threshold and save filtered crowns
    threshold_value <- 1.3
    li_crowns_joined_filtered <- li_crowns_joined %>%
      filter(tree_score > threshold_value)
    
    st_write(li_crowns_joined_filtered, li_crowns_file_path, delete_layer = TRUE)
    
    # 8. Add metrics to the summary dataframe
    li_df <- rbind(li_df, data.frame(
      dt1 = dt1[i],
      R = radius_li[j],
      li_time = as.numeric(li_time_taken) # Ensure it's numeric
    ))
  }
}

# Save the summary dataframe
write.csv(li_df, file.path("F:/Thesis/Data/Toronto/Sunnybrook/ITS/2019/run1/results", "li_results.csv"), row.names = FALSE)

# Save the summary dataframe
write.csv(li_df, file.path("F:/Thesis/Data/Toronto/Sunnybrook/ITS/2019/run1/results", "li_results.csv"), row.names = FALSE)


# # Apply general and regional height-DBH equations
# li_crowns <- li_crowns %>%
#   mutate(
#     crown_area = st_area(geometry),  # Area in square meters
#     dbh_combined = 0.05 * (Z^1.3) * (as.numeric(crown_area)^0.5),  # General equation with height and area
#     dbh_general = 0.12 * (Z^1.4),        # General temperate forest
#     dbh_urban = 0.14 * (Z^1.3),         # Urban adjustment
#     dbh_mixed = 0.11 * (Z^1.35),        # Southern Ontario mixed
#     dbh_hardwoods = 0.13 * (Z^1.4)      # Southern Ontario hardwoods
#   )


#Validation

tree_validates <- st_read("F:/Thesis/Data/Toronto/Sunnybrook/Field_Survey/tree_2018.geojson")
tree_validates <- st_transform(tree_validates, 26917)
tree_validates <- tree_validates %>%
  mutate(geometry = st_centroid(geometry))

tree_validates <- st_buffer(tree_validates, 0.5)

plot(tree_validates[1])


# tree_valid <- st_read("F:/Thesis/AnnexTrees/Data/Survey/Annex_Survey.shp")
# tree_valid <- st_transform(tree_valid, 26917)
# tree_validates <- tree_valid[12:13  
tree_validates$ID <- 1:nrow(tree_validates)

# Initialize the results dataframe
silva_df <- data.frame()

# Define the directory and list the shapefiles
silva_crown <- paste0(silva_dir, "crown/")
silva_list <- list.files(silva_crown, full.names = TRUE, pattern = ".shp")

# Assign names to the silva_crowns list based on the file names (without extensions)
silva_names <- tools::file_path_sans_ext(basename(silva_list))
silva_crowns <- setNames(lapply(silva_list, function(x) st_read(x)), silva_names)

length(silva_crowns)

# Process each crown and store results
for (i in 1:length(silva_crowns)) {
  tryCatch({
    s_names <- silva_names[i]
    crowns_silva <- silva_crowns[[i]]
    crowns_silva$Area <- st_area(crowns_silva)
    crowns_silva <- crowns_silva[crowns_silva$Area >= set_units(3, "m^2"), ]
    
    silva_test <- st_intersection(tree_validates, crowns_silva)
    
    if (nrow(silva_test) == 0) {
      # If no intersections found, skip to next iteration safely
      silva_df <- rbind(silva_df, data.frame(
        run_name = s_names,
        r2 = NA,
        p_value = NA,
        rmse = NA,
        num_trees = 0,
        status = "No intersections"
      ))
    } else {
      # Drop geometry and select only needed columns
      test_nums <- silva_test %>%
        st_drop_geometry() %>%
        select(ID, Z, HTOTAL)
      
      # Handle multiple matches: keep only the closest one
      test_nums <- test_nums %>%
        group_by(ID) %>%
        slice_max(HTOTAL) %>%
        ungroup()
      
      # Check again: if no valid matches after slice_min
      if (nrow(test_nums) == 0) {
        silva_df <- rbind(silva_df, data.frame(
          run_name = s_names,
          r2 = NA,
          p_value = NA,
          rmse = NA,
          num_trees = 0,
          status = "No valid matches after filtering"
        ))
      } else {
        # Now perform regression
        model <- lm(Z ~ HTOTAL, data = test_nums)
        r_2 <- summary(model)$r.squared
        individual_p_value <- summary(model)$coefficients["HTOTAL", "Pr(>|t|)"]
        rmse_value <- sqrt(mean((test_nums$Z - test_nums$HTOTAL)^2))
        
        # Store results
        silva_df <- rbind(silva_df, data.frame(
          run_name = s_names,
          r2 = r_2,
          p_value = individual_p_value,
          rmse = rmse_value,
          num_trees = nrow(test_nums),
          status = "Regression completed"
        ))
      }
    }
  }, error = function(e) {
    message(paste("Error in iteration for",  s_names, ":", e$message))
    silva_df <- rbind(silva_df, data.frame(
      run_name = s_names,
      r2 = NA,
      p_value = NA,
      rmse = NA,
      num_trees = NA,
      status = paste("Error:", e$message)
    ))
  })
}

# View the results
print(silva_df)

## Dalponte outputs

# Initialize the results dataframe for dalponte
dalponte_df <- data.frame()

# Define the directory and list the shapefiles
dalponte_crown <- paste0(dalponte_dir, "/crown")
dalponte_list <- list.files(dalponte_crown, full.names = TRUE, pattern = ".shp")

# Assign names to the dalponte_crowns list based on the file names (without extensions)
dalponte_names <- tools::file_path_sans_ext(basename(dalponte_list))
dalponte_crowns <- setNames(lapply(dalponte_list, function(x) st_read(x)), dalponte_names)

length(dalponte_crowns)

# Process each crown and store results
for (i in 1:length(dalponte_crowns)) {
  tryCatch({
    
    d_names <- dalponte_names[i]
    crowns_dalponte <- dalponte_crowns[[i]]
    crowns_dalponte$Area <- st_area(crowns_dalponte)
    crowns_dalponte <- crowns_dalponte[crowns_dalponte$Area >= set_units(3, "m^2"), ]
    
    dalponte_test <- st_intersection(tree_validates, crowns_dalponte)
    num_trees <- nrow(dalponte_test)  # Count intersected trees
    
    # Drop geometry and select only needed columns
    test_nums <- dalponte_test %>%
      st_drop_geometry() %>%
      select(ID, Z, HTOTAL)
    
    # Handle multiple matches: keep only the closest one
    test_nums <- test_nums %>%
      group_by(ID) %>%
      slice_max(HTOTAL) %>%
      ungroup()
    
    # Perform regression analysis
    model <- lm(Z ~ HTOTAL, data = test_nums)
    r_2 <- summary(model)$r.squared
    individual_p_value <- summary(model)$coefficients["HTOTAL", "Pr(>|t|)"]
    rmse_value <- sqrt(mean((test_nums$Z - test_nums$HTOTAL)^2))
    
    # Store results in the dataframe
    dalponte_df <- rbind(dalponte_df, data.frame(
      run_name = d_names,
      r2 = r_2,
      p_value = individual_p_value,
      rmse = rmse_value,
      num_trees = num_trees,
      status = "Regression completed"
    ))
  }, error = function(e) {
    message(paste("Error in iteration for", d_names, ":", e$message))
    dalponte_df <- rbind(dalponte_df, data.frame(
      run_name = d_names,
      r2 = NA,
      p_value = NA,
      rmse = NA,
      num_trees = NA,
      status = paste("Error:", e$message)
    ))
  })
}

# View the results
print(dalponte_df)

## Watershed outputs

# Initialize the results dataframe for watershed
watershed_df <- data.frame()

# Define the directory and list the shapefiles
watershed_crown <- paste0(ws_dir,"crown/")
watershed_list <- list.files(watershed_crown, full.names = TRUE, pattern = ".shp")

# Assign names to the watershed_crowns list based on the file names (without extensions)
watershed_names <- tools::file_path_sans_ext(basename(watershed_list))
watershed_crowns <- setNames(lapply(watershed_list, function(x) st_read(x)), watershed_names)

length(watershed_crowns)


for (i in 1:length(watershed_crowns)) {
  tryCatch({
    # Perform intersection
    w_names <- watershed_names[i]
    crowns_ws <- watershed_crowns[[i]]
    crowns_ws$Area <- st_area(crowns_ws)
    crowns_ws <-  crowns_ws[crowns_ws$Area >= set_units(3, "m^2"), ]
    
    # Perform intersection with tree_validates
    watershed_inter <- st_intersection(tree_validates, crowns_ws)
    
    # Count the number of trees in the intersection
    num_trees <- nrow(watershed_inter)
    
    # Drop geometry and select only needed columns
    test_nums <- watershed_inter %>%
      st_drop_geometry() %>%
      select(ID, Z, HTOTAL)
    
    # Handle multiple matches: keep only the closest one
    test_nums <- test_nums %>%
      group_by(ID) %>%
      slice_max(HTOTAL) %>%
      ungroup()
    
    # Perform regression analysis
    model <- lm(Z ~ HTOTAL, data = test_nums)
    r_2 <- summary(model)$r.squared
    individual_p_value <- summary(model)$coefficients["HTOTAL", "Pr(>|t|)"]
    rmse_value <- sqrt(mean((test_nums$Z - test_nums$HTOTAL)^2))
    
    # Store results in the dataframe
    watershed_df <- rbind(watershed_df, data.frame(
      run_name = w_names,
      r2 = r_2,
      p_value = individual_p_value,
      rmse = rmse_value,
      num_trees = num_trees,
      status = "Regression completed"
    ))
  }, error = function(e) {
    message(paste("Error in iteration for", w_names, ":", e$message))
    # Record the error in the dataframe
    watershed_df <- rbind(watershed_df, data.frame(
      run_name = w_names,
      r2 = NA,
      p_value = NA,
      rmse = NA,
      num_trees = NA,
      status = paste("Error:", e$message)
    ))
  })
}

# View the results
print(watershed_df)

# Crowns Dirs

# Initialize the results dataframe
li_df <- data.frame()

# Define the directory and list the shapefiles
li_crown <- paste0(li_dir,"crown")
li_list <- list.files(li_crown, full.names = TRUE, pattern = ".shp")

# Assign names to the li_crowns list based on the file names (without extensions)
li_names <- tools::file_path_sans_ext(basename(li_list))
li_crowns <- setNames(lapply(li_list, function(x) st_read(x)), li_names)

length(li_crowns)

# Process each crown and store results
for (name in names(li_crowns)) {
  tryCatch({
    crowns_li <- li_crowns[[name]]
    crowns_li$Area <- st_area(crowns_li)
    crowns_li <- crowns_li[crowns_li$Area >= set_units(3, "m^2"), ]
    
    li_test <- st_intersection(tree_validates, crowns_li)
    num_trees <- nrow(li_test)  # Count intersected trees
    
    test_nums <- st_drop_geometry(li_test)
    
    # Drop geometry and select only needed columns
    test_nums <- li_test %>%
      st_drop_geometry() %>%
      select(ID, Z, HTOTAL)
    
    # Handle multiple matches: keep only the closest one
    test_nums <- test_nums %>%
      group_by(ID) %>%
      slice_max(HTOTAL) %>%
      ungroup()
    
    # Skip if insufficient data for regression
    if (nrow(test_nums) < 2 || sd(test_nums$HTOTAL) == 0 || sd(test_nums$Z) == 0) {
      message(paste("Skipping regression for", name, "due to insufficient variability"))
      li_df <- rbind(li_df, data.frame(
        run_name = name,
        r2 = NA,
        p_value = NA,
        rmse = NA,
        num_trees = num_trees,
        status = "Removed: Insufficient variability"
      ))
      next
    }
    
    # Perform regression analysis
    model <- lm(Z ~ HTOTAL, data = test_nums)
    r_2 <- summary(model)$r.squared
    individual_p_value <- summary(model)$coefficients["HTOTAL", "Pr(>|t|)"]
    rmse_value <- sqrt(mean((test_nums$Z - test_nums$HTOTAL)^2))
    
    # Store results in the dataframe
    li_df <- rbind(li_df, data.frame(
      run_name = name,
      r2 = r_2,
      p_value = individual_p_value,
      rmse = rmse_value,
      num_trees = num_trees,
      status = "Regression completed"
    ))
  }, error = function(e) {
    message(paste("Error in iteration for", name, ":", e$message))
    li_df <- rbind(li_df, data.frame(
      run_name = name,
      r2 = NA,
      p_value = NA,
      rmse = NA,
      num_trees = NA,
      status = paste("Error:", e$message)
    ))
  })
}

# View the results
print(li_df)

silva_df$model <- "Silva"
dalponte_df$model <- "Dalponte"
watershed_df$model <- "Watershed"
li_df$model <- "Li"


# write.csv(silva_df, paste0(results_dir, "silva.csv"))
# write.csv(dalponte_df, paste0(results_dir,"dalponte.csv"))
# write.csv(watershed_df, paste0(results_dir,"watershed.csv"))
# write.csv(li_df, paste0(results_dir,"li.csv"))

names(silva_df)
names(watershed_df)

valid_df <- do.call(rbind, list(silva_df, watershed_df, dalponte_df, li_df))
# write.csv(valid_df, paste0(results_dir,"valid_df.csv"))

highest_r2_per_model <- valid_df %>%
  group_by(model) %>%
  slice_max(order_by = r2, n = 1, with_ties = FALSE) %>%
  ungroup()

best_rmse_per_model <- valid_df %>%
  group_by(model) %>%
  slice_min(order_by = rmse, n = 1, with_ties = FALSE)

highest_r2_per_model
best_rmse_per_model
