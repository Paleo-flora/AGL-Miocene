

## Fig Peruvian Biome MAT vs MAP [Ochoa et al., xx]  ## Set working directory
setwd("/Users/dianaochoa/Library/CloudStorage/OneDrive-UniversidaddeSalamanca/Peruvian coastal MAT MAP data")

# # # # Load Diana's data
library(readxl)
library (vegan)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(maps)
library(cowplot) 
library(sp)
library(raster)
library(geodata)
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(lattice)

ecoR <- st_read("/Users/dianaochoa/Library/CloudStorage/OneDrive-UniversidaddeSalamanca/Peruvian coastal MAT MAP data/Ecoregiones Sernanp/Ecoregiones_SERNANP.shp")
str(ecoR)#geometry - 168926 obs
 
ecoR_valid <- st_make_valid(ecoR) ##841 elemnts - geometry 168926 obs
str(ecoR_valid)
ECO_COD_df <- data.frame(ECO_COD = ecoR_valid$ECO_COD)  ## info of 841 polygons
ECO_COD <-as.vector(ECO_COD_df)
ECO_NAME_df <- data.frame(ECO_COD = ecoR_valid$ECO_NOMBRE)  ## info of 841 polygons


# Create a mapping table
mapping_table <- data.frame(ECO_COD = unique(ECO_COD_df$ECO_COD),
                            ECO_NAME = unique(ECO_NAME_df$ECO_COD))
print(mapping_table)

raster_resolution <- 0.01  # Set the desired raster resolution
raster_extent <- extent(ecoR)  # Set the extent based on the shapefile extent
rasterized <- raster(ncol = ncol(ecoR), nrow = nrow(ecoR), ext = raster_extent)
rasterized <- rasterize(ecoR, field = "ECO_COD", resolution = raster_resolution, extent = raster_extent)
str(rasterized)
# Print the rasterized object
print(rasterized)

wc1 <- worldclim_global(var = "bio", res = 2.5, path = getwd())
sel1 <- wc1[[1]]#MAT
sel12 <- wc1[[12]]#MAP
# Convert sel12 SpatRaster to raster object
sel1_raster <- as(sel1, "Raster")
sel12_raster <- as(sel12, "Raster")

# Extract coordinates for each polygon
polygon_coordinates <- st_coordinates(ecoR_valid$geometry)
poly_coordinates<-as.data.frame(polygon_coordinates)
# Access the x and y coordinates
x_coordinates <- polygon_coordinates[, "X"]
y_coordinates <- polygon_coordinates[, "Y"]
polygon <- polygon_coordinates[, "L2"]

#### Coord per polygon
polygon_objects <- list()
# Iterate over the polygon numbers (1 to 841)
for (i in 1:841) {
  # Subset the X and Y coordinates based on the polygon number
  x_subset <- x_coordinates[polygon == i]
  y_subset <- y_coordinates[polygon == i]
  
  coordinates <- cbind(x_subset, y_subset)
  
  # Create an object for each group of coordinates
  object_name <- paste0("Pol_", i)
  assign(object_name, coordinates, envir = .GlobalEnv)
  
  # Store the object in the polygon_objects list
  polygon_objects[[i]] <- get(object_name)
}
ls()


create_raster_files <- function() {
  crs <- CRS("+proj=longlat +datum=WGS84")
  
  for (i in 1:841) {
    pol_name <- paste0("Pol_", i)
    pol_i <- get(pol_name, envir = .GlobalEnv)  # Retrieve Pol_i matrix from the global environment
    
    # Creating SpatialPolygon object
    pol_i_polygon <- SpatialPolygons(list(Polygons(list(Polygon(pol_i)), pol_name)))
    proj4string(pol_i_polygon) <- crs
    
    # Creating rasterized polygon
    extent <- extent(pol_i_polygon)
    raster_empty <- raster(extent, ncol = 100, nrow = 100)
    rasterized_pol_i <- rasterize(pol_i_polygon, raster_empty, field = 1)
    projection(rasterized_pol_i) <- crs
    
    # Assigning raster to the global environment
    raster_name <- paste0("RastPol_", i)
    assign(raster_name, rasterized_pol_i, envir = .GlobalEnv)
  }
}

create_raster_files()
ls()


MAT.Peru <- list()
MAP.Peru <- list()

# Loop through the raster objects and extract values for each
for (i in 1:841) {
  raster_name <- paste0("RastPol_", i)
  raster_obj <- get(raster_name, envir = .GlobalEnv)
  
  # Extract sel1 values for current raster object
  sel1_values <- extract(sel1_raster, as(raster_obj, "SpatialPoints"))
  
  # Extract sel12 values for current raster object
  sel12_values <- extract(sel12_raster, as(raster_obj, "SpatialPoints"))
  
  # Save the extracted values to the respective lists
  MAT.Peru[[i]] <- sel1_values
  MAP.Peru[[i]] <- sel12_values
}

save(MAT.Peru, file = "MAT.Peru.RData")
save(MAP.Peru, file = "MAP.Peru.RData")


# Create a data frame with MAT.Peru, MAP.Peru, and ECO_COD_vector
ECO_COD_vector <- ECO_COD$ECO_COD
plot_data <- data.frame(MAT.Peru = unlist(MAT.Peru),
                        MAP.Peru = unlist(MAP.Peru),
                        ECO_COD = rep(ECO_COD_vector, sapply(MAT.Peru, length)))

# Filter the data based on the desired biomes/conditions
filtered_data <- subset(plot_data, MAP.Peru <= 1000 & MAT.Peru >= 15 & MAT.Peru <= 25)
Desert_data <- subset(plot_data, ECO_COD == "L0" & MAT.Peru <= 25)
DryForest_data <- subset(plot_data, ECO_COD == "G0" & MAT.Peru <= 25)

Desert_data <- subset(plot_data, ECO_COD == "L0" )
DryForest_data <- subset(plot_data, ECO_COD == "G0" )

min(Desert_data$MAT.Peru)
max(Desert_data$MAT.Peru)
mean(Desert_data$MAT.Peru)
min(DryForest_data$MAT.Peru)
max(DryForest_data$MAT.Peru)
mean(DryForest_data$MAT.Peru)

min(Desert_data$MAP.Peru)
max(Desert_data$MAP.Peru)
mean(Desert_data$MAP.Peru)
min(DryForest_data$MAP.Peru)
max(DryForest_data$MAP.Peru)
mean(DryForest_data$MAP.Peru)

DF.MAP.q25 <- quantile(DryForest_data $MAP.Peru, probs = 0.25)  # 25th percentile (first quartile)
DF.MAP.q75 <- quantile(DryForest_data $MAP.Peru, probs = 0.75)  # 75th percentile (third quartile)

CD.MAP.q25 <- quantile(Desert_data$MAP.Peru, probs = 0.25)  # 25th percentile (first quartile)
CD.MAP.q75 <- quantile(Desert_data$MAP.Peru, probs = 0.75)  # 75th percentile (third quartile)

DF.MAT.q25 <- quantile(DryForest_data $MAT.Peru, probs = 0.25)  # 25th percentile (first quartile)
DF.MAT.q75 <- quantile(DryForest_data $MAT.Peru, probs = 0.75)  # 75th percentile (third quartile)

CD.MAT.q25 <- quantile(Desert_data$MAT.Peru, probs = 0.25)  # 25th percentile (first quartile)
CD.MAT.q75 <- quantile(Desert_data$MAT.Peru, probs = 0.75)  # 75th percentile (third quartile)


unique_colors <- rainbow(nlevels(as.factor(ECO_COD_vector)))

Biome_plot <- ggplot(plot_data, aes(x = MAP.Peru, y = MAT.Peru, color = ECO_COD)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = unique_colors) +
  labs(x = "MAP (mm)", y = "MAT (ºC)", title = "Peruvian Biomes") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))+ 
  geom_vline(xintercept = c(38, 72), linetype = "dashed", color = "gold4") + ## Leaves MAP
  geom_hline(yintercept = c(19, 20, 22, 23), linetype = "dashed", color = c("red", "red", "gold4", "gold4")) + # today MAT/miocene MAT
  geom_vline(xintercept = c(585, 798), linetype = "dashed", color = "navy blue") + ## holocene pollen-based MAP
  geom_vline(xintercept = c(668, 850), linetype = "dashed", color = "gold4") + ## miocene pollen-based MAP
  geom_vline(xintercept = 6, linetype = "dashed", color = "red")

ggsave("Biomes_MATMAP.jpg", plot = Biome_plot, width = 8, height = 6, units = "in", dpi = 150)


############################
## Grid Analysis  - Fig 4B
############################

y.grid <-MAT.grid<-seq(15, 25, by = 1) ###
x.grid <-MAP.grid<-seq(0, 100, by = 12.5) ##


# Function to find the closest grid point
find_closest_smallest_grid_point <- function(val, grid) {
  filtered_grid <- grid[grid <= val]
  idx <- which.min(val - filtered_grid)  # Find the index of the smallest difference
  return(idx)
}


# Group IDs assignment based on the intervals
x_intervals <- sapply(filtered_data$MAP.Peru, function(x_val) find_closest_smallest_grid_point(x_val, x.grid))
y_intervals <- sapply(filtered_data$MAT.Peru, function(y_val) find_closest_smallest_grid_point(y_val, y.grid))
group_ids2 <- (x_intervals - 1) * length(y.grid) + y_intervals+1

y.lab<-c("A","B","C","D","E","F")
y.label <- y.lab[y_intervals]
group_ids<-cbind(y.label, x_intervals)

# Data frame of group_ids with xy_matrix
Gp_df <- data.frame(
  X = filtered_data$MAP.Peru,
  Y = filtered_data$MAT.Peru,
  Group_ID = group_ids
)
Gp_df$gp_ID <- paste(Gp_df$Group_ID.y.label, Gp_df$Group_ID.x_intervals, sep = "")
str(Gp_df)
head(Gp_df)

summary_table <- table(Gp_df$gp_ID) # count occurrences per interval defined by group_ids
str(summary_table)
head(summary_table)
sum(summary_table)
Gp_df $ECO_data <- filtered_data$ECO_COD # Add ECO_data 

eco_table <- table(Gp_df$gp_ID, Gp_df$ECO_data) # count occurrences of each ECO_data per interval defined by group_ids
str(eco_table) 
head(eco_table)
#print(eco_table)
proportion_df <- prop.table(eco_table, margin = 1)  # Calculate the proportions of ECO_data categories per dataframe

# Export the data frame to a txt file
write.table(proportion_df, file = "Biome proportion.txt", sep = "\t", row.names = TRUE)

# Create a data frame with the positions and group_ids values
positions <- expand.grid(x = seq(0, num_columns - 1), y = seq(15, 15 + num_rows - 1))
positions$group_ids <- positions$x * num_rows + positions$y - 14

# Create the xy plot using ggplot2
ggplot(positions, aes(x, y, label = group_ids)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(size = 3) +
  scale_x_continuous(
    breaks = seq(0, num_columns - 1, by = 6),
    labels = x.grid[seq(1, num_columns, by = 6)],
    sec.axis = dup_axis()
  ) +
  scale_y_continuous(
    breaks = seq(15, 15 + num_rows - 1, by = 1),
    labels = y.grid,
    sec.axis = dup_axis()
  ) +
  labs(x = "X Position", y = "Y Position", title = "Index (group_id) of Elements") +
  theme_minimal()


#### MAT density plot of desert and Dryforet
x_limits <- c(15, 25)
y_limits <- c(0, 1)

# Create the density plots for Desert_data
Var1 <- ggplot(Desert_data, aes(x = MAT.Peru, color = ECO_COD)) +
  geom_density() +
  scale_color_manual(values = c("blue")) +
  labs(x = "Mean Annual Temperature", y = "Density", title = "Desert_data") +
  ylim(y_limits)+
  xlim(x_limits)

# Create the density plots for DryForest_data
Var2 <- ggplot(DryForest_data, aes(x = MAT.Peru, color = ECO_COD)) +
  geom_density() +
  scale_color_manual(values = c("light green")) +
  labs(x = "Mean Annual Temperature", y = "Density", title = "DryForest_data") +
  ylim(y_limits)+
  xlim(x_limits)

# Arrange the plots side by side
grid.arrange(Var1, Var2, ncol = 2)

#### MAP density plot of desert and Dryforet
x_limits <- c(0, 1000)
y_limits <- c(0, 0.075)

# Create the density plots for Desert_data
Var1 <- ggplot(Desert_data, aes(x = MAP.Peru, color = ECO_COD)) +
  geom_density() +
  scale_color_manual(values = c("blue")) +
  labs(x = "MAP", y = "Density", title = "Desert_data") +
  ylim(y_limits)+
  xlim(x_limits)

# Create the density plots for DryForest_data
Var2 <- ggplot(DryForest_data, aes(x = MAP.Peru, color = ECO_COD)) +
  geom_density() +
  scale_color_manual(values = c("light green")) +
  labs(x = "MAP", y = "Density", title = "DryForest_data") +
  ylim(y_limits)+
  xlim(x_limits)

# Arrange the plots side by side
grid.arrange(Var1, Var2, ncol = 2)

#### #### #### ####  that´s all!

