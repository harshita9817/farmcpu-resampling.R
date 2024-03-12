library(stringr)
library(dplyr)
library(ggplot2)
library(wesanderson)

# Set your working directory
# setwd("~/Documents/Sorghumseedimage/351_biallelic/GWAS_biallelic/")

# Get the file for the single trait (modify the pattern to match your file name)
  # files <- list.files(pattern = "Zbluessignals.csv")
  # 
  # # Read the data
  # df <- read.csv(files)
  # 
  # # Extract the trait name from the file name
  # trait <- str_extract(files, ".*(?=\\.csv$)")
  # 
  # # Remove the row names column
  # row.names(df) <- NULL
  # 
  # # Add the trait column to the dataframe
  # df$trait <- trait
  # 
  # # Subset the data for support greater than or equal to 0.1
  # df <- df[df$support >= 0.1, ]
  # 
  # # Write the processed data to a CSV file
  # write.csv(df, "FTsorghum.csv", row.names = FALSE)
  # 
  # # Plotting code
  # svg("FTsor_manual.svg", width = 10, height = 7)

df=read.csv("Zbluessignals.csv")
data_cum <- df %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHROM, bp_add)

gwas_data <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)
axis_set <- gwas_data %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum))

ggplot(gwas_data, aes(x = bp_cum, y = support)) + 
  geom_vline(xintercept = data_cum$bp_add, color = "darkgrey", linetype = "dashed", lwd=0.2)+
  geom_point(alpha = 5, size=3) +
  geom_hline(yintercept = 0.1, color = "black", linetype = "dashed", lwd=1)+
  scale_x_continuous(breaks = axis_set$center, n.breaks = 2000,
                     label=c("Chr01", "Chr02", "Chr03","Chr04", "Chr05", "Chr06","Chr07", "Chr08", "Chr09", "Chr10"))
