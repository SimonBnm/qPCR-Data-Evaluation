
########################################################################### 
#### Set parameter ####

# Set the housekeeper gene
housekeeper <- "RP-49"

# Set Control name(s)
cont_name <- "Cherry|Chery"

# Set foldername of data
foldername <- "Example data"


# Set number of target genes (primers) including the housekeeper
target_number <- 13

# Set Control name
con_name <- "46F-Cherry"
exp_name <- "46F-Drgx"

# leave out target genes
target_remove <- c("")

# Set stat. methode
stat_method <- "wilcox.test"

# Adjust the names displayed in the plot (alphabetical order)

rename <- c("AdamTS-B","bbg","CG4476","CG5687","Cht2","daw","Indy","Mdr49","Mthl10","Orct2","Slow","stl")


###########################################################################
# Check wd
getwd()

#### Packages ####

library("dplyr")
library("writexl")
library("ggplot2")
library("ggpubr")
library("Hmisc")
library("openxlsx")
library("readxl")
library ("svglite")

#### Functions ####

# unnest dataframes

unnest_dataframes <- function(x) {
  
  y <- do.call(data.frame, x)
  
  if("data.frame" %in% sapply(y, class)) unnest_dataframes(y)
  
  y
  
}
############################################################################
#### Script ####


# Get a list of all files in the folder
file_list <- list.files(foldername, pattern = ".xlsx", full.names = TRUE) 
file_list

# Filter out files that begin with '#'
filtered_files <- file_list[!grepl("^#", basename(file_list)) &
                              !file.info(file_list)$isdir]
filtered_files

# Initialize an empty data frame to store the results
agg_data_all <- data.frame()
agg_data_all

# Loop through each file
for (i in seq_along(filtered_files)) {
  
  i = i
  # Read in the data
  data <- as.data.frame(read_excel(filtered_files[i]))
  #data <- read.csv(file_list[i]) for csv files
  #data <- read.xlsx(filtered_files[i], sheetName  = 1, header=TRUE)
  data


  # Add Column for Group (Experimental or Control)
  #data$Group<- ifelse(grepl("Drgx", data$Sample), "Experimental", data$Sample)
  #data$Group <- ifelse(grepl("Cherry|Chery", data$Sample), "Control", data$Sample)
  
  data$Group <- ifelse(grepl(cont_name, data$Sample), "Control", "Experimental")
  
  # Aggregate the Cq values
  agg_data <- aggregate(Cq ~ Target + Group + Sample, data = data, FUN = function(x) sort(x))
  # Calculate the mean and standard deviation of the Cq values by Target and Sample
  Cq_stats <- aggregate(Cq ~ Target + Group + Sample, data = data, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  # Merge the aggregated data with the Cq stats
  agg_data <- merge(agg_data, Cq_stats, by = c("Target", "Group", "Sample"))
  # Unnest data
  agg_data2 <- unnest_dataframes(agg_data)
  ## Extract Cq mean values for housekeeper
  # subset the dataframe for the Target "RP-49"
  df_housekeeper <- subset(agg_data2, Target == housekeeper)
  # extract the mean_cq value for the Group "Control"
  housekeeper_cq_control <- df_housekeeper[df_housekeeper$Group == "Control", "Cq.y.mean"]
  # extract the mean_cq value for the Group"Experimental"
  housekeeper_cq_experimental <- df_housekeeper[df_housekeeper$Group == "Experimental", "Cq.y.mean"]
  ### Calculate delta Cq values (∆Cq = Cq (gene of interest) – Cq (housekeeping gene))
  agg_data2$delta_cq <- ifelse(agg_data2$Group == 
  "Control", agg_data2$Cq.y.mean - housekeeper_cq_control, agg_data2$Cq.y.mean - housekeeper_cq_experimental)
  
  ### Create new df with only important columns
  #cut_df <- agg_data2[, c("Target","Sample", "Cq.y.mean", "Cq.y.sd", "delta_cq")]
  cut_df <- agg_data2[, c("Target","Group","Sample","Cq.y.mean", "delta_cq")]
  # Reorder data according to Sample name
  cut_df <- arrange(cut_df, Group)
  # Bind the results to the overall data frame
  agg_data_all <- rbind(agg_data_all, cut_df)
}

# View the final results
agg_data_all

# Sort data by Group and Target name
delta_Cq_data <- arrange(agg_data_all, Group, Target)
delta_Cq_data

## Calculate the mean delta Cq (∆Cq) values for each Target gene
mean_d_cq_data <- aggregate(delta_cq ~ Target + Group, data = delta_Cq_data, FUN = function(x) mean = mean(x))
mean_d_cq_data

# Merge data
delta_Cq_data<- merge(delta_Cq_data, mean_d_cq_data, by = c("Target", "Group"))
delta_Cq_data


## Get dataframe with Delta Cq control average values (∆Ct Control average)

# subset the data for the "Control" Group
Control_data <- subset(mean_d_cq_data, Group == "Control")

# Change colname
colnames(Control_data)[colnames(Control_data) == "delta_cq"] = "mean_control_cq"
Control_data

# Select only the columns "Target" and "delta_cq.y"
Control_delta_cq <- Control_data[, c("Target", "mean_control_cq")]
Control_delta_cq

# Merge data
data_cont_mean<- merge(delta_Cq_data,Control_delta_cq, by = "Target")
data_cont_mean


# Delete delta cq mean
data_cont_mean$delta_cq.y <- NULL
data_cont_mean

# Calculate delta delta Cq (∆∆Cq)
data_cont_mean$dd_cq <- data_cont_mean$delta_cq.x - data_cont_mean$mean_control_cq
data_cont_mean


# Calculate fold expression (Fold gene expression = 2^-(∆∆Cq))
data_cont_mean$fold_expr <- 2^(-data_cont_mean$dd_cq)
data_cont_mean

# Calculate mean fold expression (mean (2^-(∆∆Cq))and SD for each target gene
mean_d_cq_data <- aggregate(fold_expr ~ Target + Group, data = data_cont_mean,
                            FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_d_cq_data

# Merge data
final_data1<- merge(data_cont_mean, mean_d_cq_data, by = c("Target", "Group"))

# Unnest data
final_data1 <- unnest_dataframes(final_data1)
final_data1

# rename data
colnames(final_data1) <- c("Target","Group", "Replicate","Mean_Cq","∆Cq",
                           "mean_∆Cq_Control","∆∆Cq","2^-(∆∆Cq)","Mean_Fold_Expression","SE_Mean_Fold")
fin_data <- arrange(final_data1, Target)


fin_data   # Data for analysis

###############################################################################
#### Quality check ####
###############################################################################

# Check Target gene (Primer) names
unique_names <- unique(fin_data$Target)
num_unique_names <- length(unique_names)
num_unique_names # number of genes
unique_names     # names of genes


# Get the Replikate number for each Target
grouped_test <- aggregate(fin_data$Replicate, by = list(fin_data$Target), FUN = paste)

# Rename the columns
colnames(grouped_test) <- c("Target", "Replikate")

grouped_test

#### Stop code and print error message if there is something wrong with the data ####

# Check if 'test' has a specific value
if (target_number!=num_unique_names) {
  print(grouped_test)
  # Print an error message
  print("Error: Wrong number of target genes")
  
  # Stop the code
  stop() } 

# Leave out specific Targets for further evaluation

filtered_df <- fin_data[!fin_data$Target %in% target_remove, ]



##############################################################################
#### Statistical analysis####
##############################################################################


## Clean up data

# Delete housekeeper data
new_df <- subset(filtered_df, Target != housekeeper)
new_df

# New df
stat_data <-new_df
stat_data

# Rename the column 
colnames(stat_data )[8] <- "Fold_Expr"
stat_data 

# Cut data
stat_cut  <- stat_data[, c("Target", "Group","Replicate","Fold_Expr")]
stat_cut

## Compare groups (Check if significant differences between Control and Exp for all Targets)

# Reorder and clean data
stat_cut


# Perform test
compare_means(Fold_Expr ~ Group, data = stat_cut, paired = FALSE,method = stat_method,group.by = "Target")

# Save Test results as dataframe
stats <- compare_means(Fold_Expr ~ Group, data = stat_cut, paired = FALSE,method = stat_method,group.by = "Target")
stats2 <- as.data.frame(stats)
stats2

###############################################################################
#### Generate plot ####
###############################################################################

### Find highest values with sd for significant targets
stat_data

# Cut data only Mean Fold Expression and Std 
df1 <- stat_data[, c("Target","Group","Mean_Fold_Expression", "SE_Mean_Fold")]

# Remove duplicate rows
df2 <- distinct(df1, Target, Group, .keep_all = TRUE)

# Create a new column that contains the sum of Mean_Fold_Expression and Std_Mean_Fold
df2$sum <- df2$Mean_Fold_Expression + df2$SE_Mean_Fold

## Change to same order as in plot (alphabetically)
stats2_sorted<-stats2[order(stats2$Target),]

# Order data to get x_pos of each Target
stats_ranked <- stats2_sorted %>% mutate(x_pos = row_number())
stats_ranked

# Filter only signifikant Targets
filt_stats <- stats_ranked %>% filter(p.signif != "ns")

# Remove non significant targets
df2_subset <- df2[df2$Target %in% filt_stats$Target, ]

# Find the row index with the highest sum value
index <- which.max(df2_subset$sum)

# Dubset the row with the highest sum value
highest_row <- df2_subset[index, ]

# Extract x value for line position
ypos <- highest_row$sum 

#set y value for sig stars and lines in plot
yp <- ypos+0.2


### Plot with significance stars

rename <- c("AdamTS-B","bbg","CG4476","CG5687","Cht2","daw","Indy","Mdr49","Mthl10","Orct2","Slow","stl")


final_plot1 <- ggplot(stat_cut, aes(x = Target, y = Fold_Expr, fill = Group)) +
  
       # Add barplot with mean Fold Expression
       stat_summary(fun = mean, geom = "bar", position = "dodge",colour="black", linewidth=0.8, width=0.8, alpha=0.8)+
       # Add errorbars with standard deviation
       stat_summary(fun.data = "mean_se", geom = "errorbar", 
               position = position_dodge(width = 0.8), width = 0.2) +
      
       # Further adjustements of the plot
       scale_fill_manual(values=c("grey10","grey80"), labels=c(con_name ,exp_name))  +
       scale_y_continuous(expand = c(0, 0))+
       expand_limits(y = c(0,yp))+ 
       theme_bw()+
       theme(axis.text.x = element_text(angle = 45, hjust = 1))+
       labs(x="", y="Relative mRNA expression (fold change)")+
       theme(axis.text=element_text(size=20),
             axis.title=element_text(size=22,face="bold"),
             legend.title = element_blank(),
             axis.text.x = element_text(colour = "black", size=20),
             legend.text = element_text(size=20))+
      #rename
      scale_x_discrete(labels = rename)
  

#final_plot1 




### Ad lines and significance stars

#Add lines for comparison between Groups of same Target Gene (only if significant p vales)

## Change to same order as in plot (alphabetically)
stats2_sorted<-stats2[order(stats2$Target),]

# Order data to get x_pos of each Target
stats_ranked <- stats2_sorted %>% mutate(x_pos = row_number())
stats_ranked


# Add horizontal lines using the subsetted data frame

#set distance to error bar
space <- 0.03

final_plot2 <- final_plot1 + geom_segment(data = filt_stats,
                    aes(x = x_pos - 0.3,
                        xend = x_pos + 0.3,
                        y = ypos+space, yend = ypos+space),
                    linewidth = 1,
                    inherit.aes = FALSE) +
  
                # Add Significance stars manually
              geom_text(data = filt_stats,aes(x = x_pos, y = ypos+ space+ 0.01, 
                                              label = p.signif),
              size = 9, inherit.aes = FALSE, check_overlap = T)


final_plot2
###############################################################################
#### Save Data: Save Evaluation Data and GGPlot ####
###############################################################################

### generate Evaluation folder ###

# Set name of Evaluation folder
evalfolder <- paste(foldername,"Evaluation")

# Generate folder
dir.create(file.path(dirname(foldername),evalfolder))

file.path(foldername)
#### Generate Excel file with statistical evaluations ####

# Set name of test
stat_test_name <- stats2$method[1]

# Create Excel file and add data frames as sheets
wb <- createWorkbook()

# Add Delta Delta cq calculation data
addWorksheet(wb, "∆∆Cq method calculations")
writeData(wb, sheet = 1, fin_data, startRow = 1, startCol = 1)

# Add Statistical evaulation data
addWorksheet(wb, stat_test_name)
writeData(wb, sheet = 2, stats2, startRow = 1, startCol = 1)

# Set name of file
evalname  <- paste(evalfolder,"/",foldername,"_Evaluation.xlsx",sep="")

# Save Excel file
saveWorkbook(wb,evalname, overwrite = TRUE)


### Save GGPLOT ###

#set plot name
plotname <- paste0(foldername, "_Relative mRNA expression plot - PDF.pdf")

#save plot as pdf
ggsave(filename=plotname,
       path= evalfolder,
       plot = final_plot2, 
       device = cairo_pdf, 
       width = 297, 
       height = 210, 
       units = "mm")


#save plot as .tif
plotnametif <- paste0(foldername, "_Relative mRNA expression plot - tiff.tiff")

ggsave(filename=plotnametif,
       path= evalfolder,
       plot = final_plot2, 
       device = "tiff", 
       width = 297, 
       height = 210, 
       units = "mm")

#save plot as .svg
plotnametif <- paste0(foldername, "_Relative mRNA expression plot - svg.svg")

ggsave(filename=plotnametif,
       path= evalfolder,
       plot = final_plot2, 
       device = "svg", 
       width = 297, 
       height = 210, 
       units = "mm")

