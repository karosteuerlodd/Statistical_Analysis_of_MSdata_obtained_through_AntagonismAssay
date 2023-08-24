
# Analyzing LC-MS/MS data between single cells and interaction

install.packages("tidyverse")
install.packages("factoextra")
install.packages("svglite") 
install.packages("ggsci")
install.packages("cowplot")
install.packages("FSA")
install.packages("vegan")
install.packages("matrixStats")
install.packages("BiocManager")

library(tidyverse)
library(factoextra)
library(svglite)
library(vegan)
library (ggsci)
library(cowplot)
library(FSA)
library(matrixStats)
library(BiocManager)

#set directory
setwd("")

# Get the feature and metadata table -----------------------------------------------------------------------
file_names <- list.files('.') #list all the files in the working directory (mentioned by 'dot symbol')
file_names

print(file_names)
# change the respective numbers after seeing the headers!!!
ft <- read.csv(file_names[#"column number"], header = T, check.names = F, sep = ';') # ; for german excel as csv.file
md <- read.csv(file_names[#"column number"], header = T, check.names = F, sep = '\t') #seperator is mentioned in case of txt or tsv files
an <- read.csv(file_names[#"column number"], header = T, check.names = F, sep = '\t') 

head(ft) #returns the first 6 rows as default
dim(ft) 

head(md) #returns the first 6 rows as default
dim(md) 

head(an) #returns the first 6 rows as default
dim(an) 

# Exploring metadata
InsideLevels <- function(metatable){
  LEVELS <- c() #creating empty vector to store information
  typ <-c()
  COUNT <- c()
  for(i in 1:ncol(metatable)){ # for each metadata column
    temp <- as.data.frame(table(metatable[,i])) #table function gives the category in each column and the count of each category
    x <- temp$Var1 #getting the name of each category in every column
    if(is.double(metatable[,i])==T){ # for numeric columns in metadata table, round the category values
      x=round(as.double(x),2)
    } 
    LEVELS <- rbind(LEVELS,toString(x)) # adding all the category values in a row
    COUNT <- rbind(COUNT,toString(temp$Freq)) # getting the frequency of each level in every column
    typ <- rbind(typ,class(metatable[,i])) # getting the class of each column
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)), #creating an output dataframe with 1st column as INDEX
                    ATTRIBUTES=colnames(metatable), #2nd column ATTRIBUTES will be the column name of metadata table
                    LEVELS, #3rd column LEVELS will give the different categories in each ATTRIBUTE
                    COUNT, #4th column COUNT will give the number of files present with each category
                    'ATTRIBUTE_CLASS'=typ,row.names=NULL) #Final column indicating the Class or datatype of each ATTTRIBUTE
  return(out)
}

InsideLevels(md[,2:5])

#Arranging metadata and feature table---------------------------------------------------------

new_ft <- ft #storing the files under different names to preserve the original files
new_md <- md

colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft)) #removing Peak area extensions from the column names of ft

new_ft <- new_ft[order(new_ft$`row ID`),,drop=F] #arranging the rows of ft file by  by ascending order of row ID
new_ft <- new_ft[,colSums(is.na(new_ft))<nrow(new_ft)] #removing if any NA columns present in the ft file
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #removing if any NA columns present in the md file

#Add annotation by GNPS
identical(class(ft$`row ID`),class(an$`#Scan#`)) #check whether class types are the same
ft_an <- merge(ft, an, by.x="row ID", by.y="#Scan#", all.x= TRUE) #merging ft and an by the respective column names; all.x=T refers to keeping all rows of 'ft'
dim(ft_an) # checking the dimension: here, ft_an has 14709 rows and 194 cols:(149+45)

#remove the (front & tail) spaces, if any present, from the filenames of md
new_md$filename <- trimws(new_md$filename, which = c("both")) 

#should return TRUE if you have annotation file
if(exists("ft_an")){identical(ft_an$`row ID`,new_ft$`row ID`)} 

#Changing the row names of the files into the combined name as "XID_mz_RT":
rownames(new_ft) <- paste(paste0("X",new_ft$`row ID`),
                          round(new_ft$`row m/z`,digits = 3),
                          round(new_ft$`row retention time`,digits = 3),
                          if(exists("ft_an")){ft_an$Compound_Name}, 
                          sep = '_') 

rownames(new_ft) <- sub("_$", "", rownames(new_ft)) #to remove the trailing underscore at rownames

#Picking only the files with column names containing 'mzML'.
new_ft <- new_ft[,grep('\\.mzML$',colnames(new_ft))] 

new_ft <- new_ft[,order(colnames(new_ft)),drop=F] #ordering the ft by its column names
new_md <- new_md[order(new_md$filename),,drop=F] #ordering the md by the 1st column filename

# how many files in the metadata are also present in the feature table
table(new_md$filename %in% colnames(new_ft))

identical(new_md$filename, colnames(new_ft))

# which file names in the metadata are not in the feature table?
#setdiff(new_md$filename,colnames(new_ft))

print(colnames(new_ft))

#checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

#Transpose new_ft to be combined with new_md
ft_t <- as.data.frame(t(new_ft)) #transposing the imputed table
ft_t <- ft_t[match(new_md$filename,rownames(ft_t)),] # ordering the rows in the Imp_t table according to the order of fiilenames column of md_Samples
identical(rownames(ft_t),new_md$filename) #should return TRUE

#Merging ft with md
Data <- merge(new_md,ft_t,by.x="filename",by.y="row.names") #merging metadata with sclaed data
write.csv(Data, paste0(Sys.Date(),"_ft_metadata_merged.csv"),row.names =T) #writing the file as csv

############-----------------------------------------------------------------------------Blank removal
# Blank Removal
#---------------------
#split ft into blanks and samples using metadata
InsideLevels(new_md[2:ncol(new_md)])

#separate blank files
md_Blank <- new_md %>% filter(`ATTRIBUTE_SampleType` == "blank") #filtering the rows from metadata with the condition = blank
Blank <- ft_t[which(rownames(ft_t) %in% (md_Blank$`filename`)),,drop=F] #getting the corresponding rows from ft_t

head(Blank,n=2)
dim(Blank) 

#separtate sample files
md_Samples <- new_md %>% filter(`ATTRIBUTE_SampleType` == "sample") #filtering the rows from metadata with the condition = sample
Samples <- ft_t[which(rownames(ft_t) %in% (md_Samples$`filename`)),,drop=F] #getting the corresponding rows from ft_t

head(Samples, n=2)
dim(Samples)

####define cutoff, 30% advisable
#When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3

#perform blank removal
#Getting mean for every feature in blank and Samples in a data frame named 'Avg_ft'
Avg_ft <- data.frame(Avg_blank=colMeans(Blank, na.rm= F)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ft$`Avg_samples` <- colMeans(Samples, na.rm= F) # adding another column 'Avg_samples' for feature means of samples

#Getting the ratio of blank vs Sample
Avg_ft$`Ratio_blank_Sample` <- (Avg_ft$`Avg_blank`+1)/(Avg_ft$`Avg_samples`+1)

# Creating a bin with 1s when the ratio>Cutoff, else put 0s
Avg_ft$`Bg_bin` <- ifelse(Avg_ft$`Ratio_blank_Sample` > Cutoff, 1, 0 )

#Calculating the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ft)))
print(paste("No.of Background or noise features:",sum(Avg_ft$`Bg_bin` ==1,na.rm = T)))
print(paste("No.of features after excluding noise:",(ncol(Samples) - sum(Avg_ft$`Bg_bin` ==1,na.rm = T))))

blk_rem <- merge(as.data.frame(t(Samples)), Avg_ft, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) %>% #removing the last 4 columns
  column_to_rownames(var="Row.names")

#################### review blank removed table
head(blk_rem, 2)
dim(blk_rem)
        # metadata without the blanks info
head(md_Samples, 2)
dim(md_Samples)

#write a new file with blank removed features
write.csv(blk_rem, paste0(Sys.Date(),'_Blanks_Removed_with_cutoff_',Cutoff,'.csv'),row.names =TRUE)

###############bring "blk_rem" file into same structre as "Data" ffor further processing

        #Transpose blk_rem to be combined with blk_t
blk_t <- as.data.frame(t(blk_rem)) # transposing the imputed table
blk_t_new <- blk_t[match(md_Samples$filename,rownames(blk_t)),] # ordering the rows in the Imp_t table according to the order of filenames column of md_Samples
identical(rownames(blk_t_new),md_Samples$filename) #should return TRUE

#Merging ft with md
Data_blkrem <- merge(md_Samples, blk_t_new, by.x="filename", by.y="row.names") #merging metadata with scaled data
write.csv(Data_blkrem, paste0(Sys.Date(),"_ft_metadata_blankremoved_merged.csv"),row.names =T) #writing the file as csv

Data <- Data_blkrem

####################################################-----------------------
# Statistical Analysis of single cultures against the tested interaction
#--------------------------------------------------------------------------------------------------------------------
#grab all the rows from interaction/ cut "n" and 
#put them into folder n with the single analysis from propagator and indicator (LB) 
#grab each interaction and single culture and create a table for each

names <- unique(Data$ATTRIBUTE_Group) # list all the characters which are in the column
names

# cut names are defined as thinteractions between propagator-strain a and indicator-strain b
cut_names <- names[1:18] # group all the cuts into one vector
cut_names
length(cut_names) #show how many characters are inside


#filter the dataframe 'Data' for the selected cut names in the respective column:---------------

# Now let's create a loop to analyze all the interactions---------------------------------

interaction <- list() # create an empty list

for (i in 1:length(cut_names)){ #from 1 to length of interaction list
  
  interaction[[i]] <- Data %>% filter(ATTRIBUTE_Group == cut_names[i])
  
  bac_prop <- interaction[[i]]$ATTRIBUTE_Propagator[1]
  bac_ind <- interaction[[i]]$ATTRIBUTE_Indicator[1]
  
  add1 <- Data %>%
    filter(ATTRIBUTE_Group == paste(bac_prop, "_LB", sep="")|
             ATTRIBUTE_Group == paste(bac_ind, "_LB",sep = ""))
  
  interaction[[i]] <- rbind(interaction[[i]], add1)

}

names(interaction) <- cut_names

####################################################################
#grouping rows where stat. analysis should be performed
row_int <- c(1:3)
row_propcell <- c(4:5)
row_prosus <- c(6)
row_indcell <- c(7:8)
row_indsus <- c(9)

colstart_calc <- 8 #(7 metadata columns + 1)
colend_calc <- ncol(interaction[[1]]) # no.of columns of interaction element

collength_calculation <- c(colstart_calc:colend_calc) #peak area columns

#-----------------------------------------------------------------------------
#creating loop for calculating mean of every interactions and single cultures

# creating empty lists

log2_intsingle_list <- list()

for (g in 1:length(interaction)) { #length(interaction)
  
  #getting the corresponding rows:
  int_rows <- interaction[[g]][row_int, collength_calculation]
  propcell_rows <- interaction[[g]][row_propcell, collength_calculation]
  indcell_rows <- interaction[[g]][row_indcell, collength_calculation]
  
  #calculation of means of the grouped rows
  intmean <- colMeans(int_rows)
  propcellmean <- colMeans(propcell_rows) #meancalc. between row 4:5 because they are the replicates
  indcellmean <- colMeans(indcell_rows)
  
  #calculation of standard deviation of the grouped rows
  intsd <- apply(int_rows, 2, sd)
  propcellsd <- apply(propcell_rows, 2, sd)
  indcellsd <- apply(indcell_rows, 2, sd)
 
  #combine these individual dataframes into temporary dataframe:
  temp_meanSD <- as.data.frame(rbind(intmean, intsd, propcellmean, propcellsd, indcellmean, indcellsd))
  rownames(temp_meanSD) <- c('intmean','intsd', 'propcellmean', 'propcellsd', 'indcellmean','indcellsd')
  
####################################################################################
  #removed imputation here
  
  #calculation of quotient between propagator/indicator and interaction 
  quotcell_propint <- (temp_meanSD['propcellmean',]) / #row 3 single propagator 
                      (temp_meanSD['intmean',])   #and row 1 interaction zone
  
  quotsus_propint <- (interaction[[g]][6,-c(1:7)])/ #Without first 7 columns
                     (temp_meanSD['intmean',]) 

  quotcell_indint <- (temp_meanSD['indcellmean',])/ #row 5 single indicator
                     (temp_meanSD['intmean',])  
  
  quotsus_indint <- (interaction[[g]][9,-c(1:7)]) / 
                    (temp_meanSD['intmean',])

  #calculating log2 of means for x-axis
  pos1 <- log(quotcell_propint, base=2)
  pos2 <- log(quotsus_propint, base=2)
  pos3 <- log(quotcell_indint, base=2)
  pos4 <- log(quotsus_indint, base=2)
  #print(length(pos4))
 
  #putting the calculated log2 values in a dataframe and list  
  log2_intsingle_list[[g]] <- as.data.frame(rbind(pos1, pos2, pos3, pos4)) 
  rownames(log2_intsingle_list[[g]]) <- c('log2_propagator_interaction_cell', 
                                          'log2_propagator_interaction_suspension', 
                                          'log2_indicator_interaction_cell', 
                                          'log2_indicator_interaction_suspension')
  
  #get the name of the respective interaction/ cut_name
  names(log2_intsingle_list)[g] <- names(interaction)[g] 
  print(g)
}

###########################################################################
#doing ttest analysis for the each interaction list in a for loop---------------------------------------------------------

propcell_list <- list() #create an empty list to put all the t-test dataframes in it (one dataframe per interaction )
Indcell_list  <- list()

for (k in 1:length(interaction)){   #from one to n number of interactions
  
  int_short  <- interaction[[k]] #each current interaction-dataframe is taken under int-short
  rownames(int_short) <- int_short$filename #naming the rows with the first column
  int_short  <- int_short[,-c(1:7),drop=F] #crop the dataframe with just relevant information
  print(ncol(int_short)) #print the amount of columns or no.of features
  
  orga_propcell <- data.frame() #create empty dataframes to store t-test results
  orga_indcell <- data.frame()
  
  print(k)
  #new for-loop to perform t-test for every feature (column) with the respective replicates/groups
  for (t in 1:ncol(int_short)) { 
    
    # grouping: Extract the values for the current feature and the respective group
    ft_int <- int_short[row_int,t]
    
    ft_propcell <- int_short[row_propcell,t]
    #ft_propsus <- int_short[row_prosus,t]
    
    ft_indcell <- int_short[row_indcell,t]
    #ft_indsus <- int_short[row_indsus,t]
    
  #perform t.test on the selected rows of propagator----
    model1 <- broom::tidy(t.test(ft_int, ft_propcell))
    orga_propcell <- bind_rows(orga_propcell, model1) #combine outputs from all features
    
  #perform t.test on the selected rows of indicator------
    model2 <- broom::tidy(t.test(ft_int, ft_indcell))
    orga_indcell <- bind_rows(orga_indcell, model2) #combine outputs from all features
    
  }
  rownames(orga_propcell) <- colnames(int_short) #naming the rows of these t-test dataframes with colnames of int-short
  rownames(orga_indcell) <- colnames(int_short)
  
propcell_list[[k]] <- orga_propcell #put them as a list element
Indcell_list[[k]] <- orga_indcell
}

#naming the new list elements as same as interaction list
names(propcell_list) <- names(interaction)
names(Indcell_list) <- names(interaction)

##################################
# creating volcano plots for each interaction in propcell list:-------------------------------------------------
propcell_plots <- list()
for (v in 1:length(propcell_list)){
  
  my_propcell <- propcell_list[[v]]
  temp_log <- data.frame(t(log2_intsingle_list[[v]])) #get the log2 values as temp_log
  identical(rownames(orga_propcell),rownames(temp_log)) #should return TRUE
  
  #add the log2 values of propagator-interaction as another column in my_propcell
  my_propcell$log2value <- temp_log$log2_propagator_interaction_cell 
  
  # If you want more values to be considered as 'significant', 
  # allow for more features to pass through here: p_bonferroni<0.05
  my_propcell <- arrange(my_propcell, p.value) # arranging my_propcell by p values
  my_propcell["p_bonferroni"] <- p.adjust(my_propcell$p.value,method="bonferroni") #perform bonferroni correction on the p values
  my_propcell["significant"] <- ifelse(my_propcell$p_bonferroni<0.05,"Significant","Nonsignificant") #assign signifance for the bonferroni-corected p value
  my_propcell <- filter(my_propcell, significant!="NA") #remove the features with NA values in significant column
  
# Replace infinite values in log2value with NAs and add it as a column 'log2value' to my_propcell
  my_propcell$log2value <- replace(my_propcell$log2value, 
                                   is.infinite(my_propcell$log2value), 
                                   NA)
  
  propcell_plots[[v]] <- ggplot(my_propcell,
                   aes(x = log2value, #choosing log2value column (Quotient) as x axis
                       y = -log(p.value,base=10), #here, I do -log10 on the original p values
                       color = significant)) + # reminder: significance is assigned on bonferroni corrected p-values
    geom_point() +
    theme_classic() +
    scale_color_manual(values=c("black", "firebrick")) +
    geom_vline(xintercept = 0, linetype="dashed") +
    xlab("log2(Quotient)") +
    ylab("-log(p)") +
    theme(legend.title = element_blank()) +
    ggrepel::geom_text_repel(data = head(my_propcell), 
                             aes(label = rownames(my_propcell)[1:6]), 
                             size=3, 
                             show.legend = F, max.overlaps = 100) +
    labs(title = "Volcano Plot of Interaction vs. Single-culture Propagator", x = "Log2(Quotient)", y = "-Log10 (p-value)")
  
  temp_name <- names(propcell_list)[v]
  ggsave(paste("volcano_plot_propagator_",temp_name,".pdf"), propcell_plots[[v]])
}

################################---------------------------
##### vulcano plots for indicator strains -----------------------------------------------
indcell_plots <- list()
for (v in 1:length(Indcell_list)){
  
  my_indcell <- Indcell_list[[v]]
  temp_log <- data.frame(t(log2_intsingle_list[[v]])) #get the log2 values as temp_log
  identical(rownames(orga_indcell),rownames(temp_log)) #should return TRUE
  
  #add the log2 values of indicator-interaction as another column in my_indcell
  my_indcell$log2value <- temp_log$log2_indicator_interaction_cell 
  
  my_indcell <- arrange(my_indcell, p.value) # arranging my_propcell by p values
  my_indcell["p_bonferroni"] <- p.adjust(my_indcell$p.value,method="bonferroni") #perform bonferroni correction on the p values
  my_indcell["significant"] <- ifelse(my_indcell$p_bonferroni<0.05,"Significant","Nonsignificant") #assign signifance for the bonferroni-corected p value
  my_indcell <- filter(my_indcell, significant!="NA") #remove the features with NA values in significant column
  
  # Replace infinite values in log2value with NAs and add it as a column 'log2value' to my_indcell
  my_indcell$log2value <- replace(my_indcell$log2value, 
                                   is.infinite(my_indcell$log2value), 
                                   NA)
  
  indcell_plots[[v]] <- ggplot(my_indcell,
                                aes(x = log2value, #choosing log2value column (Quotient) as x axis
                                    y = -log(p.value,base=10), #here, I do -log10 on the original p values
                                    color = significant)) + # reminder: significance is assigned on bonferroni corrected p-values
    geom_point() +
    theme_classic() +
    scale_color_manual(values=c("lightgrey", "indianred")) +
    geom_vline(xintercept = 0, linetype="dashed") +
    xlab("log2(Quotient)") +
    ylab("-log(p)") +
    theme(legend.title = element_blank()) +
    ggrepel::geom_text_repel(data = head(my_indcell), 
                             aes(label = rownames(my_indcell)[1:6]), 
                             size=3, 
                             show.legend = F, max.overlaps = 100) +
    labs(title = "Volcano Plot of Interaction vs. Single-culture Indicator", x = "Log2(Quotient)", y = "-Log10 (p-value)")
  
  temp_name <- names(Indcell_list)[v]
  ggsave(paste("volcano_plot_indicator_",temp_name,".pdf"), indcell_plots[[v]])
}


# if didn't give dev.off here, because we store each plot as a list element and save that as a pdf. 
# It worked well without dev.off

