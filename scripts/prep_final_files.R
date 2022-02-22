### prep final files
# daphnia_output.txt
#head -1 /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/daphnia_output_xaa_20190201.txt > daphnia_output_all_20190624.txt; 
#tail -n +2 -q /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/daphnia_output_x*_20190201.txt >> daphnia_output_all_20190624.txt

# shape_output.txt
#head -1 /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/shape_output_xaa_20190201.txt > shape_output_all_20190624.txt; 
#tail -n +2 -q /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/shape_output_x*_20190201.txt >> shape_output_all_20190624.txt


# use output in R 

library(data.table)
### daphnia_output
daphnia_output_all <- fread('daphnia_output_all_20190624.txt')

daphnia_output_missing_a <- fread('daphnia_output_missing_20190524.txt')
daphnia_output_missing_b <- fread('daphnia_output_missing_NEW_20190524.txt')
daphnia_output_procrustes_outlier <- fread('daphnia_output_outlier_20190524.txt') 
daphnia_output_check_again <- fread('daphnia_output_check_again_20190524.txt')
daphnia_output_later <- fread('daphnia_output_later_NEW_20190524.txt')
daphnia_output_final_outlier <- fread('daphnia_output_final1_20190524.txt')

# all outliers
rm_filebase <- rbind(daphnia_output_missing_a, 
                     daphnia_output_missing_b,
                     daphnia_output_procrustes_outlier,
                     daphnia_output_check_again,
                     daphnia_output_later,
                     daphnia_output_final_outlier)

# rm outlier from data set
daphnia_output_all_removed <- daphnia_output_all[!filebase %in% rm_filebase$filebase]

# daphnia_output_final_outlier has 2 duplications 
issue_filebase <- c("full_110908_D8_663_juju_1D_RigB_20171011T124221.bmp","full_100850_AD8_14_juju3_1B_RigA_20170613T143455.bmp")
daphnia_output_final_outlier2 <- fread('daphnia_output_final2_20190524.txt')

# modified outlier data - to be used
missing_a <- daphnia_output_missing_a[!c(modified == 0 & modification_notes %like% 'broken')]
missing_b <- daphnia_output_missing_b[!c(modified == 0 & modification_notes %like% 'broken')][!filebase %in% missing_a$filebase]
procrustes_outlier <- daphnia_output_procrustes_outlier[!c(modified == 0 & modification_notes %like% 'broken')]
check_again <- daphnia_output_check_again[!c(modified == 0 & modification_notes %like% 'broken')][!filebase %in% procrustes_outlier$filebase]
later <- daphnia_output_later[modified == 1][!modification_notes %like% 'broken'][!filebase %in% procrustes_outlier$filebase]
final_outlier <- daphnia_output_final_outlier[!c(modified == 0 & modification_notes %like% 'broken')][!filebase %in% issue_filebase][!filebase %in% procrustes_outlier$filebase][!filebase %in% check_again$filebase][!filebase %in% later$filebase]
final_outlier2 <- daphnia_output_final_outlier2[!c(modified == 0 & modification_notes %like% 'broken')][!filebase %in% procrustes_outlier$filebase][!filebase %in% check_again$filebase][!filebase %in% later$filebase]   

daphnia_output_modified <- rbind(missing_a,
                                 missing_b,
                                 procrustes_outlier,
                                 check_again, 
                                 later, 
                                 final_outlier,
                                 final_outlier2)

# rm 'broken' in modification notes
daphnia_output_modified[,modification_notes := gsub("broken", "", modification_notes)]

daphnia_output_all_final <- rbind(daphnia_output_all_removed, daphnia_output_modified)
write.table(daphnia_output_all_final, 'daphnia_output_all_final_20190626.txt', row.names = F)


### shape_output
shape_output_all <- fread('shape_output_all_20190624.txt')

shape_output_missing_a <- fread('shape_output_missing_20190524.txt')
shape_output_missing_b <- fread('shape_output_missing_NEW_20190524.txt')
shape_output_procrustes_outlier <- fread('shape_output_outlier_20190524.txt') 
shape_output_check_again <- fread('shape_output_check_again_20190524.txt')
shape_output_later <- fread('shape_output_later_NEW_20190524.txt')
shape_output_final_outlier <- fread('shape_output_final1_20190524.txt')
shape_output_final_outlier2 <- fread('shape_output_final2_20190524.txt')

# rm outlier from data set
shape_output_all_removed <- shape_output_all[!filebase %in% rm_filebase$filebase]

# modified outlier data - to be used
shape_missing_a <- shape_output_missing_a[filebase %in% missing_a$filebase]
shape_missing_b <- shape_output_missing_b[filebase %in% missing_b$filebase]
shape_procrustes_outlier <- shape_output_procrustes_outlier[filebase %in% procrustes_outlier$filebase]
shape_check_again <- shape_output_check_again[filebase %in% check_again$filebase]
shape_later <- shape_output_later[filebase %in% later$filebase]
shape_final_outlier <- shape_output_final_outlier[filebase %in% final_outlier$filebase]
shape_final_outlier2 <- shape_output_final_outlier2[filebase %in% final_outlier2$filebase]

shape_output_modified <- rbind(shape_missing_a,
                               shape_missing_b,
                               shape_procrustes_outlier,
                               shape_check_again, 
                               shape_later, 
                               shape_final_outlier,
                               shape_final_outlier2)


shape_output_all_final <- rbind(shape_output_all_removed, shape_output_modified)
write.table(shape_output_all_final, 'shape_output_all_final_20190626.txt', row.names = F)




## corrected shape data
# use '*_output_all_final_20190626.txt' as input, rm filebases in shape and daphnia outputs, and replace data with modified data
# first, merge head only, head and tail, and tail only files

# files:
# daphnia_output_finalhead_20190712.txt
# daphnia_output_finalheadandtail_20190712.txt
# daphnia_output_finaltailONLY1_20190712.txt
# daphnia_output_finaltailONLY2_20190712.txt
# daphnia_output_finaltailONLY3_20190712.txt
# daphnia_output_finaltailONLY4_20190712.txt
# daphnia_output_finaltailONLY5_20190712.txt
# daphnia_output_finaltailONLY6_20190712.txt
# daphnia_output_finaltailONLY7_20190712.txt
# daphnia_output_finaltailONLY8_20190712.txt
# daphnia_output_finaltailONLY9_20190712.txt


#head -1 /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/daphnia_output_finalhead_20190712.txt > daphnia_output_adjshape_all_20190715.txt; 
#tail -n +2 -q /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/daphnia_output_*_20190712.txt >> daphnia_output_adjshape_all_20190715.txt


#head -1 /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/shape_output_finalhead_20190712.txt > shape_output_adjshape_all_20190715.txt; 
#tail -n +2 -q /mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/shape_output_*_20190712.txt >> shape_output_adjshape_all_20190715.txt



# use output in R 
library(data.table)
library(plyr)

### daphnia_output
daphnia_output_all <- fread('daphnia_output_all_final_20190626.txt')  # 7137 files

### modified data
daphnia_adjshape <- fread('daphnia_output_adjshape_all_20190715.txt')  # 563 files

# rm shape outlier from data set
daphnia_output_all_removed <- daphnia_output_all[!filebase %in% daphnia_adjshape$filebase]

# modified shape outlier data - to be used
daphnia_shape_outlier <- daphnia_adjshape[!modification_notes %like% 'exclude'][!modification_notes %like% 'no_checkpoints']


# merge modification notes from two data sets
daphnia_output_all_small <- daphnia_output_all[filebase %in% daphnia_shape_outlier$filebase]

setkey(daphnia_output_all_small, filebase)
setkey(daphnia_shape_outlier, filebase)

daphnia_output_modified <- daphnia_shape_outlier
daphnia_output_modified[, modification_notes := paste(daphnia_output_all_small$modification_notes, daphnia_shape_outlier$modification_notes, sep = ";")]

daphnia_output_all_final <- rbind.fill(daphnia_output_all_removed, daphnia_output_modified)
write.table(daphnia_output_all_final, 'daphnia_output_all_final_20190715.txt', row.names = F)



### daphnia_shape
shape_output_all <- fread('shape_output_all_final_20190626.txt')  # 7137 files

### modified data
shape_adjshape <- fread('shape_output_adjshape_all_20190715.txt')  # 563 files

# rm shape outlier from data set
shape_output_all_removed <- shape_output_all[!filebase %in% daphnia_adjshape$filebase]  

# modified shape outlier data - to be used
shape_shape_outlier <- shape_adjshape[filebase %in% daphnia_shape_outlier$filebase]  # 557 files

shape_output_all_final <- rbind(shape_output_all_removed, shape_shape_outlier)
write.table(shape_output_all_final, 'shape_output_all_final_20190715.txt', row.names = F)




