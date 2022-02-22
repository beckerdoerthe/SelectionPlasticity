rm(list=ls()) 

#################
### libraries ###
#################
    library(ggplot2)
    library(cowplot)
    library(ggbeeswarm)
    library(gridExtra)
    library(data.table)
    library(mclust)
    library(geomorph)
    library(foreach)
    library(abind)
    library(doMC)
    registerDoMC(20)
    library(lme4)
    #library(lmerTest)
    library(viridis)
    library(dplyr)
    library(lattice)
    library(plotrix)
    library(broom)



#################
### data prep ###
#################
    ### NOTE: SEASON > 'spring_1' in '2017' == April2017; 'spring_2' in '2017' == May2017
    ### functions
        load_output_data <- function(dir="/mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/",
                         dateStamp="20190715",
                         trimColumns=T) {

                            ### load & parse summary data
                                daphnia_output <- fread(paste(dir, "daphnia_output_all_final_", dateStamp, ".txt", sep=""), header=T)

                            ### fix problematic inductiondate
                                daphnia_output[,inductiondate.orig := inductiondate]
                                daphnia_output[,inductiondate := gsub("2018", "2017", inductiondate)]
                             
                            ### flag summary data
                                daphnia_output[!modification_notes%like%'broken' &
                                !modification_notes%like%'deformed' & 
                                !modification_notes%like%'tilted' &
                                !modification_notes%like%'male' &
                                !daphnia_output$manual_PF_reason %like%'male' &
                                !daphnia_output$manual_PF_reason %like%'duplicate' & 
                                !daphnia_output$cloneid %like%'Male',
                                use.mod:=T]
                               daphnia_output[is.na(use.mod), use.mod:=F]

                           # get day (days)
                               getDay <- function(induction, picture) {

                                 ind <- strptime(substr(induction, 0, 8), format="%Y%m%d")
                                 pic <- strptime(substr(picture, 0, 8), format="%Y%m%d")

                                 as.numeric(pic - ind)
                               }

                               daphnia_output[, day := getDay(induction=inductiondate, picture=datetime)]
                               daphnia_output[, first_day :=  min(day), by = barcode]  # get first day of pic for each barcode
                               daphnia_output[, age := 1 + (as.numeric(day) - first_day)]  # add extra col for AGE


                           # get hours - CHECK THIS AGAIN
                               getHours <- function(induction, picture) {

                                 ind <- strptime(substr(induction, 0,15), format="%Y%m%dT%H%M%S")
                                 pic <- strptime(substr(picture, 0, 15), format="%Y%m%dT%H%M%S")

                                 as.numeric(pic - ind)
                               }

                               daphnia_output[, hours := getHours(induction=inductiondate, pic=datetime)*24]
                               daphnia_output[, first_hours :=  min(hours), by = barcode]  # get first day of pic for each barcode
                               daphnia_output[, age_hours := 1 + (as.numeric(hours) - first_hours)]  # add extra col for AGE

                           # modify data (pixel to mm conversion)
                               daphnia_output[, pixel_to_mm := as.numeric(as.character(pixel_to_mm))]

                               daphnia_output[, animal_dorsal_area := as.numeric(as.character(animal_dorsal_area))]
                               daphnia_output[, animal_dorsal_area_mm := animal_dorsal_area/(pixel_to_mm^2)]

                               daphnia_output[, eye_area := as.numeric(as.character(eye_area))]
                               daphnia_output[, eye_area_mm := eye_area/(pixel_to_mm^2)]

                               daphnia_output[, animal_length := as.numeric(as.character(animal_length))]
                               daphnia_output[, animal_length_mm :=  animal_length/pixel_to_mm]

                               daphnia_output[, tail_spine_length := as.numeric(as.character(tail_spine_length))]
                               daphnia_output[, tail_spine_length_mm :=  tail_spine_length/pixel_to_mm]

                               daphnia_output[, pedestal_max_height := as.numeric(as.character(pedestal_max_height))]
                               daphnia_output[, pedestal_max_height_mm :=  pedestal_max_height/pixel_to_mm]

                               daphnia_output[, pedestal_area := as.numeric(as.character(pedestal_area))]
                               daphnia_output[, pedestal_area_mm :=  pedestal_area/(pixel_to_mm^2)]

                               ## add pedestal_area_corrected and pedestal_height_corrected
                               daphnia_output[, pedestal_height_mm_cor := pedestal_max_height_mm / animal_length_mm, by = filebase]
                               daphnia_output[, pedestal_area_mm_cor := pedestal_area_mm / animal_dorsal_area_mm, by = filebase]
                                
                               ## add 'deme' (i.e., independence level of offspring)
							   daphnia_output[, deme := substring(replicate, 1, 1)]
                        
                           ### add neckteeth data
                               # NOTE - neckteeth data based on 2 rounds w/ UGs; RP checked 'mis-classified' data (i.e., shifted values between 1st and 2nd run) 
                               neckteeth <- fread("/mnt/spicy_4/daphnia/analysis/analysis_results/neckteeth_data_1to3_final_wide.csv")
                               neckteeth[, filebase_cor := substring(neckteeth$filebase, 7)]

                               # add filebase_cor
                               daphnia_output[, filebase_cor := substring(daphnia_output$filebase, 6)]

                               # merge neckteeth data to pheno data
                               daphnia_output <- merge(daphnia_output, neckteeth[, c("filebase_cor", "nteeth")], by = 'filebase_cor', all.x=T)
                               daphnia_output[, nteeth := as.numeric(nteeth)]

						  ### add supercluster assignment
						  		superclones <- fread("/mnt/spicy_4/daphnia/analysis/analysis_results/superclones_Karen_May2019_updated.csv")
						       	#superclones <- fread("/mnt/spicy_4/daphnia/analysis/analysis_results/superclones_Karen_April2020.csv")

							   # merge superclone data to pheno data
                               daphnia_output <- merge(daphnia_output, superclones[, c("cloneid_geno", "cloneid_geno_adj", "cloneid", "SC", "medrd")], by = 'cloneid', all.x=T)
                        
                           ### subset?
                               if(trimColumns==T) {
                                   daphnia_output <- daphnia_output[,c("filebase", "inductiondate", "pixel_to_mm", "rig", "season", "pond",
													"cloneid", "cloneid_geno", "cloneid_geno_adj", "replicate", "deme", "treatment", "day", "first_day", "barcode",
													"age", "age_hours", "use.mod", "SC", "experimenter", "inducer", "modifier", 
													"animal_dorsal_area_mm", "eye_area_mm", "animal_length_mm", "tail_spine_length_mm",
													"pedestal_max_height_mm", "pedestal_area_mm", "nteeth", "modification_notes", "medrd")]
                               }

                           ### return
                                return(daphnia_output)
                        }

        
        load_procrust_data <- function(dir="/mnt/spicy_4/daphnia/analysis/analysis_results/GUI_Outputs/Outputs/",
             dateStamp="20190715") {

                ### load shape data
                    shape_output <- fread(paste(dir, "shape_output_all_final_", dateStamp, ".txt", sep=""), header=T)
               
                ### summary stats per filebase
                    shape_output.ag <- shape_output[,list(n=length(qi[!is.na(qi)])), list(filebase)]

                ### flag those with >700
                    setkey(shape_output, filebase)
                    setkey(shape_output.ag, filebase)

                    shape_output[,n:=shape_output[shape_output.ag]$n]
                    shape_output[,use:=(i<=699 & n==700)]

                ### return
                    shape_output

             }



        assign_instar <- function(dat, quantile.trim=.01, trt=c(0,.1,.25,.5,1), ageset=c(1,2,3), useMod=T) {
         # NOTE_DB: mixture model now run on all treatments (instead of only using ctrl and 0.5 juju)
        
            ### trim on age & pheno  & mod
                dat <- dat[treatment%in%trt][age%in%ageset][use.mod==useMod]

            ### normalize
                dat[,animal_length_mm.norm := (animal_length_mm - mean(animal_length_mm, na.rm=T))/sd(animal_length_mm, na.rm=T)]
                dat[,animal_dorsal_area_mm.norm := (animal_dorsal_area_mm - mean(animal_dorsal_area_mm, na.rm=T))/sd(animal_dorsal_area_mm, na.rm=T)]

            ### trim on quantile
                dat <- dat[animal_length_mm.norm > quantile(animal_length_mm.norm, quantile.trim) &
                        animal_length_mm.norm < quantile(animal_length_mm.norm, 1 - quantile.trim)]

                dat <- dat[animal_dorsal_area_mm.norm > quantile(animal_dorsal_area_mm.norm, quantile.trim) &
                        animal_dorsal_area_mm.norm < quantile(animal_dorsal_area_mm.norm, 1 - quantile.trim)]

          
            ### Fit mixture model using Mclust; using length and dorsal area gives sensible results; adding eye area seems to mess things up
                #mod_bic <- MclustBIC(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], prior = priorControl())
            
            ### run mixture model w/ prior: model with 3 components is good (and justifiable)
            	set.seed(0)
            	mod <- Mclust(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], G=3, prior = priorControl())
                print(summary(mod$BIC), parameters=T)

            ### perform dimensionality reduction
                drmod <- MclustDR(mod, lambda = 1)
                plot(drmod, what="classification")

            ### conduct bootstrap LRT
                # LRT <- mclustBootstrapLRT(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], modelName = mod$modelName)

            ### tack in classification to dat
                dat[,instar_old:=mod$classification]

            ### tack in 'uncertainty' (1-probability of most likely classification) to dat
                dat[,uncertainty_old:=mod$uncertainty]
				
			
	    ### run mixture model AGAIN - w/ prior and adding noise/outliers (> poisson distribution of noise)
		dat_est <- apply(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], 2, range) 
		nNoise <- 2000
				
		set.seed(0)
		poissonNoise <- apply(dat_est, 2, function(x, n)
                runif(n, min = min(x)-.1, max = max(x)+.1), n = nNoise)
                datNdata <- rbind(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], poissonNoise)
                
                set.seed(0)
                datNoiseInit <- sample(c(TRUE,FALSE),size=nrow(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F]) + nNoise, replace=TRUE)				
				mod_noise <- Mclust(datNdata[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F], G=3, prior = priorControl(), 
								initialization = list(noise = datNoiseInit))
				print(summary(mod_noise$BIC), parameters=T)

				plot(dat[,c("animal_length_mm", "animal_dorsal_area_mm"), with=F])
				points(poissonNoise, pch = 20, cex = 0.3, col = "lightgrey") 
				mclust2Dplot(mod_noise$data, classification=mod_noise$classification, 
								parameters=mod_noise$parameters) 
				
		### tack in classification to dat
                dat[,instar:=mod_noise$classification[1:nrow(dat)]]

            	### tack in 'uncertainty' (1-probability of most likely classification) to dat
                dat[,uncertainty:=mod_noise$uncertainty[1:nrow(dat)]]


            ### make return object
                ret <- dat[,c("filebase", "instar_old", "uncertainty_old", "instar", "uncertainty"), with=F]
                setkey(ret, filebase)
                return(ret)
        }



    ### load data
        dap.sum <- na.omit(load_output_data())     
        dim(dap.sum)		#  6863   31    
        
        dap.shape <- load_procrust_data()     
        dim(dap.shape)		#4966810       9
        
    
    ### assign instar
        setkey(dap.sum, filebase)
        instar <- assign_instar(dat=dap.sum, quantile.trim=.005)
        dap.sum <- instar[J(dap.sum)]
        
        instar_old.age <- na.omit(dap.sum[,list(meanAge=mean(age)), list(instar_old)])
        dap.sum[,instar_old := as.numeric(factor(instar_old, levels=instar_old.age[order(meanAge)]$instar_old))]

        instar.age <- na.omit(dap.sum[!instar == 0][,list(meanAge=mean(age)), list(instar)])
        dap.sum[,instar := as.numeric(factor(instar, levels=instar.age[order(meanAge)]$instar))]

        ### is our head screwed on straight?
            ### basic info
                table(dap.sum$instar_old, dap.sum$age)
                table(dap.sum$instar, dap.sum$age)

                range(dap.sum$uncertainty_old, na.rm=T) 
                quantile(dap.sum$uncertainty_old, na.rm=T)
                
                range(dap.sum$uncertainty, na.rm=T) 
                quantile(dap.sum$uncertainty, na.rm=T)
                
        
    	### basic plots
    	# instar assignment
		instar_old_plot <- ggplot(data=dap.sum[!is.na(instar_old)], aes(x=animal_length_mm, y=animal_dorsal_area_mm, color=as.factor(instar_old), fill=as.factor(instar_old), group=instar_old)) +
					geom_point() + 
					theme(legend.position = "bottom")

		instar_plot <- ggplot(data=dap.sum[!is.na(instar)], aes(x=animal_length_mm, y=animal_dorsal_area_mm, color=as.factor(instar), fill=as.factor(instar), group=instar)) +
					geom_point() + 
					theme(legend.position = "bottom")		
		
		plot_grid(instar_old_plot, instar_plot)
        # NOTE: use "instar" (and not "instar_old") for downstream analyses (i.e., output from noise-added model)

	# animal length
        animal_length_plot <- ggplot(data=dap.sum[!is.na(instar)], aes(animal_length_mm, color=as.factor(instar), fill=as.factor(instar), group=instar)) +
					geom_freqpoly(bins=20) + facet_grid(treatment~age) +
					geom_hline(yintercept=0, colour="white", size=1) +
					theme(legend.position = "bottom")

        animal_dorsalArea_plot <- ggplot(data=dap.sum[!is.na(instar)], aes(animal_dorsal_area_mm, color=as.factor(instar), fill=as.factor(instar), group=instar)) +
					geom_freqpoly(bins=20) + facet_grid(treatment~age) +
					geom_hline(yintercept=0, colour="white", size=1) +
					theme(legend.position = "bottom")

        plot_grid(animal_length_plot, animal_dorsalArea_plot)

              
        ### are reversals common? - clunky code / modify (!)
			per.dap <- dap.sum[,list(n=length(instar)), list(barcode)]
			setkey(dap.sum, barcode)
				
			o <- foreach(i=list(c(2,1), c(3,1), c(3,2)))%do%{
							dap.sum[J(per.dap[n>1]$barcode)][,list(
								instar.delta=instar[age==i[1]] - instar[age==i[2]],
                                age.delta=i[1]-i[2],
                                t1=instar[age==i[1]],
                                t0=instar[age==i[2]],
                                a1=age[age==i[1]],
                                a0=age[age==i[2]],
                                length.t1=animal_length_mm[age==i[1]],
                                length.t0=animal_length_mm[age==i[2]]),
                                list(barcode)]
                            }

            o <- na.omit(rbindlist(o))
                 
		### summarize
       		 	o[t1 < t0]  # 56
       		 	table(o$t1, o$t0)
				prop.table(table(o$instar.delta>=0))
				# good. ~1% of individuals seem to get younger. exclude from downstream analyses. 
					
				# plot
				ggplot(data=o, aes(x=length.t1, y=length.t0, color=as.factor(t1-t0)))+ geom_abline(intercept = 0, slope = 1) + geom_point() + facet_grid(~age.delta)
					

		### cleaned data to use for downstream analyses (i.e., lmer(), GWAS, etc)
			dap.sum_cor <- dap.sum[!barcode %in% o[t1 < t0]$barcode][pond %in% c('DBunk', 'D8', 'D10')][!is.na(instar)]
			
		### add Julian day for LMER() 
			tmp_julian <- as.POSIXlt(dap.sum_cor$inductiondate, format="%Y%m%dT%H%M%S")
			dap.sum_cor[, inductiondate_Julian := tmp_julian$yday]
						
			#add 'batch'
			# code is clunky, but works (!)
			dap.sum_cor[, batch :=  as.factor(ifelse((dap.sum_cor$inductiondate_Julian == '142'), 'batch1', 
									ifelse(dap.sum_cor$inductiondate_Julian == '148' | dap.sum_cor$inductiondate_Julian == '149', 'batch2', 
								 	ifelse(dap.sum_cor$inductiondate_Julian == '154' | dap.sum_cor$inductiondate_Julian == '155', 'batch3', 
									ifelse(dap.sum_cor$inductiondate_Julian == '161' | dap.sum_cor$inductiondate_Julian == '162', 'batch4', 
									ifelse(dap.sum_cor$inductiondate_Julian == '168' | dap.sum_cor$inductiondate_Julian == '169' | dap.sum_cor$inductiondate_Julian == '170', 'batch5', 
									ifelse(dap.sum_cor$inductiondate_Julian == '189' | dap.sum_cor$inductiondate_Julian == '190', 'batch6', 
									ifelse(dap.sum_cor$inductiondate_Julian == '196' | dap.sum_cor$inductiondate_Julian == '197', 'batch7', 
									ifelse(dap.sum_cor$inductiondate_Julian == '203' | dap.sum_cor$inductiondate_Julian == '204', 'batch8', 
									ifelse(dap.sum_cor$inductiondate_Julian == '210' | dap.sum_cor$inductiondate_Julian == '211', 'batch9', 
									ifelse(dap.sum_cor$inductiondate_Julian == '217' | dap.sum_cor$inductiondate_Julian == '218', 'batch10', 
									ifelse(dap.sum_cor$inductiondate_Julian == '224' | dap.sum_cor$inductiondate_Julian == '225', 'batch11', 
									ifelse(dap.sum_cor$inductiondate_Julian == '259' | dap.sum_cor$inductiondate_Julian == '260', 'batch12', 
									ifelse(dap.sum_cor$inductiondate_Julian == '274', 'batch13', 
									ifelse(dap.sum_cor$inductiondate_Julian == '280' | dap.sum_cor$inductiondate_Julian == '281', 'batch14', 
									ifelse(dap.sum_cor$inductiondate_Julian == '294' | dap.sum_cor$inductiondate_Julian == '295', 'batch15', 
									ifelse(dap.sum_cor$inductiondate_Julian == '315' | dap.sum_cor$inductiondate_Julian == '316', 'batch16', 
									ifelse(dap.sum_cor$inductiondate_Julian == '329', 'batch17', 
									ifelse(dap.sum_cor$inductiondate_Julian == '336' | dap.sum_cor$inductiondate_Julian == '337', 'batch18', 
									ifelse(dap.sum_cor$inductiondate_Julian == '343' | dap.sum_cor$inductiondate_Julian == '344', 'batch19', 
									ifelse(dap.sum_cor$inductiondate_Julian == '350', 'batch20', 'batchNA')))))))))))))))))))))]


### add unique supercluster level for "Os"
### code is pretty clunky right now - needs to be tidied up & generalized (!)
	table(unique(dap.sum_cor[, list(SC, cloneid)])$SC) 
	#    A  AD  AE  AF  AO  AP  AQ   C   D   E   F   H   K   L   M  OO poW   R   W 
 	#    62   2   1   1   1   2   2   9  21  17   8   9  10   5   4  39   1   4   2 
    
	dap.sum_cor[, SC_adj := SC]
	dap.sum_cor[, SC_adj := gsub("AE", "OO", SC_adj)]
	dap.sum_cor[, SC_adj := gsub("AF", "OO", SC_adj)]
	dap.sum_cor[, SC_adj := gsub("AO", "OO", SC_adj)]
	dap.sum_cor[, SC_adj := gsub("poW", "OO", SC_adj)]
		
	table(unique(dap.sum_cor[, list(SC_adj, cloneid)])$SC_adj)
	
	# A AD AP AQ  C  D  E  F  H  K  L  M OO  R  W 
    # 62  2  2  2  9 21 17  8  9 10  5  4 43  4  2 

	supercluster <- unique(dap.sum_cor[, list(SC_adj, cloneid)])
	setkey(supercluster, SC_adj, cloneid)
	
    supercluster[, SC_unique := ifelse((supercluster$SC_adj %like% "OO"), paste0("OO", .I), SC_adj)]
    # modify if Os are supposed to be O1:O42 (instead of O152:O193)
    
    # merge unique cluster ID to dap.sum_cor
	setkey(dap.sum_cor, cloneid)
	setkey(supercluster, cloneid)
		
    dap.sum_cor <- supercluster[, c('cloneid', 'SC_unique')][J(dap.sum_cor)]


### interpolation for shape data.
     setkey(dap.shape, filebase, i)
      
      tmp_shape <- dap.shape[use == 'TRUE'][, c("filebase","qi","q")][filebase %in% dap.sum_cor$filebase]  ## qi and q from Austins calculations
      colnames(tmp_shape) <- c('filebase', 'x', 'y')
      
      ## force x data 'really' into 0:1 space (i.e., divide x by x.max (per sample))
      tmp_shape_adj <- foreach(filebase.i = unique(tmp_shape$filebase), .errorhandling="remove", .combine="rbind") %do% {
       
      						tmp_data <- tmp_shape[filebase == filebase.i]
      						tmp_data[, x.max := max(tmp_data$x)]
      						tmp_data[, x.adj := x/x.max]
      			      			
      					}
      
      ## interpolate x data to get evenly spaced coordinates
      tmp_shape_adj_interpol <- tmp_shape_adj[, c("filebase","x.adj","y")]
      colnames(tmp_shape_adj_interpol) <- c('filebase', 'x', 'y')

      setkey(tmp_shape_adj_interpol, filebase)


      approxData_use <- foreach(sample = unique(tmp_shape_adj_interpol$filebase)) %do% {
          data.frame(filebase = sample, 
                   with(tmp_shape_adj_interpol[sample], 
                        approx(x, y, xout = seq(min(tmp_shape_adj_interpol$x), max(tmp_shape_adj_interpol$x), by = (1/699)), method = "linear")
                   ), 
                   method = "approx()"
          )
        }
    
      approxData_use_unlist <- rbindlist(lapply(approxData_use, `[`, c('filebase', 'x', 'y')))
      # NAs in data all due to missing values at x==0

      # add 'i' to data - 1:700 for each filebase
      setkey(approxData_use_unlist, filebase)
      
      approxData_use_unlist[,i:= rep(1:700, times=6141, length.out=4298700)]
      

      ### tack instar assignment (and other info) into shape data
      setkey(approxData_use_unlist, filebase)
      setkey(dap.sum_cor, filebase)

      dap.shape_sum <- merge(approxData_use_unlist, dap.sum_cor)
      
      # normalize x and y so that data is in mm
      dap.shape_sum[,height_mm:=y/pixel_to_mm]
      dap.shape_sum[, height_mm_norm := height_mm / animal_length_mm]  
      dap.shape_sum[,length_mm:=x*animal_length_mm]
      
      
      ## final shape data 
      shape_use <- dap.shape_sum[, c("filebase","barcode","cloneid_geno","cloneid_geno_adj","i","instar","treatment","SC","SC_unique","batch","pond","season","deme","replicate","pixel_to_mm",
                               "uncertainty","height_mm","height_mm_norm","length_mm","animal_length_mm","animal_dorsal_area_mm","eye_area_mm","nteeth","tail_spine_length_mm",
                               "pedestal_max_height_mm","pedestal_area_mm","medrd","modification_notes","inducer")][!modification_notes %like% 'noise']
    
      
      shape_use.ag <- shape_use[,list(height_mm=mean(height_mm, na.rm=T), 
      								  height_mm_norm=mean(height_mm_norm, na.rm=T), 
      								  length_mm=mean(length_mm, na.rm=T),
      								  animal_length_mm=mean(animal_length_mm, na.rm=T),
      								  eye_area_mm=mean(eye_area_mm, na.rm=T),
      								  nteeth=median(nteeth, na.rm=T)),       								  
                          list(cloneid_geno, barcode, SC_unique, treatment, instar, i, medrd, pond, season, batch, replicate, deme, uncertainty)]

      shape_use.ag_filtered <- shape_use.ag[treatment %in% c(0,0.5)][instar %in% c(1,2)][i>1]
      
      
      save(shape_use, file='/mnt/spicy_4/daphnia/analysis/analysis_results/dataANDresults/shape_use_12May2020.Rdata')
      save(shape_use.ag, file='/mnt/spicy_4/daphnia/analysis/analysis_results/dataANDresults/shape_use.ag_12May2020.Rdata')
      save(shape_use.ag_filtered, file='/mnt/spicy_4/daphnia/analysis/analysis_results/dataANDresults/shape_use.ag_filtered_12May2020.Rdata')
