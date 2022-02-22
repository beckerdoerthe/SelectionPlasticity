#module load bcftools intel/18.0 intelmpi/18.0 R/3.6.0; R


#################
### libraries ###
#################
  library(SeqArray)
  library(data.table)
  library(foreach)

# modify gds/genofile so that there is one sample.id per each pheno sample (i.e. barcode)
# open genofile & loop through barcodes and extract gds for each barcode
  load(file="/project/berglandlab/doerthe/GDS_for_AOB/data_use_O.medrd.RData")  ## pheno data for O cluster (w/ highest median read depth, 40ish clones)

### open this once
  genofile <- seqOpen("/project/berglandlab/doerthe/GDS_for_AOB/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds", allow.duplicate=T)

### output directory
  outdir <- "/project/berglandlab/alan/forDorthe/"

### open filter file
  snps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/snpsvarpulexpresentinhalf_table_wPulicaria_20200401")
  setkey(snps, chr)
  snps.ag <- snps[,list(.N, use=T), chr]
  snps <- snps[J(snps.ag[order(N, decreasing=T)][1:12])]


## extract out individual VCFs for each barcode
  sc <- unique(data_use_O.medrd$SC_unique)

  foreach(sc.i = sc)%do%{
    #sc.i<-sc[2]
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=data_use_O.medrd[SC_unique==sc.i]$Geno[1], variant.id=snps$variant.ids)

    ### export master VCF file per superclone
    seqGDS2VCF(genofile,
              vcf.fn=paste(outdir, sc.i, ".vcf.gz", sep=""),
              info.var=character(0),
              fmt.var=character(0),
              use_Rsamtools=F,
              verbose=TRUE)

    system(paste("gunzip", paste(outdir, sc.i, ".vcf.gz", sep=""), sep=" "))

    ### copy and modify
    tmp <- data_use_O.medrd[SC_unique==sc.i]

    foreach(i=1:dim(tmp)[1])%do%{

      #i<-1

      barcodePlus.i <- tmp$cloneID_barcodePLUS[i]
      geno.i <- tmp$Geno[i]

      message(paste(sc.i, barcode.i, sep=" / "))
      system(paste("cp",
                   paste(outdir, sc.i, ".vcf", sep=""),
                   paste(outdir, barcodePlus.i, ".vcf", sep=""),
                   sep=" "))

      write.table(data.table(old=geno.i, new=barcodePlus.i), sep="\t", quote=F, row.names=F, col.names=F, file=
                  paste(outdir, barcodePlus.i, ".rename", sep=""))

      system(paste("bcftools reheader -s",
                   paste(outdir, barcodePlus.i, ".rename", sep=""),
                   paste(outdir, barcodePlus.i, ".vcf ", sep=""),
                   "|",
                   "bcftools view -O b >",
                   paste(outdir, barcodePlus.i, ".reheader.vcf.gz", sep=""), sep=" "))

      system(paste("bcftools index", paste(outdir, barcodePlus.i, ".reheader.vcf.gz", sep=""), sep=" "))

    }

     ### clean up
      system(paste("rm",
                   paste(outdir, sc.i, ".vcf", sep=""),
                   paste(outdir, barcode.i, ".rename", sep=""),
                   sep=" "))
   }

### merge vcf files
  system(paste("bcftools merge ",
               "-o ", outdir, "samps.bcf ",
               "-O b ",
               outdir, "*.reheader.vcf.gz ", sep=""))
