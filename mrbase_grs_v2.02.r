##########################################################################################
#Script name: mrbase_grs_v2.01.r
#Project: GRS_SEP
#Script author: Sean Harrison
#Script purpose: Derive genetic risk scores using summary data from GWAS
#Last edited: 15/03/2019
#Date created: 01/05/2018
#Notes: This script has extensive documentation in "mrbase_grs documentation.docx"

#Update 2.02: MAF updated to 0.42 to fit with MR Base
#Update 2.01: Updated to work with new UK Biobank genetic data (still BC3 though)
#Update 1.04: Changed the bgen module used on BC3 to 1.1.4
##########################################################################################

mrbase_grs = function(output=NULL,category=NULL,subcategory=NULL,trait=NULL,population=NULL,sex=NULL,mr=NULL,samplesize=0,
                      notstudies=NULL,studies=NULL,p=5e-08,r2=0.8, 
                      exposure_dat=NULL,exposure_file=NULL,maf=0.42,ipd=FALSE,
                      proxies=TRUE,clump=TRUE,clumped=TRUE,gwas="biggest",suffix="",
                      snpstats_file="SNPstats.txt",ld_file="1000_genomes_ld.csv.gz",keep_files=FALSE,plink_grs=FALSE,
                      bgen_folder="/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen") {
   if(is.null(exposure_dat)){
    return_dat = TRUE
  } else {
    return_dat = FALSE
  }
  "%ni%" = Negate("%in%")
  #package install
  print("Installing packages if necessary")
  if("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
  library(devtools)
  if("TwoSampleMR" %in% rownames(installed.packages()) == FALSE) {install_github("MRCIEU/TwoSampleMR")}
  if("MRInstruments" %in% rownames(installed.packages()) == FALSE) {devtools::install_github("MRCIEU/MRInstruments")}
  if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
  if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
  if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}
  
  library(TwoSampleMR)
  library(MRInstruments)
  library(data.table)
  library(plyr); library(dplyr)
  
  #Warnings and errors
  if(is.null(output)) {
    warning("Please enter a value for 'output'")
  }
  
  ###################Subcategories, Traits, Studies##################################    
  
  else if(tolower(output)=="subcategories" | tolower(output)=="traits" | tolower(output)=="studies") {
    
    ##############Preprocessing##############################
    
    #Grab the available outcomes file and do some preprocessing
    print(paste("Output = ",output))
    mrbase_grs_ao = available_outcomes()
    mrbase_grs_ao$sample_size[is.na(mrbase_grs_ao$sample_size)] = 0
    mrbase_grs_ao$consortium[is.na(mrbase_grs_ao$consortium)] = "NA"
    
    #Replace any null values with ALL the values in the available outcomes
    if(is.null(category)){
      category = unique(mrbase_grs_ao$category)
    }
    if(is.null(population)){
      population = unique(mrbase_grs_ao$population)
    }
    if(is.null(subcategory)){
      subcategory = unique(mrbase_grs_ao$subcategory)
    }
    if(is.null(trait)){
      trait = unique(mrbase_grs_ao$trait)
    }
    if(is.null(sex)){
      sex = unique(mrbase_grs_ao$sex)
    }
    if(is.null(mr)){
      mr = unique(mrbase_grs_ao$mr)
    }
    if(is.null(studies)){
      studies = unique(mrbase_grs_ao$id)
    }
    
    #Restrict the outcomes to only those of interest
    mrbase_grs_ao = mrbase_grs_ao[which(mrbase_grs_ao$"category" %in% category & mrbase_grs_ao$"subcategory" %in% subcategory &
                                          mrbase_grs_ao$"population" %in% population & mrbase_grs_ao$"sample_size" >= samplesize & 
                                          mrbase_grs_ao$"id" %ni% notstudies & mrbase_grs_ao$"trait" %in% trait & 
                                          mrbase_grs_ao$"sex" %in% sex & mrbase_grs_ao$"mr" %in% mr & mrbase_grs_ao$"id" %in% studies),]
    print("Outcomes processed")
    
    ###################Subcategories##################################
    
    if (tolower(output)=="subcategories"){
      print("Processing subcategories")
      mrbase_grs_ao = mrbase_grs_ao[c("id","subcategory","trait")]
      mrbase_grs_subcategories = unique(mrbase_grs_ao$"subcategory")
      
      #Create new summary dataframe
      mrbase_grs_subcatsum = data.frame(subcategory=mrbase_grs_subcategories)
      for(subcat in mrbase_grs_subcategories) {
        mrbase_grs_subcatsum$"studies"[mrbase_grs_subcatsum$"subcategory" == subcat] = length(mrbase_grs_ao$"id"[mrbase_grs_ao$"subcategory" == subcat])
        mrbase_grs_subcatsum$"traits"[mrbase_grs_subcatsum$"subcategory" == subcat] = length(unique(mrbase_grs_ao$"trait"[mrbase_grs_ao$"subcategory" == subcat]))
      }
      print("Processing complete")
      return(mrbase_grs_subcatsum)
    }
    
    ###################Traits##################################    
    
    else if(tolower(output)=="traits"){
      print("Processing traits")
      #Return the number of studies & number of SNPs
      mrbase_grs_ao = mrbase_grs_ao[c("id","trait","subcategory")]
      mrbase_grs_traits = unique(mrbase_grs_ao$"trait")
      if(!exists("mrbase_grs_traits")){
        warning("No traits found in MR BASE using criteria supplied, make sure your options are correctly specified")
      }
      else {
        #Create new summary dataframe
        mrbase_grs_traitsum = data.frame("trait"=mrbase_grs_traits)
        #Many traits could be passed, cycle through
        for(trait_list in mrbase_grs_traits) {
          #Number of studies looking at the trait
          mrbase_grs_traitsum$"studies"[mrbase_grs_traitsum$"trait" == trait_list] = length(mrbase_grs_ao$"id"[mrbase_grs_ao$"trait" == trait_list])
          mrbase_grs_subcategories = mrbase_grs_ao$"subcategory"[mrbase_grs_ao$"trait" == trait_list]
          mrbase_grs_traitsum$"subcategory"[mrbase_grs_traitsum$"trait" == trait_list] = mrbase_grs_subcategories[1]
        }
        mrbase_grs_traitsum = mrbase_grs_traitsum[c("subcategory","trait","studies")]
        print("Processing complete")
        return(mrbase_grs_traitsum)
      }
    }  
    
    ###################Studies##################################    
    
    else if(tolower(output)=="studies"){
      #Here I need to list the individual studies, and everything interesting about them
      #Subcat, trait, lots of things in the outcomes file
      print("Processing studies")
      #mrbase_grs_ao = mrbase_grs_ao[c("author","consortium","id", "mr", "ncase", "ncontrol","note","nsnp","pmid","population","priority",
      #                                "sample_size","sd","subcategory","trait","unit","year")]
      return(mrbase_grs_ao)
    }
  }
  
  ###################SNPS, Code, GRS##################################    
  
  else if(tolower(output)=="snps" | tolower(output)=="code") {
    
    ###################SNPs##################################  
    
    if(tolower(output)=="snps"){
      #Here is where I extract all SNP data
      #Include other things like study, trait, units etc.
      if(is.null(studies)){
        #Grab the available outcomes file and do some preprocessing
        print(paste("Output = ",output))
        mrbase_grs_ao = available_outcomes()
        mrbase_grs_ao$sample_size[is.na(mrbase_grs_ao$sample_size)] = 0
        mrbase_grs_ao$consortium[is.na(mrbase_grs_ao$consortium)] = "NA"
        
        #Replace any null values with ALL the values in the available outcomes
        if(is.null(category)){
          category = unique(mrbase_grs_ao$category)
        }
        if(is.null(population)){
          population = unique(mrbase_grs_ao$population)
        }
        if(is.null(subcategory)){
          subcategory = unique(mrbase_grs_ao$subcategory)
        }
        if(is.null(trait)){
          trait = unique(mrbase_grs_ao$trait)
        }
        if(is.null(sex)){
          sex = unique(mrbase_grs_ao$sex)
        }
        if(is.null(mr)){
          mr = unique(mrbase_grs_ao$mr)
        }
        if(is.null(studies)){
          studies = unique(mrbase_grs_ao$studies)
        }
        
        #Restrict the outcomes to only those of interest
        mrbase_grs_ao = mrbase_grs_ao[which(mrbase_grs_ao$"category" %in% category & mrbase_grs_ao$"subcategory" %in% subcategory &
                                              mrbase_grs_ao$"population" %in% population & mrbase_grs_ao$"sample_size" >= samplesize & 
                                              mrbase_grs_ao$"id" %ni% notstudies & mrbase_grs_ao$"trait" %in% trait & 
                                              mrbase_grs_ao$"sex" %in% sex & mrbase_grs_ao$"mr" %in% mr),]
        if(gwas == "biggest"){
          min = min(mrbase_grs_ao$priority)
          mrbase_grs_ao = mrbase_grs_ao[mrbase_grs_ao$priority == min,]
        }
        mrbase_grs_studies = mrbase_grs_ao$id
        print("Studies acquired")
        
      }
      else{
         mrbase_grs_studies = studies
      }
      #Run through studies (if existant)
      j = length(mrbase_grs_studies)
      if(j > 0){
        mrbase_grs_instruments = extract_instruments(outcomes=mrbase_grs_studies, p1 = p, "clump" = clump)
        print("Processing complete")
        return(mrbase_grs_instruments)  
      }
      else{
        print("No studies found")
      }
    }
    
    ###################Code##################################  
    
    else if(tolower(output)=="code"){
      print(paste("Output = ",output))
      #This is the code that will create a GRS
      #Input types: 
      # 1) any arguments that mean SNPs need to be downloaded
      # 2) a list of SNPs with effect estimates etc.
      # 3) SNP list, e.g. all SNPs in UK Biobank plus effect alleles, EAF etc.
      
      #####
      #Step 0
      #####
      
      #Step 0 is to find SNPs if not specified
      if(is.null(exposure_file) & is.null(exposure_dat)){
        if(is.null(studies)){
          ##############Preprocessing##############################
          
          #Grab the available outcomes file and do some preprocessing
          print("Finding studies")
          mrbase_grs_ao = available_outcomes()
          mrbase_grs_ao$sample_size[is.na(mrbase_grs_ao$sample_size)] = 0
          mrbase_grs_ao$consortium[is.na(mrbase_grs_ao$consortium)] = "NA"
          
          #Replace any null values with ALL the values in the available outcomes
          if(is.null(category)){
            category = unique(mrbase_grs_ao$category)
          }
          if(is.null(population)){
            population = unique(mrbase_grs_ao$population)
          }
          if(is.null(subcategory)){
            subcategory = unique(mrbase_grs_ao$subcategory)
          }
          if(is.null(trait)){
            trait = unique(mrbase_grs_ao$trait)
          }
          if(is.null(sex)){
            sex = unique(mrbase_grs_ao$sex)
          }
          if(is.null(mr)){
            mr = unique(mrbase_grs_ao$mr)
          }
          if(is.null(studies)){
            studies = unique(mrbase_grs_ao$studies)
          }
          
          #Restrict the outcomes to only those of interest (and specify mr == 1)
          mrbase_grs_ao = mrbase_grs_ao[which(mrbase_grs_ao$"category" %in% category & mrbase_grs_ao$"subcategory" %in% subcategory &
                                                mrbase_grs_ao$"population" %in% population & mrbase_grs_ao$"sample_size" >= samplesize & 
                                                mrbase_grs_ao$"id" %ni% notstudies & mrbase_grs_ao$"trait" %in% trait & 
                                                mrbase_grs_ao$"sex" %in% sex & mrbase_grs_ao$"mr" == 1),]
          #The gwas option comes in here - does the user want all SNPs, or SNPs from the biggest study, or a meta-analysis?
          trait = unique(mrbase_grs_ao$"trait")
          if(gwas == "biggest"){
            mrbase_grs_ao$keep = 0
            mrbase_grs_studies = c()
            for(t in trait) {
              mrbase_grs_ao$rank[mrbase_grs_ao$"trait" == t] = rank(-mrbase_grs_ao$sample_size[mrbase_grs_ao$"trait" == t],ties.method = "first")
              min = min(mrbase_grs_ao$rank,na.rm = TRUE)
              mrbase_grs_ao$keep[mrbase_grs_ao$rank == min & mrbase_grs_ao$"trait" == t] = 1
              mrbase_grs_ao$rank = NULL
            }
            mrbase_grs_ao = mrbase_grs_ao[which(mrbase_grs_ao$keep == 1),]
          }
          mrbase_grs_studies = mrbase_grs_ao$id
          print("Studies:")
          print(mrbase_grs_studies)
          print("Studies acquired")
          
        } else{
          mrbase_grs_studies = studies
          mrbase_grs_ao = available_outcomes()
          mrbase_grs_ao$sample_size[is.na(mrbase_grs_ao$sample_size)] = 0
          mrbase_grs_ao$consortium[is.na(mrbase_grs_ao$consortium)] = "NA"
          mrbase_grs_ao = mrbase_grs_ao[which(mrbase_grs_ao$"id" %in% studies & mrbase_grs_ao$"mr" == 1),]
          trait = unique(mrbase_grs_ao$"trait")
        }
        #Run through studies (if existant)
        j = length(mrbase_grs_studies)
        clumped = TRUE
        if(j > 0){
          exposure_dat = data.frame()
          for(study in mrbase_grs_studies){
            print(paste("Extracting data for study",study))
            temp = try(extract_instruments(outcomes=study, p1 = p, "clump" = TRUE),silent=TRUE)
            if(class(temp) == "try-error"){
              print(paste("Too many SNPs to clump using MR-Base, downloading unclumped SNPs for study",study))
              print(paste("If downloading takes too long, consider downloading from the GWAS website and coercing the download into an MR-Base format CSV"))
              temp = extract_instruments(outcomes=study, p1 = p, "clump" = FALSE)
              temp$clumped = FALSE
              clumped = FALSE
            } else {
              temp$clumped = TRUE
            }
            if(class(temp) == "data.table" | class(temp) == "data.frame"){
              exposure_dat = rbind(exposure_dat,temp)
            }
            rm(temp)
          }

          #The all option doesn't bother with meta-analysis, just takes all SNPs except lower priority duplicates
          if(gwas == "all"){
            #merge
            merge = mrbase_grs_ao[c("id","priority","trait")]
            exposure_dat = merge(exposure_dat, merge, by.x = "id.exposure", by.y = "id")
            exposure_dat$rank = 1
            for(t in trait){
              #Check to see how many studies are available for the trait
              study_count = length(unique(exposure_dat$id[exposure_dat$"trait" == t]))

              #If >1 study in the trait, need to remove duplicate SNPs and set clumped to FALSE (both in exposure_dat and for indicator value)
              if(study_count > 1){
                
                clumped = FALSE
                exposure_dat$clumped[exposure_dat$"trait" == t] = FALSE
                
                #For each SNP ID, rank by priority
                snps = unique(exposure_dat$SNP[exposure_dat$"trait" == t])
                
                for(snp in snps){
                  #This code ranks SNPs depending on their information, prioritising those with beta.exposure, effect_allele.exposure and eaf.exposure
                  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t] = 1
                  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t & !is.na(exposure_dat$beta.exposure)] =  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t]+5
                  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t & !is.na(exposure_dat$effect_allele.exposure)] =  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t]+3
                  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t & !is.na(exposure_dat$eaf.exposure)] =  exposure_dat$rank1[exposure_dat$SNP == snp & exposure_dat$"trait" == t]+1
                  exposure_dat$rank1x = rank(-exposure_dat$rank1,ties.method = "min",na.last="keep")
                  
                  #This code ranks SNPs from different studies by the sample size of the study
                  exposure_dat$samplesize2[exposure_dat$SNP == snp & exposure_dat$"trait" == t] = exposure_dat$samplesize.exposure[exposure_dat$SNP == snp & exposure_dat$"trait" == t]
                  exposure_dat$samplesize2[exposure_dat$SNP == snp & exposure_dat$"trait" == t & exposure_dat$rank1x > 1] = 0
                  exposure_dat$rank2 = rank(-exposure_dat$samplesize2,ties.method = "last",na.last="keep")
                  exposure_dat$rank[exposure_dat$SNP == snp & exposure_dat$"trait" == t] = exposure_dat$rank2[exposure_dat$SNP == snp & exposure_dat$"trait" == t]
                  exposure_dat$samplesize2 = NULL
                  exposure_dat$rank2 = NULL
                  exposure_dat$rank1x = NULL
                  exposure_dat$rank1 = NULL
                }
              }
            }
            exposure_dat = exposure_dat[exposure_dat$rank == 1,]
            exposure_dat$rank = NULL
          } else {
            merge = mrbase_grs_ao[c("id","priority","trait")]
            exposure_dat = merge(exposure_dat, merge, by.x = "id.exposure", by.y = "id")
          }
          
          print("Processing complete")
        } else{
          print("No studies found matching your criteria")
        }
      }
      
      #At this point, exposure_dat should exist unless no studies match the criteria
      #Need to create several files:
      #exposure_dat.csv - a CSV of all SNP info from MR-Base
      #harmonise.R - for harmonising SNPs between MR-Base and LD file (1000 genomes)/UK Biobank      
      #proxies.R - code to grab proxies from LD file (1000 genomes) - if proxies=FALSE, doesn't code
      #code.R - the code to write the code for the GRS (done after harmonisation)
      #run.sh - shell script to run everything on BC
      
      #Also need access to UK Biobank data (including SNP stats) & scratch (to copy over LD file - 1000 genomes)
      
      ##############################################################################################################
      
      #exposure_dat.csv
      if(!is.null(exposure_file)){
        exposure_dat_file = exposure_file
      } else {
        exposure_dat_file = paste("exposure_dat",suffix,".csv",sep="")
      }
      exposure_dat_harmonised_file = paste("exposure_dat_harmonised",suffix,".csv",sep="")
      
      snp_list_file = paste("snp_list",suffix,".txt",sep="")
      snp_score_list_file = paste("snp_score_list",suffix,".txt",sep="")
      snp_list_out_file = paste("snp_list_out",suffix,".txt",sep="")
      snp_list_proxies_file = paste("snp_list_proxies",suffix,".txt",sep="")
      snp_ipd_file = paste("snp_ipd",suffix,".csv",sep="")
      snp_ld_file = paste("snp_ld",suffix,".csv",sep="")
      snp_ld_test_file = paste("snp_ld_test",suffix,".csv",sep="")
      
      included_SNPs_file = paste("included_SNPs",suffix,".csv",sep="")
      included_proxy_SNPs_file = paste("included_proxy_SNPs",suffix,".csv",sep="")
      temp_geno_prefix = paste("temp_genos",suffix,sep="")
      
      instruments_file = paste("instruments",suffix,sep="")
      grs_file = paste("grs",suffix,".csv",sep="")
      script_file = paste("script",suffix,".R",sep="")
      run_file = paste("run",suffix,".sh",sep="")
      
      if(is.null(exposure_file)){
        write.csv(exposure_dat, exposure_dat_file)
      }
      
      ##############################################################################################################
      
      #script.R
      file = script_file
      write('#!/cm/shared/languages/R-3.5-ATLAS/bin/Rscript',file,append=FALSE)
      write('#PBS -l nodes=1:ppn=16',file,append=TRUE)
      write('#PBS -l walltime=00:12:00:00',file,append=TRUE)
      write('#PBS -N R_conversion',file,append=TRUE)
      write('library(plyr); library(dplyr); library(data.table)',file,append=TRUE)
      
      write(paste('exposure_dat = read.csv("',exposure_dat_file,'", stringsAsFactors = FALSE, strip.white=TRUE)',sep=""),file,append = TRUE)
      write('#Check whether the exposure_dat file has a "trait" variable - if so, create a trait list, otherwise trait = trait
if("trait" %in% names(exposure_dat) == TRUE){
            exposure_dat$trait[is.na(exposure_dat$trait)] = "missing"
            trait = unique(exposure_dat$trait)
    } else if("Trait" %in% names(exposure_dat) == TRUE){
            exposure_dat$trait[is.na(exposure_dat$Trait)] = "missing"
            trait = unique(exposure_dat$Trait)
            exposure_dat = rename(exposure_dat,trait = Trait)
    } else if("exposure.trait" %in% names(exposure_dat) == TRUE){
            exposure_dat$exposure.trait[is.na(exposure_dat$exposure.trait)] = "missing"
            trait = unique(exposure_dat$exposure.trait)
            exposure_dat=rename(exposure_dat,trait = exposure.trait)  
    } else if("exposure.Trait" %in% names(exposure_dat) == TRUE){
            exposure_dat=rename(exposure_dat,trait = exposure.Trait) 
            exposure_dat$exposure.trait[is.na(exposure_dat$exposure.trait)] = "missing"
            trait = unique(exposure_dat$exposure.trait)
    } else {
            trait = "trait"
            exposure_dat$trait = "trait"
    }',file,append=TRUE)
      
      write('exposure_dat$trait = gsub(" ","_",exposure_dat$trait)',file,append=TRUE)
      write(paste('exposure_dat$trait = gsub("',"'",'","",exposure_dat$trait)',sep=""),file,append=TRUE)
      write(paste("exposure_dat$trait = gsub('",'"',"','',exposure_dat$trait)",sep=""),file,append=TRUE)
      write('trait = gsub(" ","_",trait)',file,append=TRUE)
      write(paste('trait = gsub("',"'",'","",trait)',sep=""),file,append=TRUE)
      write(paste("trait = gsub('",'"',"','',trait)",sep=""),file,append=TRUE)
      
      #If any traits need clumping, do that here
      if(clumped == FALSE | "clumped" %in% names(exposure_dat)){
        write('#Deal with the clumped traits',file,append=TRUE)
        
        #Check for a clumped variable
        write('if("clumped" %in% names(exposure_dat) == FALSE){
  exposure_dat$clumped = FALSE
}',file,append=TRUE)
        
        #Cycle through traits and use Plink2 to clump - exposure_dat$clumped=FALSE needs dealing with
        write('exposure_dat$clump_num[exposure_dat$clumped == TRUE] = 1',file,append=TRUE)
        write('exposure_dat$clump_num[exposure_dat$clumped == FALSE] = 0',file,append=TRUE)
        write(paste('for(t in trait){
  #The whole trait should either be TRUE or FALSE
  clumped_mean = mean(as.numeric(exposure_dat$clump_num[exposure_dat$trait == t]))
  if(clumped_mean == 0){
    unclumped_snps = exposure_dat[exposure_dat$trait == t,c("SNP","pval.exposure")]
    unclumped_snps = rename(unclumped_snps,P = pval.exposure)
    print(paste("Clumping SNPs for trait: ",t,sep=""))
    write.table(unclumped_snps,"unclumped_snps.txt",row.names=FALSE,quote=FALSE)
    system(\'module rm apps/plink-2.00
module add apps/plink-1.90b4.1
plink -bfile eur --clump unclumped_snps.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --out clumped_snps
tr -s [:blank:] < clumped_snps.clumped | cut -f 4 -d " " > clumped_snps.txt
module rm apps/plink-1.90b4.1
module add apps/plink-2.00\')

  #Then bring the SNP list back into R & remove pruned SNPs
  clumped_snps = read.table("clumped_snps.txt",stringsAsFactors = FALSE, strip.white=TRUE)
  clumped_snps = clumped_snps$V1
  exposure_dat = exposure_dat[which(exposure_dat$SNP %in% clumped_snps | exposure_dat$trait != t),]
  }
}',sep=""),file,append=TRUE)
      }
      
      write(paste('#Check which SNPs are in the SNPstats file'),file,append = TRUE)
      write("exposure_dat$effect_allele.exposure[exposure_dat$effect_allele.exposure == TRUE]='T'",file,append=TRUE)
      write("exposure_dat$other_allele.exposure[exposure_dat$other_allele.exposure == TRUE]='T'",file,append=TRUE)
      write('snp_list = unique(exposure_dat$SNP)',file,append = TRUE)
      write(paste('write(as.character(snp_list),"',snp_list_file,'", append = FALSE)',sep=""),file,append=TRUE)
      
      l1 = "system('folder=/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats"
      l2 = paste('if [ -f "',snpstats_file,'" ]
then
  echo SNPstats file present, proceeding
else
  echo Concatenating UK Biobank SNP Stats in home directory
  cmd=""
  for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do
    tail -n +17 ${folder}/data.chr${chrom}.snp-stats > data.chr${chrom}.txt
    cut -f 2- data.chr${chrom}.txt > chr${chrom}.txt
    cmd="${cmd} chr${chrom}.txt"
  done
  cat ${cmd} > ',snpstats_file,'
  rm data.chr*
  rm chr*
fi
echo Searching SNPstats for SNPs from MR-Base
zgrep -w -F -f ',snp_list_file,' ',snpstats_file,' > ',included_SNPs_file,'
echo Search complete, ',included_SNPs_file,' file created, proceeding in R\')',sep="")

      write(paste(l1,l2,sep="\n"),file,append = TRUE)
      write(paste('snp_list_in = read.csv("',included_SNPs_file,'",header=FALSE,sep="",stringsAsFactors = FALSE, strip.white=TRUE)',sep=""),file,append = TRUE)
      #Remove triallelics (duplicates in the SNPstats file) from exposure_dat
      write('"%ni%" = Negate("%in%")',file,append=TRUE)
      write("x = length(exposure_dat$SNP)",file,append=TRUE)
      write("exposure_dat = exposure_dat[which(exposure_dat$SNP %ni% snp_list_in$V1[duplicated(snp_list_in$V1)]),]",file,append = TRUE)
      write(paste("y = length(exposure_dat$SNP)","if(x!=y){","  print(paste(x-y,' triallelic SNPs dropped',sep=''))","}",sep='\n'),file,append=TRUE)
      write("snp_list_in$V4[snp_list_in$V4==TRUE]='T'",file,append = TRUE)
      write("snp_list_in$V5[snp_list_in$V5==TRUE]='T'",file,append = TRUE)
      
      if(proxies == TRUE) { 
        write("snp_list_in_list = unique(snp_list_in$V1)",file,append = TRUE)
        write('snp_list_out = setdiff(snp_list,snp_list_in_list) #SNPs not included in SNPstats file',file,append = TRUE)
        write(paste('','if(length(snp_list_out)!=0){',sep="\n"),file,append = TRUE)
        write(paste('  write(as.character(snp_list_out),"',snp_list_out_file,'", append = FALSE)',sep=""),file,append=TRUE)
        
        system_line_0 = paste('  print(paste("Extracting ",length(snp_list_out)," proxy SNPs in LD from ',ld_file,'",sep=""))',sep="")
        system_line_1 = paste('system("zgrep -w -F -f ',snp_list_out_file,' ',ld_file,' > ',snp_ld_file,'")',sep="")

        write(paste(system_line_0,system_line_1,'  print("Extraction complete")',sep="\n"),file,append = TRUE)
        
        #Note: If proxies are needed, but none are found, then the snp_ld.csv file will be empty.
        #To get around this, make the following code conditional on there being data in the snp_ld.csv file
        write(paste("system('wc -l ",snp_ld_file," > ",snp_ld_test_file,"')",sep=""),file,append = TRUE)
        write(paste("ld_x = read.delim('",snp_ld_test_file,"',header=FALSE,stringsAsFactors = FALSE, strip.white=TRUE,sep=' ')",sep=""),file,append = TRUE)
        write("if(ld_x$V1[1] != 0)  {",file,append = TRUE)
        write(paste("  ld_snps = read.csv('",snp_ld_file,"',header=FALSE,stringsAsFactors = FALSE, strip.white=TRUE)",sep = ""),file,append = TRUE)
        write("  #Remove SNPs with an R2 of less than specified amount",file,append = TRUE)
        write(paste("  ld_snps = ld_snps[which(ld_snps$V5 >= ",r2,"),]",sep=""),file,append = TRUE)
        write("  #Rank the LD SNPs by R2 and remove palindromic with moderate MAF",file,append = TRUE)
        write("ld_snps$V6[ld_snps$V6==TRUE]='T'",file,append=TRUE)
        write("ld_snps$V7[ld_snps$V7==TRUE]='T'",file,append=TRUE)
        write("ld_snps$V8[ld_snps$V8==TRUE]='T'",file,append=TRUE)
        write("ld_snps$V9[ld_snps$V9==TRUE]='T'",file,append=TRUE)
        write(paste("  ","  #Remove proxy SNPs not in SNPstats","  proxies_snps_list = unique(ld_snps$V2)",sep = "\n"),file,append = TRUE)
        write(paste('  write(as.character(proxies_snps_list),"',snp_list_proxies_file,'", append = FALSE)',sep=""),file,append=TRUE)

        l1 = "system('folder=/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats"
        l2 = paste('if [ -f "',snpstats_file,'" ]
then
  echo SNPstats file present, proceeding
else
  echo Concatenating UK Biobank SNP Stats in home directory
  cmd=""
  for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do
    tail -n +17 ${folder}/data.chr${chrom}.snp-stats > data.chr${chrom}.txt
    cut -f 2- data.chr${chrom}.txt > chr${chrom}.txt
    cmd="${cmd} chr${chrom}.txt"
  done
  cat ${cmd} > ',snpstats_file,'
  rm data.chr*
  rm chr*
fi
echo Searching UK Biobank SNP Stats for proxy SNPs from ',ld_file,'
zgrep -w -F -f ',snp_list_proxies_file,' ',snpstats_file,' > ',included_proxy_SNPs_file,'
echo Search complete, ',included_proxy_SNPs_file,' file created, proceeding in R\')',sep="")

        write(paste(l1,l2,sep="\n"),file,append = TRUE)
        write(paste("  proxy_snps_in = read.csv('",included_proxy_SNPs_file,"',header=FALSE,sep='', stringsAsFactors = FALSE, strip.white=TRUE)",sep=""),file,append = TRUE)
        write("  proxy_snps_in$V4[proxy_snps_in$V4==TRUE]='T'",file,append = TRUE)
        write("  proxy_snps_in$V5[proxy_snps_in$V5==TRUE]='T'",file,append = TRUE)
        write("  proxy_snps_in_list = unique(proxy_snps_in$V1)",file,append = TRUE)
        write(paste("  for(snp in unique(ld_snps$V2)){","    if(snp %in% proxy_snps_in_list){",sep="\n"),file,append = TRUE)
        write(paste("      ld_snps$in_snpstat[ld_snps$V2 == snp] = 1","    }  else {","      ld_snps$in_snpstat[ld_snps$V2 == snp] = 0","    }","  }",sep = "\n"),file,append = TRUE)
        write(paste("  ld_snps = ld_snps[which(ld_snps$in_snpstat == 1),]","  for(snp in unique(ld_snps$V1)){",sep="\n"),file,append = TRUE)
        write(paste("    ld_snps$rank[ld_snps$V1 == snp] = rank(-ld_snps$V5[ld_snps$V1 == snp],ties.method = 'first')","  }","  ld_snps = ld_snps[which(ld_snps$rank == 1),]",sep="\n"),file,append = TRUE)
        
        write(paste("","  #Rank SNPs","  exposure_dat$proxy=0","  for(snp in unique(ld_snps$V1)){",sep = "\n"),file,append = TRUE)
        write(paste("    exposure_dat$proxy[exposure_dat$SNP == snp] = 1","  }",sep = "\n"),file,append = TRUE)
        write(paste('  ld_snps = ld_snps[c("V1","V2","V3","V4","V5","V6","V7","V8","V9")]',sep = "\n"),file,append = TRUE)
        write(paste("  exposure_dat=merge(exposure_dat, ld_snps, by.x = 'SNP', by.y = 'V1', all = TRUE)",sep = "\n"),file,append = TRUE)
        write("  exposure_dat=exposure_dat[!is.na(exposure_dat$proxy),]",file,append = TRUE)
        write(paste("  exposure_dat=rename(exposure_dat,proxy_SNP=V2,eaf.outcome=V3,proxy_r2=V5,effect_allele.outcome=V6,other_allele.outcome=V7)",sep = "\n"),file,append = TRUE)
        write(paste(  '#Harmonise exposure_dat to LD file
        
  exposure_dat$effect_allele.exposure = toupper(exposure_dat$effect_allele.exposure)
  exposure_dat$other_allele.exposure = toupper(exposure_dat$other_allele.exposure)
  exposure_dat$effect_allele.outcome = toupper(exposure_dat$effect_allele.outcome)
  exposure_dat$other_allele.outcome = toupper(exposure_dat$other_allele.outcome)
  #And not blank
  exposure_dat$effect_allele.exposure[exposure_dat$effect_allele.exposure == ""] = NA
  exposure_dat$other_allele.exposure[exposure_dat$other_allele.exposure == ""] = NA
  exposure_dat$effect_allele.outcome[exposure_dat$effect_allele.outcome == ""] = NA
  exposure_dat$other_allele.outcome[exposure_dat$other_allele.outcome == ""] = NA
  
  exposure_dat$exposure_alleles = paste(exposure_dat$effect_allele.exposure,exposure_dat$other_allele.exposure, sep = "")
  exposure_dat$outcome_alleles = paste(exposure_dat$effect_allele.outcome,exposure_dat$other_allele.outcome, sep = "")
  
  #Mark palindromes
  exposure_dat$palindrome[(exposure_dat$exposure_alleles == "AT" | exposure_dat$exposure_alleles == "TA" | 
                             exposure_dat$exposure_alleles == "CG" | exposure_dat$exposure_alleles == "GC"|
                             exposure_dat$outcome_alleles == "AT" | exposure_dat$outcome_alleles == "TA" | 
                             exposure_dat$outcome_alleles == "CG" | exposure_dat$outcome_alleles == "GC")] = 1
  exposure_dat$palindrome[is.na(exposure_dat$palindrome)] = 0
  
  #Harmonise(OUTCOME)
  
  exposure_dat$harmonise = NA
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$effect_allele.outcome & 
                           exposure_dat$other_allele.exposure == exposure_dat$other_allele.outcome] = 0
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$other_allele.outcome & 
                           exposure_dat$other_allele.exposure == exposure_dat$effect_allele.outcome] = 1
  #Some missing values for other_allele - for non-palindromic SNPs this is not an issue
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$effect_allele.outcome & 
                           is.na(exposure_dat$other_allele.exposure) & is.na(exposure_dat$palindrome)] = 0
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$other_allele.outcome & 
                           is.na(exposure_dat$other_allele.exposure) & is.na(exposure_dat$palindrome)] = 1
  #If neither "0" nor "1", must be reverse strand, or entire strand is missing (i.e. does not need a proxy)
  exposure_dat$reverse[is.na(exposure_dat$harmonise)] = 1
  exposure_dat$reverse[is.na(exposure_dat$reverse)] = 0
  #Cycle through studies - if the study is all reverse, excluding the palindromes, then reverse whole study
  x = unique(exposure_dat$id.exposure)
  for(study in x){
    continue = length(exposure_dat$reverse[which(exposure_dat$id.exposure == study & exposure_dat$outcome_alleles != "NANA")])
    if(continue > 0){
      mean = mean(as.numeric(exposure_dat$reverse[which(exposure_dat$id.exposure == study & exposure_dat$palindrome != 1 & exposure_dat$outcome_alleles != "NANA")]))
      if(is.na(mean)==TRUE){
        mean = 0.5
      }
      if(mean == 1) {
        exposure_dat$reverse[exposure_dat$id.exposure == study] = 1
      } else if(mean == 0) {
        exposure_dat$reverse[exposure_dat$id.exposure == study] = 0
      } else {
        #Drop any palindromic SNPs with MAF >= maf IF reverse =/= 1 or 0 for all non-palindromic SNPs
        x = length(exposure_dat$SNP)
        exposure_dat = exposure_dat[which(((exposure_dat$palindrome == 0 | exposure_dat$eaf.exposure < ',maf,') & exposure_dat$id.exposure == study) | exposure_dat$id.exposure != study),]
        #Reverse palindromes with a low MAF that need reversing
        #Either the effect alleles are equal but the EAFs are opposite
        exposure_dat$reverse[exposure_dat$palindrome == 1 & exposure_dat$id.exposure == study & 
                               exposure_dat$effect_allele.exposure == exposure_dat$effect_allele.outcome & exposure_dat$eaf.outcome > 0.5] = 1
        #Or the effect alleles are opposit but the EAFs are equal
        exposure_dat$reverse[exposure_dat$palindrome == 1 & exposure_dat$id.exposure == study & 
                               exposure_dat$effect_allele.exposure != exposure_dat$effect_allele.outcome & exposure_dat$eaf.outcome < 0.5] = 1
        y = length(exposure_dat$SNP)
        if(x != y){
          print(paste("Removed ",x-y," palindromic SNPs as strand could not be ascertained",sep=""))
        }
      }
    } else {
      exposure_dat$reverse[exposure_dat$id.exposure == study] = 0
    }
  }
  #Reverse strand for 1000 genomes main SNPs
  exposure_dat$effect_allele.outcome_new = exposure_dat$effect_allele.outcome
  exposure_dat$other_allele.outcome_new = exposure_dat$other_allele.outcome
  exposure_dat$effect_allele.outcome_new[exposure_dat$effect_allele.outcome == "T" & exposure_dat$reverse == 1] = "A"
  exposure_dat$effect_allele.outcome_new[exposure_dat$effect_allele.outcome == "G" & exposure_dat$reverse == 1] = "C"
  exposure_dat$effect_allele.outcome_new[exposure_dat$effect_allele.outcome == "C" & exposure_dat$reverse == 1] = "G"
  exposure_dat$effect_allele.outcome_new[exposure_dat$effect_allele.outcome == "A" & exposure_dat$reverse == 1] = "T"
  exposure_dat$other_allele.outcome_new[exposure_dat$other_allele.outcome == "T" & exposure_dat$reverse == 1] = "A"
  exposure_dat$other_allele.outcome_new[exposure_dat$other_allele.outcome == "G" & exposure_dat$reverse == 1] = "C"
  exposure_dat$other_allele.outcome_new[exposure_dat$other_allele.outcome == "C" & exposure_dat$reverse == 1] = "G"
  exposure_dat$other_allele.outcome_new[exposure_dat$other_allele.outcome == "A" & exposure_dat$reverse == 1] = "T"
  exposure_dat$effect_allele.outcome = exposure_dat$effect_allele.outcome_new
  exposure_dat$other_allele.outcome = exposure_dat$other_allele.outcome_new
  exposure_dat$effect_allele.outcome_new = NULL
  exposure_dat$other_allele.outcome_new = NULL
  
  #Reverse strand for 1000 genomes proxy SNPs
  exposure_dat$V8_new = exposure_dat$V8
  exposure_dat$V9_new = exposure_dat$V9
  exposure_dat$V8_new[exposure_dat$V8 == "T" & exposure_dat$reverse == 1] = "A"
  exposure_dat$V8_new[exposure_dat$V8 == "G" & exposure_dat$reverse == 1] = "C"
  exposure_dat$V8_new[exposure_dat$V8 == "C" & exposure_dat$reverse == 1] = "G"
  exposure_dat$V8_new[exposure_dat$V8 == "A" & exposure_dat$reverse == 1] = "T"
  exposure_dat$V9_new[exposure_dat$V9 == "T" & exposure_dat$reverse == 1] = "A"
  exposure_dat$V9_new[exposure_dat$V9 == "G" & exposure_dat$reverse == 1] = "C"
  exposure_dat$V9_new[exposure_dat$V9 == "C" & exposure_dat$reverse == 1] = "G"
  exposure_dat$V9_new[exposure_dat$V9 == "A" & exposure_dat$reverse == 1] = "T"
  exposure_dat$V8 = exposure_dat$V8_new
  exposure_dat$V9 = exposure_dat$V9_new
  exposure_dat$V8_new = NULL
  exposure_dat$V9_new = NULL
  
  exposure_dat$outcome_alleles = paste(exposure_dat$effect_allele.outcome,exposure_dat$other_allele.outcome, sep = "")
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$effect_allele.outcome & exposure_dat$other_allele.exposure == exposure_dat$other_allele.outcome] = 0
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$other_allele.outcome & exposure_dat$other_allele.exposure == exposure_dat$effect_allele.outcome] = 1
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$effect_allele.outcome & is.na(exposure_dat$other_allele.exposure) & is.na(exposure_dat$palindrome)] = 0
  exposure_dat$harmonise[exposure_dat$effect_allele.exposure == exposure_dat$other_allele.outcome & is.na(exposure_dat$other_allele.exposure) & is.na(exposure_dat$palindrome)] = 1
  #Deal with the "1"s - the effect and other alleles are the wrong way round
  #Main SNPs
  exposure_dat$harmonise[is.na(exposure_dat$harmonise)] = 0
  exposure_dat$temp = exposure_dat$other_allele.outcome
  exposure_dat$other_allele.outcome[exposure_dat$harmonise == 1] = exposure_dat$effect_allele.outcome[exposure_dat$harmonise == 1]
  exposure_dat$effect_allele.outcome[exposure_dat$harmonise == 1] = exposure_dat$temp[exposure_dat$harmonise == 1]
  exposure_dat$eaf.outcome[exposure_dat$harmonise == 1] = 1-exposure_dat$eaf.outcome[exposure_dat$harmonise == 1]
  
  #Proxy SNPs
  exposure_dat$temp = exposure_dat$V9
  exposure_dat$V9[exposure_dat$harmonise == 1] = exposure_dat$V8[exposure_dat$harmonise == 1]
  exposure_dat$V8[exposure_dat$harmonise == 1] = exposure_dat$temp[exposure_dat$harmonise == 1]
  exposure_dat$V4[exposure_dat$harmonise == 1] = 1-exposure_dat$V4[exposure_dat$harmonise == 1]
  
  #Remove useless stuff
  exposure_dat$temp = NULL
  exposure_dat$harmonise = NULL
  exposure_dat$exposure_alleles = NULL
  exposure_dat$outcome_alleles = NULL
  exposure_dat$palindrome = NULL
  exposure_dat$reverse = NULL
  exposure_dat$study_id = NULL',sep=""),file,append = TRUE)
        
        write("  exposure_dat <- data.frame(lapply(exposure_dat, as.character), stringsAsFactors=FALSE)",file,append = TRUE)
        write("  exposure_dat$beta.exposure = as.numeric(exposure_dat$beta.exposure)",file,append = TRUE)
        write("  exposure_dat$eaf.exposure = as.numeric(exposure_dat$eaf.exposure)",file,append = TRUE)
        write("  exposure_dat$V4 = as.numeric(exposure_dat$V4)",file,append = TRUE)
        write("  #Keep relevant proxy info and swap the missing SNPs & alleles for the proxies",file,append = TRUE)
        write("  #Delete anything without a beta - some proxies might not link, this will get rid of them",file,append = TRUE)
        write("  exposure_dat = exposure_dat[which(!is.na(exposure_dat$beta.exposure)),]",file,append = TRUE)
        write("  #Change over SNP names",file,append = TRUE)
        write("  exposure_dat$temp[exposure_dat$proxy == 1] = exposure_dat$SNP[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$SNP[exposure_dat$proxy == 1] = exposure_dat$proxy_SNP[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$SNP.original[exposure_dat$proxy == 1] = exposure_dat$temp[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  #Change over effect alleles",file,append = TRUE)
        write("  exposure_dat$temp[exposure_dat$proxy == 1] = exposure_dat$effect_allele.exposure[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$effect_allele.exposure[exposure_dat$proxy == 1] = exposure_dat$V8[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$effect_allele.original[exposure_dat$proxy == 1] = exposure_dat$temp[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  #Change over other alleles",file,append = TRUE)
        write("  exposure_dat$temp[exposure_dat$proxy == 1] = exposure_dat$other_allele.exposure[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$other_allele.exposure[exposure_dat$proxy == 1] = exposure_dat$V9[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$other_allele.original[exposure_dat$proxy == 1] = exposure_dat$temp[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  #Change over effect allele frequencies (should be similar)",file,append = TRUE)
        write("  exposure_dat$temp[exposure_dat$proxy == 1] = exposure_dat$eaf.exposure[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$eaf.exposure[exposure_dat$proxy == 1] = exposure_dat$V4[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  exposure_dat$eaf.original[exposure_dat$proxy == 1] = exposure_dat$temp[exposure_dat$proxy == 1]",file,append = TRUE)
        write("  #Remove variables",file,append = TRUE)
        write(paste("  exposure_dat$eaf.outcome = NULL","  exposure_dat$effect_allele.outcome = NULL","  exposure_dat$other_allele.outcome = NULL",sep="\n"),file,append = TRUE)
        write("  #Overwrite SNP list with new SNPs",file,append = TRUE)
        write("  snp_list = unique(exposure_dat$SNP)",file,append = TRUE)
        write(paste('  write.table(snp_list,"',snp_list_file,'",row.names=FALSE,quote=FALSE,col.names = FALSE)',sep=""),file,append = TRUE)
        write(paste('system("zgrep -w -F -f ',snp_list_file,' ',snpstats_file,' > ',included_SNPs_file,'")',sep=''),file,append=TRUE)
        write(paste('snp_list_in = read.csv("',included_SNPs_file,'",header=FALSE,sep="",stringsAsFactors = FALSE, strip.white=TRUE)',sep=""),file,append = TRUE)
        write('  print("All possible proxies acquired")',file,append = TRUE)
        write("} else {",file,append=TRUE)
        write("print('No proxies found, continuing')",file,append=TRUE)
        write("}",file,append=TRUE)
        write(paste('} else {',"  print('All SNPs in UK Biobank, no proxies necessary')","}",sep = "\n"),file,append = TRUE)
      }
      write('print("Harmonising MR-Base and UK-Biobank")',file,append = TRUE)
      write("x = length(exposure_dat$SNP)",file,append=TRUE)
      write("exposure_dat = exposure_dat[which(exposure_dat$SNP %ni% snp_list_in$V1[duplicated(snp_list_in$V1)]),]",file,append = TRUE)
      write(paste("y = length(exposure_dat$SNP)","if(x!=y){","  print(paste(x-y,' triallelic proxy SNPs dropped',sep=''))","}",sep='\n'),file,append=TRUE)
      write("snp_list_in = snp_list_in[snp_list_in$V1 %ni% snp_list_in$V1[duplicated(snp_list_in$V1)],]",file,append = TRUE)
      #UKBiobank SNPstats have allele B as the effect allele when passed through Plink2
      write('snp_list_in = snp_list_in[c("V1","V12","V5","V4")]',file,append = TRUE)
      write('snp_list_in = rename(snp_list_in,eaf.outcome = V12, effect_allele.outcome = V5, other_allele.outcome = V4, SNP = V1)',file,append = TRUE)
      #Merge with exposure_dat
      write("exposure_dat_harmonised = merge(exposure_dat, snp_list_in, by.x = 'SNP', by.y = 'SNP', all = TRUE)",file,append = TRUE)
      write("exposure_dat_harmonised = exposure_dat_harmonised[which(!is.na(exposure_dat_harmonised$beta.exposure)),]",file,append = TRUE)
      #Harmonise
      write('source("harmonise.R")',file,append = TRUE)
      write("exposure_dat_harmonised = harmonise(dat = exposure_dat_harmonised)",file,append = TRUE)
      write(paste("exposure_dat_harmonised$effect_allele.outcome = NULL","exposure_dat_harmonised$other_allele.outcome = NULL",
                  "exposure_dat_harmonised$eaf.outcome = NULL","exposure_dat_harmonised$V4 = NULL","exposure_dat_harmonised$V8 = NULL",
                  "exposure_dat_harmonised$V9 = NULL","exposure_dat_harmonised$rank = NULL","exposure_dat_harmonised$proxy_SNP = NULL",
                  sep = "\n"),file,append = TRUE)
      write(paste('write.csv(exposure_dat_harmonised,"',exposure_dat_harmonised_file,'")',sep=""),file,append = TRUE)
      
      #GRS code creation code
      #Best to start fresh, so drop everything and re-import exposure_dat_harmonised as dat
      write(paste("",'remove(list=ls())',sep="\n"),file,append = TRUE)
      write(paste("dat = read.csv('",exposure_dat_harmonised_file,"', stringsAsFactors = FALSE, strip.white=TRUE)",sep=""),file,append = TRUE)
      write("dat$effect_allele.exposure[dat$effect_allele.exposure == TRUE]='T'",file,append=TRUE)
      write("dat$other_allele.exposure[dat$other_allele.exposure == TRUE]='T'",file,append=TRUE)
      write("snp_list = unique(dat$SNP)",file,append = TRUE)
      write(paste('write.table(snp_list,"',snp_list_file,'",row.names=FALSE,quote=FALSE,col.names = FALSE)',sep=""),file,append = TRUE)
      
      #Perform Bgenix and Plink2 commands using system()
      #The majority of cases will create the GRS in R
      
      if(plink_grs == FALSE){
        write("system('",file,append=TRUE)
        write('# Bgen patterns',file,append=TRUE)
        write(paste('bgen_pattern=',bgen_folder,'/data.chrCHROM.bgen',sep=""),file,append=TRUE)
        write(paste('bgen_index_pattern=',bgen_folder,'/data.chrCHROM.bgen.bgi',sep=""),file,append=TRUE)
        write('# Args',file,append=TRUE)
        write(paste('snp_list=',snp_list_file,sep=""),file,append=TRUE)
        write('echo Extracting variants using bgenix',file,append=TRUE)
        write(paste('temp_geno_prefix=',temp_geno_prefix,sep=""),file,append=TRUE)
        write('for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do',file,append=TRUE)
        write('  inbgen=${bgen_pattern/CHROM/$chrom}',file,append=TRUE)
        write('  inbgenidx=${bgen_index_pattern/CHROM/$chrom}',file,append=TRUE)
        write('  bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list > ${temp_geno_prefix}.${chrom}.bgen',file,append=TRUE)
        write('done',file,append=TRUE)
        write('cmd=""',file,append=TRUE)
        write('for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do',file,append=TRUE)
        write('  cmd="${cmd} ${temp_geno_prefix}.${chrom}.bgen"',file,append=TRUE)
        write('done',file,append=TRUE)
        write(paste('cat-bgen -g ${cmd} -og ',instruments_file,'.bgen',sep=""),file,append=TRUE)
        write('# Remove temp genos',file,append=TRUE)
        write('rm $temp_geno_prefix*',file,append=TRUE)
        write('echo Creating dosages using Plink2',file,append=TRUE)
        write(paste('plink2 --bgen ',instruments_file,'.bgen --hard-call-threshold 0.4999 --export A --out ',instruments_file,sep=""),file,append=TRUE)
        write(paste('rm ',instruments_file,'.bgen',sep=""),file,append=TRUE)
        write("')",file,append=TRUE)
        
        #instruments.raw has now been created, so proceed with generating GRS
  
        if(ipd == TRUE) {
          write(paste("print('",instruments_file,".raw created, reading into R and outputting ",snp_ipd_file,"')",sep=""),file,append=TRUE)
        } else {
          write(paste("print('",instruments_file,".raw created, reading into R')",sep=""),file,append=TRUE)
        }
        
        write(paste('snp_ipd = as.data.frame(fread("',instruments_file,'.raw", header = T, sep = "\t"))',sep=""),file,append=TRUE)
        write("print('Instruments read into R, creating GRS')",file,append=TRUE)
        if(ipd==TRUE) {
          write(paste('write.csv(snp_ipd,"',snp_ipd_file,'",row.names=FALSE)',sep=""),file,append=TRUE)
        }
      
        write('#Check whether the dat file has a "trait" variable - if so, create a trait list, otherwise trait = trait
if("trait" %in% names(dat) == TRUE){
  dat$trait[is.na(dat$trait)] = "missing"
  trait = unique(dat$trait)
} else if("Trait" %in% names(dat) == TRUE){
  dat$trait[is.na(dat$Trait)] = "missing"
  trait = unique(dat$Trait)
  dat = rename(dat,trait = Trait)
} else if("exposure.trait" %in% names(dat) == TRUE){
  dat$exposure.trait[is.na(dat$exposure.trait)] = "missing"
  trait = unique(dat$exposure.trait)
  dat=rename(dat,trait = exposure.trait)  
} else if("exposure.Trait" %in% names(dat) == TRUE){
  dat=rename(dat,trait = exposure.Trait) 
  dat$exposure.trait[is.na(dat$exposure.trait)] = "missing"
  trait = unique(dat$exposure.trait)
} else {
  trait = "trait"
  dat$trait = "trait"
}',file,append=TRUE)
        
        write('trait = gsub(" ","_",trait)',file,append=TRUE)
        write(paste('trait = gsub("',"'",'","",trait)',sep=""),file,append=TRUE)
        write(paste("trait = gsub('",'"',"','',trait)",sep=""),file,append=TRUE)
        write('dat$trait = gsub(" ","_",dat$trait)',file,append=TRUE)
        write(paste('dat$trait = gsub("',"'",'","",dat$trait)',sep=""),file,append=TRUE)
        write(paste("dat$trait = gsub('",'"',"','',dat$trait)",sep=""),file,append=TRUE)
        write("grs = data.frame(id = snp_ipd$FID)",file,append=TRUE)
        write("dat$included = 0",file,append=TRUE)
        write('for(t in trait){
  #Take all SNPs remaining in dat
  snp_list = as.vector(dat$SNP[!is.na(dat$beta.exposure) & dat$"trait" == t])
  #Take all betas from the EXPOSURE file
  effect_list = as.vector(dat$beta.exposure[!is.na(dat$beta.exposure) & dat$"trait" == t])
  #Also need the effect alleles
  effect_allele_list = as.vector(dat$effect_allele.exposure[!is.na(dat$effect_allele.exposure) & dat$"trait" == t])
  
  grs[,t] = 0
  i = 1
  for(snp in snp_list){
    snp2 = paste(snp,"_",effect_allele_list[i],sep="")
    if(snp2 %in% names(snp_ipd) == TRUE) {
      dat$included[dat$SNP == snp & dat$trait == t] = 1
      grs[,t] = grs[,t] + snp_ipd[,snp2]*as.numeric(effect_list[i])
    }
    i=i+1
  }
}',file,append=TRUE)
        write(paste('write.csv(dat,"',exposure_dat_harmonised_file,'",row.names=FALSE)',sep=""),file,append=TRUE)
        write(paste('write.csv(grs,"',grs_file,'",row.names=FALSE)',sep=""),file,append=TRUE)
      } else {
        
        #But in some cases, it will be easier to create the GRS using plink
        #The SNP list has been created, add alleles codes and weights to the next 2 columns
        write("snp_score_list = dat[,c('SNP','effect_allele.exposure','beta.exposure')]",file,append = TRUE)
        write(paste('write.table(snp_score_list,"',snp_score_list_file,'",row.names=FALSE,quote=FALSE,col.names = FALSE)',sep=""),file,append = TRUE)
        
        #Now need to use plink 2.00 (for reasons I don't understand...)
        #Use plink 2.00 to create a GRS per chromosome, then combine at the end
        write("system('#Bgen patterns",file,append=TRUE)
        write(paste('bgen_pattern=',bgen_folder,'/data.chrCHROM.bgen',sep=""),file,append=TRUE)
        write(paste('bgen_index_pattern=',bgen_folder,'/data.chrCHROM.bgen.bgi',sep=""),file,append=TRUE)
        write('# Args',file,append=TRUE)
        write(paste('snp_list=',snp_list_file,sep=""),file,append=TRUE)
        write('echo Extracting variants using bgenix',file,append=TRUE)
        write(paste('temp_geno_prefix=',temp_geno_prefix,sep=""),file,append=TRUE)
        write('for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do',file,append=TRUE)
        write('  chrom_padd=$(printf "%0*d\\n" 2 $chrom)',file,append=TRUE)
        write('  inbgen=${bgen_pattern/CHROM/$chrom_padd}',file,append=TRUE)
        write('  inbgenidx=${bgen_index_pattern/CHROM/$chrom_padd}',file,append=TRUE)
        write('  echo "#FID	IID	NMISS_ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG" > ${temp_geno_prefix}.${chrom_padd}.sscore',file,append=TRUE)
        write('  bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list > $temp_geno_prefix.$chrom_padd.bgen',file,append=TRUE)
        write(paste('  plink --bgen $temp_geno_prefix.$chrom_padd.bgen --score ',snp_score_list_file,' --out $temp_geno_prefix.$chrom_padd',sep=""),file,append=TRUE)
        write('  rm ${temp_geno_prefix}*.bgen',file,append=TRUE)
        write("done')",file,append=TRUE)
        #Note: this may screw up here if there are 0 SNPs on the 1st chromosome. Probably unlikely, but if grs.csv is created with no IDs, probably this.
        write("grs = data.frame()",file,append=TRUE)
        write(paste('temp = read.delim("',temp_geno_prefix,'.01.sscore",sep="\t")',sep=""),file,append=TRUE)
        write("grs = data.frame(id=temp[,c('IID')])",file,append=TRUE)
        write("grs$grs = temp[,3]*temp[,c('SCORE1_AVG')]",file,append=TRUE)
        write("for(i in c(01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,X)) {",file,append=TRUE)
        write(paste('  temp = read.delim(paste("',temp_geno_prefix,'.",i,".sscore",sep=""),sep="\t")',sep=""),file,append=TRUE)
        write("  if(length(unique(temp$IID))>0) {",file,append=TRUE)
        write("    merge = data.frame(id=temp[,c('IID')])",file,append=TRUE)
        write("    merge$grs_temp = temp[,3]*temp[,c('SCORE1_AVG')]",file,append=TRUE)
        write("    grs= merge(grs, merge, by.x='id', by.y='id')",file,append=TRUE)
        write("    grs$grs = grs$grs+grs$grs_temp",file,append=TRUE)
        write("    grs$grs_temp = NULL",file,append=TRUE)
        write("  }",file,append=TRUE)
        write("}",file,append=TRUE)
        write(paste('write.csv(grs,"',grs_file,'",row.names=FALSE)',sep=""),file,append=TRUE)
      }
      write(paste("print('Genetic risk scores created, see ",grs_file,"')",sep=""),file,append=TRUE)
      write("",file,append=TRUE)
      if(keep_files==FALSE){
        instruments_log_file = paste(instruments_file,".log",sep="")
        instruments_file = paste(instruments_file,".raw",sep="")
        write(paste("system('rm",snp_list_file,included_SNPs_file,"')",sep=" "),file,append=TRUE)
        write(paste('system(\'if [ -f "',snp_list_out_file,'" ]
then
  rm ',snp_list_out_file,'
fi
if [ -f "',snp_list_proxies_file,'" ]
then
  rm ',snp_list_proxies_file,'
fi
if [ -f "',snp_ld_test_file,'" ]
then
rm ',snp_ld_test_file,'
fi
if [ -f "',snp_score_list_file,'" ]
then
  rm ',snp_score_list_file,'
fi
if [ -f "',instruments_file,'" ]
then
  rm ',instruments_file,'
fi
if [ -f "',snp_ld_file,'" ]
then
  rm ',snp_ld_file,'
fi
if [ -f "',included_proxy_SNPs_file,'" ]
then
  rm ',included_proxy_SNPs_file,'
fi
if [ -f "clumped_snps.txt" ]
then
  rm  *clumped_snps*
fi
if [ -f "',temp_geno_prefix,'.01.sscore" ]
then
  rm  *sscore
fi
if [ -f "',temp_geno_prefix,'.01.log" ]
then
  rm  *log
fi
if [ -f "',instruments_log_file,'" ]
then
  rm ',instruments_log_file,'
fi\')',sep=""),file,append=TRUE)        
      }

      write("print('Remember all IDs are included, so remove any participants you do not want in the analysis')",file,append=TRUE)
      
      
      ##############################################################################################################
      
      #run.sh
      file = run_file
      write("#!/bin/bash",file,append = FALSE)
      write("#PBS -l nodes=1:ppn=16",file,append=TRUE)
      write("#PBS -l walltime=00:12:00:00",file,append=TRUE)
      write("#PBS -N GRS_creation",file,append=TRUE)
      write("cd $PBS_O_WORKDIR",file,append=TRUE)
      write("",file,append=TRUE)
      write("module load apps/bgen-1.1.4",file,append=TRUE)
      write("module load apps/plink-2.00",file,append=TRUE)
      write("module load languages/R-3.5-ATLAS-gcc-7.1.0",file,append=TRUE)
      write("",file,append=TRUE)
      write(paste("./",script_file,sep=""),file,append=TRUE)

      ##############################################################################################################
      
      #harmonise.R
      #harmonise.R [Create whenever code is written]
      file = "harmonise.R"
      write("#Harmonisation of SNPs",file,append = FALSE)
      
      write(paste('harmonise = function(dat = "dat", exposure_dat = "exposure_dat", outcome_dat = "outcome_dat", study_id = "id.exposure", merge = FALSE, maf = ',maf,') {
  #Merge if specified
  if(merge == "TRUE") {
    print("Merging dataframes")
    dat = merge(exposure_dat, outcome_dat, by.x = "SNP", by.y = "SNP")
  }
  #Make sure all alleles are capitalised
  dat$effect_allele.exposure = toupper(dat$effect_allele.exposure)
  dat$other_allele.exposure = toupper(dat$other_allele.exposure)
  dat$effect_allele.outcome = toupper(dat$effect_allele.outcome)
  dat$other_allele.outcome = toupper(dat$other_allele.outcome)
  #And not blank
  dat$effect_allele.exposure[dat$effect_allele.exposure == ""] = NA
  dat$other_allele.exposure[dat$other_allele.exposure == ""] = NA
  dat$effect_allele.outcome[dat$effect_allele.outcome == ""] = NA
  dat$other_allele.outcome[dat$other_allele.outcome == ""] = NA

  dat$exposure_alleles = paste(dat$effect_allele.exposure,dat$other_allele.exposure, sep = "")
  dat$outcome_alleles = paste(dat$effect_allele.outcome,dat$other_allele.outcome, sep = "")
  
  #Mark palindromes
  dat$palindrome[(dat$exposure_alleles == "AT" | dat$exposure_alleles == "TA" | dat$exposure_alleles == "CG" | dat$exposure_alleles == "GC"|
  dat$outcome_alleles == "AT" | dat$outcome_alleles == "TA" | dat$outcome_alleles == "CG" | dat$outcome_alleles == "GC")] = 1
  dat$palindrome[is.na(dat$palindrome)] = 0
  
  #Harmonise
  dat$harmonise = NA
  dat$harmonise[dat$effect_allele.exposure == dat$effect_allele.outcome & dat$other_allele.exposure == dat$other_allele.outcome] = 0
  dat$harmonise[dat$effect_allele.exposure == dat$other_allele.outcome & dat$other_allele.exposure == dat$effect_allele.outcome] = 1
  #Some missing values for other_allele - for non-palindromic SNPs this is not an issue
  dat$harmonise[dat$effect_allele.exposure == dat$effect_allele.outcome & is.na(dat$other_allele.exposure) & is.na(dat$palindrome)] = 0
  dat$harmonise[dat$effect_allele.exposure == dat$other_allele.outcome & is.na(dat$other_allele.exposure) & is.na(dat$palindrome)] = 1
  #If neither "0" nor "1", must be reverse strand
  dat$reverse[is.na(dat$harmonise)] = 1
  dat$reverse[is.na(dat$reverse)] = 0
  #Cycle through studies - if the study is all reverse, excluding the palindromes, then reverse whole study
  dat$"study_id" = dat[,as.character(study_id)]
  x = unique(dat$"study_id")
  for(study in x){
    mean = mean(as.numeric(dat$reverse[which(dat$study_id == study & dat$palindrome != 1 & dat$outcome_alleles != "NANA")]))
    if(is.na(mean) == TRUE){
      mean = 0.5
    }
    if(mean == 1) {
      dat$reverse[dat$study_id == study] = 1
    } else if(mean == 0) {
      dat$reverse[dat$study_id == study] = 0
    } else {
      #Drop any palindromic SNPs with MAF >= maf IF reverse =/= 1 or 0 for all non-palindromic SNPs
      x = length(dat$SNP)
      dat = dat[which(((dat$palindrome == 0 | dat$eaf.exposure < maf | 1-dat$eaf.exposure < maf) & dat$study_id == study)| dat$study_id != study),]
      #Reverse palindromes with a low MAF that need reversing
      #Either the effect alleles are equal but the EAFs are opposite
      dat$reverse[dat$palindrome == 1 & dat$study_id == study & dat$effect_allele.exposure == dat$effect_allele.outcome & dat$eaf.outcome > 0.5] = 1
      #Or the effect alleles are opposite but the EAFs are equal
      dat$reverse[dat$palindrome == 1 & dat$study_id == study & dat$effect_allele.exposure != dat$effect_allele.outcome & dat$eaf.outcome < 0.5] = 1
      y = length(dat$SNP)
      if(x != y){
        print(paste("Removed ",x-y," palindromic SNPs as strand could not be ascertained",sep=""))
      }
    }
  }
  dat$effect_allele.exposure_new = dat$effect_allele.exposure
  dat$other_allele.exposure_new = dat$other_allele.exposure
  dat$effect_allele.exposure_new[dat$effect_allele.exposure == "T" & dat$reverse == 1] = "A"
  dat$effect_allele.exposure_new[dat$effect_allele.exposure == "G" & dat$reverse == 1] = "C"
  dat$effect_allele.exposure_new[dat$effect_allele.exposure == "C" & dat$reverse == 1] = "G"
  dat$effect_allele.exposure_new[dat$effect_allele.exposure == "A" & dat$reverse == 1] = "T"
  dat$other_allele.exposure_new[dat$other_allele.exposure == "T" & dat$reverse == 1] = "A"
  dat$other_allele.exposure_new[dat$other_allele.exposure == "G" & dat$reverse == 1] = "C"
  dat$other_allele.exposure_new[dat$other_allele.exposure == "C" & dat$reverse == 1] = "G"
  dat$other_allele.exposure_new[dat$other_allele.exposure == "A" & dat$reverse == 1] = "T"
  dat$effect_allele.exposure = dat$effect_allele.exposure_new
  dat$other_allele.exposure = dat$other_allele.exposure_new
  dat$effect_allele.exposure_new = NULL
  dat$other_allele.exposure_new = NULL
  #Retry this
  dat$exposure_alleles = paste(dat$effect_allele.exposure,dat$other_allele.exposure, sep = "")
  dat$harmonise[dat$effect_allele.exposure == dat$effect_allele.outcome & dat$other_allele.exposure == dat$other_allele.outcome] = 0
  dat$harmonise[dat$effect_allele.exposure == dat$other_allele.outcome & dat$other_allele.exposure == dat$effect_allele.outcome] = 1
  dat$harmonise[dat$effect_allele.exposure == dat$effect_allele.outcome & is.na(dat$other_allele.exposure) & is.na(dat$palindrome)] = 0
  dat$harmonise[dat$effect_allele.exposure == dat$other_allele.outcome & is.na(dat$other_allele.exposure) & is.na(dat$palindrome)] = 1
  #Deal with the "1"s - the effect and other alleles are the wrong way round
  dat$harmonise[is.na(dat$harmonise)] = 0
  dat$temp = dat$other_allele.exposure
  dat$other_allele.exposure[dat$harmonise == 1] = dat$effect_allele.exposure[dat$harmonise == 1]
  dat$effect_allele.exposure[dat$harmonise == 1] = dat$temp[dat$harmonise == 1]
  dat$beta.exposure[dat$harmonise == 1] = -dat$beta.exposure[dat$harmonise == 1]
  dat$eaf.exposure[dat$harmonise == 1] = 1-dat$eaf.exposure[dat$harmonise == 1]
  dat$temp = NULL
  dat$harmonise = NULL
  dat$exposure_alleles = NULL
  dat$outcome_alleles = NULL
  dat$palindrome = NULL
  dat$reverse = NULL
  dat$study_id = NULL
  return(dat)
    }',sep=""),file,append=TRUE)
      
      ##############################################################################################################
      
      print("All code and files produced, proceed to Blue Crystal")
      print(paste("Remember to set the ",run_file," and ",script_file," files to executable (chmod +x ...)",sep=""))
      print(paste("Also remember to convert the ",run_file," and ",script_file," files to unix (dos2unix ...)",sep=""))
      
      #Also, return the dat (if exposure_dat NOT specified)
      if(return_dat==TRUE){
        return(exposure_dat)
      }
    }
 
  }
  
  else {
    warning("The value for 'output' you entered is not recognised.
            Please enter a value of: subcategories; traits; studies; snps; code; grs")
  }
  
}
