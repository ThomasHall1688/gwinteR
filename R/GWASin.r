#' A Function for GWAS integration
#' This package allows you to integrate a gene list of interest with large GWAS datasets. 
#' @param GWAS_file. GWAS .csv file containing single nucleotide polymorphisms. First five rows must be chromosome, posiition of SNP, SNP ID, P values and Q values.
#' @param gene_lists. Genes of interest that will be intergrated with the GWAS dataset in .csv format. First three rows must be chromosome, start position and end position of gene.
#' @param output_dir. Output directory to create SNP and GWAS integration output directories. Defaults to current directory
#' @param start_pos. Up and downstream maximum starting distance to extract proximal SNPs. Defaults to 10000. 
#' @param interval. Amount the search distance will increase per jump. Defaults to 10000
#' @param interval_jumps. Amount of times the distance will increase by interval value. Defaults to 5
#' @param permutations. Amount of sets of random SNPs equal the observed set that will be used to create the null distribution. Defaults to 1000
#' @keywords GWAS_file, gene_lists, output_dir, start_pos, interval, interval_jumps, permutations
#' @description In summary, this package determines the significance of the distribution of SNP nominal P values within and up to (start_pos + (interval X interval_jumps))kb up and downstream to each gene of interest suspected to have an effect on GWAS phenotype. A total of (permutations) random samplings (with replacement) from the a GWAS P value data set  representing the size of each of selected SNP subsets will be generated. The q values for each SNP P value subset and all its permuted equivalents are calculated using the qvalue library in R, a package incorporated into the GWASin library. The subsequent significance level (Pperm) assigned to each of the SNP subsets is equivalent to the proportion of permutations in which at least the same number of q values < 0.05/0.1 as the SNP subset were obtained, i.e. by chance.
#' @return description The package will create two directories and output multiple datasets : Extracted_SNPS which will contain all the SNPs extracted from each gene in single column text files. The second directory will be Integration_output,and will contain a text file containing all extracted SNPs with the orginal pvalues and qvalues along side the new permuted qvalues. It will also contain histograms plotting the distributions of SNPs with qvalues < 0.05/0.1.Finally, a summary text file will contain information about each SNP file analysed, containing the file name, total SNP amount,amount of SNPs with an original qvalue < 0.05, amount of SNPs with a new qvalue < 0.05, the probability these newqvalues are by chance (0 - 100%), average number of SNPs per (permutations) set that had a new qvalue < 0.05, amount of SNPs with an original qvalue < 0.1, amount of SNPs with a new qvalue < 0.1, the probability these newqvalues are by chance (0 - 100%), average number of SNPs per (permutations) set that had a new qvalue < 0.1. 
#' @imports qvalue
#' @export
#' @examples 
#' GWAS_integration(GWAS_file = "/home/workspace/GWAS.csv", gene_lists = "/home/workspace/gene_lists/", output_dir = "/home/workspace/output/", start_pos = 1000, interval = 2000, interval_jumps = 3, permutations = 100)	
#' GWAS_integration(GWAS_file = "./GWAS.csv", gene_lists = "./")

GWAS_integration <- function(GWAS_file, gene_lists, output_dir, start_pos, interval, interval_jumps, permutations){
#defaults
if(missing(output_dir)){
	output_dir <-  get(wd)
	}

if(missing(start_pos)){
	start_pos <-  10000
	}
	
if(missing(interval)){
	interval <-  10000
	}
	
if(missing(interval_jumps)){
	interval_jumps <-  5
	}

if(missing(permutations)){
	permutations <-  1000
	}


#making directories 
SNPpath <- paste0(output_dir, "/Extracted_SNPS")
if (!file.exists(SNPpath)){
  dir.create(SNPpath)
}

GWASPath <- paste0(output_dir, "/Integration_output")
if (!file.exists(GWASPath)){
  dir.create(GWASPath)
}



sumRes <- paste0(GWASPath, "/Summary_results.txt")
M <- paste("Reading GWAS dataset", GWAS_file, "\n\n")
cat(M)
SNP_list = read.csv(GWAS_file)
SNP_list[, 2] <- as.integer(SNP_list[, 2]) #Sometimes the position will import as a factor. This will cause errors. 
SNP_list$X <- NULL #read.csv can import row.names as col1. Will work on a fix
all_qval <- SNP_list[, c(3,4,5)] #import the relevant columns by which to compare new permuted qvalues (position, pvalue and qvalue)
colnames(all_qval) <- c("SNP", "P_value", "Original qval") #label as such
Gene_SNPs <-  data.frame(NULL) #create empty frame to populated with SNPs
distance = start_pos #This value will increase for ever f in interval_jump value by interval value
i = 1
f = 1

#A user might have the GWAS in the same folder as the gene lists. This will take the name of the GWAS file and exlude it on line XXX
#It also accounts for the fact the path may be on a windows machine. 
if(.Platform$OS.type == "windows"){ 
GWAS_file_parts <- strsplit(GWAS_file, split="\\\\")[[1]]
} else {
GWAS_file_parts <- strsplit(GWAS_file, split="/")[[1]]
}
list.files(gene_lists)
gene_list <- list.files(gene_lists)
File_GWAS <- gene_list[grepl(gene_list, pattern=GWAS_file_parts[length(GWAS_file_parts)]) == TRUE]
gene_list <- gene_list[grepl(gene_list, pattern=GWAS_file_parts[length(GWAS_file_parts)]) == FALSE]
gene_list <- gene_list[!grepl(gene_list, pattern="Extracted_SNPS")] 
gene_list <- gene_list[!grepl(gene_list, pattern="Integration_output")]
setwd(gene_lists)

N_status <- paste("Now exctracting SNPs, this may take some time", "\n\n")
cat(N_status)
	
for(k in gene_list){
    Genes = read.csv(k) #read file from current directory
    Genes[,2] <- as.numeric(as.character(Genes[,2]))
    Genes[,3] <- as.numeric(as.character(Genes[,3]))
    for (f in 1:interval_jumps){
		Extraction_status = paste("Extracting SNPs", distance, "bases from gene list", k, "\n")
		cat(Extraction_status)
        for(i in 1:nrow(Genes)){
            chr = toString(Genes[i,1]) #list and create chromosome vector from gene file
            common = SNP_list[which(SNP_list[, 1] == chr),] #compare chromsome in SNP list to gene file
            snps = common[which(common[, 2] >= (Genes[i,2] - distance)),] #extract SNPS from SNP_list that lie within the gene and 'distance' upstream
            snp = snps[which(snps[, 2] <= (Genes[i,3] + distance)),] #extract SNPS from SNP_list that lie within the gene and 'distance' downstream
            Gene_SNPs = rbind(Gene_SNPs, snp)
            i = i + 1
        }
        Gene_SNPs = unique(Gene_SNPs) #remove duplicates
        Gene_SNPs = Gene_SNPs[order(Gene_SNPs[2]),] #order SNPs by position
		name=paste(c(k, "_", File_GWAS, "_SNP_set_", distance, "_bases.txt"), collapse="") 
		write.table(Gene_SNPs[,3], file = file.path(SNPpath, name), quote = F, row.names = F, col.names = F)
		Gene_SNPs = data.frame(NULL)
		distance = distance + interval #increase search space 
		f = f + 1
		#setwd(working_dir)
	}
distance = start_pos #reset starting point for search space with new gene file. 
}

setwd(output_dir)
	
cat("\nSNP extraction has finished. New qvalues will now be calculated\n \n")

#install.packages("fdrtool")
#library(fdrtool)


library("qvalue")

list.files(SNPpath)
SNP_lists=list.files(SNPpath)

header=paste(c("File ","SNPS_total ","Orig_q<0.05 ","new_q<0.05 ","%prob1000x ","Av_num_qval0.05_per_samp ","Orig_q<0.1 ","new_q<0.1 ","%prob1000x ","Av_num_qval0.1_per_samp"), collapse="")
write.table(header, file=sumRes, sep="\t", quote = F, append = TRUE)

for(i in SNP_lists){
  setwd(SNPpath)
  data=read.table(i, header=FALSE, sep="") #This is the file "gene_snps_100000_ChIP". This is so it will parse through multiple SNP files you make. 
  colnames(data)=c("SNP")
  number_snps=length(data$SNP)
  data2=na.omit(data)
  number_snps2=length(data2$SNP)

  data3=merge(data2,all_qval, by="SNP", all.x=TRUE)
  final_number_snps=length(data2$SNP)
  summary(data2)
  final_number_snps
  name=paste(c(i," n=",final_number_snps,"_original_pval_qval.txt"), collapse="")


  write.table(data3, file = file.path(GWASPath, name), sep="\t", row.names=FALSE, quote = F)


  qval_sig_0.05=subset(data3,data3$`Original qval`<0.05)
  num_sig_0.05_orig=length(qval_sig_0.05$SNP)
  qval_sig_0.1=subset(data3,data3$`Original qval`<0.1)
  num_sig_0.1_orig=length(qval_sig_0.1$SNP)

  data4=merge(data,all_qval, by="SNP")
  number_snps_left=length(data4$SNP)


  pvalue=data4$P_value

  fdr=qvalue(pvalue,fdr.level = NULL, pi0 = 1)

  qval=as.numeric(fdr$qvalues)


  qval2=qval

  num_qval_0.05=length(subset(qval2,qval2<0.05))
  num_qval_0.1=length(subset(qval2,qval2<0.1))

  data5=cbind(data4, qval)

  name=paste(c(i," n=",number_snps_left,"_new_qval.txt"), collapse="")

  write.table(data5, file = file.path(GWASPath, name), sep="\t", row.names=FALSE, quote = F)
  
####### permutations #######

  data_p=all_qval$P_value

  y=replicate(permutations,(sample(data_p, number_snps_left, replace=TRUE)))

  results<-vector("list", permutations)



  for(a in 1:permutations){
    results[[a]]<-qvalue(y[,a], fdr.level = NULL, pi0 = 1)["qval"]
  }





##############################################################################################

## count the number of qvals per sampling that are <0.05

  results2<-vector("list", permutations)
  results3<-vector("list", permutations)

  for(b in 1:permutations){
    results2[[b]]<-length(subset(do.call(c, unlist(results[b], recursive=FALSE)), do.call(c, unlist(results[b], recursive=FALSE))<0.05))
  }


  qvals_0.05=do.call(c, results2)

  for(d in 1:permutations){
    results3[[d]]<-length(subset(do.call(c, unlist(results[d], recursive=FALSE)), do.call(c, unlist(results[d], recursive=FALSE))<0.1))
  }


  qvals_0.1=do.call(c, results3)

  greater_num_qval_0.05=subset(qvals_0.05,qvals_0.05>=num_qval_0.05)
  greater_num_qval_0.1=subset(qvals_0.1,qvals_0.1>=num_qval_0.1)

  number_samplings_0.05_greater=length(greater_num_qval_0.05)
  number_samplings_0.1_greater=length(greater_num_qval_0.1)

  percent_0.05=number_samplings_0.05_greater/10
  percent_0.1=number_samplings_0.1_greater/10

  total_num_qval_0.05=sum(qvals_0.05)
  total_num_qval_0.1=sum(qvals_0.1)

  average_num_qvals0.05_sampling=total_num_qval_0.05/permutations
  average_num_qvals0.1_sampling=total_num_qval_0.1/permutations

  name_0.05=paste(c(i,"_n=", number_snps_left,"_snps_qvals_0.05=",total_num_qval_0.05,"_1000x_",percent_0.05,"pc_",num_qval_0.05,"_or_more.pdf"), collapse="")
  main_0.05=paste(c(i," n=", number_snps_left," 1000x"), collapse="")

  xlab_0.05="# q values per sampling <0.05"

  name_0.1=paste(c(i,"_n=", number_snps_left,"_snps_qvals_0.1=",total_num_qval_0.1,"_1000x_",percent_0.1,"pc_",num_qval_0.1,"_or_more.pdf"), collapse="")
  main_0.1=paste(c(i," n=", number_snps_left," 1000x"), collapse="")

  xlab_0.1="# q values per sampling <0.1"

#  pdf(file=file.path(GWASPath, name_0.05))
#  hist(qvals_0.05, main=main_0.05, labels=TRUE,xlab=xlab_0.05, breaks=20)
#  dev.off()

#  pdf(file=file.path(GWASPath, name_0.1))
#  hist(qvals_0.1, main=main_0.1, labels=TRUE,xlab=xlab_0.1, breaks=20)
#  dev.off()


#######################

####################### 

# To flatten the list

#  flat=do.call(c, unlist(results, recursive=FALSE))

#  name=paste(c(i,"_Distrib_q_vals 1000x n=",number_snps_left,".pdf"),collapse="")
#  pdf(file=file.path(GWASPath, name))
#  hist(flat, main=i, xlab="q values")
#  dev.off()

  summary=paste(c(i," ",number_snps," ",num_sig_0.05_orig," ",num_qval_0.05," ",percent_0.05," ",average_num_qvals0.05_sampling," ",num_sig_0.1_orig," ",num_qval_0.1," ",percent_0.1," ",average_num_qvals0.1_sampling), collapse="")
  write(summary, file=sumRes, sep="\t", append = TRUE)
}
}
