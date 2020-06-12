# Bioinformatics-assignment3
Bioinformatics assignment 3

#Assignmnet Part 2

library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")


R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl","-parse_seqids")


#Question 2
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile= "sample.fa")
x <- read.fasta("sample.fa") 
head(x)
str(x)

o <- x [[31]] 

seqinr::getLength(o)

seqinr::GC(o)


#Question 3:  

myblastn_tab <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(myseq,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 6 -evalue 0.05 -ungapped >"
                       , mytmpfile2))
  
  res <- NULL
  if (file.info(mytmpfile2)$size > 0 ) {
    res <- read.csv(mytmpfile2,sep="\t",header=FALSE)
    colnames(res) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                       "qstart","qend","sstart","send","evalue","bitscore")
  }
  unlink(c(mytmpfile1,mytmpfile2))
  if (!is.null(res)  ) {
    res <- res[order(-res$bitscore),]
  }
  res
}
mutator <- function(myseq,nmut) {
  myseq_mod <- myseq
  mypos<-sample(seq_along(myseq),nmut)
  myseq_mod[mypos] <- sample(c("a","c","g","t"),length(mypos),replace = TRUE)
  return(myseq_mod)
}

res <- myblastn_tab(myseq = o, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(res)
str(res)

mysequence <- as.character(res$sseqid[1:1047])


View(res)
res


#Question 4
z_mut <-mutator (myseq = o, 100) 

z_mut

myblastn_tab(z_mut,db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
res <- myblastn_tab(z_mut,db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(res)




z_mut_1 <- DNAString(c2s(z_mut))

x_1 <- DNAString(c2s(o)) 

aln <- pairwiseAlignment(x_1,z_mut_1) 


pid(aln)

nmismatch(aln)



