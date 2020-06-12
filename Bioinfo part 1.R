#Assignment Part 1
#Download the data
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile="gene_expression.tsv")

#Question1: 
x<- read.table("gene_expression.tsv")
x <-read.table("gene_expression.tsv", header= TRUE, row.names = 1)

head(x)
str(x)
#Quesn 2
x$mean <- rowMeans(x)
head(x)
 #Quesn 3
order(-x$mean)
x[order(-x$mean), ]
x_sorted <- x[order(-x$mean), ]
head(x_sorted,10)

#Question 4
b <- subset(x, mean<10) 
head(b)
str(b)

#Question 5
hist(x$mean)
hist(x$mean, breaks = 20)
