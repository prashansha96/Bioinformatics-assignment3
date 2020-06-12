Question 6
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv", 
              destfile = "growth_data.csv")
a<- read.csv("growth_data.csv", header = TRUE)
head(a)
str(a)


#Question 7:
mean(a$Circumf_2004_cm)
subset(a, Site=="northeast")
head(a)
str(a)
subset(a, Site=="southwest")
head(a)
str(a)
ne <- subset(a, Site=="northeast")
head(ne)
str(ne)
sw <- subset(a, Site=="southwest")
head(sw)
str(sw)
mean(ne$Circumf_2004_cm)
mean(sw$Circumf_2004_cm)
mean(ne$Circumf_2019_cm)
mean(sw$Circumf_2019_cm)
sd(sw$Circumf_2019_cm)
sd(ne$Circumf_2019_cm)
sd(sw$Circumf_2004_cm)
sd(ne$Circumf_2004_cm)
head(ne)
head(sw)

# Question 8: Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm)
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm,sw$Circumf_2004_cm,sw$Circumf_2019_cm,names= c("ne2004","ne2019","sw2004","sw2019"), ylab="Circumference (cm)", xlab="site and years", main="Growth at plantation site")
#Question 9
mean(ne$Circumf_2009_cm)
c <- mean(ne$Circumf_2009_cm)
mean(ne$Circumf_2019_cm)
d <- mean(ne$Circumf_2019_cm)
x <- sum(c,d)
mean <- x/2
mean(sw$Circumf_2009_cm)
p <- mean(sw$Circumf_2009_cm)
mean(sw$Circumf_2019_cm)
q <- mean(sw$Circumf_2019_cm)
y <- sum(p,q)
mean <- y/2
nemean <- ne$Circumf_2019_cm-ne$Circumf_2009_cm
ne$growth <- ne$Circumf_2019_cm-ne$Circumf_2009_cm
mean(ne$growth)
head(ne$growth)
head(ne)
str(ne)
swmean <- sw$Circumf_2019_cm-sw$Circumf_2009_cm
sw$growth <- sw$Circumf_2019_cm-sw$Circumf_2009_cm
mean(sw$growth)
head(sw$groth)
head(sw)
str(sw)

#Question 10

t.test(ne$growth,sw$growth)
wilcox.test(ne$growth,sw$growth)

