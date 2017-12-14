#####################################
###   CORES VS SPEED


setwd('/home/roland/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11')
library(readr)
benchmark <- read_delim("~/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11/benchmark_c.csv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

results <- table(benchmark)
str(benchmark)
plot(benchmark$Cores, benchmark$`Program as a whole`, col=as.factor(benchmark$Sample))


#####################################
###   PROBE SUITABILITY   


#setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/test_oct11/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
library(readr)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
str(probes)
species <- c("GC","PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
samples <- names(probes[grep("001_identity", names(probes))])
str(samples)
str(species)
identities <- probes[grep("identity|GC|Gene_Name", names(probes))]
identities <- cbind(identities, species)
str(identities)

plot.new()
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',xlim=c(0, 60), ylim=c(0, 100))
points(c(identities[1,]), c(identities[4,]), pch=".", cex=2)

str(identities)
str(identities$`102-11C_S3_L001_identity`)
str(identities[,1])
identities[1,0]

plot(identities$`%GC_ref`, identities$`102-11C_S3_L001_identity`)
plot(identities[1,], identities[2,])
identities <- identities[ which(identities$`%GC_ref`<40),]


plot(identities$`%GC_ref`, identities$`126-11c_S9_L001_identity`, pch=".", cex=2, ylim=c(-5,115), xlab="gene GC%", ylab="Percentage of Gene covered by 3 reads", col="red")

points(identities$`%GC_ref`, identities$`512022_S1_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`102-11C_S3_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`2KR049010_S2_L001_identity`, pch=".", cex=2)

points(identities$`%GC_ref`, identities$`1ES86798_S2_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`1EV02973_S11_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`1EV02981_S8_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`2KK58589_S7_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$BL37577_S10_L001_identity, pch=".", cex=2)
points(identities$`%GC_ref`, identities$BL37590_S4_L001_identity, pch=".", cex=2)
points(identities$`%GC_ref`, identities$C00413_S5_L001_identity, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`CC50362-2_S1_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`CC50362-3_S4_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$IMIN44_S12_L001_identity, pch=".", cex=2)
points(identities$`%GC_ref`, identities$ZL10_S3_L001_identity, pch=".", cex=2)

reg1 <- lm(identities$`102-11C_S3_L001_identity` ~ identities$`%GC_ref`)
reg2 <- lm(identities$`126-11c_S9_L001_identity` ~ identities$`%GC_ref`)
reg7 <- lm(identities$`2KR049010_S2_L001_identity` ~ identities$`%GC_ref`)
reg3 <- lm(identities$`1ES86798_S2_L001_identity` ~ identities$`%GC_ref`)
reg4 <- lm(identities$`1EV02973_S11_L001_identity` ~ identities$`%GC_ref`)
reg5 <- lm(identities$`1EV02981_S8_L001_identity` ~ identities$`%GC_ref`)
reg6 <- lm(identities$`2KK58589_S7_L001_identity` ~ identities$`%GC_ref`)
reg8 <- lm(identities$`512022_S1_L001_identity` ~ identities$`%GC_ref`)
reg9 <- lm(identities$BL37577_S10_L001_identity ~ identities$`%GC_ref`)
reg10 <- lm(identities$BL37590_S4_L001_identity ~ identities$`%GC_ref`)
reg11 <- lm(identities$C00413_S5_L001_identity ~ identities$`%GC_ref`)
reg12 <- lm(identities$C00434_S6_L001_identity ~ identities$`%GC_ref`)
reg13 <- lm(identities$`CC50362-2_S1_L001_identity` ~ identities$`%GC_ref`)
reg14 <- lm(identities$`CC50362-3_S4_L001_identity` ~ identities$`%GC_ref`)
reg15 <- lm(identities$IMIN44_S12_L001_identity ~ identities$`%GC_ref`)
reg16 <- lm(identities$Undetermined_S0_L001_identity ~ identities$`%GC_ref`)
reg17 <- lm(identities$ZL10_S3_L001_identity ~ identities$`%GC_ref`)


total_gc=rep(identities$`%GC_ref`, 17)
total_id=unlist(identities[,2:18])
x <- data.frame(total_gc, total_id)
reg18 <- lm(x$total_id ~x$total_gc)

total_gc=rep(identities$`%GC_ref`, 15)
total_id=unlist(identities[,4:18])
x <- data.frame(total_gc, total_id)
reg19 <- lm(x$total_id ~x$total_gc)

str(identities$`%GC_ref`)
str(identities[,2:18])

abline(reg1, col="green")
abline(reg2, col="green")
abline(reg3, col="green")
abline(reg4, col="green")
abline(reg5, col="green")
abline(reg6, col="green")
abline(reg7, col="green")
abline(reg8, col="green")
abline(reg9, col="green")
abline(reg10, col="green")
abline(reg11, col="green")
abline(reg12, col="green")
abline(reg13, col="green")
abline(reg14, col="green")
abline(reg15, col="green")
abline(reg16, col="green")
abline(reg17, col="green")
abline(reg18, col="red")
abline(reg19, col="blue")

title(xlab="gene GC%") #Ought to look at probe GC specifically?
?Extract
str(identities[2])
str(identities$`102-11C_S3_L001_identity`)
str(c(identities[2]))

s <- summary(reg18)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('topright', legend = rp, bty = 'n')

############################
## Exon-specific GC Vs Id ##
############################

library(readr)
exons <- read_delim("~/Malaria/SequenceCapture2AV4NU/Consensus/exon_gc_vs_id.tsv", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)

str(exons)
plot(exons$X2, exons$X3, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="red")
points(exons$X2, exons$X4, pch=".", cex=2, col="red")
points(exons$X2, exons$X5, pch=".", cex=2, col="red")
points(exons$X2, exons$X9, pch=".", cex=2, col="red")
points(exons$X2, exons$X14, pch=".", cex=2, col="red")
points(exons$X2, exons$X18, pch=".", cex=2, col="red")

points(exons$X2, exons$X10, pch=".", cex=2)

points(exons$X2, exons$X6, pch=".", cex=2)
points(exons$X2, exons$X7, pch=".", cex=2)
points(exons$X2, exons$X8, pch=".", cex=2)
points(exons$X2, exons$X12, pch=".", cex=2)
points(exons$X2, exons$X11, pch=".", cex=2)
points(exons$X2, exons$X13, pch=".", cex=2)
points(exons$X2, exons$X15, pch=".", cex=2)
points(exons$X2, exons$X16, pch=".", cex=2)
points(exons$X2, exons$X17, pch=".", cex=2)
points(exons$X2, exons$X19, pch=".", cex=2)

reg1 <- lm(exons$X3 ~ exons$X2)
reg2 <- lm(exons$X4 ~ exons$X2)
reg3 <- lm(exons$X5 ~ exons$X2)
reg4 <- lm(exons$X6 ~ exons$X2)
reg5 <- lm(exons$X7 ~ exons$X2)
reg6 <- lm(exons$X8 ~ exons$X2)
reg7 <- lm(exons$X9 ~ exons$X2)
reg8 <- lm(exons$X10 ~ exons$X2)
reg9 <- lm(exons$X11 ~ exons$X2)
reg10 <- lm(exons$X12 ~ exons$X2)
reg11 <- lm(exons$X13 ~ exons$X2)
reg12 <- lm(exons$X14 ~ exons$X2)
reg13 <- lm(exons$X15 ~ exons$X2)
reg14 <- lm(exons$X16 ~ exons$X2)
reg15 <- lm(exons$X17 ~ exons$X2)
reg16 <- lm(exons$X18 ~ exons$X2)
reg17 <- lm(exons$X19 ~ exons$X2)

total_gc=rep(exons$X2, 17)
total_id=unlist(exons[,3:19])
str(exons)
str(total_gc)
str(total_id)
x <- data.frame(total_gc, total_id)
str(x)
reg18 <- lm(x$total_id ~x$total_gc)

abline(reg1, col="green")
abline(reg2, col="green")
abline(reg3, col="green")
abline(reg4, col="green")
abline(reg5, col="green")
abline(reg6, col="green")
abline(reg7, col="green")
abline(reg8, col="green")
abline(reg9, col="green")
abline(reg10, col="green")
abline(reg11, col="green")
abline(reg12, col="green")
abline(reg13, col="green")
abline(reg14, col="green")
abline(reg15, col="green")
abline(reg16, col="green")
abline(reg17, col="green")
abline(reg18, col="red")

s <- summary(reg18)
summary(reg10)
coef(reg10)

r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

plot(total_id, pch=".", cex=2)
plot(exons$X10, pch=".", cex=2)
plot(identities$`512022_S1_L001_identity`, ylim=c(-5,125))

rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('topright', legend = rp, bty = 'n')


#######################
## Probe suitability ##
#######################
## Investigate GC Vs Length 

setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
library(readr)

probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
str(probes)

lengths <- probes[grep("ength|GC|Gene_Name", names(probes))]
#Remove reads with GC above 40
lengths <- lengths[ which(lengths$`%GC_ref`<40),]

#Make dataframe of which sample is which species
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
str(species)
samples2 <- names(probes[grep("001_length", names(probes))])
str(samples2)
spec_sam2 <- data.frame(species, samples2)

unique_species <- unique(species[2:17])
str(unique_species)




#Proof of concept
#i <- "EMCIR1"
i <- "PARUS1"
col_num <- as.numeric(rownames(subset(spec_sam2, species==i)[2]))
col_i <- rowMeans(subset(lengths, select = col_num+3), na.rm = FALSE)
temp_length <- 100*col_i/lengths[,2]
test <- subset(lengths, select = col_num+3)
names <- colnames(test)
test <- cbind(test, temp_length)
test <- cbind(test, lengths$`%GC_ref`)
colnames(test)<- c(names, "means", "%GC")
str(test)


#Make dataframe of mean length per species, instead of length per sample
gc_vs_length <- lengths[,3]
for (i in unique_species){
  col_num <- as.numeric(rownames(subset(spec_sam2, species==i)[2]))
  col_num+3
  col_i <- rowMeans(subset(lengths, select = col_num+3), na.rm = FALSE)
  temp_length <- 100*col_i/lengths[,2]
  
  names <- colnames(gc_vs_length)
  
  gc_vs_length <- cbind(gc_vs_length, temp_length)
  colnames(gc_vs_length)<- c(names, i)
}
str(gc_vs_length)


#Plot all species, overlapping. Kind of a mess.
plot(gc_vs_length$`%GC_ref`, gc_vs_length$EMCIR1, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="white")
points(gc_vs_length$`%GC_ref`, gc_vs_length$EMCIR1, pch=".", cex=2, col="black")
points(gc_vs_length$`%GC_ref`, gc_vs_length$PARUS1_PHSIB1, pch=".", cex=2, col="red")
points(gc_vs_length$`%GC_ref`, gc_vs_length$PARUS1, pch=".", cex=2, col="orange")
points(gc_vs_length$`%GC_ref`, gc_vs_length$SISKIN1, pch=".", cex=2, col="blue")
points(gc_vs_length$`%GC_ref`, gc_vs_length$GRW01, pch=".", cex=2, col="green")
points(gc_vs_length$`%GC_ref`, gc_vs_length$CraneSp, pch=".", cex=2, col="cyan")
points(gc_vs_length$`%GC_ref`, gc_vs_length$WW2, pch=".", cex=2, col="blue")
points(gc_vs_length$`%GC_ref`, gc_vs_length$PHSIB1, pch=".", cex=2, col="grey")
points(gc_vs_length$`%GC_ref`, gc_vs_length$SYAT02, pch=".", cex=2, col="purple")
points(gc_vs_length$`%GC_ref`, gc_vs_length$CXPIP23, pch=".", cex=2, col="brown")


#Plot total mean coverage length vs GC Content, add linear regression. 
valid <- c(2:9,11) #Remove 'undetermined' and 'GC' columns
valid
total_means <- rowMeans(subset(gc_vs_length, select = valid), na.rm = TRUE)

################################
#### MAKE EXONS_VS_GC_GRAPH ####
################################
?dir.create
?png
?dev.print
dir.create("graphs", showWarnings = FALSE)
png(filename="graphs/exons_vs_gc.png", width=828, height=504)
  plot(gc_vs_length$`%GC_ref`, total_means, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exons covered by 3 reads")
  reg_gc_vs_length <- lm(total_means ~ gc_vs_length$`%GC_ref`)
  abline(reg_gc_vs_length, col="green")


#Add formula for linear regression, confidence value, to plot.
coef(reg_gc_vs_length)
s <- summary(reg_gc_vs_length)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_gc_vs_length)[1], digits = 2), 
                        b_val = format(coef(reg_gc_vs_length)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')
dev.off()

######################################################
## Select good genes  ##
########################

setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
library(readr)
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
identities <- genes[grep("identity|GC|Gene_Name", names(genes))]
identities <- identities[ which(identities$`%GC_ref`<40),]
str(genes)
lengths <- genes[grep("ength|GC|Gene_Name", names(genes))]
#Make dataframe of which sample is which species
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
str(species)
samples2 <- names(genes[grep("001_length", names(genes))])
str(samples2)
spec_sam2 <- data.frame(species, samples2)

unique_species <- unique(species[2:17])
str(unique_species)

#Make dataframe of mean length per species, instead of length per sample
gc_vs_length <- lengths[,3]
for (i in unique_species){
  col_num <- as.numeric(rownames(subset(spec_sam2, species==i)[2]))
  col_num+3
  col_i <- rowMeans(subset(lengths, select = col_num+3), na.rm = FALSE)
  temp_length <- 100*col_i/lengths[,2]
  
  names <- colnames(gc_vs_length)
  
  gc_vs_length <- cbind(gc_vs_length, temp_length)
  colnames(gc_vs_length)<- c(names, i)
}
str(gc_vs_length)


#Plot total mean coverage length vs GC Content, add linear regression. 
valid <- c(2:9,11) #Remove 'undetermined' and 'GC' columns
valid
gene_data <- genes[grep("Gene_Name", names(genes))]

gene_data <- as.data.frame(genes$`#Gene_Name`)
means <- rowMeans(subset(gc_vs_length, select = valid), na.rm = TRUE)
gene_data <- cbind(gene_data, means)
gene_data <- cbind(gene_data, GC=genes$`%GC_ref`)
colnames(gene_data) <- c("name", "coverage", "GC")
str(gene_data)
good_genes <- gene_data[which(gene_data[,2]>60),]
bad_genes <- gene_data[which(gene_data[,2]<60),]
top_genes <- gene_data[gene_data$coverage > quantile(gene_data$coverage,prob=1-50/100),]
bottom_genes <- gene_data[gene_data$coverage < quantile(gene_data$coverage,prob=1-50/100),]

plot(gene_data$GC, gene_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="green")
points(gene_data$GC, gene_data$coverage, pch=".", cex=2, col="red")

plot(gene_data$GC, gene_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="black")
points(bottom_genes$`GC`, bottom_genes$coverage, pch=".", cex=2, col="red")
points(top_genes$`GC`, top_genes$coverage, pch=".", cex=2, col="blue")
total_means
str(total_means)
str(genes)
str(good_genes)
write.csv(top_genes, "top_genes.csv")

######################################################
#######################
## GC Vs Identity 

setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
library(readr)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
identities <- probes[grep("identity|GC|Gene_Name", names(probes))]
identities <- identities[ which(identities$`%GC_ref`<40),]
str(probes)
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")

samples <- names(probes[grep("001_identity", names(probes))])
spec_sam <- data.frame(species, samples)
str(spec_sam)
str(samples)
str(species)
unique_species <- unique(species[2:17])
str(unique_species)

str(identities$`%GC_ref`)
str(identities$`102-11C_S3_L001_identity`)

gc_vs_id <- identities[,2]

for (i in unique_species){
  col_num <- as.numeric(rownames(subset(spec_sam, species==i)[2]))
  col_num+2
  col_i <- rowMeans(subset(identities, select = col_num+2), na.rm = FALSE)
  names <- colnames(gc_vs_id)
  col_i
  gc_vs_id <- cbind(gc_vs_id, col_i)
  colnames(gc_vs_id)<- c(names, i)
}

str(gc_vs_id)
plot(gc_vs_id$`%GC_ref`, gc_vs_id$EMCIR1, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of CDS identical to reference", col="black")
points(gc_vs_id$`%GC_ref`, gc_vs_id$PARUS1_PHSIB1, pch=".", cex=2, col="red")
points(gc_vs_id$`%GC_ref`, gc_vs_id$PARUS1, pch=".", cex=2, col="orange")
points(gc_vs_id$`%GC_ref`, gc_vs_id$SISKIN1, pch=".", cex=2, col="blue")
points(gc_vs_id$`%GC_ref`, gc_vs_id$GRW01, pch=".", cex=2, col="green")
points(gc_vs_id$`%GC_ref`, gc_vs_id$CraneSp, pch=".", cex=2, col="cyan")
points(gc_vs_id$`%GC_ref`, gc_vs_id$WW2, pch=".", cex=2, col="blue")
points(gc_vs_id$`%GC_ref`, gc_vs_id$PHSIB1, pch=".", cex=2, col="grey")
points(gc_vs_id$`%GC_ref`, gc_vs_id$SYAT02, pch=".", cex=2, col="purple")
points(gc_vs_id$`%GC_ref`, gc_vs_id$CXPIP23, pch=".", cex=2, col="brown")

str(gc_vs_id)
valid <- c(2:9,11) #Remove 'undetermined'
valid
total_means <- rowMeans(subset(gc_vs_id, select = valid), na.rm = TRUE)
plot(gc_vs_id$`%GC_ref`, total_means, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of CDS identical to reference")
reg_gc_vs_id <- lm(total_means ~ gc_vs_id$`%GC_ref`)
abline(reg_gc_vs_id, col="green")

coef(reg_gc_vs_id)
s <- summary(reg_gc_vs_id)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_gc_vs_length)[1], digits = 2), 
                        b_val = format(coef(reg_gc_vs_length)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')

######################################################
######################################################
## Select good probes by exons ##
#################################

setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
library(readr)
exons <- read_delim("~/Malaria/SequenceCapture2AV4NU/Consensus/exon_gc_vs_id.tsv", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
identities <- probes[grep("identity|GC|Gene_Name", names(probes))]
identities <- identities[ which(identities$`%GC_ref`<40),]
str(probes)
lengths <- probes[grep("ength|GC|Gene_Name", names(probes))]
#Make dataframe of which sample is which species
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
str(species)
samples2 <- names(probes[grep("001_length", names(probes))])
str(samples2)
samples2
spec_sam2 <- data.frame(species, samples2)

unique_species <- unique(species[2:17])
str(unique_species)

#Make dataframe of mean length per species, instead of length per sample
gc_vs_length <- lengths[,3]
i <- unique_species[]
for (i in unique_species){
  col_num <- as.numeric(rownames(subset(spec_sam, species==i)[2]))
  col_num+3
  col_i <- rowMeans(subset(lengths, select = col_num+3), na.rm = FALSE)
  temp_length <- 100*col_i/lengths[,2]
  
  names <- colnames(gc_vs_length)
  
  gc_vs_length <- cbind(gc_vs_length, temp_length)
  colnames(gc_vs_length)<- c(names, i)
}
str(gc_vs_length)


#Plot total mean coverage length vs GC Content, add linear regression. 
valid <- c(2:9,11) #Remove 'undetermined' and 'GC' columns
valid
gene_data <- probes[grep("Gene_Name", names(probes))]

gene_data <- as.data.frame(probes$`#Gene_Name`)
means <- rowMeans(subset(gc_vs_length, select = valid), na.rm = TRUE)
gene_data <- cbind(gene_data, means)
gene_data <- cbind(gene_data, GC=probes$`%GC_ref`)
colnames(gene_data) <- c("name", "coverage", "GC")
str(gene_data)
good_probes <- gene_data[which(gene_data[,2]>60),]
bad_probes <- gene_data[which(gene_data[,2]<60),]
top_probes <- gene_data[gene_data$coverage > quantile(gene_data$coverage,prob=1-50/100),]
bottom_probes <- gene_data[gene_data$coverage < quantile(gene_data$coverage,prob=1-50/100),]

plot(gene_data$GC, gene_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="green")
points(gene_data$GC, gene_data$coverage, pch=".", cex=2, col="red")

plot(gene_data$GC, gene_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="black")
points(bottom_probes$`GC`, bottom_probes$coverage, pch=".", cex=2, col="red")
points(top_probes$`GC`, top_probes$coverage, pch=".", cex=2, col="blue")
total_means
str(total_means)
str(probes)
str(good_probes)
write.csv(top_probes, "top_probes.csv")

######################################################
######################################################
## Select good probes by probes ##
#################################

probes <- read_delim("probes_stats.tsv", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
library(readr)
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/Consensus')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
probes_stats <- read_delim("probes_stats.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
str(probes_stats)
coverage <- probes_stats[grep("_coverage|GC|probe", names(probes_stats))]
str(coverage)

str(species)
samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
str(samples3)
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
unique_species <- unique(species[2:17])
spec_sam3 <- data.frame(species, samples3)

gc_vs_probe_coverage <- probes_stats[,2]

for (i in unique_species){
  col_num <- as.numeric(rownames(subset(spec_sam3, species==i)[2]))
  col_num+2
  col_i <- rowMeans(subset(coverage, select = col_num+2), na.rm = FALSE)
  names <- colnames(gc_vs_probe_coverage)
  col_i
  gc_vs_probe_coverage <- cbind(gc_vs_probe_coverage, col_i)
  colnames(gc_vs_probe_coverage)<- c(names, i)
}

str(gc_vs_probe_coverage)

plot(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$EMCIR1, pch=".", cex=2, ylim=c(-5,105), xlab="probe GC%", ylab="Percentage of probe covered by 3 reads", col="white")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$EMCIR1, pch=".", cex=2, col="black")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$PARUS1_PHSIB1, pch=".", cex=2, col="red")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$PARUS1, pch=".", cex=2, col="orange")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$SISKIN1, pch=".", cex=2, col="blue")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$GRW01, pch=".", cex=2, col="green")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$CraneSp, pch=".", cex=2, col="cyan")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$WW2, pch=".", cex=2, col="blue")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$PHSIB1, pch=".", cex=2, col="grey")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$SYAT02, pch=".", cex=2, col="purple")
points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$CXPIP23, pch=".", cex=2, col="brown")


valid <- c(2:9,11) #Remove 'undetermined' and 'GC' columns
valid
probe_data <- probes_stats[grep("probe", names(probes_stats))]
#probe_data <- as.data.frame(probes_stats$probe)
probe_means <- rowMeans(subset(gc_vs_probe_coverage, select = valid), na.rm = TRUE)
probe_data <- cbind(probe_data, probe_means)
probe_data <- cbind(probe_data, GC=probes_stats$GC)
colnames(probe_data) <- c("name", "coverage", "GC")
str(probe_data)

plot(probe_data$coverage, pch=".", cex=2)
top_probes <- probe_data[probe_data$coverage > quantile(probe_data$coverage,prob=1-50/100),]
bottom_probes <- probe_data[probe_data$coverage < quantile(probe_data$coverage,prob=1-50/100),]

##############################
### MAKE GC_VS_PROBE_GRAPH ###
##############################

dir.create("graphs", showWarnings = FALSE)
png(filename="graphs/probes_vs_gc.png", width=828, height=504)

plot(probe_data$GC, probe_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="probe GC%", ylab="Percentage of probe covered by 3 reads")

reg_gc_vs_probe_coverage <- lm(probe_means ~ gc_vs_probe_coverage$GC)
abline(reg_gc_vs_probe_coverage, col="green")

coef(reg_gc_vs_probe_coverage)
s <- summary(reg_gc_vs_probe_coverage)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_gc_vs_probe_coverage)[1], digits = 2), 
                        b_val = format(coef(reg_gc_vs_probe_coverage)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')
dev.off()

### Visualization of top_probes(kept) and bottom_probes(discarded)
plot(probe_data$GC, probe_data$coverage, pch=".", cex=2, ylim=c(-5,105), xlab="gene GC%", ylab="Percentage of exon covered by 3 reads", col="black")
points(bottom_probes$`GC`, bottom_probes$coverage, pch=".", cex=2, col="red")
points(top_probes$`GC`, top_probes$coverage, pch=".", cex=2, col="blue")


##################################################
###  SORT AWAY BAD PROBES 

library(readr)
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria')
probe_db <- read_delim("~/bin/SeqCapAnalyst/probe_to_exon_to_gene.list", "\t", escape_double = FALSE, col_names = FALSE, na = "NA", trim_ws = TRUE)
colnames(probe_db) <- c("probe_name","probe_pos","exon_pos","gene_name")
#View(probe_db)
str(probe_db)

# Find good genes
genes <- read_delim("/home/roland/Malaria/SequenceCapture2AV4NU/Consensus/gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(queryLength = col_integer()), trim_ws = TRUE, skip = 1)
valid <- c(2:9,11) #Remove 'undetermined' and 'GC' columns
valid
gene_data <- genes[grep("Gene_Name", names(genes))]

gene_data <- as.data.frame(genes$`#Gene_Name`)
means <- rowMeans(subset(gc_vs_length, select = valid), na.rm = TRUE)
gene_data <- cbind(gene_data, means)
gene_data <- cbind(gene_data, GC=genes$`%GC_ref`)
colnames(gene_data) <- c("name", "coverage", "GC")
str(gene_data)
good_genes <- gene_data[which(gene_data[,2]>60),]
bad_genes <- gene_data[which(gene_data[,2]<60),]
top_genes <- gene_data[gene_data$coverage > quantile(gene_data$coverage,prob=1-50/100),]
bottom_genes <- gene_data[gene_data$coverage < quantile(gene_data$coverage,prob=1-50/100),]


#Remove probes in genes with bad coverage
probes_less_bad_genes <- subset(probe_db, gene_name %in% good_genes$name)
str(probes_less_bad_genes)
probes_less_bad_genes$probe_pos
#remove probes with bad coverage even in high-coverage genes
probes_only_good <- subset(probes_less_bad_genes, probe_pos %in% top_probes$name)
str(probes_only_good)

write.csv(probes_only_good, "good_probes.csv")

#plot(probe_data$identity ~ probe_data$coverage, pch=".", cex=2)


############################################################
### BENCHMARKING: TIMES VS READ NUMBERS

library(readr)
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/test_execution_time_all_samples_oct25')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria\\execution_times')
times <- read_delim("times_with_seconds.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sizes <- read_delim("sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
matches <- read_delim("match_stats.txt", " ", skip = 1, escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(times) <- c("names", "h", "m", "s", "total_seconds")
colnames(sizes) <- c("names", "sizes")
colnames(matches) <- c("names", "total_reads", "matching_reads", "percent", "trash" )

times$names
sizes$names
matches$names
str(times)
str(sizes)
str(matches)

plot(times$total_seconds ~ matches$matching_reads, ylab="execution time (s)", xlab="Number of matching reads")

reg_time_vs_reads <- lm(times$total_seconds ~ matches$matching_reads)
abline(reg_time_vs_reads, col="green")

coef(reg_time_vs_reads)
s <- summary(reg_time_vs_reads)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_time_vs_reads)[1], digits = 2), 
                        b_val = format(coef(reg_time_vs_reads)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')

##############################################################
### BENCHMARKING: TIMES VS FILESIZE

library(readr)
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/test_execution_time_all_samples_oct25')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria\\execution_times')
times <- read_delim("times_with_seconds.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sizes <- read_delim("sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
matches <- read_delim("match_stats.txt", " ", skip = 1, escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(times) <- c("names", "h", "m", "s", "total_seconds")
colnames(sizes) <- c("names", "sizes")
colnames(matches) <- c("names", "total_reads", "matching_reads", "percent", "trash" )

times$names
sizes$names
matches$names
str(times)
str(sizes)
str(matches)

plot(times$total_seconds ~ sizes$sizes, ylab="execution time (s)", xlab="Size of fastq files (bytes)")

reg_time_vs_size <- lm(times$total_seconds ~ sizes$sizes)
abline(reg_time_vs_size, col="green")

coef(reg_time_vs_size)
s <- summary(reg_time_vs_size)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_time_vs_size)[1], digits = 2), 
                        b_val = format(coef(reg_time_vs_size)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')

##############################################
### BENCHMARKING INPUT SIZE VS OUTPUT SIZE

library(readr)
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/test_execution_time_all_samples_oct25')
#setwd('C:\\Documents and Settings\\Roland\\Dropbox\\Bioinformatics\\Malaria\\execution_times')
sizes <- read_delim("sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
folder_sizes <- read_delim("folder_sizes.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(sizes) <- c("names", "sizes")
colnames(folder_sizes) <- c("names", "sizes", "trash")

sizes$names

str(folder_sizes)
str(sizes)


plot(folder_sizes$sizes ~ sizes$sizes, ylab="output folder size (bytes)", xlab="Size of fastq files (bytes)")

reg_input_vs_size <- lm(folder_sizes$sizes ~ sizes$sizes)
abline(reg_input_vs_size, col="green")

coef(reg_input_vs_size)
s <- summary(reg_input_vs_size)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
rp[1] = substitute(expression(italic(y) == a_val + b_val %.% italic(x) ), 
                   list(a_val = format(coef(reg_input_vs_size)[1], digits = 2), 
                        b_val = format(coef(reg_input_vs_size)[2], digits = 3)))[2]
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

legend('bottomright', legend = rp, bty = 'n')



####################################################
### BENCHMARKING PIPELINE AS A WHOLE

setwd('/home/roland/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11')
library(readr)
benchmark1 <- read_delim("~/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11/benchmark_c.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
benchmark2 <- read_delim("~/Malaria/SequenceCapture2AV4NU/benchmark_test_sep17/benchmark_d.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
x <- as.data.frame(benchmark1$Sample)
x <- cbind(x, benchmark1$Cores)
x <- cbind(x, benchmark1$`Program as a whole`)
x <- cbind(x, benchmark2$`Program as a whole`)
colnames(x) <- c("Sample", "Cores", "Run1", "Run2")
x$average <- rowMeans(x[c("Run1","Run2")], na.rm = TRUE)
plot(x$Cores, x$average, col=as.factor(x$Sample), xlab="Number of Cores", pch=19, ylab="Average completion time (second)")
legend('topright', c("BL37590 (2.6GB)","BL37577 (1.4GB)","IMIN44 (0.7GB)") , lty=0, col=c("green","black","red"), bty='n', cex=1, pch=19)
#legend('topright', levels(x$Sample) , lty=1, col=as.numeric(unique(factor(x$Sample))), bty='n', cex=.75)


###############################################
### EXECUTION TIME OF DIFFERENT MODULES
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11')
library(readr)
benchmark1 <- read_delim("~/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11/benchmark_c.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
benchmark2 <- read_delim("~/Malaria/SequenceCapture2AV4NU/benchmark_test_sep17/benchmark_d.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
library(ggplot2)

module_times_1 <- as.numeric(benchmark1[15,c(3:23, 25)])
module_times_2 <- as.numeric(benchmark1[23,c(3:23, 25)])
module_times_3 <- as.numeric(benchmark1[19,c(3:23, 25)])

seconds <- c(module_times_1, module_times_2, module_times_3)

cores <-c(rep(" 1 core", 22), rep(" 6 cores", 22), rep("36 cores", 22))
time_data <-data.frame(modules, seconds, cores)

p <-ggplot(time_data, aes(modules, seconds ))
p +geom_bar(stat = "identity", aes(fill = cores), position = "dodge") + coord_flip()

#################################################
### EFFICIENCY OF SUBSAMPLER 
setwd('/home/roland/Malaria/SequenceCapture2AV4NU/benchmark_test_sep11')

IMIN_sub <- read_delim("~/Malaria/SequenceCapture2AV4NU/test_sub_sample/IMIN_sub_depth.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
IMIN_main <- read_delim("~/Malaria/SequenceCapture2AV4NU/test_sub_sample/IMIN_main_depth.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
colnames(IMIN_sub) <- c("scaffold", "pos", "reads")
colnames(IMIN_main) <- c("scaffold", "pos", "reads")

plot(IMIN_main$reads ~ IMIN_main$pos, pch=".", ylim=c(70,2200), cex=3, xlab="scaffold position", ylab="read coverage", col="#555588")
points(IMIN_sub$reads ~ IMIN_sub$pos, pch=".", cex=2, col="#55FF5588")

legend('topleft', c("Before subsampling","After subsampling") , lty=0, col=c("#555588","#55FF55"), bty='n', cex=1, pch=19)


