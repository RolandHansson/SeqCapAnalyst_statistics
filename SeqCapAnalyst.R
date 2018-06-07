#####################################
###   CORES VS SPEED

setwd('.')
library(readr)
benchmark <- read_delim("benchmark_test_sep11/benchmark_c.csv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

results <- table(benchmark)
str(benchmark)
plot(benchmark$Cores, benchmark$`Program as a whole`, col=as.factor(benchmark$Sample))


#####################################
###   PROBE SUITABILITY   

setwd('.')
library(readr)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
str(probes)
species <- c("GC","PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
samples <- names(probes[grep("001_identity", names(probes))])
str(samples)
str(species)
identities <- probes[grep("identity|GC|Gene_Name", names(probes))]
#identities <- cbind(identities, species)
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


plot(identities$`%GC_ref`, identities$`126-11c_S9_L001_identity`, pch=".", cex=2, ylim=c(-5,115), xlab="gene GC%", ylab="Percentage of Gene covered by 3 reads", col="white")

points(identities$`%GC_ref`, identities$`102-11C_S3_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`126-11c_S9_L001_identity`, pch=".", cex=2)

points(identities$`%GC_ref`, identities$`512022_S1_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`2KR049010_S2_L001_identity`, pch=".", cex=2)

points(identities$`%GC_ref`, identities$C00413_S5_L001_identity, pch=".", cex=2)

points(identities$`%GC_ref`, identities$`1ES86798_S2_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`1EV02973_S11_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`1EV02981_S8_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$`2KK58589_S7_L001_identity`, pch=".", cex=2)
points(identities$`%GC_ref`, identities$BL37577_S10_L001_identity, pch=".", cex=2)
points(identities$`%GC_ref`, identities$BL37590_S4_L001_identity, pch=".", cex=2)
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


#######################
## Probe suitability ##
#######################
## Investigate GC Vs Length 

setwd('.')
library(readr)

probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
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
#### MAKE GENES_VS_GC_GRAPH ####
################################
?dir.create
?png
?dev.print
dir.create("graphs", showWarnings = FALSE)
png(filename="graphs/genes_vs_gc.png", width=828, height=504)
  plot(gc_vs_length$`%GC_ref`, total_means, pch=".", cex=2, ylim=c(-5,105), xlab=" GC%", ylab="Percentage of genes covered by 3 reads")
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

setwd('.')
library(readr)
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
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
## Gene GC Vs Identity 

setwd('.')
library(readr)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
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
## Select good probes by genes ##
#################################

setwd('.')
library(readr)
exons <- read_delim("exon_length_and_id2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
probes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
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

setwd('.')
library(readr)
probes <- read_delim("probes_stats.tsv", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
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
png(filename="graphs/probes_vs_gc2.png", width=828, height=504)

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
setwd('.')
probe_db <- read_delim("probe_to_exon_to_gene.list", "\t", escape_double = FALSE, col_names = FALSE, na = "NA", trim_ws = TRUE)
colnames(probe_db) <- c("probe_name","probe_pos","exon_pos","gene_name")
str(probe_db)

# Find good genes
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
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
setwd('.')
times <- read_delim("test_execution_time_all_samples_oct25/times_with_seconds.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sizes <- read_delim("test_execution_time_all_samples_oct25/sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
matches <- read_delim("test_execution_time_all_samples_oct25/match_stats.txt", " ", skip = 1, escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

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
#######################################
### BENCHMARKING: TIMES VS FILESIZE

library(readr)
setwd('.')
times <- read_delim("test_execution_time_all_samples_oct25/times_with_seconds.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sizes <- read_delim("test_execution_time_all_samples_oct25/sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
matches <- read_delim("test_execution_time_all_samples_oct25/match_stats.txt", " ", skip = 1, escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

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


########################################################
##############################################
### BENCHMARKING INPUT SIZE VS OUTPUT SIZE

library(readr)
setwd('.')
sizes <- read_delim("test_execution_time_all_samples_oct25/sizes_added.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
folder_sizes <- read_delim("test_execution_time_all_samples_oct25/folder_sizes.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

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



##########################################################################
#######################################
### BENCHMARKING PIPELINE AS A WHOLE

setwd('.')
library(readr)
benchmark1 <- read_delim("benchmark_test_sep11/benchmark_c.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
benchmark2 <- read_delim("benchmark_test_sep17/benchmark_d.tsv", 
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
setwd('.')
library(readr)
benchmark1 <- read_delim("benchmark_test_sep11/benchmark_c.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
benchmark2 <- read_delim("benchmark_test_sep17/benchmark_d.tsv", 
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
setwd('.')

IMIN_sub <- read_delim("test_sub_sample/IMIN_sub_depth.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
IMIN_main <- read_delim("test_sub_sample/IMIN_main_depth.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
colnames(IMIN_sub) <- c("scaffold", "pos", "reads")
colnames(IMIN_main) <- c("scaffold", "pos", "reads")

plot(IMIN_main$reads ~ IMIN_main$pos, pch=".", ylim=c(70,2200), cex=3, xlab="scaffold position", ylab="read coverage", col="#555588")
points(IMIN_sub$reads ~ IMIN_sub$pos, pch=".", cex=2, col="#55FF5588")

legend('topleft', c("Before subsampling","After subsampling") , lty=0, col=c("#555588","#55FF55"), bty='n', cex=1, pch=19)


#########################################################
##################################################
## PROBE GC VS PROBE COVERAGE, Majoris ONLY ##
#############################################

library(readr)
setwd('.')
probes_stats <- read_delim("probes_stats.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
str(probes_stats)
coverage <- probes_stats[grep("_coverage|GC|probe", names(probes_stats))]
str(coverage)

samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
str(samples3)
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")

name_species <- data.frame(matrix(ncol = 2, nrow = 0))
names(name_species) <- c("name", "species")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1")
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02")
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2")
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1")
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1")
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01")
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1")
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2")
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23")
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?")
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp")
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""))
names(spec_sam3) <- c("species", "name")

relevant_species <- c("PARUS1", "WW2", "PHSIB1")
#spec_sam3 <- data.frame(species, samples3)

gc_vs_probe_coverage <- probes_stats[,2]

subset(spec_sam3, species %in% c("PARUS1", "WW2", "PHSIB1"))

  col_num <- as.numeric(rownames(subset(spec_sam3, species %in% relevant_species)[2]))
  col_num+2
  colnames(subset(coverage, select = col_num+2))
  col_i <- rowMeans(subset(coverage, select = col_num+2), na.rm = FALSE)
  names <- colnames(gc_vs_probe_coverage)
  col_i
  gc_vs_probe_coverage <- cbind(gc_vs_probe_coverage, col_i)
  colnames(gc_vs_probe_coverage)<- c(names, "MAJORIS")

str(gc_vs_probe_coverage)

plot(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$MAJORIS, pch=".", cex=2, ylim=c(-5,105), xlim=c(-2,85), xlab="probe GC%", ylab="Percentage of probe covered by 3 reads", col=rgb(0, 0, 1, 0.33))
#points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$PARUS1, pch=".", cex=2, col="orange")
#points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$WW2, pch=".", cex=2, col="blue")
#points(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$PHSIB1, pch=".", cex=2, col="grey")

reg_gc_vs_probe_coverage <- lm(gc_vs_probe_coverage$MAJORIS ~ gc_vs_probe_coverage$GC)
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

##########
## WITH ALL POINTS INSTEAD OF MEAN VALUE
plot(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$MAJORIS, pch=".", cex=2, ylim=c(-5,105), xlim=c(0,85), xlab="probe GC%", ylab="Percentage of probe covered by 3 reads", col=rgb(0, 0, 1, 0))
points(gc_vs_probe_coverage$GC, coverage$`1ES86798_S2_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$GC, coverage$`1EV02981_S8_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$GC, coverage$BL37577_S10_L001_coverage, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$GC, coverage$BL37590_S4_L001_coverage, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$GC, coverage$`CC50362-2_S1_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$GC, coverage$`CC50362-3_S4_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))

total <- c(coverage$`1ES86798_S2_L001_coverage`, coverage$`1EV02981_S8_L001_coverage`, coverage$BL37577_S10_L001_coverage, coverage$BL37590_S4_L001_coverage, coverage$`CC50362-2_S1_L001_coverage`, coverage$`CC50362-3_S4_L001_coverage`) 
total
total_gc <- rep(gc_vs_probe_coverage$GC, 6)
str(total)
str(total_gc)
str(gc_vs_probe_coverage$GC)
reg_gc_vs_probe_coverage <- lm(total ~ total_gc)
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

#########
## With mean values of y

means_vs_gc <- data.frame(matrix(ncol = 3, nrow = 0))
names(means_vs_gc) <- c("GC", "gc%", "mean")
 
i <- 36
for (i in 1:120){
  x <- round(100*i/120, 1)
  y <- mean(subset(gc_vs_probe_coverage, GC==x)[,2])
  means_vs_gc[nrow(means_vs_gc) + 1,] = list(i, x, y)
}
plot(means_vs_gc$`gc%`, means_vs_gc$mean, pch=".", cex=5, ylim=c(-5,105), xlim=c(0,85), xlab="probe GC%", ylab="Mean probe coverage", col=rgb(0, 0, 1, 1))

reg_gc_vs_probe_coverage <- lm(means_vs_gc$mean ~ means_vs_gc$`gc%`)
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


##############################################
## EXON GC VS EXON COVERAGE, Majoris ONLY ##
##############################################

library(readr)
setwd('.')
exons <- read_delim("exon_length_and_id2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
str(exons)

coverage <- exons[grep("coverage", names(exons))]
str(coverage)

samples3 <- names(coverage[grep("001_coverage", names(coverage))])

str(samples3)
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")

name_species <- data.frame(matrix(ncol = 3, nrow = 0))
names(name_species) <- c("name", "species", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",6.4)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",4.3)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",1.6)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",6.1)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",4.6)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",4.3)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",9.2)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",9.2)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""))
names(spec_sam3) <- c("species", "name")

relevant_species <- c("PARUS1", "WW2", "PHSIB1")
#spec_sam3 <- data.frame(species, samples3)

gc_vs_probe_coverage <- exons[,3]

subset(spec_sam3, species %in% c("PARUS1", "WW2", "PHSIB1"))

col_num <- as.numeric(rownames(subset(spec_sam3, species %in% relevant_species)[2]))
col_num+2
colnames(subset(coverage, select = col_num))
col_i <- rowMeans(subset(coverage, select = col_num), na.rm = FALSE)
names <- colnames(gc_vs_probe_coverage)
col_i
gc_vs_probe_coverage <- cbind(gc_vs_probe_coverage, col_i)
colnames(gc_vs_probe_coverage)<- c(names, "MAJORIS")

str(gc_vs_probe_coverage)

plot(gc_vs_probe_coverage$gc, gc_vs_probe_coverage$MAJORIS, pch=".", cex=2, ylim=c(-5,105), xlim=c(2,70), xlab="exon GC%", ylab="Percentage of exon covered by 3 reads", col=rgb(0, 0, 1, 0.33))

reg_gc_vs_probe_coverage <- lm(gc_vs_probe_coverage$MAJORIS ~ gc_vs_probe_coverage$GC)
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

##########
## WITH ALL POINTS INSTEAD OF MEAN VALUE
plot(gc_vs_probe_coverage$gc, gc_vs_probe_coverage$MAJORIS, pch=".", cex=2, ylim=c(-5,105), xlim=c(0,85), xlab="exon GC%", ylab="Percentage of exon covered by 3 reads", col=rgb(0, 0, 1, 0))
points(gc_vs_probe_coverage$gc, coverage$`1ES86798_S2_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$gc, coverage$`1EV02981_S8_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$gc, coverage$BL37577_S10_L001_coverage, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$gc, coverage$BL37590_S4_L001_coverage, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$gc, coverage$`CC50362-2_S1_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))
points(gc_vs_probe_coverage$gc, coverage$`CC50362-3_S4_L001_coverage`, pch=".", cex=2, col=rgb(0, 0, 1, 0.5))

total <- c(coverage$`1ES86798_S2_L001_coverage`, coverage$`1EV02981_S8_L001_coverage`, coverage$BL37577_S10_L001_coverage, coverage$BL37590_S4_L001_coverage, coverage$`CC50362-2_S1_L001_coverage`, coverage$`CC50362-3_S4_L001_coverage`) 
total
total_gc <- rep(gc_vs_probe_coverage$GC, 6)
str(total)
str(total_gc)
str(gc_vs_probe_coverage$GC)
reg_gc_vs_probe_coverage <- lm(total ~ total_gc)
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

legend('topright', legend = rp, bty = 'n')

#########
## With mean values of y

means_vs_gc <- data.frame(matrix(ncol = 3, nrow = 0))
names(means_vs_gc) <- c("GC", "gc%", "mean")
i <- 1
x <- i
z <- subset(gc_vs_probe_coverage, x>gc)
z2 <- subset(z, gc>(x-1))
for (i in 1:200){
  x <- i/2
  z1 <- subset(gc_vs_probe_coverage, x>=gc)
  z2 <- subset(z1, gc>(x-0.5))
  y <- mean(z2[,2])
  means_vs_gc[nrow(means_vs_gc) + 1,] = list(i, x, y)
}
plot(means_vs_gc$`gc%`, means_vs_gc$mean, pch=".", cex=5, ylim=c(-5,105), xlim=c(0,85), xlab="exon GC%", ylab="Mean exon coverage", col=rgb(0, 0, 1, 1))

reg_gc_vs_probe_coverage <- lm(means_vs_gc$mean ~ means_vs_gc$`gc%`)
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

legend('topright', legend = rp, bty = 'n')

#############
## With 50 groups of about 50, MAJORIS only


means_vs_gc <- data.frame(matrix(ncol = 3, nrow = 0))
names(means_vs_gc) <- c("GC", "gc%", "mean")
i <- 3
x <- i
z3 <- list()
str(gc_vs_probe_coverage$gc)
for (i in 1:50){
  x <- i/2
  z1 <- gc_vs_probe_coverage[gc_vs_probe_coverage$gc > quantile(gc_vs_probe_coverage$gc,prob=1-i/50),]
  z2 <- z1[z1$gc <= quantile(gc_vs_probe_coverage$gc,prob=1-(i-1)/50),]
  str(z2)
  x <- mean(z2[,1])
  y <- mean(z2[,2])
  z3 <- c(z3, z2$gc)
  means_vs_gc[nrow(means_vs_gc) + 1,] = list(i, x, y)
}
plot(means_vs_gc$`gc%`, means_vs_gc$mean, pch=".", cex=5, ylim=c(0,60), xlim=c(5,60), xlab="exon GC%", ylab="Mean exon coverage", col=rgb(0, 0, 1, 1))

reg_gc_vs_probe_coverage <- lm(means_vs_gc$mean ~ means_vs_gc$`gc%`)
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

legend('topright', legend = rp, bty = 'n')

#############
## With 100 groups of about 26, MAJORIS only

means_vs_gc <- data.frame(matrix(ncol = 3, nrow = 0))
names(means_vs_gc) <- c("GC", "gc%", "mean")
i <- 2
x <- i
z3 <- list()
z4 <- list()
str(gc_vs_probe_coverage$gc)
for (i in 1:100){
  x <- i/2
  z1 <- gc_vs_probe_coverage[gc_vs_probe_coverage$gc > quantile(gc_vs_probe_coverage$gc,prob=1-i/100),]
  z2 <- z1[z1$gc <= quantile(gc_vs_probe_coverage$gc,prob=1-(i-1)/100),]
  str(z2)
  l <- length(z2$gc)
  x <- mean(z2[,1])
  y <- mean(z2[,2])
  z3 <- c(z3, z2$gc)
  z4 <- c(z4, l)
  means_vs_gc[nrow(means_vs_gc) + 1,] = list(i, x, y)
}
plot(means_vs_gc$`gc%`, means_vs_gc$mean, pch=".", cex=5, ylim=c(0,60), xlim=c(5,60), xlab="exon GC%", ylab="Mean exon coverage", col=rgb(0, 0, 1, 1))

reg_gc_vs_probe_coverage <- lm(means_vs_gc$mean ~ means_vs_gc$`gc%`)
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

legend('topright', legend = rp, bty = 'n')

#############
## With groups of about 100, MAJORIS only, without messing with quantiles.

means_vs_gc <- data.frame(matrix(ncol = 3, nrow = 0))
names(means_vs_gc) <- c("GC", "gc%", "mean")
i <- 1
x <- i
str(gc_vs_probe_coverage$gc)
gc_vs_probe_coverage[0:1,]
gc_vs_probe_sorted <- gc_vs_probe_coverage[order(gc_vs_probe_coverage$gc, decreasing=FALSE),]
str(gc_vs_probe_sorted)
bitesize <- length(gc_vs_probe_coverage$gc)/100
for (i in 1:100){
  x <- i
  min <- round(bitesize*(i-1))+1
  max <- round(bitesize*i)
  z1 <- gc_vs_probe_sorted[min:max,]
  x <- mean(z1[,1])
  y <- mean(z1[,2])
  means_vs_gc[nrow(means_vs_gc) + 1,] = list(i, x, y)
}
plot(means_vs_gc$`gc%`, means_vs_gc$mean, pch=".", cex=5, ylim=c(0,60), xlim=c(5,60), xlab="Mean exon GC%", ylab="Mean exon coverage (%)", col=rgb(0, 0, 1, 1))

reg_gc_vs_probe_coverage <- lm(means_vs_gc$mean ~ means_vs_gc$`gc%`)
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

legend('topright', legend = rp, bty = 'n')


############################################################
#####################################################
## PHYLOGENETIC DISTANCE VS EXON COVERAGE ######
#############################################


library(readr)
setwd('.')
exons <- read_delim("exon_length_and_id2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
str(exons)

coverage <- exons[grep("coverage", names(exons))]
str(coverage)

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)

species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 3, nrow = 0))
names(name_species) <- c("name", "species", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",6.4)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",4.3)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",1.6)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",6.1)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",4.6)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",4.3)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",3.9)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",9.2)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",9.2)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,3])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "CXPIP23", "CraneSp")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
rel_coverage <- rel_coverage[-remove_this,]

cov_vs_dis <- data.frame(matrix(ncol = 4, nrow = 0))
names(cov_vs_dis) <- c("sample", "has100", "has50", "has40")

i <- toString(rel_coverage[3,2])
x <- coverage[,i]
for(i in rel_coverage[,2]){
  name <- toString(i)
  has30 <- sum(coverage[,name]>=30)
  has50 <- sum(coverage[,name]>=50)
  has100 <- sum(coverage[,name] >= 99.9)
  cov_vs_dis[nrow(cov_vs_dis) + 1,] = list(name, has100, has50, has30)
}

str(cov_vs_dis)
cov_vs_dis <- cbind(cov_vs_dis, rel_coverage[,3])
names(cov_vs_dis) <- c("sample", "has100", "has50", "has30", "distance")

#Genes with 100% coverage vs distance 
plot(cov_vs_dis$distance, cov_vs_dis$has100, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 100% coverage", col=rgb(0, 0, 1, 1))

#Genes with 50% coverage vs distance
plot(cov_vs_dis$distance, cov_vs_dis$has50, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 50% coverage", col=rgb(0, 0, 1, 1))

#Genes with 40% coverage vs distance
plot(cov_vs_dis$distance, cov_vs_dis$has30, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 30% coverage", col=rgb(0, 0, 1, 1))

#Genes with 100% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has100, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of genes with 100% coverage", col=rgb(0, 0, 1, 1))

#Genes with 50% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has50, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of genes with 50% coverage", col=rgb(0, 0, 1, 1))

#Genes with 40% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has30, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of genes with 30% coverage", col=rgb(0, 0, 1, 1))



#################################################################################3
#############################################################################
## PHYLOGENETIC DISTANCE VS EXON COVERAGE, NEW DISTANCE VALUES ########
#################################################################


library(readr)
setwd('.')
exons <- read_delim("exon_length_and_id2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
str(exons)

coverage <- exons[grep("coverage", names(exons))]
str(coverage)

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)

species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 3, nrow = 0))
names(name_species) <- c("name", "species", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",0.05)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",0.082)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",0.05)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",0.025)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",0.068)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",0.054)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",0.05)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",0.05)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",0.05)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",0.112)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",0.103)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,3])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "XPIP23", "CraneSp")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
rel_coverage <- rel_coverage[-remove_this,]

cov_vs_dis <- data.frame(matrix(ncol = 3, nrow = 0))
names(cov_vs_dis) <- c("sample", "has100", "has50")

i <- toString(rel_coverage[1,2])
x <- coverage[,i]
for(i in rel_coverage[,2]){
  name <- toString(i)
  has50 <- sum(coverage[,name]>=50)
  has100 <- sum(coverage[,name] >= 99.9)
  cov_vs_dis[nrow(cov_vs_dis) + 1,] = list(name, has100, has50)
}

str(cov_vs_dis)
cov_vs_dis <- cbind(cov_vs_dis, rel_coverage[,3])
names(cov_vs_dis) <- c("sample", "has100", "has50", "distance")

#Genes with 100% coverage vs distance 
plot(cov_vs_dis$distance, cov_vs_dis$has100, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 100% coverage", col=rgb(0, 0, 1, 1))

#Genes with 50% coverage vs distance
plot(cov_vs_dis$distance, cov_vs_dis$has50, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 50% coverage", col=rgb(0, 0, 1, 1))


###################################################################
## PHYLOGENETIC DISTANCE VS EXON COVERAGE, NEWER DISTANCE VALUES ##
###################################################################


library(readr)
setwd('.')
exons <- read_delim("exon_length_and_id2.tsv", "\t", escape_double = FALSE, col_types = cols(length = col_integer()), trim_ws = TRUE, skip = 1)
str(exons)

coverage <- exons[grep("coverage|length", names(exons))]

coverage <- coverage[which(coverage$length>=120),]

str(coverage)

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)

species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 4, nrow = 0))
names(name_species) <- c("name", "species", "similarity", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",92.484,7.516)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",(95.816+95.198)/2,(4.184+4.802)/2)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",97.286,2.714)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",93.933,6.067)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",95.198,4.802)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",89.77,10.23)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA,NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",89.77,10.23)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,4])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "CXPIP23", "CraneSp", "PARUS1_PHSIB1")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
#remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
#rel_coverage <- rel_coverage[-remove_this,]

cov_vs_dis <- data.frame(matrix(ncol = 4, nrow = 0))
names(cov_vs_dis) <- c("sample", "has100", "has50", "has40")

i <- toString(rel_coverage[3,2])
x <- coverage[,i]
for(i in rel_coverage[,2]){
  name <- toString(i)
  has30 <- sum(coverage[,name]>=30)
  has50 <- sum(coverage[,name]>=50)
  has100 <- sum(coverage[,name] >= 99.9)
  cov_vs_dis[nrow(cov_vs_dis) + 1,] = list(name, has100, has50, has30)
}

str(cov_vs_dis)
cov_vs_dis <- cbind(cov_vs_dis, rel_coverage[,3])
names(cov_vs_dis) <- c("sample", "has100", "has50", "has30", "distance")

#Exons with 100% coverage vs distance 
plot(jitter(cov_vs_dis$distance,3), jitter(cov_vs_dis$has100,5), pch=".", cex=5, xlab="CytB distance (%) from SISKIN1", ylab="Number of exons with 100% coverage", col=rgb(0, 0, 1, 1))

#Exons with 50% coverage vs distance
plot(cov_vs_dis$distance, cov_vs_dis$has50, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of exons with 50% coverage", col=rgb(0, 0, 1, 1))

#Exons with 40% coverage vs distance
plot(cov_vs_dis$distance, cov_vs_dis$has30, pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of exons with 30% coverage", col=rgb(0, 0, 1, 1))

#Exons with 100% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has100, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of exons with 100% coverage", col=rgb(0, 0, 1, 1))

#Exons with 50% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has50, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of exons with 50% coverage", col=rgb(0, 0, 1, 1))

#Exons with 40% coverage vs distance, logarithmic
plot(cov_vs_dis$distance, cov_vs_dis$has30, pch=".", cex=5, log="y", xlab="SISKIN1 Distance", ylab="Number of exons with 30% coverage", col=rgb(0, 0, 1, 1))

###############
### ABLINE LOGARITHMIC
#Genes with 100% coverage vs distance



y <- cov_vs_dis$has100
x <- cov_vs_dis$distance
fit <- lm(log(y)~x)
distances <- seq(0, 10, 0.1)
fit2 <- predict(fit, newdata=list(x=distances))
fit3 <- exp(fit2)

coef(fit)
s <- summary(fit)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

rp = vector('expression',3)
if(coef(fit)[2]>=0){
  rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                     list(a_val = format(coef(fit)[1], digits = 2), 
                          b_val = format(coef(fit)[2], digits = 3)))[2]
} else {
  rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                     list(a_val = format(coef(fit)[1], digits = 2), 
                          b_val = format(0 - coef(fit)[2], digits = 3)))[2]
}
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]

plot(jitter(x,5), jitter(y,5), pch=".", yaxt="n", cex=7,log="y", ylim=c(1,2000), xlab="CytB distance (%) from SISKIN1", ylab="Number of exons with 100% coverage", col=rgb(0, 0, 1, 1))
lines(distances, fit3,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
axis(2, at=c(1,2,5,10,20,50,100,200,500,1000,2000),labels=TRUE)
legend('topright', legend = rp, bty = 'n')


##########
## And with 50% coverage

y <- cov_vs_dis$has50
x <- cov_vs_dis$distance

fit <- lm(log(y)~x)
distances <- seq(0, 10, 0.1)
fit2 <- predict(fit, newdata=list(x=distances))
fit3 <- exp(fit2)


coef(fit)
s <- summary(fit)
r2 <- s$adj.r.squared
r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
p <- s$coefficients[2,4]
plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
?expression
rp = vector('expression',3)
if(coef(fit)[2]>=0){
  rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                     list(a_val = format(coef(fit)[1], digits = 2), 
                          b_val = format(coef(fit)[2], digits = 3)))[2]
} else {
rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                   list(a_val = format(coef(fit)[1], digits = 2), 
                        b_val = format(0 - coef(fit)[2], digits = 3)))[2]
}
rp[2] = substitute(expression(italic(R)^2 == r_val), 
                   list(r_val = format(r2,dig=3)))[2]
rp[3] = substitute(expression(italic(p) == p_val), 
                   list(p_val = format(p, digits = 2)))[2]
?jitter
plot(jitter(x,5), jitter(y,5), pch=".", yaxt="n", cex=7,log="y", ylim=c(10, 3000), xlab="CytB distance (%) from SISKIN1", ylab="Number of exons with 50% coverage", col=rgb(0, 0, 1, 1))
lines(distances, fit3,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
axis(2, at=c(1,2,5,10,20,50,100,200,500,1000,2000),labels=TRUE)
legend('topright', legend = rp, bty = 'n')


#################################
### CATEGORY-SPECIFIC IDENTITY ###
#################################


library(readr)
library(data.table)
#install.packages("data.table")

setwd('.')
genes <- read_delim("Consensus/gene_length_and_id.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
mips <- read_delim("genes_vs_mips2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(mips) <- c("gene","mips")
#NOTE: Add CoxI etc to gene definition
str(genes)

coverage <- genes[grep("_identity|GC|Gene_Name", names(genes))]
colnames(coverage)
str(coverage)
coverage$`#Gene_Name`
mips$gene
coverage <- merge(coverage, mips, by.x = "#Gene_Name", by.y = "gene")

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)
mips_types = c("01", "02", "10", "11", "12", "14", "16", "18", "20", "30", "32", "34", "36", "41", "42", "43", "47", "70", "73", "77")
species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 4, nrow = 0))
names(name_species) <- c("name", "species", "similarity", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",92.484,7.516)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA,NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",97.286,2.714)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",93.933,6.067)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",95.198,4.802)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",89.77,10.23)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA,NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",89.77,10.23)
name_species
paste(name_species[,1],"_identity", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_identity", sep=""), name_species[,4])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "CXPIP23", "CraneSp")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_identity")
rel_coverage <- rel_coverage[-remove_this,]

description <- data.frame(matrix(ncol = 2, nrow = 0))
names(description) <- c("class", "text")
description[nrow(description) + 1,] = list("01", "METABOLISM")
description[nrow(description) + 1,] = list("02", "ENERGY")
description[nrow(description) + 1,] = list("10", "CELL CYCLE AND DNA PROCESSING") 
description[nrow(description) + 1,] = list("11", "TRANSCRIPTION")
description[nrow(description) + 1,] = list("12", "PROTEIN SYNTHESIS")
description[nrow(description) + 1,] = list("14", "PROTEIN FATE")
description[nrow(description) + 1,] = list("16", "PROTEIN WITH BINDING FUNCTION OR COFACTOR")
description[nrow(description) + 1,] = list("18", "REGUL. METABOLISM AND PROTEIN FUNCTION")
description[nrow(description) + 1,] = list("20", "CELLULAR TRANSPORT, FACILITIES, ROUTES")
description[nrow(description) + 1,] = list("30", "COMMUNICATION/SIGNAL TRANSDUCTION")
description[nrow(description) + 1,] = list("32", "CELL RESCUE, DEFENSE AND VIRULENCE")
description[nrow(description) + 1,] = list("34", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("36", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("41", "DEVELOPMENT (Systemic)")
description[nrow(description) + 1,] = list("42", "BIOGENESIS OF CELLULAR COMPONENTS")
description[nrow(description) + 1,] = list("43", "CELL TYPE DIFFERENTIATION")
description[nrow(description) + 1,] = list("47", "ORGAN DIFFERENTIATION")
description[nrow(description) + 1,] = list("70", "SUBCELLULAR LOCALIZATION")
description[nrow(description) + 1,] = list("73", "CELL TYPE LOCALIZATION")
description[nrow(description) + 1,] = list("77", "ORGAN LOCALIZATION")
description

mips_coverage <- list()
i <- mips_types[1]
for(i in mips_types){
  x <- coverage[coverage$mips %like% i, ] 
  cov_vs_dis <- data.frame(matrix(ncol = 4, nrow = 0))
  names(cov_vs_dis) <- c("sample", "has100", "has50", "has30")
  j <- toString(rel_coverage[3,2])
  for(j in rel_coverage[,2]){
    name <- toString(j)
    has30 <- sum(x[,name]>=30)
    has50 <- sum(x[,name]>=50)
    has100 <- sum(x[,name] >= 99.9)
    cov_vs_dis[nrow(cov_vs_dis) + 1,] = list(name, has100, has50, has30)
  }
  cov_vs_dis <- cbind(cov_vs_dis, rel_coverage[,3])
  names(cov_vs_dis) <- c("sample", "has100", "has50", "has30", "distance")
  mips_coverage[[i]] <- cov_vs_dis 
}
str(mips_coverage)


#Genes with 100% identity vs distance 
i <- mips_types[1]
x <- mips_coverage[[i]]
mips_coverage
plot(jitter(x$distance,2), jitter(x$has50,2), pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 100% identity", col=rgb(0, 0, 1, 1))

### ABLINE LOGARITHMIC
#Genes with 100% coverage vs distance

### SORT BY MIN_VALUE
class_vs_number <- data.frame(matrix(ncol = 2, nrow = 0))
names(class_vs_number) <- c("class", "number")
i <- "16"
for(i in mips_types){
y1 <- mips_coverage[[i]]
y2 <- y1[c(4,13,14),]
number <- sum(y2$has100)
class_vs_number[nrow(class_vs_number) + 1,] = list(i, number)
}
class_vs_number <- class_vs_number[with(class_vs_number, order(-number)), ]


### PLOT
name <- "102-11C_S3_L001_identity"
ymax <- sum(genes[,name] >= 99.9)
ymin <- 0
xmax <- max(spec_sam3$distance, na.rm=TRUE)
xmin <- 0

plot(x, y, pch=".", cex=5,log="y", xlim=c(xmin,xmax+3), ylim=c(0.5,200), xlab="SISKIN1 Distance", ylab="Number of genes with 100% identity", col=rgb(0, 0, 0, 0))
color_vs_legend <- data.frame(matrix(ncol = 4, nrow = 0))
names(color_vs_legend) <- c("class", "color", "amount", "text")
color_vs_legend
palette(rainbow(14)) 
i <- "01"
color<-1
for(i in class_vs_number$class){
  y <- mips_coverage[[i]]$has100
  x <- mips_coverage[[i]]$distance
  number <- sum(y)
  if(number>0){
    points_and_line(x,y, color)
    descr <- subset(description, class %in% i)$text
    color_vs_legend[nrow(color_vs_legend) + 1,] = list(i, color, number, descr)
    color <- (color+1)
  }
}
color_vs_legend
legend('bottomright', bty='n', legend = color_vs_legend$class, col = "black", cex = 1.0, pch = 21, pt.bg=color_vs_legend$color)

### ABLINE LOGARITHMIC
#Genes with 50% coverage vs distance

### FUNCTION POINTS_AND_LINE

points_and_line <- function(x, y, color=rgb(0, 0, 1, 1)){
  
  logy <- log(y)
  logy <- replace(logy, is.infinite(logy), NA)
  
  fit <- lm(logy~x, na.action = na.omit)
  distances <- seq(0, 10, 0.1)
  fit2 <- predict(fit, newdata=list(x=distances))
  fit3 <- exp(fit2)
  
  coef(fit)
  s <- summary(fit)
  r2 <- s$adj.r.squared
  r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p <- s$coefficients[2,4]
  plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  
  rp = vector('expression',3)
  if(coef(fit)[2]>=0){
    rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                       list(a_val = format(coef(fit)[1], digits = 2), 
                            b_val = format(coef(fit)[2], digits = 3)))[2]
  } else {
    rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                       list(a_val = format(coef(fit)[1], digits = 2), 
                            b_val = format(0 - coef(fit)[2], digits = 3)))[2]
  }
  rp[2] = substitute(expression(italic(R)^2 == r_val), 
                     list(r_val = format(r2,dig=3)))[2]
  rp[3] = substitute(expression(italic(p) == p_val), 
                     list(p_val = format(p, digits = 2)))[2]
  
  points(jitter(x), jitter(y), pch=".", cex=5, col=color)
  lines(distances, fit3,lwd=3, col = color, xlab = "Time (s)", ylab = "Counts")
  #legend('topright', legend = rp, bty = 'n')
}

### SORT BY MIN_VALUE
class_vs_number <- data.frame(matrix(ncol = 2, nrow = 0))
names(class_vs_number) <- c("class", "number")
i <- "16"
for(i in mips_types){
y1 <- mips_coverage[[i]]
y2 <- y1[c(4,13,14),]
number <- sum(y2$has50)
class_vs_number[nrow(class_vs_number) + 1,] = list(i, number)
}
class_vs_number <- class_vs_number[with(class_vs_number, order(-number)), ]


### PLOT 
plot(x, y, pch=".", cex=5,log="y", xlim=c(xmin,xmax+3), ylim=c(0.5,500), xlab="SISKIN1 Distance", ylab="Number of genes with 50% identity", col=rgb(0, 0, 0, 0))
color_vs_legend <- data.frame(matrix(ncol = 4, nrow = 0))
names(color_vs_legend) <- c("class", "color", "amount", "text")
color_vs_legend
palette(rainbow(14)) 
i <- "34"
color<-3
mips_coverage
for(i in class_vs_number$class){
  y <- mips_coverage[[i]]$has50
  x <- mips_coverage[[i]]$distance
  number <- sum(y)
  if(number>0){
    points_and_line(x,y, color)
    descr <- subset(description, class %in% i)$text
    color_vs_legend[nrow(color_vs_legend) + 1,] = list(i, color, number, descr)
    color <- (color+1)
  }
}
color_vs_legend
legend('right', bty='n', legend = color_vs_legend$class, col = "black", cex = 1.2, pch = 21, pt.bg=color_vs_legend$color)





#################################
### CATEGORY-SPECIFIC COVERAGE ###
#################################


library(readr)
library(data.table)
#install.packages("data.table")

setwd('.')
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(Gene_Length = col_integer()), trim_ws = TRUE, skip = 1)
mips <- read_delim("genes_vs_mips2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(mips) <- c("gene","mips")
#NOTE: Add CoxI etc to gene definition
str(genes)

coverage <- genes[grep("_coverage|GC|Gene_Name|Gene_Length", names(genes))]
colnames(coverage)
str(coverage)
coverage$`#Gene_Name`
mips$gene
coverage <- merge(coverage, mips, by.x = "#Gene_Name", by.y = "gene")

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)
mips_types = c("01", "02", "10", "11", "12", "14", "16", "18", "20", "30", "32", "34", "36", "41", "42", "43", "47", "70", "73", "77", "host")
species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 4, nrow = 0))
names(name_species) <- c("name", "species", "similarity", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",92.484,7.516)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA,NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",97.286,2.714)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",93.933,6.067)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",95.198,4.802)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",89.77,10.23)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA,NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",89.77,10.23)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,4])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "CXPIP23", "CraneSp")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
#remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
#rel_coverage <- rel_coverage[-remove_this,]

description <- data.frame(matrix(ncol = 2, nrow = 0))
names(description) <- c("class", "text")
description[nrow(description) + 1,] = list("01", "METABOLISM")
description[nrow(description) + 1,] = list("02", "ENERGY")
description[nrow(description) + 1,] = list("10", "CELL CYCLE AND DNA PROCESSING") 
description[nrow(description) + 1,] = list("11", "TRANSCRIPTION")
description[nrow(description) + 1,] = list("12", "PROTEIN SYNTHESIS")
description[nrow(description) + 1,] = list("14", "PROTEIN FATE")
description[nrow(description) + 1,] = list("16", "PROTEIN WITH BINDING FUNCTION OR COFACTOR")
description[nrow(description) + 1,] = list("18", "REGUL. METABOLISM AND PROTEIN FUNCTION")
description[nrow(description) + 1,] = list("20", "CELLULAR TRANSPORT, FACILITIES, ROUTES")
description[nrow(description) + 1,] = list("30", "COMMUNICATION/SIGNAL TRANSDUCTION")
description[nrow(description) + 1,] = list("32", "CELL RESCUE, DEFENSE AND VIRULENCE")
description[nrow(description) + 1,] = list("34", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("36", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("41", "DEVELOPMENT (Systemic)")
description[nrow(description) + 1,] = list("42", "BIOGENESIS OF CELLULAR COMPONENTS")
description[nrow(description) + 1,] = list("43", "CELL TYPE DIFFERENTIATION")
description[nrow(description) + 1,] = list("47", "ORGAN DIFFERENTIATION")
description[nrow(description) + 1,] = list("70", "SUBCELLULAR LOCALIZATION")
description[nrow(description) + 1,] = list("73", "CELL TYPE LOCALIZATION")
description[nrow(description) + 1,] = list("77", "ORGAN LOCALIZATION")
description[nrow(description) + 1,] = list("host", "HOST INTERACTION")
description
host_interaction<-c('HtGene0063', 'HtGene0318', 'HtGene0562', 'HtGene0708', 'HtGene0721', 'HtGene0852', 'HtGene0856', 'HtGene0925', 'HtGene0958', 'HtGene0998', 'HtGene1040', 'HtGene1453', 'HtGene1574', 'HtGene1595', 'HtGene1657', 'HtGene1707', 'HtGene1957', 'HtGene1994', 'HtGene2049', 'HtGene2264', 'HtGene2426', 'HtGene2433', 'HtGene2524', 'HtGene3169', 'HtGene3484', 'HtGene3996', 'HtGene4101', 'HtGene4136', 'HtGene4137', 'HtGene4607', 'HtGene6386')
#sum(coverage[colnames(coverage)==i])
i <- 'HtGene0318'
for(i in host_interaction)
{
  x <- coverage[which(coverage$`#Gene_Name`==i),]$mips
  x <- paste(x, "host")
  coverage[which(coverage$`#Gene_Name`==i),]$mips <- x
}

### CALCULATE PERCENTAGE
#coverage_percentage <- coverage[grep("GC|Gene_Name|Gene_Length|mips", names(coverage))]
#i <- rel_coverage$name[1]
#for(i in rel_coverage$name){
#  x <- 100*coverage[colnames(coverage)==i]/coverage$Gene_Length
#  coverage_percentage <- cbind(coverage_percentage, x)
#}
#str(coverage_percentage)

### CALCULATE HAS100, HAS50, HAS30
mips_coverage <- list()
i <- "32"
for(i in mips_types){
  x <- coverage[coverage$mips %like% i, ]
  cov_vs_dis <- data.frame(matrix(ncol = 4, nrow = 0))
  names(cov_vs_dis) <- c("sample", "has100", "has50", "has30")
  j <- toString(rel_coverage[3,2])
  for(j in rel_coverage[,2]){
    name <- toString(j)
    has30 <- sum(x[,name]>=30)
    has50 <- sum(x[,name]>=50)
    has100 <- sum(x[,name] >= 99.9)
    cov_vs_dis[nrow(cov_vs_dis) + 1,] = list(name, has100, has50, has30)
  }
  cov_vs_dis <- cbind(cov_vs_dis, rel_coverage[,3])
  names(cov_vs_dis) <- c("sample", "has100", "has50", "has30", "distance")
  mips_coverage[[i]] <- cov_vs_dis 
}
str(mips_coverage)


#Genes with 100% coverage vs distance 
i <- mips_types[1]
x <- mips_coverage[[i]]
mips_coverage
plot(jitter(x$distance,2), jitter(x$has50,4), pch=".", cex=5, xlab="SISKIN1 Distance", ylab="Number of genes with 100% coverage", col=rgb(0, 0, 1, 1))

### ABLINE LOGARITHMIC
#Genes with 100% coverage vs distance

### SORT BY MIN_VALUE
class_vs_number <- data.frame(matrix(ncol = 2, nrow = 0))
names(class_vs_number) <- c("class", "number")
i <- "16"
for(i in mips_types){
y1 <- mips_coverage[[i]]
y2 <- y1[c(4,13,14),]
number <- sum(y2$has100)
class_vs_number[nrow(class_vs_number) + 1,] = list(i, number)
}
class_vs_number <- class_vs_number[with(class_vs_number, order(-number)), ]

### PLOTTING FUNCTION 

points_and_line <- function(x, y, color=rgb(0, 0, 1, 1)){
  
  logy <- log(y)
  logy <- replace(logy, is.infinite(logy), NA)
  
  fit <- lm(logy~x, na.action = na.omit)
  distances <- seq(0, 10, 0.1)
  fit2 <- predict(fit, newdata=list(x=distances))
  fit3 <- exp(fit2)
  
  coef(fit)
  s <- summary(fit)
  r2 <- s$adj.r.squared
  r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  p <- s$coefficients[2,4]
  plabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  
  rp = vector('expression',3)
  if(coef(fit)[2]>=0){
    rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                       list(a_val = format(coef(fit)[1], digits = 2), 
                            b_val = format(coef(fit)[2], digits = 3)))[2]
  } else {
    rp[1] = substitute(expression(italic(log(y)) == a_val - b_val %.% italic(x) ), 
                       list(a_val = format(coef(fit)[1], digits = 2), 
                            b_val = format(0 - coef(fit)[2], digits = 3)))[2]
  }
  rp[2] = substitute(expression(italic(R)^2 == r_val), 
                     list(r_val = format(r2,dig=3)))[2]
  rp[3] = substitute(expression(italic(p) == p_val), 
                     list(p_val = format(p, digits = 2)))[2]
  
  points(jitter(x), jitter(y), pch=".", cex=5, col=color)
  lines(distances, fit3,lwd=3, col = color, xlab = "Time (s)", ylab = "Counts")
  #legend('topright', legend = rp, bty = 'n')
  return(fit)
}

### PLOT

name <- "102-11C_S3_L001_coverage"
ymax <- sum(genes[,name] >= 99.9)
ymin <- 0
xmax <- max(spec_sam3$distance, na.rm=TRUE)
xmin <- 0
x <- 0
y <- 1
#?plot
plot(x, y, pch=".", cex=5,log="y", xlim=c(1,xmax+3), ylim=c(0.5,1000), xlab="SISKIN1 Distance", ylab="Number of genes with 100% coverage", col=rgb(0, 0, 0, 0))
#?axis
#axis(2, at=c(1, 10, 100))
color_vs_legend <- data.frame(matrix(ncol = 6, nrow = 0))
names(color_vs_legend) <- c("class", "color", "amount", "text", "slope","p")
color_vs_legend
palette(rainbow(14)) 
i <- "32"
color<-1
for(i in class_vs_number$class){
  y <- mips_coverage[[i]]$has100
  x <- mips_coverage[[i]]$distance
  number <- sum(y)
  if(number>0){
    fit <- points_and_line(x,y, color)
    slope <- coef(fit)[2]
    p <- summary(fit)$coefficients[2,4]
    descr <- subset(description, class %in% i)$text
    color_vs_legend[nrow(color_vs_legend) + 1,] = list(i, color, number, descr, round(slope, digits=2), round(p, digits=4))
    color <- (color+1)
  }
}
color_vs_legend
color_vs_legend$p
legend_strings <- paste(color_vs_legend$class, ' (',color_vs_legend$slope, 'x)', ' p=',color_vs_legend$p, sep='' )
legend('right', bty='n', legend = legend_strings, col = "black", cex = 1.0, pch = 21, pt.bg=color_vs_legend$color)

#Save graph
folder='graphs/'
name="mips_coverage_has100_2"
date=toString(Sys.Date())
dev.copy(png, paste(folder, name, '_', date, sep=''))
dev.off()


### ABLINE LOGARITHMIC
#Genes with 50% coverage vs distance

### SORT BY MIN_VALUE
class_vs_number <- data.frame(matrix(ncol = 2, nrow = 0))
names(class_vs_number) <- c("class", "number")
i <- "16"
for(i in mips_types){
y1 <- mips_coverage[[i]]
y2 <- y1[c(4,13,14),]
number <- sum(y2$has50)
class_vs_number[nrow(class_vs_number) + 1,] = list(i, number)
}
class_vs_number <- class_vs_number[with(class_vs_number, order(-number)), ]


### PLOT

name <- "102-11C_S3_L001_coverage"
ymax <- sum(genes[,name] >= 99.9)
ymin <- 0
xmax <- max(spec_sam3$distance, na.rm=TRUE)
xmin <- 0
?plot
x <- 0
y <- 1
plot(x, y, pch=".", cex=5,log="y", xlim=c(1,xmax+4), ylim=c(0.5,1000), xlab="SISKIN1 Distance", ylab="Number of genes with 50% coverage", col=rgb(0, 0, 0, 0))
?axis
#axis(2, at=c(1, 10, 100))
color_vs_legend <- data.frame(matrix(ncol = 6, nrow = 0))
names(color_vs_legend) <- c("class", "color", "amount", "text", "slope","p")
color_vs_legend
palette(rainbow(14)) 
i <- "host"
color<-1
for(i in class_vs_number$class){
  y <- mips_coverage[[i]]$has50
  x <- mips_coverage[[i]]$distance
  number <- sum(y)
  if(number>0){
    fit <- points_and_line(x,y, color)
    slope <- coef(fit)[2]
    p <- summary(fit)$coefficients[2,4]
    descr <- subset(description, class %in% i)$text
    color_vs_legend[nrow(color_vs_legend) + 1,] = list(i, color, number, descr, round(slope, digits=2), round(p, digits=4))
    color <- (color+1)
  }
}
color_vs_legend
color_vs_legend$p
legend_strings <- paste(color_vs_legend$class, ' (',color_vs_legend$slope, 'x)', ' p=',color_vs_legend$p, sep='' )
legend('right', bty='n', legend = legend_strings, col = "black", cex = 1.0, pch = 21, pt.bg=color_vs_legend$color)

#Save graph
folder='/home/roland/Documents/graphs/'
name="mips_coverage_has50_2"
date=toString(Sys.Date())
dev.copy(png, paste(folder, name, '_', date, sep=''))
dev.off()



##########################
### SIMPLER MIPS GRAPH ###
##########################

library(readr)
library(data.table)
#install.packages("data.table")

setwd('.')
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(Gene_Length = col_integer()), trim_ws = TRUE, skip = 1)
mips <- read_delim("genes_vs_mips2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(mips) <- c("gene","mips")
#NOTE: Add CoxI etc to gene definition
str(genes)

coverage <- genes[grep("_coverage|GC|Gene_Name|Gene_Length|mips", names(genes))]
colnames(coverage)
str(coverage)
coverage$`#Gene_Name`
mips$gene
coverage <- merge(coverage, mips, by.x = "#Gene_Name", by.y = "gene")

#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)
mips_types = c("01", "02", "10", "11", "12", "14", "16", "18", "20", "30", "32", "34", "36", "41", "42", "43", "47", "70", "73", "77", "host")
species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 4, nrow = 0))
names(name_species) <- c("name", "species", "similarity", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",92.484,7.516)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA,NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",97.286,2.714)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",93.933,6.067)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",95.198,4.802)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",89.77,10.23)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA,NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",89.77,10.23)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,4])
names(spec_sam3) <- c("species", "name", "distance")

relevant_species <- c("SISKIN1", "PARUS1")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
rel_coverage <- rel_coverage[-remove_this,]

description <- data.frame(matrix(ncol = 2, nrow = 0))
names(description) <- c("class", "text")
description[nrow(description) + 1,] = list("01", "METABOLISM")
description[nrow(description) + 1,] = list("02", "ENERGY")
description[nrow(description) + 1,] = list("10", "CELL CYCLE AND DNA PROCESSING") 
description[nrow(description) + 1,] = list("11", "TRANSCRIPTION")
description[nrow(description) + 1,] = list("12", "PROTEIN SYNTHESIS")
description[nrow(description) + 1,] = list("14", "PROTEIN FATE")
description[nrow(description) + 1,] = list("16", "PROTEIN WITH BINDING FUNCTION OR COFACTOR")
description[nrow(description) + 1,] = list("18", "REGUL. METABOLISM AND PROTEIN FUNCTION")
description[nrow(description) + 1,] = list("20", "CELLULAR TRANSPORT, FACILITIES, ROUTES")
description[nrow(description) + 1,] = list("30", "COMMUNICATION/SIGNAL TRANSDUCTION")
description[nrow(description) + 1,] = list("32", "CELL RESCUE, DEFENSE AND VIRULENCE")
description[nrow(description) + 1,] = list("34", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("36", "INTERACTION WITH THE ENVIRONMENT")
description[nrow(description) + 1,] = list("41", "DEVELOPMENT (Systemic)")
description[nrow(description) + 1,] = list("42", "BIOGENESIS OF CELLULAR COMPONENTS")
description[nrow(description) + 1,] = list("43", "CELL TYPE DIFFERENTIATION")
description[nrow(description) + 1,] = list("47", "ORGAN DIFFERENTIATION")
description[nrow(description) + 1,] = list("70", "SUBCELLULAR LOCALIZATION")
description[nrow(description) + 1,] = list("73", "CELL TYPE LOCALIZATION")
description[nrow(description) + 1,] = list("77", "ORGAN LOCALIZATION")
description[nrow(description) + 1,] = list("host", "HOST INTERACTION")
description
host_interaction<-c('HtGene0063', 'HtGene0318', 'HtGene0562', 'HtGene0708', 'HtGene0721', 'HtGene0852', 'HtGene0856', 'HtGene0925', 'HtGene0958', 'HtGene0998', 'HtGene1040', 'HtGene1453', 'HtGene1574', 'HtGene1595', 'HtGene1657', 'HtGene1707', 'HtGene1957', 'HtGene1994', 'HtGene2049', 'HtGene2264', 'HtGene2426', 'HtGene2433', 'HtGene2524', 'HtGene3169', 'HtGene3484', 'HtGene3996', 'HtGene4101', 'HtGene4136', 'HtGene4137', 'HtGene4607', 'HtGene6386')
#sum(coverage[colnames(coverage)==i])
i <- 'HtGene0318'
for(i in host_interaction)
{
  x <- coverage[which(coverage$`#Gene_Name`==i),]$mips
  x <- paste(x, "host")
  coverage[which(coverage$`#Gene_Name`==i),]$mips <- x
  }

### CALCULATE PERCENTAGE
#coverage_percentage <- coverage[grep("GC|Gene_Name|Gene_Length|mips", names(coverage))]
#i <- rel_coverage$name[6]
#for(i in c("102-11C_S3_L001_coverage", "CC50362-3_S4_L001_coverage")){
#  x <- 100*coverage[colnames(coverage)==i]/coverage$Gene_Length
#  coverage_percentage <- cbind(coverage_percentage, x)
#}
str(coverage_percentage)


### CALCULATE HAS100, HAS50, HAS30
SISKIN_vs_samples_100 <- list()
PARUS_vs_samples_100 <- list()
SISKIN_vs_samples_50 <- list()
PARUS_vs_samples_50 <- list()

i <- "32"
for(i in mips_types){
  x <- coverage[coverage$mips %like% i, ]
  SISKIN_has100 <- sum(x[, "102-11C_S3_L001_coverage"]>=99.9)
  PARUS_has100 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=99.9)
  SISKIN_has50 <- sum(x[, "102-11C_S3_L001_coverage"]>=50)
  PARUS_has50 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=50)
  SISKIN_vs_samples_100[[i]] <- SISKIN_has100 
  PARUS_vs_samples_100[[i]] <- PARUS_has100
  SISKIN_vs_samples_50[[i]] <- SISKIN_has50 
  PARUS_vs_samples_50[[i]] <- PARUS_has50
}
str(SISKIN_vs_samples_100) 

#i <- "host"
#x <- coverage[coverage$mips %like% i, ]
#SISKIN_has100 <- sum(x[, "102-11C_S3_L001_coverage"]>=99.9)
#PARUS_has100 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=99.9)
#SISKIN_has50 <- sum(x[, "102-11C_S3_L001_coverage"]>=50)
#PARUS_has50 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=50)
#SISKIN_vs_samples_100[[i]] <- SISKIN_has100 
#PARUS_vs_samples_100[[i]] <- PARUS_has100
#SISKIN_vs_samples_50[[i]] <- SISKIN_has50 
#PARUS_vs_samples_50[[i]] <- PARUS_has50

mips_numbers <- as.numeric(c(SISKIN_vs_samples_100, PARUS_vs_samples_100))
mips_samples <- c(rep("SISKIN", 21), rep("PARUS", 21) )
mips_names <- c(mips_types, mips_types)
mips_data <- data.frame(mips_names, mips_numbers, mips_samples)
mips_data

### PLOT

library(ggplot2)
#install.packages("ggplot2")
p <-ggplot(mips_data, aes(mips_names, mips_numbers )) +geom_bar(stat = "identity", aes(fill = mips_samples), position = "dodge") + coord_flip()
p

### HAS_50
mips_numbers <- as.numeric(c(SISKIN_vs_samples_50, PARUS_vs_samples_50))
mips_samples <- c(rep("SISKIN", 21), rep("PARUS", 21) )
mips_names <- c(mips_types, mips_types)
mips_data <- data.frame(mips_names, mips_numbers, mips_samples)
mips_data
p <-ggplot(mips_data, aes(mips_names, mips_numbers )) +geom_bar(stat = "identity", aes(fill = mips_samples), position = "dodge") + coord_flip()
p

mips_types_b <- c()
SISKIN_vs_samples_100 <- list()
PARUS_vs_samples_100 <- list()
SISKIN_vs_samples_50 <- list()
PARUS_vs_samples_50 <- list()
normalized_PARUS_50 <- list()

i <- mips_types[1]
for(i in mips_types){
  x <- coverage[coverage$mips %like% i, ]

  SISKIN_has100 <- sum(x[, "102-11C_S3_L001_coverage"]>=99.9)
  PARUS_has100 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=99.9)
  SISKIN_has50 <- sum(x[, "102-11C_S3_L001_coverage"]>=50)
  PARUS_has50 <- sum(x[, "CC50362-3_S4_L001_coverage"]>=50)
  if ( ! SISKIN_has50 > 0) {next}
  if ( ! PARUS_has50 > 0) {next}
  mips_types_b <- c(mips_types_b, i)
  SISKIN_vs_samples_100[[i]] <- SISKIN_has100 
  PARUS_vs_samples_100[[i]] <- PARUS_has100
  SISKIN_vs_samples_50[[i]] <- SISKIN_has50 
  PARUS_vs_samples_50[[i]] <- PARUS_has50
  normalized_PARUS_50[[i]] <- 100*PARUS_has50/SISKIN_has50
}

### NORMALIZED PARUS

mips_numbers <- as.numeric(c(rep("100", length(mips_types_b)), normalized_PARUS_50))
mips_samples <- c(rep("SISKIN", length(mips_types_b)), rep("PARUS", length(mips_types_b)) )
mips_names <- c(mips_types_b, mips_types_b)
mips_description <- subset(description, class %in% mips_types_b)$text

i<- description[1,]
for(i in 1:nrow(description)){
  txt <- description[i,]$text
  cls <- description[i,]$class
if( cls %in% mips_types_b ){mips_description <- c(mips_description, txt) }
}
mips_names <- c(mips_description, mips_description)
mips_data <- data.frame(mips_names, mips_numbers, mips_samples)

p <-ggplot(mips_data, aes(mips_names, mips_numbers )) +geom_bar(stat = "identity", aes(fill = mips_samples), position = "dodge") + coord_flip()
p


###########################################
### SAMPLE MULTIPLICITY ###
###########################

library(readr)
library(data.table)
#install.packages("data.table")

setwd('.')
genes <- read_delim("gene_length_and_id.tsv", "\t", escape_double = FALSE, col_types = cols(Gene_Length = col_integer()), trim_ws = TRUE, skip = 1)
#NOTE: Add CoxI etc to gene definition
str(genes)

coverage <- genes[grep("_coverage|GC|Gene_Name|Gene_Length", names(genes))]
colnames(coverage)
str(coverage)
coverage$`#Gene_Name`


#samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
#str(samples3)
mips_types = c("01", "02", "10", "11", "12", "14", "16", "18", "20", "30", "32", "34", "36", "41", "42", "43", "47", "70", "73", "77")
species <- c("SISKIN1", "SISKIN1", "EMCIR1","PARUS1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")
names <- c("102-11C_S3_L001", "126-11c_S9_L001", "1ES86798_S2_L001" 
           ,"1EV02973_S11_L001", "1EV02981_S8_L001", "2KK58589_S7_L001"
           ,"2KR049010_S2_L001",    "512022_S1_L001",  "BL37577_S10_L001"    
           ,"BL37590_S4_L001", "C00413_S5_L001",  "C00434_S6_L001"
           ,"CC50362-2_S1_L001",    "CC50362-3_S4_L001",    "IMIN44_S12_L001"     
           ,"Undetermined_S0_L001", "ZL10_S3_L001")

name_species <- data.frame(matrix(ncol = 4, nrow = 0))
names(name_species) <- c("name", "species", "similarity", "distance")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02",92.484,7.516)
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1",NA,NA)
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1",97.286,2.714)
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01",93.933,6.067)
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1",95.198,4.802)
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2",95.616,4.384)
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1",100,0)
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1",95.816,4.184)
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23",89.77,10.23)
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?",NA,NA)
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp",89.77,10.23)
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""), name_species[,4])
names(spec_sam3) <- c("species", "name", "distance")
spec_sam4 <- data.frame(name_species[,2], paste(name_species[,1],"_length", sep=""), name_species[,4])
names(spec_sam4) <- c("species", "name", "distance")
relevant_species <- c("EMCIR1", "SISKIN1", "PARUS1", "GRW01", "WW2", "PHSIB1", "SYAT02", "CXPIP23", "CraneSp")

rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
#rel_coverage <- rel_coverage[-remove_this,]

### CALCULATE PERCENTAGE
#str(coverage_cells )
#coverage_percentage <- coverage[grep("GC|Gene_Name|Gene_Length", names(coverage))]
#i <- rel_coverage$name[6]
#for(i in rel_coverage$name){
#  x <- 100*coverage[colnames(coverage)==i]/coverage$Gene_Length
#  coverage_percentage <- cbind(coverage_percentage, x)
#}
#str(coverage_percentage)

coverage
rel_length <- subset(spec_sam4, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage")
#rel_length <- rel_coverage[-remove_this,] 

coverage_cells <- coverage[grep("_coverage", names(coverage))]

################################
### Total number of covered bases

covered_bases <- genes[grep("_length", names(genes))]
df <- as.data.frame(covered_bases)
sample_list <- list()

i <- rel_length$name[1]
covered_bases_list <- data.frame(matrix(ncol = 2, nrow = 0))
names(covered_bases_list) <- c("sample", "total_bases")
for( i in rel_length$name)
{
  x <- sum(df[,i])
  covered_bases_list[nrow(covered_bases_list) + 1,] = list(i,x)
}

covered_bases_list <- cbind(covered_bases_list, rel_length$distance)
names(covered_bases_list) <- c("sample","total_bases","distance")
plot(covered_bases_list$total_bases~covered_bases_list$distance)
x <- covered_bases_list$distance
y <- covered_bases_list$total_bases
plot(jitter(x,3), jitter(y,5), yaxt="n", pch=".", cex=7, xlab="CytB distance (%) from SISKIN1", ylab="Total coverage (bases)", col=rgb(0, 0, 1, 1))
axis(2, at=c(0,250000,500000, 750000, 1000000,1500000,2000000),labels=TRUE)

######################################
### FOR COVERAGE X, HOW MANY GENES ARE COVERED BY AT LEAST Y SAMPLES? 


coverage_amount <- (1:100)
i<-1
minimum_coverage=90
minimum_samples<-15
min(coverage_cells[1,])
x <- genes[grep("_coverage", names(genes))]
genes_covered <- list() 
for (i in 1:100)
{
  percentage_value <- rowSums(x>=i)
  genes_covered <- cbind(genes_covered, percentage_value)
  colnames(genes_covered)[[ncol(genes_covered)]] <- i
}

sample_nr_vs_genes_covered <- list()
for (i in 1:14)
{
   s <- colSums(genes_covered>=i)
   sample_nr_vs_genes_covered <- cbind(sample_nr_vs_genes_covered, s)
   colnames(sample_nr_vs_genes_covered)[[ncol(sample_nr_vs_genes_covered)]] <- i
}

str(genes_covered)
######################################
### WHICH GENES ARE IN 10 SAMPLES COVERED BY 85% COVERAGE OR BETTER? 

coverage <- genes[grep("_coverage|GC|Gene_Name", names(genes))]
rel_coverage <- subset(spec_sam3, species %in% relevant_species)
remove_this <- which(rel_coverage$name=="C00413_S5_L001_coverage|Undetermined_S0_L001_coverage|2KK58589_S7_L001_coverage")
rel_coverage <- rel_coverage[-remove_this,]

remove_this <- grep("C00413_S5_L001|Undetermined_S0_L001|2KK58589_S7_L001", names(genes))
genes2 <- genes[-remove_this]

genes_85 <- list()

#i <- 1
for(i in 1:nrow(genes)) {
  row <- genes[i,]
  if (rowSums(row[grep("_coverage",names(genes))]>=85)>=10)
  {
    genes_85 <- rbind(genes_85, row)
  }
}
write.csv(genes_85, "genes_85_3.csv")
######################################
### WHICH GENES ARE IN 8 SAMPLES COVERED BY 90% COVERAGE OR BETTER? 
#i <- 1

genes_90 <- list()
for(i in 1:nrow(genes)) {
  row <- genes[i,]
  if (rowSums(row[grep("_coverage",names(genes))]>=90)>=8)
  {
    genes_90 <- rbind(genes_90, row)
  }
}
write.csv(genes_90, "genes_90_2.csv")

######################################
### WHICH GENES ARE IN 8 SAMPLES COVERED BY 99.9% COVERAGE OR BETTER? 
#i <- 1

genes_100 <- list()
for(i in 1:nrow(genes)) {
  row <- genes[i,]
  if (rowSums(row[grep("_coverage",names(genes))]>=99.9)>=8)
  {
    genes_100 <- rbind(genes_100, row)
  }
}
write.csv(genes_100, "genes_100.csv")

######################################
### WHICH GENES ARE IN 10 SAMPLES COVERED BY 96% COVERAGE OR BETTER? 
#i <- 1

genes_96 <- list()
for(i in 1:nrow(genes)) {
  row <- genes[i,]
  if (rowSums(row[grep("_coverage",names(genes))]>=96)>=10)
  {
    genes_96 <- rbind(genes_96, row)
  }
}
write.csv(genes_96, "genes_96.csv")


###############################
### HOW MANY BASES ARE COVERED IN HOST_INTERACTION DEPENDING ON DISTANCE? 
length <- genes[grep("Gene_Name|_length", names(genes))]


hi <- subset(rel_length, `#Gene_Name` %in% host_interaction)
hi <- hi[grep("_length", names(hi))]
hi_total <- colSums(hi)
hi_list <- cbind(spec_sam3, hi_total)
remove_this <- grep("C00413_S5_L001|Undetermined_S0_L001|2KK58589_S7_L001", rownames(hi_list))
cut_hi_list <- hi_list[-remove_this,]
write.csv(hi, "host_interaction.csv")
plot(cut_hi_list$distance, cut_hi_list$hi_total, ylab="total host-interaction bases covered", xlab="phylogenetic distance")
spec_sam3$distance


###############################################################################################
###############################################################################
### WHAT IS THE DISTRIBUTION OF GC CONTENT IN SUCCESSFUL VS FAILED PROBES?
# Note: Take strand into account. 
# Load probe data for Majoris lineage, separate out sequence and average coverage.   
library(readr)
setwd('.')
probes_stats <- read_delim("probes_stats.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
probes_seqs <- read_delim("Ht_probes_c.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
str(probes_stats)
probes_seqs$Coordinates <- substring(probes_seqs$Coordinates, 4) #Remove useless letters
str(probes_seqs)
coverage <- probes_stats[grep("_coverage|GC|probe", names(probes_stats))]
str(coverage)

samples3 <- names(probes_stats[grep("001_coverage", names(probes_stats))])
str(samples3)
species <- c("PARUS1","EMCIR1","SISKIN1","PARUS1","GRW01","PARUS1","CraneSp","WW2","SISKIN1","SISKIN1","PARUS1_PHSIB1","WW2","SISKIN1","PHSIB1","SYAT02", "?", "CXPIP23")

name_species <- data.frame(matrix(ncol = 2, nrow = 0))
names(name_species) <- c("name", "species")
name_species[nrow(name_species) + 1,] = list("102-11C_S3_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("126-11c_S9_L001", "SISKIN1")
name_species[nrow(name_species) + 1,] = list("1ES86798_S2_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("1EV02973_S11_L001", "SYAT02")
name_species[nrow(name_species) + 1,] = list("1EV02981_S8_L001", "WW2")
name_species[nrow(name_species) + 1,] = list("2KK58589_S7_L001", "PARUS1_PHSIB1")
name_species[nrow(name_species) + 1,] = list("2KR049010_S2_L001","EMCIR1")
name_species[nrow(name_species) + 1,] = list("512022_S1_L001", "GRW01")
name_species[nrow(name_species) + 1,] = list("BL37577_S10_L001","PHSIB1")
name_species[nrow(name_species) + 1,] = list("BL37590_S4_L001","WW2")
name_species[nrow(name_species) + 1,] = list("C00413_S5_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("C00434_S6_L001","SISKIN1")
name_species[nrow(name_species) + 1,] = list("CC50362-2_S1_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("CC50362-3_S4_L001","PARUS1")
name_species[nrow(name_species) + 1,] = list("IMIN44_S12_L001","CXPIP23")
name_species[nrow(name_species) + 1,] = list("Undetermined_S0_L001","?")
name_species[nrow(name_species) + 1,] = list("ZL10_S3_L001","CraneSp")
name_species
paste(name_species[,1],"_coverage", sep="")
spec_sam3 <- data.frame(name_species[,2], paste(name_species[,1],"_coverage", sep=""))
names(spec_sam3) <- c("species", "name")

relevant_species <- c("PARUS1", "WW2", "PHSIB1")
#spec_sam3 <- data.frame(species, samples3)
probe_name_and_sequence <- cbind(probes_seqs[,6],probes_seqs[,5], probes_seqs[,3] )
str(probe_name_and_sequence)
gc_vs_probe_coverage <- cbind(probes_stats[,2], probes_stats[,1])

subset(spec_sam3, species %in% c("PARUS1", "WW2", "PHSIB1"))
str(coverage)
col_num <- as.numeric(rownames(subset(spec_sam3, species %in% relevant_species)[2]))
col_num+2
colnames(subset(coverage, select = col_num+2))
col_i <- rowMeans(subset(coverage, select = col_num+2), na.rm = FALSE)
names <- colnames(gc_vs_probe_coverage)
col_i
gc_vs_probe_coverage <- cbind(gc_vs_probe_coverage, col_i)
colnames(gc_vs_probe_coverage)<- c(names, "MAJORIS_coverage")

str(gc_vs_probe_coverage)
coverage_vs_seq <- merge(gc_vs_probe_coverage, probe_name_and_sequence, by.x = "probe", by.y = "Coordinates")
str(coverage_vs_seq)
# Divide probes into 3 tiers: coverage0-50, coverage_51-90, coverage_90-100
coverage_vs_seq_0_50 <- coverage_vs_seq[coverage_vs_seq$MAJORIS_coverage<50, ]
coverage_vs_seq_51_90 <- coverage_vs_seq[(coverage_vs_seq$MAJORIS_coverage>=50 & coverage_vs_seq$MAJORIS_coverage<90), ]
coverage_vs_seq_91_100 <- coverage_vs_seq[coverage_vs_seq$MAJORIS_coverage>=90, ]
# Alternatively, instead of coverage percentage: Best 10%, worst 50% and the rest.
cutoff_50 <- quantile(coverage_vs_seq$MAJORIS_coverage,prob=1-50/100)
cutoff_90 <- quantile(coverage_vs_seq$MAJORIS_coverage,prob=1-10/100)
coverage_vs_seq_b_0_50 <- coverage_vs_seq[coverage_vs_seq$MAJORIS_coverage < cutoff_50,]
coverage_vs_seq_b_51_90 <- coverage_vs_seq[(coverage_vs_seq$MAJORIS_coverage>=cutoff_50 & coverage_vs_seq$MAJORIS_coverage<cutoff_90), ]
coverage_vs_seq_b_91_100 <- coverage_vs_seq[coverage_vs_seq$MAJORIS_coverage > cutoff_90,]

plot(gc_vs_probe_coverage$GC, gc_vs_probe_coverage$MAJORIS, pch=".", cex=2, ylim=c(-5,105), xlim=c(-2,85), xlab="probe GC%", ylab="Percentage of probe covered by 3 reads", col=rgb(0, 0, 1, 0.33))

data <- coverage_vs_seq_0_50$Sequence
i<-1

average_GC_position <- function(data){
  data$Sequence <- gsub("[GC]", "1", data$Sequence)
  data$Sequence <- gsub("[AT]", "0", data$Sequence)
  cov <- data.frame(matrix(ncol = 2, nrow = 0))
  for (i in 1:120 ){
    x <- as.numeric(substr(data$Sequence, i, i))
    cov[nrow(cov) + 1,] = list(i, 100*mean(x))  
  }  
  names(cov) <- c("position","GC")
  return(cov)
}
?aggregate
?reshape

coverage_vs_seq_0_50_GC <- average_GC_position(coverage_vs_seq_0_50)  
coverage_vs_seq_51_90_GC <- average_GC_position(coverage_vs_seq_51_90)  
coverage_vs_seq_91_100_GC <- average_GC_position(coverage_vs_seq_91_100)  

plot(coverage_vs_seq_0_50_GC$position, coverage_vs_seq_0_50_GC$GC,pch=25, bg="green")

points(coverage_vs_seq_51_90_GC$position, coverage_vs_seq_51_90_GC$GC, pch=23,bg="blue")

points(coverage_vs_seq_91_100_GC$position, coverage_vs_seq_91_100_GC$GC, pch=24,bg="red")

# Display as 20-point 6-positions-per-point graphs. 
plot(colSums(matrix(coverage_vs_seq_91_100_GC$GC, nrow=6))/6, ylim=c(28,34),pch=25, bg="green")
points(colSums(matrix(coverage_vs_seq_51_90_GC$GC, nrow=6))/6, pch=23,bg="blue")
points(colSums(matrix(coverage_vs_seq_0_50_GC$GC, nrow=6))/6, pch=24,bg="red")

# Display as 10-point 12-positions-per-point graphs. 
plot(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_91_100_GC$GC, nrow=12))/12, ylim=c(28,35),pch=25, bg="green", ylab="Average GC content over 12 probe bases", xlab="probe position")
points(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_51_90_GC$GC, nrow=12))/12, pch=23,bg="blue")
points(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_0_50_GC$GC, nrow=12))/12, pch=24,bg="red")
legend_strings <- c("coverage>90%", "coverage between 50-90%", "coverage< 50%")
color_vs_legend <- c("green","blue","red")
legend('topright', bty='n', legend = legend_strings, col = "black", cex = 1.0, pch = 21, pt.bg=color_vs_legend)

coverage_vs_seq_b_0_50_GC <- average_GC_position(coverage_vs_seq_b_0_50)  
plot(coverage_vs_seq_b_0_50_GC$position, coverage_vs_seq_b_0_50_GC$GC)

coverage_vs_seq_b_51_90_GC <- average_GC_position(coverage_vs_seq_b_51_90)  
plot(coverage_vs_seq_b_51_90_GC$position, coverage_vs_seq_b_51_90_GC$GC)

coverage_vs_seq_b_91_100_GC <- average_GC_position(coverage_vs_seq_b_91_100)  
plot(coverage_vs_seq_b_91_100_GC$position, coverage_vs_seq_b_91_100_GC$GC)

# Display as 10-point 12-positions-per-point graphs. 
plot(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_b_91_100_GC$GC, nrow=12))/12, ylim=c(28,34),pch=25, bg="green", ylab="Average GC content over 12 probe bases", xlab="probe position")
points(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_b_51_90_GC$GC, nrow=12))/12, pch=24,bg="blue")
points(c(12,24,36,48,60,72,84,96,108,120), colSums(matrix(coverage_vs_seq_b_0_50_GC$GC, nrow=12))/12, pch=23,bg="red")
legend_strings <- c("top 10%", "50-90%", "bottom 50%")
color_vs_legend <- c("green","blue","red")
legend('topright', bty='n', legend = legend_strings, col = "black", cex = 1.0, pch = 21, pt.bg=color_vs_legend)


?gsub <- coverage_vs_seq_0_50$Sequence
coverage_vs_seq_0_50$Sequence <- gsub("[GC]", "1", coverage_vs_seq_0_50$Sequence)
coverage_vs_seq_0_50$Sequence <- gsub("[AT]", "0", coverage_vs_seq_0_50$Sequence)
coverage_vs_seq_0_50_GC <- data.frame(matrix(ncol = 2, nrow = 0))
names(coverage_vs_seq_0_50_GC) <- c("position", "GC")

for (i in 1:120 ){
  x <- substr(coverage_vs_seq_0_50$Sequence, i, i)
  name_species[nrow(name_species) + 1,] = list(i,)  
}

# Possible additional divider: high-GC, mid-GC, low-GC would make for 18 groups. 


