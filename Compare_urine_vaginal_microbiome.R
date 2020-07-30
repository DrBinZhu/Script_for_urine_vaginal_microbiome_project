########### library ###########
library(stringr)
library(GUniFrac)
library(vegan)
library(ggplot2)
library(Rtsne)
library(ALDEx2)
library(tidyr)
library(DiscriMiner)
library(compositions)


############################################# prepare samples ##################################################
# get reads of vaginal samples
setwd('/Users/binzhu/secure/godel/gpfs_fs/bccl/bzhu/Chris')
metadata_v = read.table('data_vagina.txt',header = F, sep = "|")
metadata_u = read.table('data_urine.txt',header = F, sep = "|")

########################### data arrangement ###########################################
metadata_v$V3 = str_replace_all(metadata_v$V3,'AT','')
metadata_u$V3 = str_replace_all(metadata_u$V3,'AT','')
 
metadata_v_2 = as.data.frame(str_c(metadata_v$V2,metadata_v$V3,sep = " "))
metadata_v = cbind(metadata_v$V1,metadata_v_2,metadata_v$V4)

rm(metadata_v_2)

sample_name <- unique(metadata_v$`metadata_v$V1`)
species_name <- unique(metadata_v$`str_c(metadata_v$V2, metadata_v$V3, sep = " ")`)

reads_table_v <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_v) = species_name 
colnames(reads_table_v) = sample_name

for (a in 1: dim(metadata_v)[1]) {
  
  column_num <- which(sample_name == metadata_v[a,1])
  row_num <- which(species_name == metadata_v[a,2])
  
  reads_table_v[row_num,column_num] = as.numeric(as.character(metadata_v[a,3]))
  
}
reads_table_v = as.data.frame(reads_table_v)

# get reads of urine samples

metadata_u_2 = as.data.frame(str_c(metadata_u$V2,metadata_u$V3,sep = " "))
metadata_u = cbind(metadata_u$V1,metadata_u_2,metadata_u$V4)
rm(metadata_u_2)

sample_name <- unique(metadata_u$`metadata_u$V1`)
species_name <- unique(metadata_u$`str_c(metadata_u$V2, metadata_u$V3, sep = " ")`)

reads_table_u <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_u) = species_name 
colnames(reads_table_u) = sample_name

for (a in 1: dim(metadata_u)[1]) {
  
  column_num <- which(sample_name == metadata_u[a,1])
  row_num <- which(species_name == metadata_u[a,2])
  
  reads_table_u[row_num,column_num] = as.numeric(as.character(metadata_u[a,3]))
}

reads_table_u = as.data.frame(reads_table_u)



########################## total reads rarefaction curve ###########################
rare_v <- matrix(data = 0, ncol =2, nrow = 100)
m=0
for (b in seq(1000, 100000, by=1000)) {
  keep <- matrix(ncol = dim(reads_table_v)[2],)
  for (a in 1: dim(reads_table_v)[2]) {
    keep[,a] = sum(reads_table_v[,a]) >= b    
  }
  m = m+1
  rare_v[m,2] = sum(keep)
  rare_v[m,1] = b
}
rare_v <- as.data.frame(rare_v)
colnames(rare_v) <- c('Total reads threshold','Sample number')
rare_v$Site = rep('Vagina', times = nrow(rare_v))


rare_u <- matrix(data = 0, ncol =2, nrow = 100)
m=0
for (b in seq(1000, 100000, by=1000)) {
  keep <- matrix(ncol = dim(reads_table_u)[2],)
  for (a in 1: dim(reads_table_u)[2]) {
    keep[,a] = sum(reads_table_u[,a]) >= b    
  }
  m = m+1
  rare_u[m,2] = sum(keep)
  rare_u[m,1] = b
}
rare_u <- as.data.frame(rare_u)
colnames(rare_u) <- c('Total reads threshold','Sample number')
rare_u$Site = rep('Urine', times = nrow(rare_u))

rare <- rbind(rare_v,rare_u)

ggplot(data=rare, aes(x=`Total reads threshold`, y=`Sample number`, group = Site,color = Site)) +
  geom_line()+
  geom_point()

ggsave('Total_reads_threshold.pdf', width=4, height=3)


#################### sample total reads threshold ####################
keep <- matrix(ncol = dim(reads_table_v)[2],)
for (a in 1: dim(reads_table_v)[2]) {
  keep[,a] = sum(reads_table_v[,a]) >=10000    # input
}
reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

keep <- matrix(ncol = dim(reads_table_u)[2],)
for (a in 1: dim(reads_table_u)[2]) {
  keep[,a] = sum(reads_table_u[,a]) >=10000    # input
}

reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

######### prepare two reads tables #########
species_list = unique(c(row.names(reads_table_v),row.names(reads_table_u)))
species_list <- species_list[order(species_list)]
reads_table_u <- reads_table_u[order(row.names(reads_table_u)), ]
reads_table_v <- reads_table_v[order(row.names(reads_table_v)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v)


keep <- species_list %in% row.names(reads_table_v)
reads_table_v_new[keep,] = reads_table_v

keep <- species_list %in% row.names(reads_table_u)
reads_table_u_new[keep,] = reads_table_u

reads_table <- cbind(reads_table_v_new,reads_table_u_new)

###################### species threshold & present or not ###################################

# get abundance table
reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
  reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

# species threshold
keep <- matrix(, ncol = nrow(reads_table_abundance))

for (a in 1: nrow(reads_table_abundance)) {
  c = sum(reads_table_abundance[a,] >= 0.001) / ncol(reads_table_abundance) >= 0.05        # input
  d = sum(reads_table_abundance[a,] >= 0.0001) / ncol(reads_table_abundance) >= 0.15        # input
  keep[a] = c|d
}

reads_table <- reads_table[keep,]


reads_table_v_new <- as.data.frame(reads_table[,(1:(ncol(reads_table)/2))])
reads_table_u_new <- as.data.frame(reads_table[,(ncol(reads_table)/2+1) : ncol(reads_table)])

reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
  reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new))
sample_list_u = as.matrix(colnames(reads_table_u_new))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))


# find species in two sites # present using 0.01% abundance for each species in each sample
present <- as.data.frame(matrix(data = NA, nrow = nrow(reads_table_u_new), ncol = ncol(reads_table_u_new)))
colnames(present) <- colnames(reads_table_u_new)
row.names(present) <- row.names(reads_table_u_new)

for (a in 1: nrow(present)) {
  for (b in 1: ncol(present)) {
    c = reads_table_abundance[a,b] > 0.0001    # input
    d = reads_table_abundance[a,b+(ncol(reads_table_u_new))] > 0.0001    # input
    
    if (c == T & d == T) {
      present[a,b] = 'Both'
    } else if (c == T & d == F) {
      present[a,b] = 'Vagina'
    } else if (c == F & d == T) {
      present[a,b] = 'Urine'
    } else {
      present[a,b] = 'None'
    }
  }
}

present <- as.data.frame(t(present))
BVAB1 <- as.character(present$`Lachnospiraceae BVAB1`)

present_new <- as.data.frame(matrix(data =0, nrow = 4, ncol = ncol(present)))
colnames(present_new) <- colnames(present)
row.names(present_new) <- c('Both','Vagina','Urine','None')

for (a in 1: nrow(present)) {
  for (b in 1: ncol(present)) {
    if (present[a,b] == 'Both') {
      present_new[1,b] = present_new[1,b]+1
    } else if (present[a,b] == 'Vagina') {
      present_new[2,b] = present_new[2,b]+1
    } else if (present[a,b] == 'Urine') {
      present_new[3,b] = present_new[3,b]+1
    } else {
      present_new[4,b] = present_new[4,b]+1
    }
  }
}
present_new <- as.data.frame(t(present_new))
present_new <- present_new[order(present_new$Both),]
present_new_name <- row.names(present_new)

present_new <- gather(present_new)

present_name <- rep( present_new_name , 4)
present_new$Species <- present_name

present_new$Species<- as.factor(present_new$Species)

present_new$Species <- factor(present_new$Species, levels=present_new_name)

colnames(present_new)[1] = 'Present'

ggplot(data=present_new, aes(x=Species, y=value, fill = `Present`)) +
  geom_bar(stat="identity") +
  labs(x = 'Species', y = "Number of sample pairs")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 6), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16)
  ) + coord_flip() 
ggsave('present_0.01.pdf', width=7, height=7)

############################################# rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)



############################################# alpha diversity ##################################################

alpha.shannon_diversity <- data.frame(diversity(reads_table_rare))
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_rare))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_rare) != 0))

sample_list$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_rare.
sample_list$alpha.evenness <- alpha.evenness$diversity.reads_table_rare.
sample_list$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_rare.....0.

# creat images for alpha diversity
geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

ggplot(sample_list, aes(x=Site, y=alpha.shannon ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Shannon index")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('shannon.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.shannon)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Shannon index")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("shannon_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Evenness")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('evenness.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Evenness")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("evenness_2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Observed Species")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('otu.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Observed species")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("otu_2.pdf", width=5.5, height=5.5)

# Compute paired Wilcoxon test for the significance of alpha diversity
Urine <- subset(sample_list,  Site == "Urine", alpha.shannon,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.shannon,
                 drop = TRUE)
pvalue.alpha.shannon <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.shannon

Urine <- subset(sample_list,  Site == "Urine", alpha.evenness,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.evenness,
                 drop = TRUE)
pvalue.alpha.evenness <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.evenness 

Urine <- subset(sample_list,  Site == "Urine", alpha.ovserved_OTU,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.ovserved_OTU,
                 drop = TRUE)
pvalue.alpha.ovserved_OTU <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.ovserved_OTU


pvalue = rbind(pvalue.alpha.shannon,pvalue.alpha.evenness,pvalue.alpha.ovserved_OTU)

########################################### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Bray_Curtis distance between paired samples
x = nrow(sample_list)/2
pair_dis = as.data.frame(matrix(data = NA, nrow = x, ncol =3))
pair_dis$V1 = sample_list[(x+1):(x*2),1]
pair_dis$V2 = sample_list[(x+1):(x*2),2]

for (a in 1: x) {
  bray_x = which(row.names(Bray_Curtis) == as.character(pair_dis[a,1]))
  bray_y = which(colnames(Bray_Curtis) == as.character(pair_dis[a,2]))
  pair_dis[a,3] = Bray_Curtis[bray_x,bray_y]
}
pair_dis <- pair_dis[order(pair_dis$V3),]
colnames(pair_dis) <- c('Sample_name_u', 'Sample_name_v','Bray_Curtis')

pair_dis$Sample_name_u <- factor(pair_dis$Sample_name_u, levels = pair_dis$Sample_name_u)  # convert to factor to retain sorted order in plot.


ggplot(data=pair_dis, aes(x=Sample_name_u, y=Bray_Curtis)) +
  geom_bar(stat="identity") +
  labs(x = 'Paired samples', y = "Bray-Curtis distance")+ 
  theme(axis.title = element_text(size = 12), 
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))   
ggsave("paired_distance.pdf", width=3.5, height=5.5)


# distance comparison 
vagina_bray_distantc = Bray_Curtis[(1:(nrow(Bray_Curtis)/2)),(1:(ncol(Bray_Curtis)/2))]
vagina_bray_distantc <- gather(vagina_bray_distantc)

distance_comparison <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(vagina_bray_distantc)), ncol =2)
distance_comparison[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison[1:nrow(pair_dis),2] = rep('Vagina_Paired urine',nrow(pair_dis))
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),1] = vagina_bray_distantc$value
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),2] = rep('Vagina_Other vagina',nrow(vagina_bray_distantc))
distance_comparison<- as.data.frame(distance_comparison)
colnames(distance_comparison) = c('Bray-Curtis Distance', 'Type')
distance_comparison$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison$`Bray-Curtis Distance`))
distance_comparison$Type <- as.character(distance_comparison$Type)

ggplot(distance_comparison, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
  geom_boxplot() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('Bray_curtis_distance_compare.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,pvalue.wil1 )

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
  metaMDS(Bray_Curtis,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 50,
          wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- sample_list$Site
mds_data$Pair <- sample_list$Pair
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
  geom_point()
ggsave("NMDS1.pdf", width=5.5, height=5.5)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
  geom_point()+
  geom_line(aes(group = Pair), size = 0.3, alpha = 0.7, color = 'black')
ggsave("NMDS2.pdf", width=5.5, height=5.5)

# bar plot for beta diversity
vagina2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),1:(length(Bray_Curtis)/2)])
vagina2vagina = as.vector(vagina2vagina)
vagina2vagina = vagina2vagina[vagina2vagina!=0]
vagina = as.data.frame(cbind(vagina2vagina,'Vagina'))
colnames(vagina) = c('Bray_Curtis_distance','Body_site')

Urine2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),(length(Bray_Curtis)/2+1):length(Bray_Curtis)])
Urine2vagina = as.vector(Urine2vagina)
urine = as.data.frame(cbind(Urine2vagina,'Urine'))
colnames(urine) = c('Bray_Curtis_distance','Body_site')

beta_plot = rbind(vagina,urine)
beta_plot$Body_site <- as.factor(beta_plot$Body_site)
beta_plot$Bray_Curtis_distance <- as.numeric(as.character(beta_plot$Bray_Curtis_distance))

ggplot(beta_plot, aes(x=Body_site, y=Bray_Curtis_distance,fill=Body_site)) +
  geom_boxplot()  +
  labs(x = NULL, y = "Distance to Vagina samples")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))

ggsave("Bray_Curtis_bar.pdf", width=4, height=5.5)

# Compute perMANOVA test for the significance of beta diversity
# adonis2 : permutational multivariate analysis of variance using distance matrices the centroid and/or the spread of the objects is different between the groups.
pvalue.adonis.site <- adonis2(reads_table_rare ~ Site, data = sample_list, method = "bray")  
pvalue.adonis.site <- pvalue.adonis.site[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.site)
pvalue.adonis.pair <- adonis2(reads_table_rare ~ Pair, data = sample_list, method = "bray")  
pvalue.adonis.pair <- pvalue.adonis.pair[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.pair)

write.csv(pvalue,'pvalue.csv')





################# species with differential abundance #############################
reads_table_das <- as.data.frame((reads_table))

conds <- sample_list$Site

x <- aldex.clr(reads_table_das, conds, mc.samples=128, denom="all", verbose=F)

# paired Wilcoxon Rank Sum test and Welch's t-test
x.tt <- aldex.ttest(x, paired.test=TRUE)

x.effect <- aldex.effect(x)

x.all <- data.frame(cbind(x.tt,x.effect))

abundance_change_das <- x.all$diff.btw # diff.btw is a vector containing the per-feature median difference between condition A and B


x.all <- cbind(abundance_change_das,x.all)
# get abundances of species
abundance_das = rowSums(reads_table_das)
abundance_das = abundance_das / sum(abundance_das)

`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the adj-pvalue
x.all$abundance_das <- abundance_das
x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`

# draw figure
das <- x.all[(x.all$`-Log10(adj-pvalue)` >=1.301 & (x.all$abundance_change_das >=1 | x.all$abundance_change_das <=-1)),]

das$Species <- row.names(das)
das <- das[order(das$abundance_change_das),]    
das$Color <- ifelse(das$abundance_change_das < 0, "Enriched in Urine", "Enriched in Vagina")  # above / below avg flag
das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.

das$abundance_change_das[das$abundance_change_das == Inf] = 5
das$abundance_change_das[das$abundance_change_das == -Inf] = -5
das$abundance_change_das[das$abundance_change_das <= -5] = -5
das$abundance_change_das[das$abundance_change_das >= 5] = 5

theme_set(theme_bw())  

ggplot(das, aes(Species, abundance_change_das)) + 
  geom_jitter(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
  coord_flip() +          # convert x y axis
  labs(x = 'Species', y = "Median difference in clr values")+ 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))    
ggsave("das.pdf", width=8, height=3.5)




########################## draw bar figure ########################

# find top n species
reads_table_v_new_bar <- reads_table_v_new

reads_table_v_new_bar$total_reads <- rowSums(reads_table_v_new_bar)
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$total_reads,decreasing = T),]
reads_table_v_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_v_new_bar[10:nrow(reads_table_v_new_bar),])))
reads_table_v_new_bar <- reads_table_v_new_bar[-c(10:nrow(reads_table_v_new_bar)), ]
reads_table_v_new_bar <- rbind(reads_table_v_new_bar,other_reads)
row.names(reads_table_v_new_bar)[10] = 'Others'

reads_table_v_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new_bar),nrow = nrow(reads_table_v_new_bar))

for (a in 1:ncol(reads_table_v_new)) {
  reads_table_v_new_bar_abundance[,a] <- reads_table_v_new_bar[,a] / colSums(reads_table_v_new_bar)[a]
  
}
row.names(reads_table_v_new_bar_abundance) = row.names(reads_table_v_new_bar)
colnames(reads_table_v_new_bar_abundance) = colnames(reads_table_v_new_bar) 

reads_table_v_new_bar <- as.data.frame(reads_table_v_new_bar_abundance)
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar <- reads_table_v_new_bar/(colSums(reads_table_v_new_bar))

plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
plot_v$key <- factor(plot_v$key, levels = colnames(reads_table_v_new_bar))

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vaginal samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("bar_v.pdf", width=9, height=3.5)


# urine
# find top n species
reads_table_u_new_bar <- reads_table_u_new

reads_table_u_new_bar$total_reads <- rowSums(reads_table_u_new_bar)
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$total_reads,decreasing = T),]
reads_table_u_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_u_new_bar[10:nrow(reads_table_u_new_bar),])))
reads_table_u_new_bar <- reads_table_u_new_bar[-c(10:nrow(reads_table_u_new_bar)), ]
reads_table_u_new_bar <- rbind(reads_table_u_new_bar,other_reads)
row.names(reads_table_u_new_bar)[10] = 'Others'

reads_table_u_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new_bar),nrow = nrow(reads_table_u_new_bar))

for (a in 1:ncol(reads_table_u_new)) {
  reads_table_u_new_bar_abundance[,a] <- reads_table_u_new_bar[,a] / colSums(reads_table_u_new_bar)[a]
  
}
row.names(reads_table_u_new_bar_abundance) = row.names(reads_table_u_new_bar)
colnames(reads_table_u_new_bar_abundance) = colnames(reads_table_u_new_bar) 

reads_table_u_new_bar <- as.data.frame(reads_table_u_new_bar_abundance)
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar <- reads_table_u_new_bar/(colSums(reads_table_u_new_bar))

plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
plot_u$key <- factor(plot_u$key, levels = colnames(reads_table_u_new_bar))

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#FDFFBA','#55ff42',
                             '#D7A7F1','#3264B8','#F7A0A0','#3F3F3F')) +
  labs(x = 'Urine samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_u.pdf", width=9, height=3.5)


# urine with marched order
plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
order <- colnames(reads_table_v_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
  n = which(sample_list_2[,2] == order[a])
  order_urine[a] = sample_list_2[n,1]
}
plot_u$key <- factor(plot_u$key, levels = order_urine)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#FDFFBA','#55ff42',
                             '#D7A7F1','#3264B8','#F7A0A0','#3F3F3F')) +
  labs(x = 'Urine samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_u_vagina_order.pdf", width=9, height=3.5)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#FDFFBA','#55ff42',
                             '#D7A7F1','#3264B8','#F7A0A0','#3F3F3F')) +
  labs(x = 'Urine samples', y = "")+ 
  coord_polar() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

ggsave("bar_u_circle_vagina_order.pdf", width=10, height=10)

# vagina with marched order
plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
order <- colnames(reads_table_u_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
  n = which(sample_list_2[,1] == order[a])
  order_urine[a] = sample_list_2[n,2]
}
plot_v$key <- factor(plot_v$key, levels = order_urine)

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vagina samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_v_urine_order.pdf", width=9, height=3.5)

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vagina samples', y = "")+ 
  coord_polar() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

ggsave("bar_v_circle_urine_order.pdf", width=10, height=10)

################## create abundance table ##################
reads_table_u_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new),nrow = nrow(reads_table_u_new))

for (a in 1:ncol(reads_table_u_new)) {
  reads_table_u_abundance[,a] <- reads_table_u_new[,a] / colSums(reads_table_u_new)[a]
}
reads_table_u_abundance <- as.data.frame(reads_table_u_abundance)
row.names(reads_table_u_abundance) = row.names(reads_table_u_new)
colnames(reads_table_u_abundance) = colnames(reads_table_v_new)

reads_table_v_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new),nrow = nrow(reads_table_v_new))

for (a in 1:ncol(reads_table_v_new)) {
  reads_table_v_abundance[,a] <- reads_table_v_new[,a] / colSums(reads_table_v_new)[a]
}
reads_table_v_abundance <- as.data.frame(reads_table_v_abundance)
row.names(reads_table_v_abundance) = row.names(reads_table_v_new)
colnames(reads_table_v_abundance) = colnames(reads_table_v_new)

total_abundance_v = as.data.frame(rowSums(reads_table_v_abundance))
total_abundance_v <- total_abundance_v/sum(total_abundance_v$`rowSums(reads_table_v_abundance)`)


total_abundance_u = as.data.frame(rowSums(reads_table_u_abundance))
total_abundance_u <- total_abundance_u/sum(total_abundance_u$`rowSums(reads_table_u_abundance)`)



################## linear regression  ################

line_list <- matrix(data = NA, nrow = nrow(reads_table_v_abundance), ncol =3)
line_list <- as.data.frame(line_list)
colnames(line_list) = c('Species','p-value','R-value')
row.names(line_list)= row.names(reads_table_v_abundance)
line_list$Species = row.names(reads_table_v_abundance)

for (a in 1:nrow(line_list)) {
  # no treatment
  line <- as.data.frame(t(rbind(reads_table_v_abundance[a,],reads_table_u_abundance[a,])))
  
  line <- log10(line)

  colnames(line) = c('Vagina','Urine')
  
  line2=line
  line2[line2 == '-Inf'] = -6
  linearMod <- lm(Vagina ~ Urine, data=line2) 
  r_squared <- summary(linearMod)$r.squared
  pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
  
  line_list[a,2]= pvalue
  line_list[a,3]= r_squared
  
  if (line_list[a,1] %in% c('Sneathia amnii ','Lachnospiraceae BVAB1 ','TM7 OTU-H1 ','Prevotella cluster2 ')) {
    ggplot(line, aes(Vagina, Urine) )  +
      geom_point(size=1) +
      xlab(paste0("Vagina abundance")) +
      ylab(paste0("Urine abundance")) + ggtitle(paste0(line_list[a,1], '\n',"R = ", r_squared, '\n',"pvalue = ", pvalue))  + 
      geom_smooth(method='lm', se=T) + 
      theme(axis.title.x = element_text( size=16),
            axis.title.y = element_text( size=16)) 
    ggsave(paste0('linear_1_',line_list[a,1],'.pdf'), width=5.5, height=5.5)
  }

}

  
line_list$adj_p = p.adjust(line_list$`p-value`)  #  Benjamini & Hochberg 

write.csv(line_list,'linear_regression.csv')




########## threshold ###########
reads_table_u_new2 <- reads_table_u_new
reads_table_v_new2 <- reads_table_v_new

#if ratio of total abundance is larger than 10 folder, the species is removed 
abundance_ratio <- matrix(data = NA, nrow = nrow(reads_table_u_new), ncol =3)

abundance_ratio[,1] <- row.names(reads_table_u_new)
abundance_ratio[,2] <- rowSums(reads_table_u_abundance)/ rowSums(reads_table_v_abundance)
abundance_ratio[,3] <- rowSums(reads_table_u_abundance) +rowSums(reads_table_v_abundance)

abundance_ratio <- as.data.frame(abundance_ratio)
abundance_ratio$V2 <- as.numeric(as.character(abundance_ratio$V2))
abundance_ratio$V3 <- as.numeric(as.character(abundance_ratio$V3))

keep = abundance_ratio$V2 < 10 

# Remove species present (0.01% threshold) urine / vagina >= 5 
present_urine <- present == 'Urine' | present == 'Both'
present_urine <- colSums(present_urine)
present_vagina <- present == 'Vagina' | present == 'Both'
present_vagina <- colSums(present_vagina)
present_ratio <- present_urine/present_vagina

keep2 = present_ratio < 5 

keep3 = keep & keep2
reads_table_v_new2 <- reads_table_v_new2[keep3,]
reads_table_u_new2 <- reads_table_u_new2[keep3,]



########## new reads and abundance table ############
# sample total reads threshold
keep <- matrix(ncol = dim(reads_table_v_new2)[2],)
for (a in 1: dim(reads_table_v_new2)[2]) {
  keep[,a] = sum(reads_table_v_new2[,a]) >=10000    # input
}
reads_table_u_new2 <- reads_table_u_new2[,keep]
reads_table_v_new2 <- reads_table_v_new2[,keep]

keep <- matrix(ncol = dim(reads_table_u_new2)[2],)
for (a in 1: dim(reads_table_u_new2)[2]) {
  keep[,a] = sum(reads_table_u_new2[,a]) >=10000    # input
}

reads_table_u_new2 <- reads_table_u_new2[,keep]
reads_table_v_new2 <- reads_table_v_new2[,keep]

reads_table = cbind(reads_table_v_new2, reads_table_u_new2)

# get abundance table
reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))

for (a in 1:ncol(reads_table)) {
  reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
}
reads_table_abundance <- as.data.frame(reads_table_abundance)
row.names(reads_table_abundance) = row.names(reads_table)
colnames(reads_table_abundance) = colnames(reads_table)

reads_table_v_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new2),nrow = nrow(reads_table_v_new2))

for (a in 1:ncol(reads_table_v_new2)) {
  reads_table_v_abundance[,a] <- reads_table_v_new2[,a] / colSums(reads_table_v_new2)[a]
}
reads_table_v_abundance <- as.data.frame(reads_table_v_abundance)
row.names(reads_table_v_abundance) = row.names(reads_table_v_new2)
colnames(reads_table_v_abundance) = colnames(reads_table_v_new2)

reads_table_u_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new2),nrow = nrow(reads_table_u_new2))

for (a in 1:ncol(reads_table_u_new2)) {
  reads_table_u_abundance[,a] <- reads_table_u_new2[,a] / colSums(reads_table_u_new2)[a]
}
reads_table_u_abundance <- as.data.frame(reads_table_u_abundance)
row.names(reads_table_u_abundance) = row.names(reads_table_u_new2)
colnames(reads_table_u_abundance) = colnames(reads_table_v_new2)

# get samples metadata
sample_list_v = as.matrix(colnames(reads_table_v_new2))
sample_list_u = as.matrix(colnames(reads_table_u_new2))

sample_list = as.data.frame(matrix(data = NA, ncol=3, nrow= 2*nrow(sample_list_u)))
colnames(sample_list) = c('Sample','Pair','Site')
sample_list$Sample = rbind(sample_list_v,sample_list_u)
sample_list$Pair = rbind(sample_list_v,sample_list_v)
sample_list$Site = c(matrix(data = 'Vagina' , nrow = nrow(sample_list_u)),matrix(data = 'Urine' , nrow = nrow(sample_list_u)))

reads_table_u_new <- reads_table_u_new2
reads_table_v_new <- reads_table_v_new2
############################################# rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)



############################################# alpha diversity ##################################################

alpha.shannon_diversity <- data.frame(diversity(reads_table_rare))
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_rare))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_rare) != 0))

sample_list$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_rare.
sample_list$alpha.evenness <- alpha.evenness$diversity.reads_table_rare.
sample_list$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_rare.....0.
#colnames(sample_list)=c('Sample','Pair','Site','Shannon','Evenness','Observed OTU')

# creat images for alpha diversity

geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

ggplot(sample_list, aes(x=Site, y=alpha.shannon ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Shannon index")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('shannon__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.shannon)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Shannon index")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("shannon_2__2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Evenness")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('evenness__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.evenness)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Evenness")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("evenness_2__2.pdf", width=5.5, height=5.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU ,color = Site)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x = NULL, y = "Observed Species")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave('otu__2.pdf', width=3.5, height=3.5)

ggplot(sample_list, aes(x=Site, y=alpha.ovserved_OTU)) +
  geom_line(aes(group = Pair), size = 0.1, alpha = 0.7) +
  geom_point(aes(color = Site), alpha = 0.5) +
  labs(x = NULL, y = "Observed species")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))
ggsave("otu_2__2.pdf", width=5.5, height=5.5)

# Compute paired Wilcoxon test for the significance of alpha diversity
Urine <- subset(sample_list,  Site == "Urine", alpha.shannon,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.shannon,
                 drop = TRUE)
pvalue.alpha.shannon <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.shannon

Urine <- subset(sample_list,  Site == "Urine", alpha.evenness,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.evenness,
                 drop = TRUE)
pvalue.alpha.evenness <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.evenness 

Urine <- subset(sample_list,  Site == "Urine", alpha.ovserved_OTU,
                drop = TRUE)
Vagina <- subset(sample_list,  Site == "Vagina", alpha.ovserved_OTU,
                 drop = TRUE)
pvalue.alpha.ovserved_OTU <- wilcox.test(Urine, Vagina, paired = TRUE)$p.value
pvalue.alpha.ovserved_OTU


pvalue = rbind(pvalue.alpha.shannon,pvalue.alpha.evenness,pvalue.alpha.ovserved_OTU)

########################################### beta diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Bray_Curtis distance between paired samples
x = nrow(sample_list)/2
pair_dis = as.data.frame(matrix(data = NA, nrow = x, ncol =3))
pair_dis$V1 = sample_list[(x+1):(x*2),1]
pair_dis$V2 = sample_list[(x+1):(x*2),2]

for (a in 1: x) {
  bray_x = which(row.names(Bray_Curtis) == as.character(pair_dis[a,1]))
  bray_y = which(colnames(Bray_Curtis) == as.character(pair_dis[a,2]))
  pair_dis[a,3] = Bray_Curtis[bray_x,bray_y]
}
pair_dis <- pair_dis[order(pair_dis$V3),]
colnames(pair_dis) <- c('Sample_name_u', 'Sample_name_v','Bray_Curtis')

pair_dis$Sample_name_u <- factor(pair_dis$Sample_name_u, levels = pair_dis$Sample_name_u)  # convert to factor to retain sorted order in plot.


ggplot(data=pair_dis, aes(x=Sample_name_u, y=Bray_Curtis)) +
  geom_bar(stat="identity") +
  labs(x = 'Paired samples', y = "Bray-Curtis distance between paired samples")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))   

ggsave("paired_distance_2.pdf", width=3.5, height=5.5)


# distance comparison 
vagina_bray_distantc = Bray_Curtis[(1:(nrow(Bray_Curtis)/2)),(1:(ncol(Bray_Curtis)/2))]
vagina_bray_distantc <- gather(vagina_bray_distantc)

distance_comparison <- matrix(data = NA, nrow =(nrow(pair_dis)+nrow(vagina_bray_distantc)), ncol =2)
distance_comparison[1:nrow(pair_dis),1] = pair_dis$Bray_Curtis
distance_comparison[1:nrow(pair_dis),2] = rep('Vagina_Paired urine',nrow(pair_dis))
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),1] = vagina_bray_distantc$value
distance_comparison[(nrow(pair_dis)+1):nrow(distance_comparison),2] = rep('Vagina_Other vagina',nrow(vagina_bray_distantc))
distance_comparison<- as.data.frame(distance_comparison)
colnames(distance_comparison) = c('Bray-Curtis Distance', 'Type')
distance_comparison$`Bray-Curtis Distance` = as.numeric(as.character(distance_comparison$`Bray-Curtis Distance`))
distance_comparison$Type <- as.character(distance_comparison$Type)

ggplot(distance_comparison, aes(x=Type, y=`Bray-Curtis Distance`,fill = Type)) +
  geom_boxplot() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('Bray_curtis_distance_compare_2.pdf', width=5, height=5)

pvalue.wil1 = wilcox.test(`Bray-Curtis Distance` ~ Type, data=distance_comparison)$p.value 
pvalue.Bray_Curtis_comparision = pvalue.wil1
pvalue = rbind(pvalue,pvalue.wil1 )

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
  metaMDS(Bray_Curtis,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 50,
          wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- sample_list$Site
mds_data$Pair <- sample_list$Pair
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
  geom_point()+
  stat_ellipse(type = "t")
ggsave("NMDS1_2.pdf", width=5.5, height=5.5)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
  geom_point()+
  geom_line(aes(group = Pair), size = 0.3, alpha = 0.7, color = 'black')
ggsave("NMDS2_2.pdf", width=5.5, height=5.5)

# bar plot for beta diversity
vagina2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),1:(length(Bray_Curtis)/2)])
vagina2vagina = as.vector(vagina2vagina)
vagina2vagina = vagina2vagina[vagina2vagina!=0]
vagina = as.data.frame(cbind(vagina2vagina,'Vagina'))
colnames(vagina) = c('Bray_Curtis_distance','Body_site')

Urine2vagina = as.matrix(Bray_Curtis[1:(length(Bray_Curtis)/2),(length(Bray_Curtis)/2+1):length(Bray_Curtis)])
Urine2vagina = as.vector(Urine2vagina)
urine = as.data.frame(cbind(Urine2vagina,'Urine'))
colnames(urine) = c('Bray_Curtis_distance','Body_site')

beta_plot = rbind(vagina,urine)
beta_plot$Body_site <- as.factor(beta_plot$Body_site)
beta_plot$Bray_Curtis_distance <- as.numeric(as.character(beta_plot$Bray_Curtis_distance))

ggplot(beta_plot, aes(x=Body_site, y=Bray_Curtis_distance,fill=Body_site)) +
  geom_boxplot()  +
  labs(x = NULL, y = "Distance to Vagina samples")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))

ggsave("Bray_Curtis_bar_2.pdf", width=4, height=5.5)

# vagitype
vagitype_v = matrix(data = NA, ncol =1, nrow = ncol(reads_table_v_new))
colnames(vagitype_v) = 'Vagitype'

for (a in 1:ncol(reads_table_v_new)) {
  n= which(reads_table_v_new[,a] == max(reads_table_v_new[,a]), arr.ind=TRUE)
  vagitype_v[a,1]= row.names(reads_table_v_new)[n]
}

vagitype_u= matrix(data = NA, ncol =1, nrow = ncol(reads_table_u_new))
colnames(vagitype_u) = 'Vagitype'

for (a in 1:ncol(reads_table_u_new)) {
  n= which(reads_table_u_new[,a] == max(reads_table_u_new[,a]), arr.ind=TRUE)
  vagitype_u[a,1]= row.names(reads_table_u_new)[n]
}

vagitype = rbind(vagitype_v,vagitype_u)
vagitype_unique = unique(vagitype)
for (a in 1:length(vagitype_unique)) {
  if (sum(vagitype == vagitype_unique[a]) < 5) {
    vagitype[vagitype == vagitype_unique[a]] = 'Others'
  }
}

mds_data = cbind(mds_data,vagitype)
ggplot(mds_data, aes(MDS1, MDS2))  +
  geom_line(aes(group = Pair), size = 0.3, alpha = 0.7) +
  geom_point(size=2,aes(shape = Site, color = Vagitype)) +
  scale_color_manual(values=c('#A50000','#FF4040','#FFD054','#FFFB00','#FC91CB','#71D4FC','#D0D0D0')) +
  coord_fixed()+ 
  theme(
    axis.title.x = element_text( size=16),
    axis.title.y = element_text( size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) 
ggsave("vagitype2.pdf", width=8, height=5.5)

same_vagitype = matrix(data = NA, nrow= nrow(mds_data)/2, ncol =2)
same_vagitype[,1] = mds_data$Vagitype[1:(nrow(mds_data)/2)]
same_vagitype[,2] = mds_data$Vagitype[(nrow(mds_data)/2+1):nrow(mds_data)]
number_vagitype_2 = sum(same_vagitype[,1] == same_vagitype[,2])
pvalue <- rbind(pvalue,number_vagitype_2)



# Compute perMANOVA test for the significance of beta diversity
# adonis2 : permutational multivariate analysis of variance using distance matrices the centroid and/or the spread of the objects is different between the groups.
pvalue.adonis.site <- adonis2(reads_table_rare ~ Site, data = sample_list, method = "bray")  
pvalue.adonis.site <- pvalue.adonis.site[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.site)
pvalue.adonis.pair <- adonis2(reads_table_rare ~ Pair, data = sample_list, method = "bray")  
pvalue.adonis.pair <- pvalue.adonis.pair[1,5]
pvalue <- rbind(pvalue,pvalue.adonis.pair)

write.csv(pvalue,'pvalue_2.csv')
















################################# species with differential abundance #############################
reads_table_das <- as.data.frame((reads_table))

conds <- sample_list$Site

x <- aldex.clr(reads_table_das, conds, mc.samples=128, denom="all", verbose=F)

# paired Wilcoxon Rank Sum test and Welch's t-test
x.tt <- aldex.ttest(x, paired.test=TRUE)

x.effect <- aldex.effect(x)

x.all <- data.frame(cbind(x.tt,x.effect))

abundance_change_das <- x.all$diff.btw 

x.all <- cbind(abundance_change_das,x.all)
# get abundances of species
abundance_das = rowSums(reads_table_das)
abundance_das = abundance_das / sum(abundance_das)

`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the pvalue
x.all$abundance_das <- abundance_das
x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`


# draw figure
das <- x.all[(x.all$`-Log10(adj-pvalue)` >=1.301 & (x.all$abundance_change_das >=1 | x.all$abundance_change_das <=-1)),]

das$Species <- row.names(das)
das <- das[order(das$abundance_change_das),]    #13. ♦ dif.btw - median difference in clr values between S and NS groups
das$Color <- ifelse(das$abundance_change_das < 0, "Enriched in Urine", "Enriched in Vagina")  # above / below avg flag
das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.

das$abundance_change_das[das$abundance_change_das == Inf] = 5
das$abundance_change_das[das$abundance_change_das == -Inf] = -5
das$abundance_change_das[das$abundance_change_das <= -5] = -5
das$abundance_change_das[das$abundance_change_das >= 5] = 5

theme_set(theme_bw())  # pre-set the bw theme.

ggplot(das, aes(Species, abundance_change_das)) + 
  geom_jitter(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
  coord_flip() +          # convert x y axis
  labs(x = 'Species', y = "Median difference in clr values")+ 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))    
ggsave("das_2.pdf", width=6, height=3.2)



########################## draw bar figure ########################

# find top n species
reads_table_v_new_bar <- reads_table_v_new

reads_table_v_new_bar$total_reads <- rowSums(reads_table_v_new_bar)
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$total_reads,decreasing = T),]
reads_table_v_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_v_new_bar[10:nrow(reads_table_v_new_bar),])))
reads_table_v_new_bar <- reads_table_v_new_bar[-c(10:nrow(reads_table_v_new_bar)), ]
reads_table_v_new_bar <- rbind(reads_table_v_new_bar,other_reads)
row.names(reads_table_v_new_bar)[10] = 'Others'

reads_table_v_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_v_new_bar),nrow = nrow(reads_table_v_new_bar))

for (a in 1:ncol(reads_table_v_new)) {
  reads_table_v_new_bar_abundance[,a] <- reads_table_v_new_bar[,a] / colSums(reads_table_v_new_bar)[a]
  
}
row.names(reads_table_v_new_bar_abundance) = row.names(reads_table_v_new_bar)
colnames(reads_table_v_new_bar_abundance) = colnames(reads_table_v_new_bar) 

reads_table_v_new_bar <- as.data.frame(reads_table_v_new_bar_abundance)
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar = reads_table_v_new_bar[order(reads_table_v_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_v_new_bar <- as.data.frame(t(reads_table_v_new_bar))
reads_table_v_new_bar <- reads_table_v_new_bar/(colSums(reads_table_v_new_bar))

plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
plot_v$key <- factor(plot_v$key, levels = colnames(reads_table_v_new_bar))

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vaginal samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("bar_v_2.pdf", width=9, height=3.5)


# urine
# find top n species
reads_table_u_new_bar <- reads_table_u_new

reads_table_u_new_bar$total_reads <- rowSums(reads_table_u_new_bar)
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$total_reads,decreasing = T),]
reads_table_u_new_bar$total_reads <- NULL
other_reads = as.data.frame(t(colSums(reads_table_u_new_bar[10:nrow(reads_table_u_new_bar),])))
reads_table_u_new_bar <- reads_table_u_new_bar[-c(10:nrow(reads_table_u_new_bar)), ]
reads_table_u_new_bar <- rbind(reads_table_u_new_bar,other_reads)
row.names(reads_table_u_new_bar)[10] = 'Others'

reads_table_u_new_bar_abundance <- matrix(data =0, ncol = ncol(reads_table_u_new_bar),nrow = nrow(reads_table_u_new_bar))

for (a in 1:ncol(reads_table_u_new)) {
  reads_table_u_new_bar_abundance[,a] <- reads_table_u_new_bar[,a] / colSums(reads_table_u_new_bar)[a]
  
}
row.names(reads_table_u_new_bar_abundance) = row.names(reads_table_u_new_bar)
colnames(reads_table_u_new_bar_abundance) = colnames(reads_table_u_new_bar) 

reads_table_u_new_bar <- as.data.frame(reads_table_u_new_bar_abundance)
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar = reads_table_u_new_bar[order(reads_table_u_new_bar$`Lactobacillus iners`,decreasing = T),]
reads_table_u_new_bar <- as.data.frame(t(reads_table_u_new_bar))
reads_table_u_new_bar <- reads_table_u_new_bar/(colSums(reads_table_u_new_bar))

plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
plot_u$key <- factor(plot_u$key, levels = colnames(reads_table_u_new_bar))

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#FDFFBA','#55ff42',
                             '#D7A7F1','#3264B8','#F7A0A0','#3F3F3F')) +
  labs(x = 'Urine samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_u_2.pdf", width=9, height=3.5)


# urine with marched order
plot_u <- gather(reads_table_u_new_bar)
plot_u_name <- row.names(reads_table_u_new_bar)
plot_u_name <- rep( plot_u_name , ncol(reads_table_u_new))
plot_u$Species <- plot_u_name

plot_u$Species <- as.factor(plot_u$Species)
plot_u$Species <- factor(plot_u$Species, levels = plot_u$Species[1:10])

plot_u$key <- as.factor(plot_u$key)
order <- colnames(reads_table_v_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
  n = which(sample_list_2[,2] == order[a])
  order_urine[a] = sample_list_2[n,1]
}
plot_u$key <- factor(plot_u$key, levels = order_urine)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#FDFFBA','#55ff42','#D7A7F1',
                             '#3264B8','#F7A0A0','#bfbfbf','#3F3F3F')) +
  labs(x = 'Urine samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_u_vagina_order_2.pdf", width=9, height=3.5)

ggplot(data=plot_u, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#3be5f7','#FC9800','#C80B0B','#FDFFBA','#55ff42',
                             '#D7A7F1','#3264B8','#F7A0A0','#3F3F3F')) +
  labs(x = 'Urine samples', y = "")+ 
  coord_polar() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

ggsave("bar_u_circle_vagina_order_2.pdf", width=10, height=10)

# vagina with marched order
plot_v <- gather(reads_table_v_new_bar)
plot_v_name <- row.names(reads_table_v_new_bar)
plot_v_name <- rep( plot_v_name , ncol(reads_table_v_new))
plot_v$Species <- plot_v_name

plot_v$Species <- as.factor(plot_v$Species)
plot_v$Species <- factor(plot_v$Species, levels = plot_v$Species[1:10])

plot_v$key <- as.factor(plot_v$key)
order <- colnames(reads_table_u_new_bar)
order_urine <- vector(mode = "character", length = length(order))
sample_list_2 <- sample_list[(length(order)+1):nrow(sample_list),]
for (a in 1: length(order)) {
  n = which(sample_list_2[,1] == order[a])
  order_urine[a] = sample_list_2[n,2]
}
plot_v$key <- factor(plot_v$key, levels = order_urine)

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vagina samples', y = "Abundance")+ 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("bar_v_urine_order_2.pdf", width=9, height=3.5)

ggplot(data=plot_v, aes(x=key, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#A8C5D3','#FC9800','#C80B0B','#D7A7F1','#55ff42','#F7A0A0',
                             '#FDFFBA','#3264B8','#bfbfbf','#3F3F3F')) +
  labs(x = 'Vagina samples', y = "")+ 
  coord_polar() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

ggsave("bar_v_circle_urine_order_2.pdf", width=10, height=10)

################## linear regression  ################

line_list <- matrix(data = NA, nrow = nrow(reads_table_v_abundance), ncol =3)
line_list <- as.data.frame(line_list)
colnames(line_list) = c('Species','p-value','R-value')
row.names(line_list)= row.names(reads_table_v_abundance)
line_list$Species = row.names(reads_table_v_abundance)

for (a in 1:nrow(line_list)) {
  # no treatment
  line <- as.data.frame(t(rbind(reads_table_v_abundance[a,],reads_table_u_abundance[a,])))
  
  line <- log10(line)
  
  colnames(line) = c('Vagina','Urine')
  
  line2=line
  line2[line2 == '-Inf'] = -6
  linearMod <- lm(Vagina ~ Urine, data=line2) 
  r_squared <- summary(linearMod)$r.squared
  pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
  
  line_list[a,2]= pvalue
  line_list[a,3]= r_squared
  
  if (line_list[a,1] %in% c('Sneathia amnii ','Lachnospiraceae BVAB1 ','TM7 OTU-H1 ','Prevotella cluster2 ')) {
    ggplot(line, aes(Vagina, Urine) )  +
      geom_point(size=1) +
      xlab(paste0("Vagina abundance")) +
      ylab(paste0("Urine abundance")) + ggtitle(paste0(line_list[a,1], '\n',"R = ", r_squared, '\n',"pvalue = ", pvalue))  + 
      geom_smooth(method='lm', se=T) + 
      theme(axis.title.x = element_text( size=16),
            axis.title.y = element_text( size=16)) 
    ggsave(paste0('linear_2_',line_list[a,1],'.pdf'), width=5.5, height=5.5)
  }
  
}
line_list$adj_p = p.adjust(line_list$`p-value`)  #  Benjamini & Hochberg 
write.csv(line_list,'linear_regression_2.csv')










