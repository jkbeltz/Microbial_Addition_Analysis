
### TWO FILES REQUIRED FOR PCA Figures, BOTH DERIVED FROM DF WHICH IS GENERATED HERE (LINK) ###

pca_data=read.csv("raw/pca_data.A.AT.LB.AllTpt.csv")
pca_data245=read.csv("raw/pca_data.A.AT.LB.Tpt2_4_5.csv")

####PCA ALL TIMEPOINTS ####
#Just looking at Apple vs bacterial conditioned - all tpts
load('./orch2020.RData')
df = cbind(samps, t(afmat))

df = df %>% filter(treatment != 'B')
df = df %>% rename(Timepoint = tpt)

samps = df[,1:ncol(samps)]
df = df[, -c(1:ncol(samps))]

#Run pca on all treatments/time points for A, AT, and LB (even though AT and LB don't have timepoints 1 an 3)
pca_res <- prcomp(na.omit(df), scale = TRUE, center = TRUE)
pca_data <- as.data.frame(pca_res$x)
pca_data <- cbind(samps, pca_data)
pca_data$Timepoint = as.numeric(pca_data$Timepoint)
write.csv(pca_data, 'pca_data.A.AT.LB.AllTpt.csv', row.names = FALSE)
pca_data=read.csv("raw/pca_data.A.AT.LB.AllTpt.csv")
View(pca_data)

b = ggplot(pca_data, aes(x = PC1, y = PC2, colour = Timepoint, shape = treatment))+
  geom_point(size = 5) +
  scale_colour_gradientn(colours = c("red","lightgrey","blue")) +
  scale_shape_manual(name = "Cage Treatment", 
                      breaks = c("A", "AT", "LB"),
                      labels = c("Control","+AT", "+LB"),
                      values = c(15,17,2))+
  theme_bw(base_size = 15)

####PCA TIMEPOINTS 2,4,5 ####
#Just looking at Apple vs bacterial conditioned for time points in which we have AT and LB data (2, 4, 5)
load('./orch2020.RData')
df = cbind(samps, t(afmat))

df = df %>% filter(treatment != 'B' & tpt %in% c('2', '4', '5'))
df = df %>% rename(Timepoint = tpt)

samps = df[,1:ncol(samps)]
df = df[, -c(1:ncol(samps))]

#Run pca on all treatments/time points for A, AT, and LB (even though AT and LB don't have timepoints 1 an 3)
pca_res <- prcomp(na.omit(df), scale = TRUE, center = TRUE)
pca_data <- as.data.frame(pca_res$x)
pca_data <- cbind(samps, pca_data)
pca_data$Timepoint = as.numeric(pca_data$Timepoint)
write.csv(pca_data, 'pca_data.A.AT.LB.Tpt2_4_5.csv', row.names = FALSE)
pca_data245=read.csv("raw/pca_data.A.AT.LB.Tpt2_4_5.csv")

a=ggplot(pca_data245, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment,group = treatment)) +
  geom_point(size = 5) +
  scale_color_gradientn(colors = c("red","lightgrey","blue")) +
  scale_shape_manual(name = "Cage Treatment", 
                     breaks = c("A", "AT", "LB"),
                     labels = c("Control","+AT", "+LB"),
                     values = c(15,17,2))+
  theme_bw(base_size = 15)

####Combined PCA plot####
Fig3PCA = ggarrange(b,a,labels = c("A","B"), ncol=2, common.legend = TRUE)
Fig3PCA

####GENERATING FILES FOR HEATMAPS####
##Using Apple cluster sites (identified from tpt 1->5) and testing their dynamics in AT and LB cages:

#First, get mean allele frequency shifts in AT, LB, and A
load('../RData/orch2020.RData')
df = cbind(samps, t(afmat))

#mean allele frequency shifts in AT and LB Cages
for (treat in c('AT', 'LB') ){
  df.treat = df %>% filter(treatment == treat)
  samps.treat = df.treat[,c(1:ncol(samps))]
  afmat.treat = df.treat[,-c(1:ncol(samps))]
  afmat.treat = as.data.frame(t(afmat.treat))
  shifts = get_af_shifts(afmat = afmat.treat, samps = samps.treat, cage_set = NULL, 
                         comparisons = c('2_4', '4_5', '2_5'))
  write.csv(shifts, paste0('shifts.', treat, '.csv'), row.names = FALSE)
}

#mean allele frequency shifts in A Cages

for (treat in c('A') ){
  df.treat = df %>% filter(treatment == treat)
  samps.treat = df.treat[,c(1:ncol(samps))]
  afmat.treat = df.treat[,-c(1:ncol(samps))]
  afmat.treat = as.data.frame(t(afmat.treat))
  shifts = get_af_shifts(afmat = afmat.treat, samps = samps.treat, cage_set = NULL, 
                         comparisons = c('2_4', '4_5', '2_5'))
  write.csv(shifts, paste0('shifts.', treat, '.csv'), row.names = FALSE)
}

#shifts for 'M' (combined AT + LB)
df = df %>% mutate(treatment = case_when(treatment == 'AT' ~ 'M',
                                         treatment == 'LB' ~ 'M',
                                         treatment == 'A' ~ 'A',
                                         treatment == 'B' ~ 'B'))

for (treat in c('M') ){
  df.treat = df %>% filter(treatment == treat)
  samps.treat = df.treat[,c(1:ncol(samps))]
  afmat.treat = df.treat[,-c(1:ncol(samps))]
  afmat.treat = as.data.frame(t(afmat.treat))
  shifts = get_af_shifts(afmat = afmat.treat, samps = samps.treat, cage_set = NULL, 
                         comparisons = c('2_4', '4_5', '2_5'))
  write.csv(shifts, paste0('shifts.', treat, '.csv'), row.names = FALSE)
}
#Now isolate allele frequency shifts for sites w/in each apple cluster
sites.clust = read.csv('./cluster.sites.Apple.csv')
sites.clust = sites.clust %>% mutate(snp = paste0(chrom, pos))
shifts.a = read.csv('./shifts.A.csv')
shifts.a = cbind(sites, shifts.a)
shifts.a = left_join(sites.clust, shifts.a)
sign.2_4 = sign(shifts.a$dAF.2_4) ##sign of apple shift to phase other treatments by
sign.4_5 =sign(shifts.a$dAF.4_5) ##sign of apple shift to phase other treatments by
sign.2_5 = sign(shifts.a$dAF.2_5) ##sign of apple shift to phase other treatments by
Joining with `by = join_by(chrom, pos, snp)`
#Shifts in AT and LB at target cluster sites
for (treat in c('AT', 'LB') ){
  shifts = read.csv(paste0('./shifts.',treat ,'.csv'))
  shifts = cbind(sites, shifts)
  shifts = shifts %>% mutate(snp = paste0(chrom, pos))
  
  sites.clust.target = sites.clust %>% dplyr::select(chrom, pos, snp)
  sites.clust.matched = sites.clust %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched) %>%
    mutate(snp = paste0(chrom, pos))
  shifts.target = left_join(sites.clust.target, shifts)
  names(shifts.target) = c('chrom', 'pos', 'snp', 'dAF.2_4.target', 'dAF.4_5.target', 'dAF.2_5.target')
  #phase target shifts in accordance w/ Apple shifts
  shifts.target = shifts.target %>% mutate(
    dAF.2_4.target = dAF.2_4.target * sign.2_4,
    dAF.4_5.target = dAF.4_5.target * sign.4_5,
    dAF.2_5.target = dAF.2_5.target * sign.2_5)
  shifts.matched = left_join(sites.clust.matched, shifts)
  names(shifts.matched) = c('chrom', 'pos.matched', 'snp', 'dAF.2_4.matched', 'dAF.4_5.matched', 'dAF.2_5.matched')
  
  shifts.clust = cbind(sites.clust %>% dplyr::select(cluster, cluster.start, cluster.end, afShift),
                       shifts.target %>% dplyr::select(pos,'dAF.2_4.target', 'dAF.4_5.target', 'dAF.2_5.target'),
                       shifts.matched %>% dplyr::select(pos.matched, 'dAF.2_4.matched', 'dAF.4_5.matched', 'dAF.2_5.matched'))
  write.csv(shifts.clust, paste0('./AClusterSiteShifts.', treat, '.csv'), row.names = FALSE)
}

#Shifts in M at target cluster sites
for (treat in c('M') ){
  shifts = read.csv(paste0('./shifts.',treat ,'.csv'))
  shifts = cbind(sites, shifts)
  shifts = shifts %>% mutate(snp = paste0(chrom, pos))
  
  sites.clust.target = sites.clust %>% dplyr::select(chrom, pos, snp)
  sites.clust.matched = sites.clust %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched) %>%
    mutate(snp = paste0(chrom, pos))
  shifts.target = left_join(sites.clust.target, shifts)
  names(shifts.target) = c('chrom', 'pos', 'snp', 'dAF.2_4.target', 'dAF.4_5.target', 'dAF.2_5.target')
  #phase target shifts in accordance w/ Apple shifts
  shifts.target = shifts.target %>% mutate(
    dAF.2_4.target = dAF.2_4.target * sign.2_4,
    dAF.4_5.target = dAF.4_5.target * sign.4_5,
    dAF.2_5.target = dAF.2_5.target * sign.2_5)
  shifts.matched = left_join(sites.clust.matched, shifts)
  names(shifts.matched) = c('chrom', 'pos.matched', 'snp', 'dAF.2_4.matched', 'dAF.4_5.matched', 'dAF.2_5.matched')
  
  shifts.clust = cbind(sites.clust %>% dplyr::select(cluster, cluster.start, cluster.end, afShift),
                       shifts.target %>% dplyr::select(pos,'dAF.2_4.target', 'dAF.4_5.target', 'dAF.2_5.target'),
                       shifts.matched %>% dplyr::select(pos.matched, 'dAF.2_4.matched', 'dAF.4_5.matched', 'dAF.2_5.matched'))
  write.csv(shifts.clust, paste0('./AClusterSiteShifts.', treat, '.csv'), row.names = FALSE)
}
Joining with `by = join_by(chrom, pos, snp)`
Joining with `by = join_by(chrom, pos, snp)`
Joining with `by = join_by(chrom, pos, snp)`
Joining with `by = join_by(chrom, pos, snp)`
#stats for each treatment and time point comparison
for (treat in c('AT', 'LB', 'M')){
  af.shifts = read.csv(paste0('./AClusterSiteShifts.',treat,'.csv'))
  data.meta.treat = data.frame()
  for (comp in c('2_4', '4_5', '2_5')){
    data.meta = data.frame()
    for(clust in (as.character(unique(af.shifts$cluster)))){
      d = data.frame()
      df.clust = af.shifts %>% filter(clust == cluster)
      cluster = clust
      shifts = df.clust %>% dplyr::select(contains(comp))
      median.target = median(shifts[,1])
      median.matched = median(shifts[,2])
      pval = as.numeric(t.test(shifts[,1], shifts[,2])$p.value)
      d = as.data.frame(cbind(cluster, median.target, median.matched, pval))
      d$comparison = comp
      data.meta = rbind(data.meta, d)
      
    }
    data.meta = data.meta %>% mutate(pval = as.numeric(as.character(pval)))
    data.meta = data.meta %>% mutate(median.target = as.numeric(as.character(median.target)))
    data.meta = data.meta %>% mutate(FDR = p.adjust(pval, method = 'BH'))
    data.meta = data.meta %>% mutate(median.target = as.numeric(as.character(median.target)),
                                     median.matched = as.numeric(as.character(median.matched)))
    data.meta = data.meta %>% rowwise() %>% 
      mutate(significant = if_else(FDR < 0.05 & abs(as.numeric(median.target)) > abs(as.numeric(median.matched)), 'yes', 'no'))
    data.meta = data.meta %>% rowwise() %>% 
      mutate(cluster.dynamics = case_when(significant == 'yes' & median.target > 0 ~ 'positive',
                                          significant == 'yes' & median.target < 0 ~ 'negative',
                                          significant == 'no' ~ 'neutral')) 
    data.meta.treat = rbind(data.meta.treat, data.meta)
    write.csv(data.meta.treat, paste0('./AClusterShiftStats.', treat, '.csv'), row.names = FALSE)
  }
}
##look at stats:
df.at = read.csv('./AClusterShiftStats.AT.csv')
df.at %>% group_by(comparison, cluster.dynamics) %>% count()
df.lb = read.csv('./AClusterShiftStats.LB.csv')
df.lb %>% group_by(comparison,cluster.dynamics) %>% count()
df.m = read.csv('./AClusterShiftStats.M.csv')
df.m %>% group_by(comparison, cluster.dynamics) %>% count()


####Importing files####
###Plotting cluster stats/shifts shifts in AT, LB, and M cages:
df.at = read.csv('raw/AClusterShiftStats.AT.csv')
df.lb = read.csv('raw/AClusterShiftStats.LB.csv')
df.m = read.csv('raw/AClusterShiftStats.M.csv')
####AT####
df = df.at %>% rowwise() %>% mutate(clust = strsplit(as.character(cluster), '[.]')[[1]][1],
                                    chrom = strsplit(as.character(cluster), '[.]')[[1]][2])

#add index values for all clusters (irrespective of chromosome)
df$idx =  c(1:nrow(df))
df = df %>% mutate(idx = as.numeric(as.character(idx)))

#add cluster numbers for each chromosome
df = df %>% 
  arrange(chrom, clust) 
df = df %>%
  group_by(chrom) %>%
  mutate(clust.num = row_number()) %>%
  ungroup()

names(df)[names(df) == 'median.target'] <- 'Median.Shift'

AT25 = ggplot(df %>% filter(comparison == '2_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="darkblue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
AT25
ggsave('../Figures/AClusters.ATTest.2_5.pdf', p, height = 10, width = 10)

AT24 = ggplot(df %>% filter(comparison == '2_4'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="darkblue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
AT24
ggsave('../Figures/AClusters.ATTest.2_4.pdf', p, height = 10, width = 10)

AT45 = ggplot(df %>% filter(comparison == '4_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="darkblue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
AT45
ggsave('../Figures/AClusters.ATTest.4_5.pdf', p, height = 10, width = 10)



####LB####
df = df.lb %>% rowwise() %>% mutate(clust = strsplit(as.character(cluster), '[.]')[[1]][1],
                                    chrom = strsplit(as.character(cluster), '[.]')[[1]][2])


#add index values for all clusters (irrespective of chromosome)
df$idx =  c(1:nrow(df))
df = df %>% mutate(idx = as.numeric(as.character(idx)))

#add cluster numbers for each chromosome
df = df %>% 
  arrange(chrom, clust) 
df = df %>%
  group_by(chrom) %>%
  mutate(clust.num = row_number()) %>%
  ungroup()

names(df)[names(df) == 'median.target'] <- 'Median.Shift'

LB25 = ggplot(df %>% filter(comparison == '2_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
LB25
ggsave('../Figures/AClusters.LBTest.2_5.pdf', p, height = 10, width = 10)

LB24 = ggplot(df %>% filter(comparison == '2_4'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
LB24
ggsave('../Figures/AClusters.LBTest.2_4.pdf', p, height = 10, width = 10)

LB45 = ggplot(df %>% filter(comparison == '4_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 4, 
            vjust = -0.5, hjust = 0) 
LB45
ggsave('../Figures/AClusters.LBTest.4_5.pdf', p, height = 10, width = 10)



####M####
df = df.m %>% rowwise() %>% mutate(clust = strsplit(as.character(cluster), '[.]')[[1]][1],
                                   chrom = strsplit(as.character(cluster), '[.]')[[1]][2])


#add index values for all clusters (irrespective of chromosome)
df$idx =  c(1:nrow(df))
df = df %>% mutate(idx = as.numeric(as.character(idx)))

#add cluster numbers for each chromosome
df = df %>% 
  arrange(chrom, clust) 
df = df %>%
  group_by(chrom) %>%
  mutate(clust.num = row_number()) %>%
  ungroup()

names(df)[names(df) == 'median.target'] <- 'Median.Shift'

M25 = ggplot(df %>% filter(comparison == '2_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 5, 
            vjust = -0.5, hjust = 0) 
M25
ggsave('../Figures/AClusters.MTest.2_5.pdf', p, height = 10, width = 10)

M24 = ggplot(df %>% filter(comparison == '2_4'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 5, 
            vjust = -0.5, hjust = 0) 
M24
ggsave('../Figures/AClusters.MTest.2_4.pdf', p, height = 10, width = 10)

M45 = ggplot(df %>% filter(comparison == '4_5'), aes(y = chrom, x = clust.num, fill = Median.Shift)) +
  geom_tile() +
  theme_classic(base_size = 15) +
  xlab('Cluster') +
  ylab('Chromosome') +
  scale_fill_gradient2(low="dark blue", mid="white", high="red", 
                       midpoint=0, limits=range(df$Median.Shift)) +
  geom_text(aes(label = ifelse(significant == "yes", "*", "")), 
            color = "black", 
            size = 5, 
            vjust = -0.5, hjust = 0) 
M45
ggsave('../Figures/AClusters.MTest.4_5.pdf', p, height = 10, width = 10)
####combine plots####
CLMAP25 = ggarrange(AT25 + rremove("xlab") ,LB25+ rremove("xlab"),M25, labels = c("A","B","C"), ncol = 1, common.legend = TRUE  )
CLMAP25
CLMAP24 = ggarrange(AT24 + rremove("xlab") ,LB24+ rremove("xlab"),M24, labels = c("A","B","C"), ncol = 1, common.legend = TRUE  )
CLMAP24
CLMAP45 = ggarrange(AT45 + rremove("xlab") ,LB45+ rremove("xlab"),M45, labels = c("A","B","C"), ncol = 1, common.legend = TRUE  )
CLMAP45

Fig3 = ggarrange(Fig3PCA, AT25 + rremove("xlab"), LB25, labels = c("", "C", "D"), ncol = 1, heights = c(1,.75,.75))
Fig3

####Compariosn between AT and LB####
##how many clusters w/ treatment specific effects?
df.a = df.at %>% dplyr::select(cluster, comparison, cluster.dynamics) 
names(df.a)[names(df.a) == 'cluster.dynamics'] <- 'at.dynamics'
df.l = df.lb %>% dplyr::select(cluster.dynamics) 
names(df.l)[names(df.l) == 'cluster.dynamics'] <- 'lb.dynamics'
d = cbind(df.a, df.l)
d = cbind(df.a, df.l)
d = d %>% mutate(concordance = if_else(at.dynamics == lb.dynamics, 'yes', 'no'))
#opposite dynamics:
d %>% filter(at.dynamics != 'neutral' & lb.dynamics != 'neutral') %>% filter(concordance == 'no')


A data.frame: 32 × 5
cluster	comparison	at.dynamics	lb.dynamics	concordance
<fct>	<fct>	<fct>	<fct>	<chr>
  c1.X.1_5	2_4	positive	negative	no
c3.X.1_5	2_4	positive	negative	no
c4.3R.1_5	2_4	positive	negative	no
c7.X.1_5	2_4	negative	positive	no
c14.2R.1_5	2_4	positive	negative	no
c25.2R.1_5	2_4	negative	positive	no
c1.X.1_5	4_5	negative	positive	no
c3.3L.1_5	4_5	negative	positive	no
c2.X.1_5	4_5	negative	positive	no
c5.X.1_5	4_5	positive	negative	no
c4.2R.1_5	4_5	negative	positive	no
c7.2R.1_5	4_5	positive	negative	no
c8.2R.1_5	4_5	positive	negative	no
c16.2L.1_5	4_5	positive	negative	no
c7.3L.1_5	4_5	negative	positive	no
c18.2L.1_5	4_5	positive	negative	no
c20.2L.1_5	4_5	positive	negative	no
c21.2L.1_5	4_5	positive	negative	no
c15.3R.1_5	4_5	positive	negative	no
c17.3R.1_5	4_5	positive	negative	no
c28.2L.1_5	4_5	negative	positive	no
c21.3L.1_5	4_5	negative	positive	no
c16.X.1_5	4_5	negative	positive	no
c25.3R.1_5	4_5	negative	positive	no
c3.2R.1_5	2_5	positive	negative	no
c6.2R.1_5	2_5	positive	negative	no
c7.2R.1_5	2_5	positive	negative	no
c7.X.1_5	2_5	negative	positive	no
c21.2L.1_5	2_5	positive	negative	no
c23.2L.1_5	2_5	positive	negative	no
c15.3R.1_5	2_5	positive	negative	no
c12.X.1_5	2_5	positive	negative	no
