library(data.table)
library(vegan)
library(forcats)
library(NearestBalance)
library(stringr)
library(plyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggh4x)
library(clusterProfiler)
library(MicrobiomeProfiler)
library(SpiecEasi)
library(igraph)
library(compositions)
library(robCompositions)
library(outliers)

source('data/functions_general.R')
source('data/functions_local.R')

dir.create('out_microbiome')

# Здесь определяется переменная n_sim, а также загружается файл и сохраняется в переменной
n_sim <- 100
meta.points <- fread('data/for_Anya_meta_points.txt')

# загрузка данных из файла otu_table_2271.txt, пропуская первую строку. Подготовка данных с использованием функции PrepareSilvaOtuTab,
# которая создает список data.all.list, содержащий матрицу данных (mat) и список таксономии (taxa). 
# Затем извлекаются уникальные таксоны и сохраняются в переменной tax.all.
data.t <-  fread('D:/Linux/Study/ИТМО/R/bariatricSurgery1/NB_local/input/otu_table_2271.txt', skip = 1)
data.all.list <- PrepareSilvaOtuTab(otu.dt = data.t, cont = 'Wolbachia')
data.all <- data.all.list$mat
tax.all <- unique(data.all.list$taxa)

# определяются уровни таксономии в переменной scen, а затем для каждого уровня происходит агрегация данных
scen =  c('species', 'genus', 'family', 'order', 'class', 'phylum')
mat.cl <- sapply(scen, function(t.i) {
  AggregateTaxa(m = data.all, taxa.i = tax.all, agg.to.level=T, level = t.i, extract = F) 
}, simplify = F)

# фильтрует данные по количеству прочтений (n.reads = 3000) и по процентному содержанию (perc = 0.5) 
mat.cl.f <- lapply(mat.cl, function(x) FilterSamples(x, n.reads = 3000))
mat.cl.f <- sapply(mat.cl.f, function(x) FilterTaxa(mat = x[meta.points$sample,], perc = 0.5, 
                                                    num.more.perc = nrow(mat.cl.f$species)/10), simplify = F)

# объединяет все переменные из матриц mat.cl.f и сохраняет результат в переменной mat.cl.all.vars
mat.cl.all.vars <- GetAllVars(mat.cl.f, scen = names(mat.cl))

# обновляет переменную meta.points, оставляя только те строки, где значение переменной sample присутствует в именах строк переменной mat.cl.all.vars$f$species.
meta.points <- meta.points[sample %in% rownames(mat.cl.all.vars$f$species)]
# Эта часть кода проверяет, есть ли строки в переменной mat.cl.all.vars$f$species, которых нет в переменной meta.points$sample.
rownames(mat.cl.all.vars$f$species)[!(rownames(mat.cl.all.vars$f$species) %in% meta.points$sample)]

####################################################
# FIGURE 1 TIME POINTS
####################################################
lala <- as.matrix(dcast(formula = batch ~ point, meta.points[,.N, by=c('batch', 'point')], value.var = 'N')[, -1])
fisher.test(dcast(formula = batch ~ point, meta.points[,.N, by=c('batch', 'point')], value.var = 'N')[, -1])
meta.points.w <- dcast(meta.points, subj ~ point, value.var = 'sample')
meta.points.w.m <- (meta.points.w[, -1])
rownames(meta.points.w.m) <- meta.points.w$subj
meta.points.w.m <- !is.na(meta.points.w.m)
meta.points.w.m.dt <- data.table(meta.points.w.m)
meta.points.w.m.dt <- meta.points.w.m.dt[, .N, by = c('1', '2', '3')]
meta.points.w.m.dt[, or := meta.points.w.m.dt[, get('1')] + meta.points.w.m.dt[, get('2')] + meta.points.w.m.dt[, get('3')]]
setorder(meta.points.w.m.dt, 'or')
molten <- melt(as.matrix(meta.points.w.m.dt[, c(1, 2, 3), with=F])*meta.points.w.m.dt$N)

pdf('out_microbiome//time_points_ini.pdf', height = 5, width = 4)
ggplot(molten,
       aes(x = str_to_title((Var2)), 
           y = Var1,
           size = value)) +
  scale_size(range = c(1, 20)) +
  geom_point(aes(size=value),shape=21,color='#02315E', fill = '#02315E', alpha = 0.5) +#fill = add.alpha('#469597', 0.8)) +
  geom_text(aes(label = value), 
            colour = "white", 
            size = 5) +
  scale_fill_manual(values = c("darkorange"))+
  theme(legend.key=element_rect(fill="white")) +
  theme(panel.grid.major=element_line(linetype=1,color=add.alpha("grey", 0.5)),
        panel.background = element_blank()) +
  scale_y_continuous(breaks = seq(0, 7, 1)) +
  ylab("") + 
  coord_cartesian(clip = 'off') +
  xlab("Time points") 
dev.off()
meta.points$subj <- paste0('subj', meta.points$subj)

pol <- list(p12=c(1, 2), p13=c(1, 3), p23=c(2, 3))
meta.list <- list(p12 = meta.points[subj %in% meta.points[point %in% c(1, 2), .N, by='subj'][N==2]$subj & point %in% c(1, 2)],
                  p13 = meta.points[subj %in% meta.points[point %in% c(1, 3), .N, by='subj'][N==2]$subj & point %in% c(1, 3)],
                  p23 = meta.points[subj %in% meta.points[point %in% c(2, 3), .N, by='subj'][N==2]$subj & point %in% c(2, 3)])


###############################################################################
# read and filter metadata
################################################################################
meta <- fread('data/all_meta_new_new.txt')
meta.add <- fread('data/ca_ter.txt')
meta.l <- PreprocMeta(meta=meta, samples=rownames(mat.cl.all.vars$f$species), meta.add=meta.add, do.grubbs.1=T, meta.points=meta.points)
meta.change <- ParamChanges(meta.f = meta.l$meta.f, meta.list = meta.list, pol = pol, do.grubbs = T, meta.points=meta.points)

meta.f <- meta.l$meta.f
meta.f.1 <- meta.l$meta.f.1

###########################
# ILR transformation and PCOA with different colors
############################
counts.ilr <- ilr(mat.cl.all.vars$z$genus)  # ILR transformation
dis.ilr <- dist(counts.ilr)
dis.ilr.m <- as.matrix(dis.ilr)
pca.res.ilr <- ape::pcoa(dis.ilr.m)
data.pic.ilr <- mat.cl.all.vars$rar$genus / rowSums(mat.cl.all.vars$rar$genus)
colnames(data.pic.ilr) <- MakeNiceNamesSilva(colnames(data.pic.ilr))
cols.ilr <- c('orange', 'cadetblue')
col.samp.ilr <- sapply(rownames(data.pic.ilr), function(x) {
  if (x %in% meta.points[batch == 'y.2021']$sample) {
    cols.ilr[1]    
  } else if (x %in% meta.points[batch == 'y.2022']$sample) {
    cols.ilr[2]
  }
})

cols.ilr.2 <- beauty.pallete(4, contrast = TRUE)[2:4]
col.samp.ilr.2 <- sapply(rownames(data.pic.ilr), function(x) {
  if (x %in% meta.points[point == 1]$sample) {
    cols.ilr.2[1]    
  } else if (x %in% meta.points[point == 2]$sample) {
    cols.ilr.2[2]
  } else if (x %in% meta.points[point == 3]$sample) {
    cols.ilr.2[3]
  }
})

cols.ilr.3 <- c('red', 'blue')
col.samp.ilr.3 <- sapply(rownames(data.pic.ilr), function(x) {
  if (x %in% meta.f[met == 1]$sample) {
    cols.ilr.3[1]    
  } else if (x %in% meta.f[met == 0]$sample) {
    cols.ilr.3[2]
  }
})

ilr_data <- data.frame(counts.ilr)
