library(zCompositions)

GrubbsFlag <- function(x)  {
  Outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  if (pv > 0.05) {
    return(NA)
  }
  while ((pv <= 0.05) & (!is.na(pv))) {
    lOutliers <- outlier(test, logical = TRUE)
    Outliers <- c(Outliers, names(test)[lOutliers])
    test <- x[!(names(x) %in% Outliers)]
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
    pv
    if (is.null(names(pv)) | is.na(pv)) {
      return(which((names(x) %in% Outliers)))
    }
  }
  return(which((names(x) %in% Outliers)))
}

BiplotPcoa2D <- function(abund, pcoa.obj, col_sample, pch_sample = 21, thresh = 0.2,
                         alpha.bact = 1, seg.from = c(), seg.to = c(), alpha.seg = 0.1,
                         arrow.length = 0.1, n1 = 1, n2 = 2, base = "hit", ...) {
  Y <- abund
  n <- nrow(Y)
  points.stand <- scale(pcoa.obj$vectors[, c(n1, n2)])
  S <- cov(Y, points.stand)
  U <- S %*% diag((pcoa.obj$values$Eigenvalues[c(n1, n2)]/(n -
                                                             1))^(-0.5))
  colnames(U) <- colnames(pcoa.obj$vectors[, c(n1, n2)])
  perc_expl <- paste(round(100 * (pcoa.obj$values$Relative_eig)[1:10],
                           1), "%", sep = "")
  pn <- paste("PC", c(n1, n2), sep = "")
  perc_expl <- paste(pn, perc_expl[c(n1, n2)], sep = ": ")
  BiPlot2D(coords1 = pcoa.obj$vectors[, n1], coords2 = pcoa.obj$vectors[,
                                                                        n2],
           rotat1 = U[, 1], rotat2 = U[, 2], pch_sample = pch_sample,
           varsnames = rownames(U), col_sample, thresh = thresh,
           alpha.bact = alpha.bact, seg.from = seg.from, seg.to = seg.to,
           base = base, alpha.seg = alpha.seg, perc_expl = perc_expl,
           arrow.length = arrow.length, ...)
}


BiPlot2D <- function(coords1, coords2, rotat1, rotat2, varsnames, col_sample,
                     thresh = 0.2, alpha.bact = 1, pch_sample = 21, seg.from = c(),
                     seg.to = c(), alpha.seg = 0.1, perc_expl = c("PC1", "PC2"),
                     arrow.length = 0.1, base = "hit", ...) {
  data <- data.frame(x1 = coords1, x2 = coords2)
  rownames(data) <- names(coords1)
  datapc <- data.frame(varnames = varsnames, x1 = rotat1,
                       x2 = rotat2)
  mult <- min((max(data[, "x1"]) - min(data[, "x1"])/(max(datapc[,
                                                                 "x1"]) -
                                                        min(datapc[, "x1"]))),
              (max(data[, "x2"]) -
                 min(data[, "x2"])/(max(datapc[, "x2"]) -
                                      min(datapc[,
                                                 "x2"]))))
  datapc <- transform(datapc, v1 = 1 * mult * x1, v2 = 1 *
                        mult * x2)
  datapc <- datapc[(datapc$x1^2 + datapc$x2^2)^0.5 > thresh,  ]
  if (base == "hit") {
    old <- datapc$varnames
    datapc$varnames <- str_replace_all(old, ";unclassified",
                                       "_u")
    datapc$varnames <- str_replace(datapc$varnames, ".*;",
                                   "")
  }   else if (base == "gg") {
    datapc$varnames <- str_replace(datapc$varnames, ";.__$",
                                   "_u")
    datapc$varnames <- str_replace(datapc$varnames, ".*;",
                                   "")
  }   else {
    datapc$varnames <- datapc$varnames
  }
  plot(data$x1, data$x2, col = col_sample, bg = add.alpha(col_sample,
                                                          0.5),
       pch = pch_sample, xlim = c(min(data$x1, datapc$v1),
                                  
                                  max(data$x1, datapc$v1)),
       ylim = c(min(data$x2, datapc$v2),
                max(data$x2, datapc$v2)), xlab = perc_expl[1],
       ylab = perc_expl[2],
       ...)
  if (length(seg.to) > 0) {
    sapply(1:length(seg.to), function(i) {
      lines(c(data[seg.from[i], ]$x1, data[seg.to[i],
      ]$x1), c(data[seg.from[i], ]$x2, data[seg.to[i],
      ]$x2), col = add.alpha("black", alpha.seg))
    })
  }
  if (length(datapc$varnames > 0)) {
    text(x = datapc$v1, y = datapc$v2, label = datapc$varnames,
         col = add.alpha(rgb(0.3, 0.3, 0.3), alpha.bact))
    arrows(x = 0, y = 0, x1 = datapc$v1, y1 = datapc$v2,
           col = add.alpha(rgb(0.3, 0.3, 0.3), alpha.bact),
           cex = 0.1, angle = 10, length = arrow.length)
  }
}


beauty.pallete <- function (n, contrast = F) {
  if (n == 2) {
    if (contrast) {
      pal <- c("azure2", "lightcoral")
    }     else {
      pal <- c("skyblue3", "azure2")
    }
  }   else if (n > 2) {
    pal <- c("azure2", "lightcoral", "palegreen3", "darkcyan",
             "aquamarine4")
  }
  pal[1:n]
}

MakeNiceNamesSilva <- function(x, level = 's') {
  x1 <- str_replace_all(str_replace_all(x, '[gsfoc]_;', '-u;'), ';[gsfoc]_$', ';-u')
  if (level=='s') {
    x1 <- str_replace(x1, ';s_', '_')
  }
  x2 <- str_remove(x1, '.*;[gsfoc]_')
  
  sapply(1:length(x2), function(i) {
    if (str_detect(x2[i], '^[A-Z0-9\\-\\_]*$')) {
      x3 <- str_replace(x1[i], paste0(';', level, '_'), '-')
      str_remove(x3, '.*;[gsfoc]_')
    } else {
      x2[i]
    }
  })
  
}

ScatterBox <- function(list.to.draw, ylab='', xlab='', main = '',
                       xlim = c(1-jitter.num, length(list.to.draw)+jitter.num),
                       jitter.num = 0.5, col.points = 'slategray3', col.box = 'grey', ...) {
  plot(jitter(unlist(lapply(1:length(list.to.draw), function(i) {rep(i, length(list.to.draw[[i]]))})), jitter.num),
       unlist(list.to.draw),  col = col.points, pch = 21, bg=add.alpha('slategray3', 0.3),
       bty='n', xaxt='n', xlab=xlab, ylab=ylab, main = main, xlim = xlim)
  boxplot(list.to.draw, add=T, frame = F, outline=F, col = col.box, ...)
}


WriteTableCust <- function(dt, filename) {
  write.table(dt, file = filename, sep = '\t', quote = F, row.names = F)
}


MakeStar <- function (x1, y1, x2, y2, y.star, cex = 1.3, label = "*") {
  segments(c(x1, x1, x2), c(y1, y2, y2), c(x1, x2, x2), c(y2,
                                                          y2, y1))
  text((x1 + x2)/2, y.star, label, cex = cex)
}

heatmap_with_split_loc <- function(abund, metadata, formula=NULL, clr=NULL,
                               type = c("perc", "clr"), bal_list=list(),
                               num_name = "num", den_name = "den",
                               abund_limits = range(abund),
                               taxa_colors = c(),
                               sample_col = "sample"){
  pl_type <- match.arg(type)
  
  if (pl_type == "clr"){
    if (is.null(clr)) {
      clr <- t(apply(abund, 1, function(x) log(x) - mean(log(x))))
    }
    dat <- clr
  }  else {
    dat <- abund*100
  }
  
  abund_melted <- reshape2::melt(as.matrix(dat))
  abund_melted$Var1 <- as.character(abund_melted$Var1)
  abund_merged <- data.table(merge(abund_melted, metadata,
                                   by.x = "Var1", by.y = "sample"))
  
  if(length(bal_list) > 0){
    taxa_descr <- rbindlist(lapply(bal_list, function(x){
      rbind(data.frame(Var2 = x$num, taxa_gr=num_name),
            data.frame(Var2 = x$den, taxa_gr=den_name))
    }), idcol = "bal")
    taxa_descr[, taxa_gr := paste(bal, taxa_gr)]
    abund_merged   <- merge(abund_merged, taxa_descr[,.(Var2, taxa_gr)],
                            all.x=T, by = "Var2")
  }
  abund_merged[ , max_abund := max(value), by = .(Var2)]
  # abund_merged$Var2 <- factor(abund_merged$Var2, levels=taxa_order)
  abund_merged %<>% mutate(Var2 = Var2  %>% fct_reorder(max_abund))
  yticks  <- unique(abund_merged$Var2)
  
  pl <- ggplot(abund_merged, aes(Var1, Var2, fill=value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_nested(formula, drop = T, scales = "free", space = "free") +
    xlab("") + ylab("")
  
  if (length(taxa_colors)>0){
    yticks <-
      ggplot_build(pl)$layout$panel_scales_y[[1]]$get_labels()
    colors <- taxa_colors[yticks]
    pl <- pl +
      theme(axis.text.y = element_text(colour = colors))
  }
  
  if (pl_type == "perc"){
    pl <- pl +
      scale_fill_gradientn(colors = c("#0072B2", "#009E73", "#F0E442", "#D55E00"),
                           breaks=c(0.001, 1, 10, 50,100), trans="log",
                           limits=abund_limits*100) +
      labs(fill="относительная\nпредставленность,\n%")
  } else if (pl_type == "clr"){
    pl <- pl +
      scale_fill_gradientn(colors = c("#0072B2", "#009E73", "#F0E442", "#D55E00"))+
      labs(fill="CLR-components")
  }
  pl
}

PrepareSilvaOtuTab <- function(otu.dt, cont=c()) {
  if (length(cont)>0) {
    otu.dt <- otu.dt[!str_detect(taxonomy, paste0(cont, collapse = '|'))]
  }
  lala.m <- t(as.matrix(otu.dt[, !c('#OTU ID', 'taxonomy'), with=F]))
  colnames(lala.m) <- otu.dt$`#OTU ID`
  taxa = otu.dt[, c('#OTU ID', 'taxonomy'), with=F]
  taxa[, taxonomy:= str_replace(taxonomy, '_+bacterium', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+human_gut', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+gut_group', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+prokaryote', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+possible_genus[^;]*', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+prokaryote', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+rumen', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+soil', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+forest', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+marine', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, '_+sediment', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, 's_+Firmicutes', 's_')]
  taxa[, taxonomy:= str_replace(taxonomy, 's_+Eubacterium', 's_')]
  taxa[, taxonomy:= str_replace(taxonomy, 's_+Desulfovibrionaceae', 's_')]
  taxa[, taxonomy:= str_replace(taxonomy, 's_+toluene-degrading_methanogenic', 's_')]
  taxa[, taxonomy:= str_replace_all(taxonomy, '_Incertae_Sedis', '_')]
  taxa[, taxonomy:= str_replace_all(taxonomy, '_+', '_')]
  taxa[, taxonomy:= str_replace(taxonomy, 's_[^_]*ae_*$', 's_')]
  list(mat=lala.m, taxa=taxa)
}


AggregateTaxa <- function(m, taxa.i, genus=F, agg.to.level=F, level = 'genus', extract = T) {
  if (extract) {
    seqs <- str_remove(str_remove(str_extract(colnames(m), '\\([ACGT]*\\)'), '\\)'), '\\(')
  } else {
    seqs <- colnames(m)
  }
  tmp <- copy(m)
  colnames(tmp) <- seqs
  tax.i <- taxa.i[taxa.i$`#OTU ID` %in% seqs]
  tax.i$taxonomy <- str_replace_all(tax.i$taxonomy, '; *', ';')
  tax.i$taxonomy <- str_replace_all(tax.i$taxonomy, '__', '_')
  if (is.na(level) & agg.to.level==F) {
    tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ' *')
    if (genus) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';s_.*')
    }
  }
  if (agg.to.level) {
    taxa <- c('species', 'genus', 'family', 'order', 'class', 'phylum')
    if (!(level %in% taxa)) {
      stop("level must be in c('species', 'genus', 'family', 'order', 'class', 'phylum')")
    }
    if (level == taxa[1]) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ' *')
    } else if (level == taxa[2]) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';s_.*')
    } else if (level == taxa[3]) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';g_.*')
    } else if (level == taxa[4]) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';f_.*')
    } else if (level == taxa[5]) {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';o_.*')
    } else {
      tax.i$taxonomy <- str_remove_all(tax.i$taxonomy, ';c_.*')
    }
  }
  arr <- tax.i[, as.list(rowSums(tmp[, unlist(.SD$`#OTU ID`), drop=F])), by = taxonomy]
  arr.m <- t(as.matrix(arr[, !'taxonomy', with=F]))
  colnames(arr.m) <- arr$taxonomy
  arr.m
}

FilterTaxa <- function(mat, num.more.perc=10, perc = NA, num.more.count=10, count.t = NA) {
  mat.i.p <- mat/rowSums(mat)*100
  mat.i <- mat
  if (!is.na(perc)) {
    mat.i <- mat.i[, apply(mat.i.p, 2, function(x) sum(x>perc))>num.more.perc]
  }
  if (!is.na(count.t)) {
    mat.i <- mat.i[, apply(mat.i, 2, function(x) sum(x>count.t))>num.more.count]
  }
  message(paste(dim(mat.i), collapse = ' '))
  mat.i
}


GetAllVars <- function(data.all.f, scen, rar=NA) {
  if (is.na(rar)) {
    data.all.rar <- sapply(scen, function(i) Rarefy(data.all.f[[i]], min(rowSums(data.all.f[[i]])))$otu.tab.rff, simplify = F)
  } else {
    data.all.rar <- sapply(scen, function(i) Rarefy(data.all.f[[i]], rar)$otu.tab.rff, simplify = F)
  }
  data.all.z.rar <- sapply(scen, function(i) {
    z <- data.all.rar[[i]]
    z.zer <- t(apply(z, 1, function(x) {x[x==0] <- 0.5; x}))
    z.zer}, simplify = F)
  data.all.z <- sapply(scen, function(i) {
    z <- data.all.f[[i]]
    z.zer <- t(apply(z, 1, function(x) {x[x==0] <- 0.5; x}))
    z.zer}, simplify = F)
  z.zc = sapply(scen, function(i) {
    if (any(data.all.f[[i]]==0)) {
      cmultRepl(data.all.f[[i]])
    } else {
      data.all.f[[i]]
    }
  }, simplify = F)
  list(f=data.all.f[scen], z = data.all.z, rar=data.all.rar, z.crepl=z.zc, z.rar=data.all.z.rar)
}

FilterSamples <- function(mat, n.reads=3000) {
  mat.i <- mat[rowSums(mat) >= n.reads,]
  message(paste(rownames(mat)[rowSums(mat) < n.reads], collapse = ' '))
  message(paste(dim(mat.i), collapse = ' '))
  mat.i
}