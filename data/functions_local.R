ParamChanges <- function(meta.f, meta.list, pol, do.grubbs=T, meta.points) {
  sapply(names(pol), function(x) {
    message(x)
    meta.f.1 <- meta.f[, ! "длительность.заболевания", with=F]
    meta.f.2 <- meta.f[, ! "длительность.заболевания", with=F]
    meta.f.2 <- meta.f.2[sample %in% meta.points[point==pol[[x]][2]]$sample]
    meta.f.1[is.na(ремиссия)]$ремиссия <- 0
    meta.f.1 <- meta.f.1[sample %in% meta.points[point==pol[[x]][1]]$sample]
    meta.points.12.w <- reshape2::dcast(formula = subj ~ point, data = meta.list[[x]], value.var = 'sample')
    meta.f.ch <- - left_join(data.frame(sample=meta.points.12.w[, as.character(pol[[x]][1])]), meta.f.1)[, -1] + 
      left_join(data.frame(sample=meta.points.12.w[, as.character(pol[[x]][2])]), meta.f.2)[, -1]
    meta.f.ch$sample <- meta.points.12.w[, as.character(pol[[x]][1])]
    meta.f.ch <- data.table(meta.f.ch)
    if (do.grubbs) {
      meta.f.ch[, (colnames(meta.f.ch)[!colnames(meta.f.ch) %in% c('sample', 'ремиссия', 'met')]) := lapply(.SD, function(x){
        names(x) <- paste0('s', 1:length(x))
        x[GrubbsFlag(x)] <- NA
        x
      }), .SDcols=colnames(meta.f.ch)[!colnames(meta.f.ch) %in% c('sample', 'ремиссия', 'met')]]
    }
    meta.f.ch
  }, simplify = F)
}

PreprocMeta <- function(meta, samples, meta.add=NULL, do.grubbs.1=T, meta.points) {
  meta[, .N, by='Номер образца кала'][N>1]
  meta[, `Номер образца кала`:=paste0('S', `Номер образца кала`)]
  meta[duplicated(meta$`Номер образца кала`)]
  meta[duplicated(meta$`Номер образца кала`), `Номер образца кала` := paste0(`Номер образца кала`, '.new')]
  
  meta <- meta[(`Номер образца кала` %in% samples)]
  samples[!(samples %in% meta$`Номер образца кала`)]
  
  n.na <- apply(meta, 2, function(x) length(which(is.na(x))))
  meta.f <- copy(meta)
  names(n.na)[n.na>=20]
  n.na <- n.na[names(n.na)!='ремиссия']
  n.na <- n.na[!names(n.na) %in% c("длительность заболевания")]
  meta.f[, names(n.na)[n.na>=20]:=NULL]
  meta.f$sample <- meta.f$`Номер образца кала`
  
  if (!is.null(meta.add)) {
    meta.add[, met := str_detect(medication, 'Метформин|ГалвусМет|Меформин')]
    meta.add[, sample:= paste0('S', sample)]
    meta.f <- merge(meta.f, meta.add[, c('Ca/alb', 'met', 'sample')], by = 'sample')
    meta.f$met <- as.numeric(meta.f$met)
  } else {
    meta.f[, met := str_detect(`ССП терапия`, 'Метформин|метформин|ГалвусМет|Меформин')]
    meta.f$met <- as.numeric(meta.f$met)
  }
  meta.f[, c( 'Возраст',
              'Масса тела (кг)', 'Рост (см)', 'Номер образца кала',
              'Офтальмолог ретинопатия да - 1, нет - 2', 'ССП терапия'):= NULL]
  
  meta.f.1 <- meta.f[sample %in% meta.points[point==1]$sample]
  meta.f.1$ремиссия <- NULL
  
  # delete correlated factors
  combs <- combn(colnames(meta.f.1)[colnames(meta.f.1)!='sample'], 2)
  crosscor <- rbindlist(apply(combs, 2, function(x) {
    res <- cor.test(meta.f.1[, get(x[1])], meta.f.1[, get(x[2])])
    data.table(f1=x[1], f2=x[2], cor = res$estimate, pval = res$p.value)
  }))
  crosscor[abs(cor)>0.7]
  excl <- c('СКФ (EPI) мл/мин/1.73м2',
            'А/Кр  мг/ммоль',
            'инсулин, мкЕ/мл')
  meta.f[, c(excl):= NULL]
  meta.f.1[, c(excl):= NULL]
  
  # delete outliers
  if (do.grubbs.1) {
    meta.f.1[, (colnames(meta.f.1)[!colnames(meta.f.1) %in% c('sample')]) := lapply(.SD, function(x){
      names(x) <- paste0('s', 1:length(x))
      x[GrubbsFlag(x)] <- NA
      x
    }), .SDcols=colnames(meta.f.1)[!colnames(meta.f.1) %in% c('sample')]]
  }
  
  colnames(meta.f.1) <- make.names(colnames(meta.f.1))
  colnames(meta.f) <- make.names(colnames(meta.f))
  list(meta.f.1=meta.f.1, meta.f=meta.f)
}

GetDiffKegg <- function(res.nb.pval.met.i, kp, back.comp.kegg) {
  k1 <- str_remove(convert_hmdb_to_kegg(res.nb.pval.met.i$hmdb.ids), 'cpd:')
  ko <- unlist(kp[names(kp) %in% paste0('cpd:', k1)])
  kpi <- unlist(kp[names(kp) %in% paste0('cpd:', back.comp.kegg)])
  kot <- data.table(comp = names(ko), pathway = (ko))
  kot.all <- data.table(comp = names(kpi), pathway = (kpi))
  kot.n <- kot[, length(unique(comp)), by = pathway]
  kot.all.n <- kot.all[, length(unique(comp)), by = pathway][pathway %in% kot.n$pathway]
  kot.n <- merge(kot.n, kot.all.n, by = 'pathway', suffixes = c('', '.all'))
  kot.n$K <- length(unique(k1))
  kot.n$M <- length(unique(back.comp.kegg))
  kot.n[, pval := phyper(V1-1, V1.all, M-V1.all, K, lower.tail = F)]
  kot.n[, FDR := p.adjust(pval, 'BH')]
  kot.n[, name := all.h.path[str_remove(kot.n$pathway, 'path:')]]
  kot.n
}