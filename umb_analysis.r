#########################################
# UTI microbiome study (UMB)            #
# Gut-bladder axis syndrome associated  #
# with recurrent UTI in humans          #
#                                       #
# Colin Worby cworby@broadinstitute.org #
# 2021-11-16                            #
#########################################

# Prepare data; define functions
#homedir <- "/path/to/dir/"
source(paste0(homedir,"code/umb_prep.r"))

#################################################
# FIGURE 2: Overview of gut microbiome profiles #
#################################################
# Taxonomic level
lvl <- 5
# Order participants by group and median RA of butyrate producers
ord <- order(ruti_pats==1,sapply(1:31,function(x)median(-scfa_ra[[1]][which(patient==x)])), decreasing=TRUE)
# Color for each cohort
cohortcolz <- c("black", "black")
# Color palette for main taxa
colz <- brewer.pal(8,"Set2")

layout(1:4, heights = c(3,3,3,0.7))
par(mar=c(1,4,0,1))

# Stacked bars for mean relative abundances
plot(NULL, xlim=c(0,31)+0.5, ylim=c(0,100),bty="n", xaxt="n", ylab="Relative abundance (%)",
     xlab="", las=1)
for (i in 1:length(ord)) {
  cml <- 0
  for (j in 1:7) {
    rect(i-0.3-ruti_pats[ord[i]],cml,i+0.3-ruti_pats[ord[i]],cml+mean(tax_abundances[[lvl]][j,which(patient==ord[i])]), col=colz[j])
    cml <- cml + mean(tax_abundances[[lvl]][j,which(patient==ord[i])])
  }
  rect(i-0.3-ruti_pats[ord[i]],cml,i+0.3-ruti_pats[ord[i]],100, col="lightgrey")
}
legend("bottom", legend = c(tax_names[[lvl]][1:7], "Other"), ncol=4, 
       pch=22, pt.bg = c(colz[1:7], "lightgrey"), pt.cex=2, bg="white")

# Plot range of butyrate relative abundances for each participant
plot(NULL, xlim=c(0,31)+0.5, ylim=c(0,50), xaxt="n", ylab="Butyrate producers (%)", 
     las=1, bty="n", xlab="")
boxer(lapply(1:31,function(x)scfa_ra[[1]][which(patient==x)])[ord],1:31-ruti_pats[ord], 
      width = 0.3, border = cohortcolz[ruti_pats[ord]+1], lwd=2)
lines(c(0,15)-0.5,rep(mean(sapply(which(ruti_pats==1),function(x)mean(scfa_ra[[1]][which(patient==x)]))),2), lwd=3,
      col=adjustcolor(cohortcolz[1],0.5))
lines(c(16,32)-0.5,rep(mean(sapply(which(ruti_pats==0),function(x)mean(scfa_ra[[1]][which(patient==x)]))),2), lwd=3,
      col=adjustcolor(cohortcolz[1],0.5))

# Plot range of microbial richness for each participant
plot(NULL, xlim=c(0,31)+0.5, ylim=c(0,80), xaxt="n", ylab="Microbial richness", 
     las=1, bty="n", xlab="")
boxer(lapply(1:31,function(x)c1[which(patient==x)])[ord],1:31-ruti_pats[ord], 
      width = 0.3, border = cohortcolz[ruti_pats[ord]+1], lwd=2)
lines(c(0,15)-0.5,rep(mean(sapply(which(ruti_pats==1),function(x)mean(c1[which(patient==x)]))),2), lwd=3,
      col=adjustcolor(cohortcolz[1],0.5))
lines(c(16,32)-0.5,rep(mean(sapply(which(ruti_pats==0),function(x)mean(c1[which(patient==x)]))),2), lwd=3,
      col=adjustcolor(cohortcolz[2],0.5))

# Plot abx and UTI data for each participant
plot(NULL, xlim=c(0,31)+0.5, ylim=c(0,2)+0.5, xaxt="n", yaxt="n", ylab="", xlab="", las=1, 
     bty="n", yaxs="i")
axis(2, at=1:2, labels=c("UTIs", "Antibiotics"), las=1)
points(1:31-ruti_pats[ord], rep(1,31), pch=16, cex=umb_participants$UTIs[ord]+1, col="red")
text(1:31-ruti_pats[ord], rep(1,31),ifelse(umb_participants$UTIs[ord]>0,umb_participants$UTIs[ord],""), col="white", font=2)
points(1:31-ruti_pats[ord], rep(2,31), pch=16, cex=named_abx_use[ord]+1, col="orange")
text(1:31-ruti_pats[ord], rep(2,31),ifelse(named_abx_use[ord]>0,named_abx_use[ord],""), col="white", font=2)

dev.off()

#################################################
# Fit models to identify differential abundance #
#################################################

pval_adj <- list()
for (tx in 2:7) {
  n_tax <- length(tax_names[[tx]])
  pval_tab <- matrix(NA,n_tax,6)
  res_tab <- matrix(NA,n_tax,6)
  for (j in 1:n_tax) {
    # patient level mixed effects
    ta <- tax_abundances[[tx]][j,]/100
    if (sum(ta==0)/n_samps<0.9) {
      # arcsin sqrt transformation
      ta_trans <- asin(sqrt(ta))
      # No correction for abx/race
      M1 <- lmer(ta_trans~(1|patient)+ruti_bin, REML=FALSE) 
      pval_tab[j,1] <- Anova(M1)[,3]  
      res_tab[j,1] <- summary(M1)$coef[2,1]
      # adjust for abx
      M2 <- lmer(ta_trans~(1|patient)+ruti_bin+umb_stool$antibiotics, REML=FALSE) 
      pval_tab[j,2] <- Anova(M2)[1,3]
      res_tab[j,2] <- summary(M2)$coef[2,1]
      # adjust for race
      M3 <- lmer(ta_trans~(1|patient)+ruti_bin+race2, REML=FALSE) 
      pval_tab[j,3] <- Anova(M3)[1,3]
      res_tab[j,3] <- summary(M3)$coef[2,1]
      # adjust for abx & race
      M4 <- lmer(ta_trans~(1|patient)+ruti_bin+race2+umb_stool$antibiotics, REML=FALSE) 
      pval_tab[j,4:6] <- Anova(M4)[1:3,3]
      res_tab[j,4:6] <- summary(M4)$coef[2:4,1]  
    }
    summary_matrix <- matrix(0,n_tax,12)
    summary_matrix[,(1:6)*2-1] <- res_tab
    summary_matrix[,(1:6)*2] <- apply(pval_tab,2,p.adjust, method="BH")
  }
  pval_adj[[tx]] <- data.frame(taxa=tax_names[[tx]],summary_matrix)
  names(pval_adj[[tx]])[-1] <- paste(paste0(rep(c("","p_"),6),rep(c("ruti_raw", "abx_adjust", "race_adjust", 
                                                                    "ruti_race_abx_adjust", "race_race_abx_adjust", "abx_race_abx_adjust"),each=2)), sep="_")
}

#####################################
# FIGURE 2B: Differential abundance #
#####################################
# List of ancestors to display in plot
root <- c("Bacteroidetes", "Firmicutes", "Methanobacteriaceae", "Verrucomicrobiaceae", "Desulfovibrionaceae", "Sutterellaceae")
root <- tax_names[[5]]
taxon_levels <- c("kingdom","phylum","class","order","family","genus","species")
taxon_abbr <- substr(taxon_levels,1,1)
pval_lim <- 0.5
plot(NULL, ylim=c(5E-3,100), xlim=c(0,0.4), bty="l", las=1, log="y", xlab="FDR", 
     ylab="Mean relative abundance (%)", yaxt="n")
axis(2, at=10^(-3:2), labels=10^(-3:2), las=1)
for (k in 1:length(root)) {
  root_level <- which(sapply(strsplit(grep(root[k], taxIDs, value=TRUE)[1], "\\|"),
                             breakup,sep="__",elements=2)==root[k])[1]
  txs <- root[k]
  lvs <- root_level
  for (v in (root_level+1):7) {
    txs <- c(txs, unique(breakup(grep(root[k],taxIDs,value=TRUE), c("\\|","__"), c(v,2), seq=TRUE)))
    lvs <- c(lvs, rep(v,length(unique(breakup(grep(root[k],taxIDs,value=TRUE), c("\\|","__"), c(v,2), seq=TRUE)))))
  }
  # Plot links
  for (i in 2:length(txs)) {
    txra1 <- tax_abundances[[lvs[i]]][which(tax_names[[lvs[[i]]]]==txs[i]),]
    pp1 <- pval_adj[[lvs[i]]][which(pval_adj[[lvs[i]]]$taxa==txs[i]),]
    adj_level1 <- mean(txra1)
    ancest <- breakup(unique(breakup(grep(paste(taxon_abbr[lvs[i]],txs[i],sep="__"),taxIDs, value=TRUE), "\\|", 
                                     1:(lvs[i]-1), join="|")), "__", fromend=1)
    txra2 <- tax_abundances[[lvs[i]-1]][which(tax_names[[lvs[i]-1]]==ancest),]
    pp2 <- pval_adj[[lvs[i]-1]][which(pval_adj[[lvs[i]-1]]$taxa==ancest),]
    adj_level2 <- mean(txra2)
    if (!is.na(pp1$p_ruti_race_abx_adjust) && !is.na(pp2$p_ruti_race_abx_adjust)) {
      if (pp1$p_ruti_race_abx_adjust<pval_lim && pp2$p_ruti_race_abx_adjust<pval_lim) {
        lines(c(pp1$p_ruti_race_abx_adjust,pp2$p_ruti_race_abx_adjust),
              c(adj_level1,adj_level2),col="lightgrey")
      }
    }
  }
  # Plot points
  for (i in 1:length(txs)) {
    txra1 <- tax_abundances[[lvs[i]]][which(tax_names[[lvs[[i]]]]==txs[i]),]
    pp1 <- pval_adj[[lvs[i]]][which(pval_adj[[lvs[i]]]$taxa==txs[i]),]
    adj_level1 <- mean(txra1)
    if (length(pp1$ruti_race_abx_adjust)>0) {
      if (!is.na(pp1$ruti_race_abx_adjust)) {
        points(pp1$p_ruti_race_abx_adjust,adj_level1,
               bg=ifelse(pp1$ruti_race_abx_adjust<0,"black","white"),
               pch=ifelse(pp1$ruti_race_abx_adjust<0,25,24),
               cex=(50*abs(pval_adj[[lvs[i]]]$ruti_race_abx_adjust[which(pval_adj[[lvs[i]]]$taxa==txs[i])]))^0.4)
        text(pp1$p_ruti_race_abx_adjust, adj_level1,sp_short(txs[i]),
             pos=ifelse(pp1$ruti_race_abx_adjust>0|pp1$p_ruti_race_abx_adjust<0.1,sample(c(1,3:4),1),sample(1:4,1)))
      }
    }
  }
}
legend("topright", legend=c("Elevated in rUTI", "Reduced in rUTI"), pch=24:25, pt.bg=c("white", "black"))

##########################################
# Figure 2e: E. coli relative abundances #
##########################################
# E. coli relative abundance
ecoli <- tax_abundances[[7]][grep("Escherichia_coli", tax_names[[7]]),]
# By cohort
ecoli_cohort <- list(ecoli[which(ruti_bin==0)],ecoli[which(ruti_bin==1)])
# Separate out zeros
ecoli_cohort_nonzero <- list(ecoli[which(ruti_bin==0 & ecoli>0)],ecoli[which(ruti_bin==1 & ecoli>0)])
# Colors for cohorts
colz <- c("blue", "red")

plot(NULL, xlim=c(8E-4,30), ylim=c(0,2)+0.5, log="x", las=1, bty="l", xlab="E. coli relative abundance", 
     yaxt="n", ylab="", xaxt="n")
boxer(data=ecoli_cohort_nonzero, x=1:2, width = 0.4, horiz=TRUE, border=colz, lwd=2, flicks = 1)
points(rep(1E-3,2), 1:2, col=colz, cex=sapply(ecoli_cohort,function(x)sum(x==0)/length(x))*10, lwd=2)
text(rep(1E-3,2), 1:2, paste0(sapply(ecoli_cohort,function(x)round(100*sum(x==0)/length(x),0)), "%"), cex=0.7)
points(sapply(1:31,function(x)median(ecoli[which(patient==x)])), ruti_pats+1, col=colz[ruti_pats+1], pch=5)
axis(2, at=1:2, labels=c("Control", "rUTI"), las=1)
axis(1, at=10^(-3:1), labels=paste0(c(0,10^(-2:1)),"%"))
axis(1, at=10^-2.5, labels="//", tick=FALSE, line=-1.5)

###########################################
# Figure 3: Strain dynamics and phylogeny #
###########################################
# Phylogroup to strain mapping
clades <- unique(ecoliclades[,2])
# Colors for clades
cladecols <- rep("grey", length(clades))
cladecols[1:8] <- c("chartreuse4", "lightblue", "red", "yellow", "orange", "purple", "salmon", "lightgreen")
# Abx classes
abx_classes <- names(sort(table(abx_table$abx_class),decreasing = TRUE))
abx_abbr <- c("N","X","B","F","M","S","T")
abx_cols <- c(brewer.pal(7,"Set1"), rep("grey",10))
# Gut above bladder, or vice versa?
flip <- 0
# Sample types
alltypes <- unique(strainge_table$type)
allcats <- unique(strainge_table$cat)
# Point styles for different sample types
pt_type <- c(16,16,18)
pt_type2 <- c(1,1,5)
pt_cex <- c(2.5,2.5,3.2)
# Subgroups of interest (UTI sufferers, all rUTI, all control)
patgroups <- list(utis=which(umb_participants$UTIs>0), ruti=which(ruti_pats==1), control_pats=which(ruti_pats==0))

for (pp in 1:length(patgroups)) {
  # participant IDs considered; ordered by no. UTIs
  ptlist <- patgroups[[pp]][-which(patgroups[[pp]]%in%c(14,16,30))]
  #ptlist <- ptlist[order(umb_participants$UTIs[ptlist])]
  # All strains detected
  alluqstrains <- unique(strainge_table$strain[which(strainge_table$pat%in%ptlist & strainge_table$duplicate==0)])
  # E. coli tree of this strain subset
  ecolitree2 <- drop.tip(ecolitree_umb, ecolitree_umb$tip.label[-which(ecolitree_umb$tip.label%in%alluqstrains)])
  # Codify strains with numeral
  strainidentifier <- diag(sapply(ecoliclades[match(ecolitree2$tip.label, ecoliclades[,1]),2],
                                  function(x)cumsum(ecoliclades[match(ecolitree2$tip.label, ecoliclades[,1]),2]==x)))
  par(mar=c(1,1,0,1),oma=c(3,1,1,1))
  #layout(cbind(c(1,rep(2,5)),c(1,rep(2,5)),c(3,rep(4,5)),c(3,rep(4,5)),c(3,rep(4,5))))
  layout(cbind(c(1,rep(3,5)),c(1,rep(3,5)),c(2,rep(3,5)),c(2,rep(3,5)),c(2,rep(3,5))))
  treeclades <- ecoliclades[match(ecolitree_umb$tip.label,ecoliclades[,1]),2]
  treeclades2 <- ecoliclades[match(ecolitree2$tip.label,ecoliclades[,1]),2]
  treecols <- edge.color(ecolitree_umb, groups=sapply(1:8,function(x)which(match(treeclades,clades)==x)), col=cladecols, bgcol="lightgrey")
  # Plot phylogeny
  plot(ecolitree_umb, direction="d", show.tip.label = FALSE, edge.color=treecols, edge.width=2)
  add.scale.bar("bottomright", length=0.01,lwd = 2, lcol = "black")

  plot(NULL, xlim=c(0,1), ylim=c(0,1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  legend("topleft", legend=sort(unique(treeclades)), pch=rep(15,length(unique(treeclades))), 
         col=cladecols[match(sort(unique(treeclades)),clades)], pt.cex=2, cex=1.2, bty="n", title="E. coli phylogroup")
  
  # Strain dynamics panel
  plot(NULL, xlim=c(-10,400),ylim=c(0,2.5*length(ptlist)+1), yaxt="n", xlab="", xaxs="i",
       ylab="", bty="n", yaxs="i")
  #abline(h=(0:length(ptlist))*2)
  tt <- 0
  for (pt in ptlist) {
    pat_strains <- strainge_table[which(strainge_table$pat==pt & strainge_table$duplicate==0),]
    uqstraintab <- sort(table(pat_strains$strain),decreasing=TRUE)
    uqstrains <- names(uqstraintab)
    uqstrainclades <- pat_strains$clade[match(uqstrains,pat_strains$strain)]
    uqstraincladeno <- diag(sapply(uqstrainclades,function(x)cumsum(uqstrainclades==x)))
    if (flip==1) {
      strainlevels <- tt+0.9*as.numeric(pat_strains$cat=="s")+0.4*as.numeric(pat_strains$cat=="u")
    } else {
      strainlevels <- tt+1.3*as.numeric(pat_strains$cat=="u")
    }
    uq_samples <- unique(pat_strains$samp_name)
    strain_adj <- numeric(length(strainlevels))
    for (i in 1:length(uq_samples)) {
      strainspresent <- pat_strains$strain[which(pat_strains$samp_name==uq_samples[i])]
      strainrank <- match(strainspresent,uqstrains)
      if (min(strainrank)>1 & length(grep("u",uq_samples[i]))==0) {
        strain_adj[which(pat_strains$samp_name==uq_samples[i])] <- rank(strainrank)+1
      } else {
        strain_adj[which(pat_strains$samp_name==uq_samples[i])] <- rank(strainrank)   
      }
    }
    strainlevels <- strainlevels+strain_adj/4
    
    rect(-10,tt+flip,400,tt+1+flip, col=adjustcolor("lightgreen",0.2), border=NA)
    rect(-10,tt+1-flip,400,tt+2-flip, col=adjustcolor("pink",0.4), border=NA)

    dayone <- min(umb_samples$date[which(as.numeric(breakup(umb_samples$ID,"B",2))==pt)])
    B <- abx_table[which(abx_table$patient==pt),]
    if (sum(abx_table$patient==pt)>0) {
      B$adj_day <- difftime(as.Date(B$date_ended,"%m/%d/%y"),dayone,units="days")
      B <- B[which(B$adj_day>0),]
      if (nrow(B)>0) {
        segments(B$adj_day[which(B$adj_day>0)], tt,B$adj_day, tt+2,col=abx_cols[match(B$abx_class,abx_classes)], lty=2)
        text(B$adj_day,rep(tt+0.2+1.6*(1-flip),nrow(B)),abx_abbr[match(B$abx_class,abx_classes)],
             pos=4, font=2, col=abx_cols[match(B$abx_class,abx_classes)])
      }
    }
    segments(difftime(as.Date(uti_table$date[which(uti_table$pat==pt)]),dayone,units="days"), rep(tt,sum(uti_table$pat==pt)),
             difftime(as.Date(uti_table$date[which(uti_table$pat==pt)]),dayone,units="days"), rep(tt+2,sum(uti_table$pat==pt)),
             col="black", lwd=2)
    urinestrains <- unique(pat_strains$strain[which(pat_strains$cat=="u")])
    persistentstrains <- uqstrains[which(uqstraintab>=2)]
    trackstrains <- union(urinestrains,persistentstrains)
    for (i in 1:length(trackstrains)) {
      lines(difftime(as.Date(pat_strains$date[which(pat_strains$strain==trackstrains[i])]),dayone,units="days"),
            strainlevels[which(pat_strains$strain==trackstrains[i])], 
            col=cladecols[which(clades==uqstrainclades[which(uqstrains==trackstrains[i])])])    
    }
    for (i in 1:sum(as.numeric(breakup(umb_stool$ID,"B",2))==pt)) {
      if (sum(pat_strains$samp_name==umb_stool$sample[which(as.numeric(breakup(umb_stool$ID,"B",2))==pt)[i]])==0) {
        points(difftime(as.Date(umb_stool$date[which(as.numeric(breakup(umb_stool$ID,"B",2))==pt)[i]]),dayone, units = "days"),
               tt+0.2+1.6*flip, col="grey", cex=2.5)
      }
    }
    points(difftime(as.Date(pat_strains$date),dayone,units="days"),strainlevels,
           cex=pt_cex[match(pat_strains$cat,allcats)], pch=pt_type[match(pat_strains$cat,allcats)], 
           col=cladecols[match(pat_strains$clade,clades)])
    text(difftime(as.Date(pat_strains$date),dayone,units="days"),strainlevels,
         strainidentifier[match(pat_strains$strain,ecolitree2$tip.label)], font=2, col="white")
    
    B <- umb_samples[intersect(which(as.numeric(breakup(umb_samples$ID,"B",2))==pt),grep("rectal",umb_samples$type)),]
    B <- B[which(!breakup(B$sample,"r",1)%in%pat_strains$samp[which(pat_strains$cat=="r")]),]
    if (nrow(B)>0) {
      points(difftime(B$date,dayone,units="days"),rep(tt+0.25,nrow(B)), pch=5, col="grey", cex=2.5)
    }
    B <- umb_samples[intersect(which(as.numeric(breakup(umb_samples$ID,"B",2))==pt),grep("urine",umb_samples$type)),]
    B <- B[which(!breakup(B$sample,"u",1)%in%pat_strains$samp[which(pat_strains$cat=="u")]),]
    if (nrow(B)>0) {
      points(difftime(B$date,dayone,units="days"),rep(tt+1.25,nrow(B)), pch=1, col="grey", cex=2.5)
    }
    text(-10, tt+1, pt, cex=2.5, col="black")
    rect(-10,tt,400,tt+2)
    tt<- tt+2.5
  }
}

#####################################################
# Figure 4: E. coli species tree & observed strains #
#####################################################

uq_patstrain <- unique(strainGST_stool$patstrain)
uq_patstrain_strain <- sapply(uq_patstrain,function(x)sapply(strsplit(x,"_"),function(y)paste(y[-length(y)],collapse="_")))
uq_patstrain_pat <- breakup(uq_patstrain, "_",fromend=1)
gutstrain_table <- data.frame(strain=uq_patstrain_strain,pat=uq_patstrain_pat,
                              freq=sapply(uq_patstrain,function(x)sum(strainGST_stool$patstrain==x))/sapply(uq_patstrain_pat,function(x)sum(patient==x)),
                              clade=strainGST_stool$clade[match(uq_patstrain_strain,strainGST_stool$strain)])
gutstrain_table$persister <- as.numeric(gutstrain_table$freq>=0.25)
gutstrain_table$persister <- as.numeric(gutstrain_table$freq>=0.25)
gutstrain_table$ruti_pat <- as.numeric(gutstrain_table$pat%in%which(ruti_pats==1))

clades <- c("A", "B1", "B2", "D", "E", "F", "G")
cladecols <- c("chartreuse4", "lightblue", "red", "orange", "purple", "salmon", "lightgreen")
baseclades <- ecoliclades[match(ecolitree_full$tip.label,ecoliclades[,1]),2]
ecolitree_full <- drop.tip(ecolitree_full, ecolitree_full$tip.label[which(!baseclades%in%clades)])
ecolitree_full <- drop.tip(ecolitree_full, c(ecolitree_full$tip.label[c(1,134)]))

layout(t(1:4), widths = c(4,1,1,1.25,1))
par(mar=c(2,0.5,1,0.5),oma=c(3,1,1,1))
treeclades <- ecoliclades[match(ecolitree_full$tip.label,ecoliclades[,1]),2]
treecols <- edge.color(ecolitree_full, groups=sapply(1:8,function(x)which(match(treeclades,clades)==x)), col=cladecols, bgcol="lightgrey")
# phylogeny
plot(ecolitree_full, edge.color=treecols, edge.width=2, show.tip.label = FALSE)
add.scale.bar()

uq_strains <- unique(strainGST_urine$strain[which(strainGST_urine$duplicate==0)])
# Urine samps
plot(NULL, xlim=c(0,1), ylim=c(1,length(ecolitree_full$tip.label)), bty="n", 
     yaxt="n", ylab="", xlab="", xaxt="n")
abline(v=0, col="lightgrey")
for (i in 1:length(uq_strains)) {
  B <- strainGST_urine[which(strainGST_urine$strain==uq_strains[i] & strainGST_urine$duplicate==0),]
  B <- B[order(ruti_pats[B$pat], decreasing = TRUE),]
  urinepts <- unique(B$pat)
  for (j in 1:length(urinepts)) {
    points(j*0.1,which(ecolitree_full$tip.label==uq_strains[i]),col=c("blue","red")[as.numeric(urinepts[j]%in%which(ruti_pats==1))+1],
           pch=ifelse(urinepts[j]%in%uti_table$pat[union(grep(gsub("Esch_coli_","",uq_strains[i]),uti_table$straingst_out),
                                                         grep(gsub("Esch_coli_","",uq_strains[i]),uti_table$straingst_raw))],16,1))
  }
}
legend("right", legend=c("Control urine", "rUTI healthy", "UTI"), col=c("blue", "red", "red"), pch=c(1,1,16))

uq_strains <- unique(gutstrain_table$strain)
colz <- list(c("pink", "red"),c("lightblue","darkblue"))
colz <- list(c("grey", "black"),c("grey", "black"))
# Observation frequency
for (grp in 0:1) {
  plot(NULL, xlim=c(0,4+grp), ylim=c(1,length(ecolitree_full$tip.label)), #bty="n", 
       yaxt="n", ylab="", xlab="", xaxs="i")
  for (i in 1:length(uq_strains)) {
    posn <- which(ecolitree_full$tip.label==uq_strains[i])
    if (length(posn)>0) {
      rect(0,posn-0.5,sum(gutstrain_table$strain==uq_strains[i] & gutstrain_table$ruti_pat==1-grp),posn+0.5, col=colz[[grp+1]][1], border=NA)
      rect(0,posn-0.5,sum(gutstrain_table$strain==uq_strains[i] & gutstrain_table$ruti_pat==1-grp & gutstrain_table$persister==1),
           posn+0.5, col=colz[[grp+1]][2], border=NA)      
    }
  }
}
dev.off()


####################################
# Figure S1a

samptimediffs <- matrix(0,n_samps,n_samps)
for (i in 1:n_samps) {
  samptimediffs[i,] <- difftime(as.Date(umb_stool$date),as.Date(umb_stool$date[i]),units="days")
}
braydists <- list()
braydists_time <- list()
for (i in 1:31) {
  bm <- bray_mat[which(patient==i),which(patient==i)]
  braydists[[i]] <- bm[lower.tri(bm)]
  
  bm_time <- bray_mat[which(patient==i),which(patient==i)]/samptimediffs[which(patient==i),which(patient==i)]
  braydists_time[[i]] <- -bm_time[lower.tri(bm_time)]  
}
braymeans <- sapply(braydists,mean)
ord <- order(ruti_pats,braymeans)
ord <- ord[-which(ord%in%c(14,16,30))]

plot(NULL, xlim=c(0,28)+0.5, ylim=c(0,1), bty="l", xaxt="n", xlab="", 
     ylab="Pairwise Bray Curtis dissimilarity", las=1)
#axis(1,at=1:length(ord),labels = ord, las=2)
for (i in 1:length(ord)) {
  points(rnorm(length(braydists[[ord[i]]]),i,0.1), braydists[[ord[i]]],
         col=adjustcolor(ifelse(studygroup[ord[i]]=="rUTI","red","blue"),0.2))
  boxer(braydists[[ord[i]]],i,width=0.5, border=ifelse(studygroup[ord[i]]=="rUTI","red","blue"),
        lwd=2, col=rgb(1,1,1,0.5))
}
lines(c(0.5,14.5),rep(mean(braymeans[ord[which(ord%in%which(ruti_pats==0))]]),2), lwd=5, col=rgb(0,0,1,0.5))
lines(c(14.5,28.5),rep(mean(braymeans[ord[which(ord%in%which(ruti_pats==1))]]),2), lwd=5, col=rgb(1,0,0,0.5))

braymeans_time <- sapply(braydists_time,mean)
ord <- order(ruti_pats,braymeans_time)
ord <- ord[-which(ord%in%c(14,16,30))]

plot(NULL, xlim=c(0,28)+0.5, ylim=c(0.0005,0.05), log="y", bty="l", xaxt="n", xlab="", 
     ylab="Pairwise Bray Curtis dissimilarity", las=1)
#axis(1,at=1:length(ord),labels = ord, las=2)
for (i in 1:length(ord)) {
  points(rnorm(length(braydists_time[[ord[i]]]),i,0.1), braydists_time[[ord[i]]],
         col=adjustcolor(ifelse(studygroup[ord[i]]=="rUTI","red","blue"),0.2))
  boxer(braydists_time[[ord[i]]],i,width=0.5, border=ifelse(studygroup[ord[i]]=="rUTI","red","blue"),
        lwd=2, col=rgb(1,1,1,0.5))
}


# Figure S1b; S5, S6
metric <- list(utiB14w,uti0w)
wdws <- c("pre-UTI", "time of UTI")
bounds <- cbind(c(-0.2,-35,-0.7),c(0.4,25,0.5))

nms <- c("Relative perturbation", "Change in Chao1", "Change in Shannon Diversity")

for (j in 1:2) {
  uti_samps <- sampIDs[which(metric[[j]]==1)]
  n_events <- sum(metric[[j]],na.rm=TRUE)
  bray_over <- numeric(n_events)
  bray_quantile <- list(n_events)
  bray_healthy <- list(n_events)
  c1_over <- numeric(n_events)
  #dv_over <- numeric(n_events)
  for (i in 1:length(uti_samps)) {
    pt <- patient[which(sampIDs==uti_samps[i])]
    healthy_samps <- which(patient==pt & metric[[j]]!=1)
    bray_dv <- bray_mat[healthy_samps,healthy_samps]
    bray_dv <- bray_dv[which(lower.tri(bray_dv))]
    mn_pt <- mean(bray_dv)
    bray_over[i] <- mean(bray_mat[which(sampIDs==uti_samps[i]),
                                  which(patient==pt & sampIDs!=uti_samps[i])])-mn_pt
    bray_quantile[[i]] <- bray_mat[which(sampIDs==uti_samps[i]),
                                   which(patient==pt & sampIDs!=uti_samps[i])]
    bray_healthy[[i]] <- bray_dv
    c1_over[i] <- c1[which(sampIDs==uti_samps[i])]-mean(c1[which(patient==pt & sampIDs!=uti_samps[i])])
    #dv_over[i] <- dv[which(sampIDs==uti_samps[i])]-mean(dv[which(patient==pt & sampIDs!=uti_samps[i])])
  }
  met_list <- list(BC=bray_over, Chao1=c1_over)#, Shannon=dv_over)
  
  plot(NULL, xlim=c(0,n_events)+0.5, ylim=c(0,1), bty="n", 
       ylab="Bray-Curtis dissimilarity", xlab="", xaxt="n", las=1)
  ord <- order(sapply(bray_healthy,median))
  boxer(bray_healthy[ord],1:n_events, width=0.6, col="grey", border="darkgrey", lwd=4)
  boxer(bray_quantile[ord],1:n_events, width=0.3, col="red", lwd=2)  
  axis(1, at=1:n_events, labels = gsub("UMB","",uti_samps[ord]), las=2)
  legend("topleft", legend=c("UTI vs. healthy", "healthy vs. healthy"), pch=22,
         pt.bg=c("red", "grey"), col=c("black", "darkgrey"), pt.cex=2, pt.lwd=c(1,2))
  
  for (i in 1:2) {
    ord <- order(met_list[[i]])
    plot(NULL, xlim=c(0,n_events)+0.5, ylim=bounds[i,], bty="n", las=1,
         ylab=nms[i], xlab="", xaxt="n", main=paste(nms[i], wdws[j], sep="; "))
    abline(h=0)
    rect(1:n_events-0.3,0,1:n_events+0.3,met_list[[i]][ord],
         col=ifelse(met_list[[i]][ord]<0,"green","red"))
    axis(1, at=1:sum(metric[[j]]), labels=gsub("UMB","",uti_samps[ord]), las=2, lwd=0)
    text(1:sum(metric[[j]]),rep(0,sum(metric[[j]])), gsub("UMB","",uti_samps[ord]), 
         srt=90, pos=ifelse(met_list[[i]][ord]<0,4,2))
  }
}

# Figure S2

ruti_cats <- ruti_pats
ruti_cats <- sapply(1:31,function(x)ifelse(umb_participants$UTIs[x]>0,2,ifelse(ruti_pats[x]==1,1,0)))

layout(rbind(1:2,3:4))
colz <- c("blue", "orange", "red")
for (k in 1:length(scfa_ra)) {
  plot(NULL, xlim=c(1,34), ylim=c(0,quantile(scfa_ra[[k]],1)),bty="l",las=1, xaxt="n", 
       xlab="", ylab="Relative abundance (%)", main=paste(names(scfa_ra)[k], "producers"), yaxs="i")
  #axis(1, at=c(8,26), labels = c("rUTI", "Control"), tick = FALSE)
  pt_lvls <- sapply(1:31,function(x)scfa_ra[[k]][which(patient==x)])
  pt_mns <- sapply(pt_lvls,mean)
  ord <- order(ruti_cats, pt_mns)
  #axis(1, at=1:31+(1-ruti_pats[ord]), labels=ord, las=2)
  boxer(pt_lvls[ord],ruti_cats[ord]+1:length(pt_mns),width = 0.4, col=colz[ruti_cats[ord]+1],
        quants=c(0,0.25,0.5,0.75,1))
  for (j in 0:max(ruti_cats)) {
    lines(range(which(ruti_cats[ord]==j))+j,rep(mean(pt_mns[ord][which(ruti_cats[ord]==j)]),2), lwd=4, 
          col=adjustcolor(colz[j+1],0.5))
  }
  if (k==1) {
    legend("topright", legend=c("Control", "rUTI (0 UTIs)", "rUTI (1+ UTIs)"), pch=22, 
           pt.bg=colz, pt.cex=1.5, lwd=1, bty="n")
  }
}
colz_sp <- c(brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))
#layout(1:2)
for (i in 1:length(scfa_ra)) {
  plot(NULL, xlim=c(0,32)+1, ylim=c(0,quantile(scfa_ra[[i]],1)),bty="l",las=1, xaxt="n", 
       xlab="", ylab="Relative abundance (%)", main=paste(names(scfa_ra)[i], "producers"), yaxs="i")
  pt_lvls <- sapply(1:31,function(x)scfa_ra[[i]][which(patient==x)])
  pt_mns <- sapply(pt_lvls,mean)
  ord <- order(ruti_cats, pt_mns)
  axis(1, at=1:31+ruti_cats[ord], labels=ord, las=2)
  for (j in 1:31) {
    cml <- 0
    for (k in 1:length(scfa_producers[[i]])) {
      taxadd <- mean(tax_abundances[[7]][scfa_idx[[i]][k],which(patient==ord[j])])
      rect(ruti_cats[ord[j]]+j-0.3,cml,ruti_cats[ord[j]]+j+0.3,cml+taxadd, col=colz_sp[which(allscfa==scfa_producers[[i]][k])])
      cml <- cml+taxadd
    }
  }
  for (j in 0:max(ruti_cats)) {
    lines(range(which(ruti_cats[ord]==j))+j,rep(mean(pt_mns[ord][which(ruti_cats[ord]==j)]),2), lwd=4, 
          col=adjustcolor(colz[j+1],0.5))
  }
  legend("top", legend=sp_short(scfa_producers[[i]]), ncol=3, pch=22, cex=0.8, 
         pt.cex=1.6, pt.bg=colz_sp[match(scfa_producers[[i]],allscfa)])
}


ord <- order(ruti_cats,sapply(1:31,function(x)mean(c1[which(patient==x)])))
colz <- c("blue", "orange", "red")

plot(NULL, xlim=c(0,31)+0.5, ylim=c(0,80), ylab="Microbial richness", 
     las=1, bty="l", xlab="Participant ID", xaxt="n")
axis(1, at=1:31+ruti_cats[ord], labels=ord, las=2)
boxer(lapply(1:31,function(x)c1[which(patient==x)])[ord],1:31+ruti_cats[ord], 
      width = 0.3, border = colz[ruti_cats[ord]+1], lwd=3)
for (j in 0:max(ruti_cats)) {
  lines(range(which(ruti_cats[ord]==j))+j,rep(mean(sapply(which(ruti_cats==j),function(x)mean(c1[which(patient==x)]))),2), lwd=3,
        col=adjustcolor(colz[j+1],0.5))
}
legend("topright", legend=c("Control", "rUTI (0 UTIs)", "rUTI (1+ UTIs)"), pch=0, 
       col=colz, pt.cex=1.5, lwd=2, bty="n")







#Figure S5

ecoli <- tax_abundances[[7]][grep("Escherichia_coli",tax_names[[7]]),]
ecoli[which(ecoli==0)] <- 1E-3

pt_nonuti_mean <- numeric(31)
pt_nonuti_median <- numeric(31)

timept <- uti0w_ecoli
for (i in 1:31) {
  pt_nonuti_mean[i] <- mean(ecoli[which(patient==i & timept==0)])
  pt_nonuti_median[i] <- median(ecoli[which(patient==i & timept==0)])  
}
target_samps <- which(timept==1)

medianchange <- numeric(sum(timept))
prevchange <- numeric(sum(timept))
for (i in 1:sum(timept)) {
  medianchange[i] <- log(ecoli[target_samps[i]])-log(pt_nonuti_median[patient[target_samps[i]]])
  if (patient[target_samps[i]-1]==patient[target_samps[i]]) {
    prevchange[i] <- log(ecoli[target_samps[i]])-log(ecoli[target_samps[i]-1])
  } else {
    prevchange[i] <- NA
  }
}
ecoli_changes <- list(median=medianchange, previous=prevchange)

layout(1:2)
ord <- order(ecoli_changes[[1]])
for (k in 1:2) {
  #ord <- order(ecoli_changes[[k]])
  plot(NULL, xlim=c(0,sum(timept))+0.5, ylim=c(-8,8), bty="l", las=1, xaxt="n",# yaxt="n",
       ylab=paste0("E. coli log fold change"),xlab="")# (vs. ", names(ecoli_changes)[k],  ")"), xlab="")
  abline(h=0, lwd=2)
  #axis(2,at=log(c(10^c(-3:3))),labels=10^c(-3:3), las=1)
  axis(1, at=1:sum(timept), labels=gsub("UMB","",sampIDs[target_samps][ord]), las=2)
  for (i in 1:sum(timept)) {
    rect(i-0.3,0,i+0.3,ecoli_changes[[k]][ord[i]], lwd=2,
         col=ifelse(ecoli_changes[[k]][ord[i]]<0,"green","red"))
    if (is.na(ecoli_changes[[k]][ord[i]])) {
      text(i,0,"X",pos=3, font=2)
    }
  }
}

# Figure S8a

plot(NULL, xlim=c(0,2)+0.5, ylim=c(0,105), las=1, bty="l", xlab="",
     ylab="Samples with E. coli strains detected (%)", yaxs="i", xaxt="n")
axis(1, at=1:4, labels=c("rUTI stool", "Control stool", "Urine (raw)", "Urine (outgrowth)"), tick=FALSE, line=2)
axis(1, at=1:4-0.375+rep((0:2+0.5)/4,4), labels=paste0(rep(1:3,4),"+"))
ruti_samples <- sampIDs[which(studygroup[patient]=="rUTI")]
for (i in 1:3) {
  rect(1-0.375+(i-1)/4,0,1-0.375+i/4,100*sum(strainGST_stool$pat%in%which(ruti_pats==1) & strainGST_stool$i==i-1)/length(ruti_samples), col="red")
  text(1-0.375+(i-0.5)/4,100*sum(strainGST_stool$pat%in%which(ruti_pats==1) & strainGST_stool$i==i-1)/length(ruti_samples),
       sum(strainGST_stool$pat%in%which(ruti_pats==1) & strainGST_stool$i==i-1), pos=3)
}
text(1,105,paste0("n=",length(ruti_samples)),pos=1)
for (i in 1:3) {
  rect(2-0.375+(i-1)/4,0,2-0.375+i/4,100*sum(strainGST_stool$pat%in%which(ruti_pats==0) & strainGST_stool$i==i-1)/(n_samps-length(ruti_samples)), col="blue")
  text(2-0.375+(i-0.5)/4,100*sum(strainGST_stool$pat%in%which(ruti_pats==0) & strainGST_stool$i==i-1)/(n_samps-length(ruti_samples)),
       sum(strainGST_stool$pat%in%which(ruti_pats==0) & strainGST_stool$i==i-1), pos=3)
}
text(2,105,paste0("n=",n_samps-length(ruti_samples)),pos=1)

# Figure S8b

strainGST_stool$patstrainclade <- paste(strainGST_stool$pat, strainGST_stool$strain, strainGST_stool$clade, sep="_")
ruti_strains <- unique(strainGST_stool$patstrainclade[which(strainGST_stool$pat%in%which(ruti_pats==1))])
control_strains <- unique(strainGST_stool$patstrainclade[which(strainGST_stool$pat%in%which(ruti_pats==0))])

clds <- c("A", "B1", "B2", "D", "E", "F", "G")
plot(NULL, xlim=c(0,7)+0.5, ylim=c(0,0.5), las=1, bty="l", xlab="Phylogroup", 
     ylab="Proportion of unique strains", xaxt="n", yaxs="i")
axis(1, at=1:length(clds), labels=clds)
rect(1:length(clds),rep(0,length(clds)),
     1:length(clds)+0.3,sapply(clds,function(x)sum(breakup(control_strains,"_",fromend=1)==x)/length(control_strains)), col="blue")
rect(1:length(clds)-0.3,rep(0,length(clds)),
     1:length(clds),sapply(clds,function(x)sum(breakup(ruti_strains,"_",fromend=1)==x)/length(ruti_strains)), col="red")
legend("topright", legend=c("rUTI", "Control"), pch=22, pt.bg = c("red", "blue"), pt.cex=2)

