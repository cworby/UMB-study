library(phyloch)
library(ape)
library(vegan)
library(RColorBrewer)

source(paste0(homedir,"code/umb_functions.r"))
datadir <- paste0(homedir, "data/")
outdir <- paste0(homedir, "outputs/")

umb_samples <- read.table(paste0(datadir, "umb_samples.tsv"), header=TRUE, sep="\t")
umb_stool <- read.table(paste0(datadir, "umb_stool.tsv"), header=TRUE, sep="\t")
umb_participants <- read.table(paste0(datadir, "umb_participants.tsv"), header=TRUE, sep="\t")
uti_table <- read.table(paste0(datadir, "umb_utis.tsv"), header=TRUE, sep="\t")
abx_table <- read.table(paste0(datadir, "abx_table.tsv"), header=TRUE, sep="\t")
abundances_mp <- read.table(paste0(outdir, "umb_metaphlan2.tsv"), header=TRUE, sep="\t")
strainge_table <- read.table(paste0(outdir, "strainge_table.tsv"), header=TRUE, sep="\t")
ecoliclades <- read.table(paste0(outdir, "ecoli_strain_clades.txt"), sep="\t", header=TRUE,
                          stringsAsFactors = FALSE)
ecolitree_umb <- read.tree(paste0(outdir, "ecoli.scc.umb_only.nwk"))
ecolitree_full <- read.tree(paste0(outdir, "ecoli.scc.midpointroot.nwk"))

patientID <- umb_participants$ID
patient <- as.numeric(breakup(umb_samples$ID[which(umb_samples$type=="stool")],"B",2))
studygroup <- umb_participants$studygroup
ruti_pats <- as.numeric(studygroup=="rUTI")
ruti_bin <- as.numeric(studygroup[patient]=="rUTI")
race <- umb_participants$race[patient]
race2 <- race
race2[which(race2==3)] <- 0
race2[which(race2==2)] <- 1
sampIDs <- umb_stool$sample
n_samps <- nrow(umb_stool)

taxIDs <- abundances_mp[,1]
# taxonomic levels for MP output
tax_abundances <- list()
tax_names <- list()
taxon_levels <- c("kingdom","phylum","class","order","family","genus","species")
for (i in 2:7) {
  txnames <- breakup(taxIDs, c("\\|","__"), c(i,2), seq=TRUE)
  uq_txnames <- unique(txnames)
  txab <- NULL
  for (j in 1:length(uq_txnames)) {
    txab <- rbind(txab, apply(abundances_mp[which(txnames==uq_txnames[j]),-1],2,sum))
  }
  ord <- order(apply(txab,1,sum), decreasing=TRUE)
  tax_abundances[[i]] <- txab[ord,]
  tax_names[[i]] <- uq_txnames[ord]
}

bray_mat <- as.matrix(vegdist(t(tax_abundances[[7]])), method="bray")

butyrate_producers <- c("Faecalibacterium_prausnitzii", "Subdoligranulum_variabile", "Eubacterium_biforme", "Eubacterium_rectale", 
                        "Roseburia_inulinivorans", "Roseburia_intestinalis", "Eubacterium_hallii", "Anaerostipes_hadrus", 
                        "Coprococcus_eutactus", "Coprococcus_catus")
propionate_producers <- c("Bacteroides_uniformis", "Bacteroides_vulgatus", "Prevotella_copri", "Alistipes_putredinis", 
                          "Roseburia_inulinivorans", "Eubacterium_hallii", "Ruminococcus_obeum", "Coprococcus_catus", 
                          "Dialister_invisus", "Phascolarctobacterium_succinatutens", "Akkermansia_muciniphila")
scfa_producers <- list(Butyrate=butyrate_producers, Propionate=propionate_producers)#, Propionate_Pdu=propionate_producers_pdu)

allscfa <- unique(unlist(scfa_producers))
txord <- order(sapply(allscfa,function(x)mean(tax_abundances[[7]][which(tax_names[[7]]==x),])), decreasing = TRUE)
allscfa <- allscfa[txord]
scfa_producers <- lapply(scfa_producers, function(x)x[order(match(x,allscfa))])

scfa_idx <- lapply(scfa_producers,function(x)sapply(x,function(y)grep(y,tax_names[[7]])))
scfa_ra <- lapply(scfa_idx,function(x)apply(tax_abundances[[7]][x,],2,sum))

c1 <- umb_samples$c1[which(umb_samples$type=="stool")]


utiB14 <- numeric(n_samps)
uti0w <- numeric(n_samps)
uti14 <- numeric(n_samps)
utiB14_ecoli <- numeric(n_samps)
uti0w_ecoli <- numeric(n_samps)
uti14_ecoli <- numeric(n_samps)
utiB14_diag <- numeric(n_samps)
uti0w_diag <- numeric(n_samps)
uti14_diag <- numeric(n_samps)

for (i in 1:n_samps) {
  if (patient[i]%in%uti_table$pat) {
    stool_date <- as.Date(umb_stool$date[i])
    uti_dates <- as.Date(uti_table$date[which(uti_table$pat==patient[i])])
    t_diff <- as.numeric(difftime(stool_date,uti_dates, units="days"))
    if (sum((-14:0)%in%t_diff)>0) {
      utiB14[i] <- 1
      bt <- which(t_diff%in%(-14:0))
      if (uti_table$ecoli_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        utiB14_ecoli[i] <- 1
      }
      if (uti_table$diagnosed_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        utiB14_diag[i] <- 1
      }
    }
    if (sum((-3:3)%in%t_diff)>0) {
      uti0w[i] <- 1
      bt <- which(t_diff%in%(-3:3))
      if (uti_table$ecoli_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        uti0w_ecoli[i] <- 1
      }
      if (uti_table$diagnosed_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        uti0w_diag[i] <- 1
      }
    }
    if (sum((0:14)%in%t_diff)>0) {
      uti14[i] <- 1
      bt <- which(t_diff%in%(0:14))
      if (uti_table$ecoli_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        uti14_ecoli[i] <- 1
      }
      if (uti_table$diagnosed_utis[which(uti_table$pat==patient[i])[bt]]==1) {
        uti14_diag[i] <- 1
      }
    }
  }
}
utiB14w <- as.numeric(utiB14==1 & uti0w==0)
uti14w <- as.numeric(uti14==1 & uti0w==0)
utiB14w_ecoli <- as.numeric(utiB14_ecoli==1 & uti0w_ecoli==0)
uti14w_ecoli <- as.numeric(uti14_ecoli==1 & uti0w_ecoli==0)
utiB14w_diag <- as.numeric(utiB14_diag==1 & uti0w_ecoli==0)
uti14w_diag <- as.numeric(uti14_diag==1 & uti0w_diag==0)

strainge_table$patstrain <- paste(strainge_table$strain, strainge_table$pat, sep="_")
strainGST_stool <- strainge_table[which(strainge_table$type=="stool"),]
strainGST_urine <- strainge_table[which(strainge_table$cat=="u"),]

# Number of anitbiotic courses reported (& named)
named_abx_use <- sapply(1:31,function(x)sum(abx_table$abx_class[which(abx_table$patient==x)]!="n/a"))

# Chao1 microbial richness
c1 <- umb_stool$c1

  