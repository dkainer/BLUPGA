library(BLUPGA)
library(dplyr)
library(magrittr)
library(cpgen)
library(SNPRelate)
library(rrBLUP)
library(foreach)
library(doSNOW)


############  Load Pheno  ##########################################################################
setwd("C:/Users/u5309568/OneDrive/DNA Analysis/")

missing   =  c("DK028", "DK044", "DK045", "DK055", "DK065", "DK066", "DK115", "DK204", "DK228", "DK364", "DK194", "DK195")

FT2       <- read.table('FT2_GWAS_pheno.csv', header=TRUE, sep=",", stringsAsFactors = F)
FT2.uniq  <- filter(FT2, ! ID %in% missing )
FT2.uniq$logapinene <- log(FT2.uniq$aPinene)
FT2.uniq$logapinene[is.na(FT2.uniq$logapinene)] <- min(FT2.uniq$logapinene, na.rm=T)
FT2.uniq$Family <- as.factor(FT2.uniq$Family)
FT2.uniq$postblock <- as.factor(FT2.uniq$postblock)
FT2.uniq$FamOrder <- as.factor(FT2.uniq$FamOrder)
FT2.uniq$Prov <- as.factor(FT2.uniq$Prov)
FT2.uniq$Terrain <- as.factor(FT2.uniq$Terrain)
FT2.uniq$Leafiness <- as.factor(FT2.uniq$Leafiness)
blockmap <- c("A"=1, "B"=2, "C"=3, "D"=4)
FT2.uniq$postblock <- blockmap[as.character(FT2.uniq$postblock)]
FT2.uniq$OC.rn <- GenABEL::rntransform(FT2.uniq$OilConc.adj)
FT2.uniq$sqrtSesqui  <- sqrt(FT2.uniq$Sesqui.adj)
FT2.uniq$logLimonene  <- log(FT2.uniq$limonene)
FT2.uniq$logmsratio  <- log(FT2.uniq$MSratio)
FT2.uniq$rnPercCin <- GenABEL::rntransform(FT2.uniq$percCin.adj)

FT1     <- read.table('FT1_GWAS_pheno.txt', header=TRUE, sep="\t", stringsAsFactors = F)
FT1.uniq  <- filter(FT1, ! ID %in% missing )
FT1.uniq$logapinene <- log(FT1.uniq$a.pinene)
FT1.uniq$logapinene[is.na(FT1.uniq$logapinene)] <- min(FT1.uniq$logapinene, na.rm=T)
FT1.uniq %<>% mutate(postblock = FT2.uniq$postblock, FamOrder = FT2.uniq$FamOrder, Prov = FT2.uniq$Prov, Terrain = as.character(FT2.uniq$Terrain))

FT2.uniq$GRW.CA      <- FT2.uniq$DiaDia/FT1.uniq$DiaDia
FT2.uniq$dCA      <- FT2.uniq$DiaDia - FT1.uniq$DiaDia

### Fixed Effects (postblock) design matrix
#cuts          <- cut(as.numeric(FT2.uniq$FamOrder), breaks=c(0,8,16,24,32,40), labels=seq(1,5) )
#FT2.uniq$cuts <- cut(as.numeric(FT2.uniq$FamOrder), breaks=c(0,8,16,24,32,40), labels=seq(1,5) )
X3            <- model.matrix(ID ~ 1 + as.factor(postblock), data=FT2.uniq)

# constants
A.ped     <- make_ped_matrix(FT2.uniq)
oilCandidates <- c("Eucgr.J00973","Eucgr.J00970","Eucgr.J00971","Eucgr.J00972","Eucgr.J00977","Eucgr.J00976",
                   "Eucgr.B01785","Eucgr.G02301","Eucgr.G00843","Eucgr.B03157","Eucgr.E00317","Eucgr.C03570",
                   "Eucgr.F00384","Eucgr.J01589","Eucgr.D00731","Eucgr.K03541","Eucgr.B03623","Eucgr.E00840",
                   "Eucgr.F02616","Eucgr.J02222","Eucgr.J01534","Eucgr.F03370","Eucgr.E01295","Eucgr.C01567",
                   "Eucgr.E02010","Eucgr.G01944","Eucgr.TPS058","Eucgr.A01085","Eucgr.K00218","Eucgr.F00203",
                   "Eucgr.A01301","Eucgr.E02010","Eucgr.A01108","Eucgr.G01947","Eucgr.C00724","Eucgr.K03539","Eucgr.F02313",
                   "Eucgr.B00860","Eucgr.A00917","Eucgr.D00732","Eucgr.H00220","Eucgr.B03598","Eucgr.F03295",
                   "Eucgr.F03296","EgranTPS023","EgranTPS061","EgranTPS062","EgranTPS063","Eucgr.D00872",
                   "EgranTPS076","Eucgr.K00877","EgranTPS091","Eucgr.D00677","Eucgr.E03598","Eucgr.E03995",
                   "Eucgr.E02451","Eucgr.H01236","Eucgr.H01235","Eucgr.I02200","Eucgr.F01818","Eucgr.F03606",
                   "Eucgr.J00991","Eucgr.C01793","Eucgr.B00316","Eucgr.B02495","Eucgr.C00644","Eucgr.D00872")

annot.GENES <- read.table("C:/devwork/NGS/Epoly/WGS/GS/SSGP/thunder.filtered.whole.culled.MAF02.VIF50.468.annotGENES", sep=' ', header = T, stringsAsFactors = F )
annot.GENES$CAND[ annot.GENES$ANNOT %in% oilCandidates  ] <- "OIL"
featnote <- read.table("c:/devwork/NGS/Epoly/WGS/gotcloud/thunder.filtered.whole.culled.MAF02.VIF50.468.FEATnote.consolidated.txt", sep='\t', header = T, stringsAsFactors = F)
featnote %<>% mutate( SNPID = paste(Chromosome,Position, sep='_') )
setwd('c:/devwork/NGS/Epoly/WGS/gotcloud/')
whole.gds <- snpgdsOpen("thunder.filtered.whole.culled.MAF02.VIF50.468.gds", readonly=T, allow.fork=T)


pheno_list <- c("OilConc.adj","Mono.adj","sqrtSesqui","cineole","percCin.adj","logapinene","Height","GRW.Height")
pheno_abbr <- c("OC","MONO","SESQ","CIN","PCIN","APIN","HT","dHT")


#===============  make random cross-val sets ===============================

for(i in 1:100)
{
  write.table( paste(sample(1:nrow(pheno), 70), sep='\n'),
             paste("c:/devwork/NGS/Epoly/WGS/GS/GBLUP/valsets/RANDOM.samples",6, sep = '_'),
             quote=F, row.names = F, col.names = F)
}


#================================= 10K =====================================

# set up full GRM
setwd("C:/devwork/NGS/Epoly/WGS/GS/GBLUP/")
pruned.id <- read.table("10k.MAF0.05.nested.snpid", header = F, stringsAsFactors = F)
pruned.id <- pruned.id$V1
geno.012      <- snpgdsGetGeno(whole.gds, snpfirstdim = F, snp.id = unlist(pruned.id), with.id = T)  # get 100,000 random SNPs
geno.012$genotype  <- geno.012$genotype - 1
colnames(geno.012$genotype) <- geno.012$snp.id
rownames(geno.012$genotype) <- geno.012$sample.id
G1 <- make_GRM(geno.012$genotype)
gc()

# set up CAND GRMs
GRMs <- annot.GENES %>% filter(CAND %in% c('OIL'), SNPID %in% pruned.id) %>%
  group_by(ANNOT) %>%
  do(grm = single_gene_GRM(., pruned.id))


# set up POS GRM
pos     <- which( geno.012$snp.id %in% filter(featnote, SNP_Type=='NONSENSE'| SNP_Type=='MISSENSE')$SNPID )
W       <- rep(1,length(pos))
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='NONSENSE')$SNPID) ]  <- 100
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='MISSENSE')$SNPID) ]  <- 10
W       <- W * (length(W)) / sum((W))
pos.genomat   <- geno.012$genotype[,pos]
S             <- make_weighted_GRM(pos.genomat, W)


# run stuff!
df <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                 VAL=character(), W=numeric(), COR=numeric())
for (i in seq_along(pheno_list))
#for (i in 8:8)
{
  pheno <- select_(FT2.uniq, ID="ID", FE="postblock", y=pheno_list[i])
  df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='10K', genomat=geno.012$genotype, TRAIT=pheno_abbr[i], REL='RANDOM', start=1, end=100))
  #df <- rbind(df, run_one_trait_EFF(pheno, A.ped, G1, S, NSNP='10K', genomat=geno.012$genotype, TRAIT=pheno_abbr[i], FE.mat=cbind(X3,myPCAir$Q1), REL='RANDOM', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='10K', TRAIT=pheno_abbr[i], REL='REL', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='10K', TRAIT=pheno_abbr[i], REL='UNREL', start=1, end=100))
}

#================================= 100K ========================================

# set up full GRM
setwd("C:/devwork/NGS/Epoly/WGS/GS/GBLUP/")
pruned.id <- read.table("100k.MAF0.05.nested.snpid", header = F, stringsAsFactors = F)
pruned.id <- pruned.id$V1
geno.012      <- snpgdsGetGeno(whole.gds, snpfirstdim = F, snp.id = unlist(pruned.id), with.id = T)
geno.012$genotype  <- geno.012$genotype - 1
colnames(geno.012$genotype) <- geno.012$snp.id
rownames(geno.012$genotype) <- geno.012$sample.id
G1 <- make_GRM(geno.012$genotype)
gc()

# set up CAND GRMs
GRMs <- annot.GENES %>% filter(CAND %in% c('OIL'), SNPID %in% pruned.id) %>%
  group_by(ANNOT) %>%
  do(grm = single_gene_GRM(., pruned.id))


# set up POS GRM
pos     <- which( geno.012$snp.id %in% filter(featnote, SNP_Type=='NONSENSE'| SNP_Type=='MISSENSE')$SNPID )
W       <- rep(1,length(pos))
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='NONSENSE')$SNPID) ]  <- 100
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='MISSENSE')$SNPID) ]  <- 10
W       <- W * (length(W)) / sum((W))
pos.genomat   <- geno.012$genotype[,pos]
S             <- make_weighted_GRM(pos.genomat, W)


# run stuff!
df <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                 VAL=character(), W=numeric(), COR=numeric())
for (i in seq_along(pheno_list))
{
  pheno <- select_(FT2.uniq, ID="ID", FE="postblock", y=pheno_list[i])
  df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='100K', genomat=geno.012$genotype, TRAIT=pheno_abbr[i], REL='RANDOM', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='100K',TRAIT=pheno_abbr[i], REL='REL', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='100K',TRAIT=pheno_abbr[i], REL='UNREL', start=1, end=100))
}


#================================= 500K ========================================

# set up full GRM
setwd("C:/devwork/NGS/Epoly/WGS/GS/GBLUP/")
pruned.id <- read.table("500k.MAF0.05.nested.snpid", header = F, stringsAsFactors = F)
pruned.id <- pruned.id$V1
geno.012      <- snpgdsGetGeno(whole.gds, snpfirstdim = F, snp.id = unlist(pruned.id), with.id = T)
geno.012$genotype  <- geno.012$genotype - 1
colnames(geno.012$genotype) <- geno.012$snp.id
rownames(geno.012$genotype) <- geno.012$sample.id
G1 <- make_GRM(geno.012$genotype)
gc()

# set up CAND GRMs
GRMs <- annot.GENES %>% filter(CAND %in% c('OIL'), SNPID %in% pruned.id) %>%
  group_by(ANNOT) %>%
  do(grm = single_gene_GRM(., pruned.id))


# set up POS GRM
pos     <- which( geno.012$snp.id %in% filter(featnote, SNP_Type=='NONSENSE'| SNP_Type=='MISSENSE')$SNPID )
W       <- rep(1,length(pos))
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='NONSENSE')$SNPID) ]  <- 100
W[ which( geno.012$snp.id[pos] %in% filter(featnote, SNP_Type=='MISSENSE')$SNPID) ]  <- 10
W       <- W * (length(W)) / sum((W))
pos.genomat   <- geno.012$genotype[,pos]
S             <- make_weighted_GRM(pos.genomat, W)


# run stuff!
df <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                 VAL=character(), W=numeric(), COR=numeric())
for (i in seq_along(pheno_list))
{
  pheno <- select_(FT2.uniq, ID="ID", FE="postblock", y=pheno_list[i])
  df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='500K',genomat=geno.012$genotype, TRAIT=pheno_abbr[i], REL='RANDOM', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='500K',TRAIT=pheno_abbr[i], REL='REL', start=1, end=100))
  #df <- rbind(df, run_one_trait(pheno, A.ped, G1, S, NSNP='500K',TRAIT=pheno_abbr[i], REL='UNREL', start=1, end=100))
}


#======================================== 100K RAND =======================
# set up full GRM
setwd("C:/devwork/NGS/Epoly/WGS/GS/GBLUP/")
pruned.id <- read.table("100k.MAF0.05.nested.snpid", header = F, stringsAsFactors = F)
pruned.id <- pruned.id$V1
geno.012      <- snpgdsGetGeno(whole.gds, snpfirstdim = F, snp.id = unlist(pruned.id), with.id = T)
geno.012$genotype  <- geno.012$genotype - 1
colnames(geno.012$genotype) <- geno.012$snp.id
rownames(geno.012$genotype) <- geno.012$sample.id
G1 <- make_GRM(geno.012$genotype)
gc()

# run 100 different sets of random genes in for each cross-val
df <- data.frame(GENESET=numeric(), TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                 VAL=character(), W=numeric(), COR=numeric())

system.time( for (geneset in 1:100)
{
  # set up CAND GRMs
  annot.GENES$CAND <- NA
  rand <- sample(levels(as.factor(annot.GENES$ANNOT)), size = 110, replace = F)
  annot.GENES$CAND[ annot.GENES$ANNOT %in% rand ] <- "RAND"

  cat("GENESET",geneset," containing", count(annot.GENES, CAND=='RAND')$n[[1]],"SNPs\n============================\n")
  GRMs <- annot.GENES %>% filter(CAND %in% c('RAND'), SNPID %in% pruned.id) %>%
    group_by(ANNOT) %>%
    do(grm = single_gene_GRM(., pruned.id))

  # run stuff in parallel!

  pheno <- select_(FT2.uniq, ID="ID", FE="postblock", y=pheno_list[1])
  df <- rbind(df, run_one_trait_RAND_para(geneset, pheno, A.ped, G1, S, GRMs, NSNP='100K',TRAIT=pheno_abbr[1], REL='RANDOM', start=1, end=100))
  pheno <- select_(FT2.uniq, ID="ID", FE="postblock", y=pheno_list[7])
  df <- rbind(df, run_one_trait_RAND_para(geneset, pheno, A.ped, G1, S, GRMs, NSNP='100K',TRAIT=pheno_abbr[7], REL='RANDOM', start=1, end=100))
  #df <- rbind(df, run_one_trait_RAND_para(geneset, pheno, A.ped, G1, S, GRMs, NSNP='100K',TRAIT=pheno_abbr[7], REL='REL', start=1, end=100))
  #df <- rbind(df, run_one_trait_RAND_para(geneset, pheno, A.ped, G1, S, GRMs, NSNP='100K',TRAIT=pheno_abbr[7], REL='UNREL', start=1, end=100))
} )


run_one_trait_RAND_para <- function(GENESET, pheno, A.ped, G1, S, GRMs, NSNP, TRAIT, REL, start=1, end=100)
{
  out <- data.frame(GENESET=numeric(), TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                    VAL=character(), W=numeric(), COR=numeric())

  cl<-makeCluster(7, type="SOCK", outfile="")
  registerDoSNOW(cl)

  out <- foreach(run = start:end, .packages=c("rrBLUP","BLUPGA","dplyr","magrittr"),
                 .verbose=TRUE, .combine=rbind, .inorder = FALSE) %dopar%
                 {
                   cat("\nRUN",run," for ",TRAIT," with ",REL," validation \n")
                   cat("==============================================================\n")
                   if(REL=='REL')        val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
                   else if(REL=='UNREL') val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
                   else                  val <- read.table(paste("valsets/RANDOM.samples",run, sep="_"), header = F)
                   val <- val$V1

                   # BLUP|GA with weighted candidate gene GRMs
                   results <-  blupga_CAND(G1, pheno, val, GRMs)
                   results %<>% mutate(GENESET=GENESET, TRAIT=TRAIT, MODEL='RAND', RUN=run, NSNP=NSNP, VAL=REL)
                   results
                 }
  stopCluster(cl)
  closeAllConnections()
  return(out)
}


run_one_trait_RAND <- function(GENESET, pheno, A.ped, G1, S, GRMs, NSNP, TRAIT, REL, start=1, end=100)
{
  out <- data.frame(GENESET=numeric(), TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                    VAL=character(), W=numeric(), COR=numeric())

  for (run in start:end)
  {
     cat("\nRUN",run," for ",TRAIT," \n")
     cat("=================================\n")
     if(REL=='REL')  val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
     else            val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
     val <- val$V1

     # BLUP|GA with weighted candidate gene GRMs
     results <-  blupga_CAND(G1, pheno, val, GRMs)
     results %<>% mutate(GENESET=GENESET, TRAIT=TRAIT, MODEL='RAND', RUN=run, NSNP=NSNP, VAL=REL)
     out <- rbind(out, results)
   }
  return(out)
}

run_one_trait_para <- function(pheno, A.ped, G1, S, genomat, NSNP, TRAIT, REL, start=1, end=100)
{

  out <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                    VAL=character(), W=numeric(), COR=numeric())
  # make design matrix for fixed effs. Only done for ABLUP
  FE.mat  <- model.matrix(ID ~ 1 + as.factor(FE), data=pheno)

  cl<-makeCluster(7, type="SOCK", outfile="")
  registerDoSNOW(cl)

  out <- foreach(run = start:end, .packages=c("rrBLUP","BLUPGA","dplyr","magrittr"),
                 .verbose=TRUE, .combine=rbind, .inorder = FALSE) %dopar%
                 {
                   cat("\nRUN",run," for ",TRAIT," with ",REL," validation \n")
                   cat("==============================================================\n")
                   if(REL=='REL')        val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
                   else if(REL=='UNREL') val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
                   else                  val <- read.table(paste("valsets/RANDOM.samples",run, sep="_"), header = F)
                   val <- val$V1

                   # standard ABLUP
                   cat("ABLUP\n")
                   pheno.blup        <- pheno
                   pheno.blup$y[ val ] <- NA

                   ABLUP   <- mixed.solve(pheno.blup$y, K=A.ped, X=FE.mat)
                   y.adj   <- pheno$y - as.vector( FE.mat %*% as.numeric(ABLUP$beta) )
                   cp      <- cor(ABLUP$u[val], y.adj[val], use = "complete.obs")
                   results <- data.frame(W=0, COR=cp)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='ABLUP', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)

                   # standard GBLUP
                   results <- blupga_GBLUP(G1, pheno, val)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='GBLUP', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)

                   # BLUP|GA with SNPs weighted by genic position
                   results <-  blupga(G1, S, pheno, val)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='POS', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)


                   # BLUP|GA with SNPs weighted by squared GWAS effect size
                   bsq <- est_SNPeffects(pheno, genomat, val, FE.mat)
                   results <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.001, flank=TRUE)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.1', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)
                   results <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.005, flank=TRUE)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.5', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)
                   results <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.01, flank=TRUE)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='EFF1.0', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)

                   # BLUP|GA with weighted candidate gene GRMs
                   results <-  blupga_CAND(G1, pheno, val, GRMs)
                   results %<>% mutate(TRAIT=TRAIT, MODEL='CAND', RUN=run, NSNP=NSNP, VAL=REL)
                   out <- rbind(out, results)

                   gc()
                   write.table(x = out, file="c:/devwork/NGS/Epoly/WGS/GS/GBLUP/Rlog.txt", sep = '\t', append = TRUE, row.names = F,col.names=F,quote=F)
                 }

  stopCluster(cl)
  closeAllConnections()
  return(out)

}

run_one_trait_EFF <- function(pheno, A.ped, G1, Smat, NSNP, genomat, TRAIT, FE.mat=NULL, REL='RANDOM', start=1, end=100)
{
  out <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                    VAL=character(), W=numeric(), COR=numeric())
  #FE.mat  <- model.matrix(ID ~ 1 + as.factor(FE), data=pheno)

  for (run in start:end)
  {
    cat("\nRUN",run," for ",TRAIT," with ",REL," validation \n")
    cat("==============================================================\n")
    if(REL=='REL')        val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
    else if(REL=='UNREL') val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
    else                  val <- read.table(paste("valsets/RANDOM.samples",run, sep="_"), header = F)
    val <- val$V1

    # BLUP|GA with SNPs weighted by squared GWAS effect size
    bsq <- est_SNPeffects(pheno, genomat, val, FE.mat, method = "emRR")
    results <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.001, flank=FALSE)
    results %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.1', RUN=run, NSNP=NSNP, VAL=REL)
    res <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.005, flank=FALSE)
    res %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.5', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)
    res <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.01, flank=FALSE)
    res %<>% mutate(TRAIT=TRAIT, MODEL='EFF1.0', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)

    gc()
    write.table(x = results, file="c:/devwork/NGS/Epoly/WGS/GS/GBLUP/Rlog.txt", sep = '\t', append = TRUE, row.names = F,col.names=F,quote=F)
    out <- rbind(out, results)
  }
  return(out)

}

run_one_trait <- function(pheno, A.ped, G1, Smat, NSNP, genomat, TRAIT, REL='RANDOM', start=1, end=100)
{
  out <- data.frame(TRAIT=character(), MODEL=character(), RUN=numeric(), NSNP=character(),
                    VAL=character(), W=numeric(), COR=numeric())

  # make design matrix for fixed effs. Only done for ABLUP
  FE.mat  <- model.matrix(ID ~ 1 + as.factor(FE), data=pheno)

  for (run in start:end)
  {
    cat("\nRUN",run," for ",TRAIT," with ",REL," validation \n")
    cat("==============================================================\n")
    if(REL=='REL')        val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
    else if(REL=='UNREL') val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
    else                  val <- read.table(paste("valsets/RANDOM.samples",run, sep="_"), header = F)
    val <- val$V1

    # standard ABLUP
    cat("ABLUP\n")
    pheno.blup        <- pheno
    pheno.blup$y[ val ] <- NA

    ABLUP   <- mixed.solve(pheno.blup$y, K=A.ped, X=FE.mat)
    y.adj   <- pheno$y - as.vector( FE.mat %*% as.numeric(ABLUP$beta) )
    cp      <- cor(ABLUP$u[val], y.adj[val], use = "complete.obs")
    results <- data.frame(W=0, COR=cp)
    results %<>% mutate(TRAIT=TRAIT, MODEL='ABLUP', RUN=run, NSNP=NSNP, VAL=REL)
    #out <- rbind(out, results)

    # standard GBLUP
    res <- blupga_GBLUP(G1, pheno, val)
    res %<>% mutate(TRAIT=TRAIT, MODEL='GBLUP', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)

    # BLUP|GA with SNPs weighted by genic position
    res <-  blupga(G1, Smat, pheno, val)
    res %<>% mutate(TRAIT=TRAIT, MODEL='POS', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)

    # BLUP|GA with SNPs weighted by squared GWAS effect size
    bsq <- est_SNPeffects(pheno, genomat, val, FE.mat, method = "emRR")
    res <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.001, flank=TRUE)
    res %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.1', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)
    res <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.005, flank=TRUE)
    res %<>% mutate(TRAIT=TRAIT, MODEL='EFF0.5', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)
    res <- blupga_EFF(G1, pheno, val, genomat, bsq, perc=0.01, flank=TRUE)
    res %<>% mutate(TRAIT=TRAIT, MODEL='EFF1.0', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)

    # BLUP|GA with weighted candidate gene GRMs
    res <-  blupga_CAND(G1, pheno, val, GRMs)
    res %<>% mutate(TRAIT=TRAIT, MODEL='CAND', RUN=run, NSNP=NSNP, VAL=REL)
    results <- rbind(results, res)

    gc()
    write.table(x = results, file="c:/devwork/NGS/Epoly/WGS/GS/GBLUP/Rlog.txt", sep = '\t', append = TRUE, row.names = F,col.names=F,quote=F)
    out <- rbind(out, results)
  }

  return(out)

}

single_gene_GRM <- function(snps_df, snpids)
{
  cat("making GRM from", length(snps_df$ANNOT)," SNPS in gene ", unique(snps_df$ANNOT),"\n")
  top <- which( snpids %in% snps_df$SNPID )
  top.M     <- as.matrix(geno.012$genotype[,top])

  W <- rep(1,length(top))
  S <- cgrm.A(top.M)
  colnames(S) <- geno.012$sample.id
  rownames(S) <- geno.012$sample.id
  return(S)
}


make_ped_matrix <- function(df_ped)
{
  # create a half-sib pedigree A-matrix
  A <- diag(468)
  rownames(A) <- df_ped$Family
  colnames(A) <- df_ped$Family
  for(f in rownames(A))
  {
    A[which(rownames(A)==f),which(colnames(A)==f)] <- 0.25
  }
  diag(A) <- 1
  return(A)
}

pheno_list <- c("OilConc.adj","Mono.adj","sqrtSesqui","cineole","percCin.adj","logapinene","Height","GRW.Height")
pheno <- select_( FT2.uniq, "ID", y="cineole",
                 postblock="postblock", FamOrder="FamOrder") %>%
                 mutate(cuts = cut(as.numeric(FT2.uniq$FamOrder), breaks=c(0,8,16,24,32,40), labels=seq(1,5)) )


genomic_h2 <- function(pheno, G, REL='REL', fixeff=FALSE)
{

  cl<-makeCluster(7, type="SOCK", outfile="")
  registerDoSNOW(cl)
  out <- foreach(run = 1:100, .packages=c("rrBLUP"),
                 .verbose=TRUE, .combine=c, .inorder = FALSE) %dopar%
  {
    if(REL=='REL')  val <- read.table(paste("valsets/REL.samples",  run, sep="_"), header = F)
    else            val <- read.table(paste("valsets/UNREL.samples",run, sep="_"), header = F)
    val <- val$V1
    pheno.blup <- pheno
    pheno.blup$y[ val ] <- NA
    if(fixeff==TRUE)
    {
      GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=G, PEV=F, fixed = c("postblock","cuts"))
    }
    else
    {
      GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=G, PEV=F)
    }
    h = GBLUP$Vg/(GBLUP$Vg+GBLUP$Ve)
    #cat("h2 = ", h, "\n")
    h
  }
  stopCluster(cl)
  closeAllConnections()
  return(out)
}

pairwise_rel <- function(G)
{
  xy <- t(combn(colnames(G), 2))
  return( data.frame(xy, dist=G1[xy]) )
}

make_val_sets(G, pheno, n=10)
{
  pairs <- pairwise_rel(G)
  pairs %<>% left_join(., select(FT2.uniq, Family, ID), by=c("X1"="ID"))
  pairs %<>% left_join(., select(FT2.uniq, Family, ID), by=c("X2"="ID"))
  pairs %<>% mutate(infam = (Family.x==Family.y))

}

#=================== cluster individuals (Saatchi 2011) =======================
# make a genetic distance matrix that accounts for inbreeding
dist <- matrix(nrow = nrow(G1), ncol = ncol(G1))
for(i in 1:nrow(G1)) {
  for(j in 1:ncol(G1)){
    dist[i,j] <- ( 1 - ( G1[i,j]/(sqrt(G1[i,i]*G1[j,j])) ) )
  }
}

# cluster individuals by genetic distance
nclust = 40
clusters    <- kmeans(dist, centers = nclust, nstart = 500, iter.max=1000) # we chose 8 clusters based on an earlier plot.

# how big are the clusters?
clusters$size

G2 <- `diag<-`(G1, 0) # set diagonal of relationship matrix to zero

# get the max and mean pairwise additive relationship per individual with those in its own cluster and with those in each other cluster
Gmax <- data.frame(ID=character(), Fam=factor(), clusterA=numeric(), clusterB=numeric(), Gmax=numeric())
for(i in 1:nclust){
  for(j in 1:nclust){
      Gmax %<>% bind_rows(., data.frame(ID=FT2.uniq$ID[which(clusters$cluster==i)], Fam=FT2.uniq$Family[which(clusters$cluster==i)], clusterA=i, clusterB=j,
                                        Gmax = apply((G2[which(clusters$cluster==i),which(clusters$cluster==j)]), MARGIN = 1, FUN = max),
                                        Gmean = apply((G2[which(clusters$cluster==i),which(clusters$cluster==j)]), MARGIN = 1, FUN = mean)))
  }
}

for(i in 1:nclust){
  for(j in 1:nclust){
    Gmax %<>% bind_rows(., data.frame(ID=FT2.uniq$ID[which(clusters$cluster==i)], Fam=FT2.uniq$Family[which(clusters$cluster==i)], clusterA=i, clusterB=j,
                                      Gmax = apply((dist[which(clusters$cluster==i),which(clusters$cluster==j)]), MARGIN = 1, FUN = max),
                                      Gmean = apply((dist[which(clusters$cluster==i),which(clusters$cluster==j)]), MARGIN = 1, FUN = mean)))
  }
}

# have a look at how the Gmax is distributed within clusters
ggplot(data=Gmax %>% filter(clusterA==clusterB) %>% group_by(clusterA), aes(color=as.factor(clusterA))) +
  geom_density(aes(Gmax)) +
  facet_wrap(~clusterA)
# have a look at how the Gmax is distributed between clusters
ggplot(data=Gmax %>% filter(clusterA!=clusterB) %>% group_by(clusterA), aes(color=as.factor(clusterA))) +
  geom_density(aes(Gmax)) +
  facet_wrap(~clusterA)

# have a look at how the Gmean is distributed within clusters
ggplot(data=Gmax %>% filter(clusterA==clusterB) %>% group_by(clusterA), aes(color=as.factor(clusterA))) +
  geom_density(aes(Gmean)) +
  facet_wrap(~clusterA)
# have a look at how the Gmean is distributed between clusters
ggplot(data=Gmax %>% filter(clusterA!=clusterB) %>% group_by(clusterA), aes(color=as.factor(clusterA))) +
  geom_density(aes(Gmean)) +
  facet_wrap(~clusterA)
