library("QuGAcomp")
library("corrplot")

genome.length.file <- file.path(system.file(package="QuGAcomp"), "data", "mm9.info")

oct4.bed.file <- file.path(
  system.file(package="QuGAcomp"),
  "data",
  "GSM288346_ES_Oct4.mm9.header.bed"
)
sox2.bed.file <- file.path(
  system.file(package="QuGAcomp"),
  "data",
  "GSM288347_ES_Sox2.mm9.header.bed"
)
nanog.bed.file <- file.path(
  system.file(package="QuGAcomp"),
  "data",
  "GSM288345_ES_Nanog.mm9.header.bed"
)

oct4.gr  <- loadBedFile(oct4.bed.file,  genome.length.file)
sox2.gr  <- loadBedFile(sox2.bed.file,  genome.length.file)
nanog.gr <- loadBedFile(nanog.bed.file, genome.length.file)

oct4.fat  <- fat(oct4.gr,  200)
sox2.fat  <- fat(sox2.gr,  200)
nanog.fat <- fat(nanog.gr, 200)

oct4.unistd  <- unifyStrand(oct4.fat)
sox2.unistd  <- unifyStrand(sox2.fat)
nanog.unistd <- unifyStrand(nanog.fat)

oct4.cov  <- coverage(oct4.unistd)
sox2.cov  <- coverage(sox2.unistd)
nanog.cov <- coverage(nanog.unistd)

oct4.bin500  <- lapply( oct4.cov,  function(x) rleBinning(x, 500) )
sox2.bin500  <- lapply( sox2.cov,  function(x) rleBinning(x, 500) )
nanog.bin500 <- lapply( nanog.cov, function(x) rleBinning(x, 500) )

oct4.bin500  <- flatRleList(oct4.bin500)
sox2.bin500  <- flatRleList(sox2.bin500)
nanog.bin500 <- flatRleList(nanog.bin500)

quga.oct4.sox2  <- qugacomp(oct4.bin500,  sox2.bin500)
quga.oct4.nanog <- qugacomp(oct4.bin500, nanog.bin500)
quga.sox2.nanog <- qugacomp(sox2.bin500, nanog.bin500)

num <- 3
mat.cor <- matrix(0, nrow=num, ncol=num)
rownames(mat.cor) <- c("Oct4", "Sox2", "Nanog")
colnames(mat.cor) <- c("Oct4", "Sox2", "Nanog")

diag(mat.cor) <- rep(1, num)

mat.cor[1,2] <- mat.cor[2,1] <- pearsonCoef(quga.oct4.sox2)
mat.cor[1,3] <- mat.cor[3,1] <- pearsonCoef(quga.oct4.nanog)
mat.cor[2,3] <- mat.cor[3,2] <- pearsonCoef(quga.sox2.nanog)

mat.cor.max <- max( mat.cor[upper.tri(mat.cor, diag=F)] )
mat.cor.min <- min( mat.cor[upper.tri(mat.cor, diag=F)] )

pdf("corrplot.pdf")
corrplot(
  mat.cor, method="circle", type="full", diag=FALSE, outline=FALSE,
  addcolorlabel="bottom", cl.lim=c(mat.cor.min, mat.cor.max), cl.ratio=0.2
)
dev.off()
