nIter  <- 12000
burnIn <- 2000

phenotype.file <- '../../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv'

AB       <- list()
AB[[1]]  <- '../../../5.E/output/EVD_E.rda'    # E kernel
AB[[2]]  <- '../../../7.L/output/EVD_L.rda' 
AB[[3]] <- '../../../4.G/output/EVD.rda'    # G kernel  
AB[[3]]  <- '../../../6.P/output/EVD_P_orig.rda'      # P kernel
type     <- c('RKHS', 'RKHS', 'RKHS')

colVAR   <- 2       # Genotype column
colSEA   <- 1       # Season column
colPhen  <- 3       # First trait column (X1) - will loop over all
ESC      <- FALSE
set.seed(1)