## Take a Genepop File and the Metadata file and calculate pairwise PCA distance - using loci in HWE

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  genpop <- args[1]
  metadata <- args[2]
  outDir <- args[3]
  BOOTS <- as.integer(args[4])
  P_CUTOFF <- as.numeric(args[5])
  max_missing <- as.numeric(args[6])
  
} else {
  
  if(Sys.info()['sysname'] == 'Linux'){
    genpop <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir/MiSeq_CHYA_chyaKonlyHaplo.2.1.Fltr19.popmap.2.1.haps.genepop'
    metadata <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/individual_metadata.shp'
  } else {
    genpop <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir/MiSeq_CHYA_chyaKonlyHaplo.2.1.Fltr19.popmap.2.1.haps.genepop'
    metadata <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/individual_metadata.shp'
    outDir <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir'
    
    BOOTS <- 0
    P_CUTOFF <- 0.001
    max_missing <- 0.75
  }
}

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#### Libraries ####
library(sf)
library(tidyverse)
library(adegenet)
library(ade4)

#### Functions ####
read_genepop <- function (file, ncode = 2L, quiet = FALSE){
  if (!toupper(.readExt(file)) %in% c("GEN", "GENEPOP")) 
    stop("File extension .gen expected")
  if (!quiet) 
    cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
  prevcall <- match.call()
  txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
  if (!quiet) 
    cat("\nFile description: ", txt[1], "\n")
  txt <- txt[-1]
  txt <- gsub("\t", " ", txt)
  NA.char <- paste(rep("0", ncode), collapse = "")
  locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
  locinfo <- txt[locinfo.idx]
  locinfo <- paste(locinfo, collapse = ",")
  loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
  loc.names <- trimws(loc.names)
  nloc <- length(loc.names)
  txt <- txt[-locinfo.idx]
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
  npop <- length(pop.idx)
  nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
  splited <- nocomma[which(!nocomma %in% pop.idx)]
  if (length(splited) > 0) {
    for (i in sort(splited, decreasing = TRUE)) {
      txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
    }
    txt <- txt[-splited]
  }
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
  txt[length(txt) + 1] <- "POP"
  nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$", 
                          toupper(txt))) - 1
  pop <- factor(rep(1:npop, nind.bypop))
  txt <- txt[-c(pop.idx, length(txt))]
  temp <- sapply(1:length(txt), function(i) strsplit(txt[i], 
                                                     ","))
  ind.names <- vapply(temp, function(e) e[1], character(1))
  ind.names <- trimws(ind.names)
  vec.genot <- vapply(temp, function(e) e[2], character(1))
  vec.genot <- trimws(vec.genot)
  X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), 
              ncol = nloc, byrow = TRUE)
  if (any(duplicated(ind.names))) {
    rownames(X) <- .genlab("", nrow(X))
    warning("Duplicate individual names detected. Coercing them to be unique.")
  }
  else {
    rownames(X) <- ind.names
  }
  colnames(X) <- loc.names
  pop.names.idx <- cumsum(table(pop))
  pop.names <- ind.names[pop.names.idx]
  levels(pop) <- pop.names
  if (!all(unique(nchar(X)) == (ncode * 2))) 
    stop(paste("some alleles are not encoded with", ncode, 
               "characters\nCheck 'ncode' argument"))
  res <- df2genind(X = X, pop = as.character(pop), ploidy = 2, 
                   ncode = ncode, NA.char = NA.char)
  res@call <- prevcall
  if (!quiet) 
    cat("\n...done.\n\n")
  return(res)
}

#### Read in Data ####
genotypes <- read_genepop(genpop, ncode = 3)

individual_data <- st_read(metadata) %>%
  filter(ID %in% rownames(genotypes@tab))

strata(genotypes) <- individual_data
strata(genotypes) <- strata(genotypes, formula = ~site/shoal, combine = FALSE)

#### Calculate HWE ####
# genotypes <- genotypes[loc = sample(3047, 10)]

hw_test <- pegas::hw.test(genotypes, B = BOOTS)
loci_hwe <- as_tibble(hw_test, rownames = 'locus') %>%
  janitor::clean_names() %>%
  # rename(pr_exact = pr_chi_2) %>%
  mutate(pr_exact = p.adjust(pr_exact, method = 'holm')) %>%
  filter(pr_exact > P_CUTOFF) %>%
  pull(locus)

message('Using ', length(loci_hwe), ' loci')

individual_use <- propTyped(genotypes[, loc = loci_hwe], by = 'ind')  %>%
  enframe(name = 'ID', value = 'prop_typed') %>%
  filter(prop_typed > max_missing) %>%
  pull(ID)

message('Using ', length(individual_use), ' individuals')

#### Calculate PCA Distance ####
full_pca <- scaleGen(genotypes[individual_use, loc = loci_hwe, drop = TRUE], center = TRUE, scale = TRUE, NA.method = "mean") %>%
  dudi.pca(., center = FALSE, scale = FALSE, scannf = FALSE, nf = min(dim(.)))


write_csv(as_tibble(full_pca$li, rownames = 'ID'), 
          str_c(outDir, '/', str_extract(genpop, '[A-Za-z0-9\\._]+$') %>% str_remove('.haps.genepop'), '_PCA_coordinates.csv'))



distance_matrix <- dist(full_pca$li) %>%
  as.matrix() %>%
  as_tibble(rownames = 'ID1') %>%
  pivot_longer(cols = -ID1,
               names_to = 'ID2',
               values_to = 'PCA_distance')

write_csv(distance_matrix, 
          str_c(outDir, '/', str_extract(genpop, '[A-Za-z0-9\\._]+$') %>% str_remove('.haps.genepop'), '_PCA_distances.csv'))
