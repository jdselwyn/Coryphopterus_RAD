## Take a Genepop File and the Metadata file and calculate pairwise Fst

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  genpop <- args[1]
  metadata <- args[2]
  outDir <- args[3]
  BOOTS <- as.integer(args[4])
  P_CUTOFF <- as.numeric(args[5])
  MIN_SHOAL <- as.integer(args[6])
  
} else {
  
  if(Sys.info()['sysname'] == 'Linux'){
    genpop <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir/MiSeq_CHYA_chyaKonlyHaplo.2.1.Fltr19.popmap.2.1.haps.genepop'
    metadata <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/individual_metadata.shp'
  } else {
    genpop <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir/MiSeq_CHYA_chyaKonlyHaplo.2.1.Fltr19.popmap.2.1.haps.genepop'
    metadata <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/individual_metadata.shp'
    outDir <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir'
  }
  BOOTS <- 0
  P_CUTOFF <- 0.001
  MIN_SHOAL <- 5
  
}

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#### Libraries ####
library(vcfR)
library(sf)
library(tidyverse)
library(adegenet)
library(hierfstat)
library(mcreplicate)

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

pairwise.fst <- function(data, boot = FALSE){
  require(adegenet)
  if(!boot){
    data <- hierfstat::genind2hierfstat(data, pop = as.character(pop(data)))
  } else {
    data <- hierfstat::genind2hierfstat(data, pop = as.character(sample(pop(data))))
  }
  as.matrix(hierfstat::genet.dist(data, method = "WC84"))
}

#### Read in Data ####
if(str_detect(genpop, 'vcf$')){
  genotypes <- read.vcfR(genpop) %>%
    vcfR2genind()
  rownames(genotypes@tab) <- str_c(rownames(genotypes@tab), '2.1', sep = '.')
  
} else {
  genotypes <- read_genepop(genpop, ncode = 3)
}


individual_data <- st_read(metadata) %>%
  filter(ID %in% rownames(genotypes@tab))

strata(genotypes) <- individual_data
strata(genotypes) <- strata(genotypes, formula = ~site/shoal, combine = FALSE)
  
#### Calculate HWE ####
hw_test <- pegas::hw.test(genotypes, B = BOOTS)
loci_hwe <- as_tibble(hw_test, rownames = 'locus') %>%
  janitor::clean_names() %>%
  # rename(pr_exact = pr_chi_2) %>%
  mutate(pr_exact = p.adjust(pr_exact, method = 'holm')) %>%
  filter(pr_exact > P_CUTOFF) %>%
  pull(locus)

message('Using ', length(loci_hwe), ' loci')

#### Pairwise Fst Site ####
setPop(genotypes) <- ~site

obs_fst_site <- pairwise.fst(genotypes[, loc = loci_hwe], FALSE)
perm_fst_site <- mc_replicate(BOOTS - 1, pairwise.fst(genotypes[, loc = loci_hwe], TRUE), 
                              mc.cores = parallel::detectCores(),
                              varlist = c('pairwise.fst', 'genotypes', 'loci_hwe'))

pairwise_fst_site <- obs_fst_site %>%
  as_tibble(rownames = 'site1') %>%
  pivot_longer(cols = -site1,
               names_to = 'site2',
               values_to = 'Fst') %>%
  mutate(p = map2_dbl(site1, site2, ~mean(c(obs_fst_site[.x,.y] < na.omit(sapply(1:(BOOTS - 1), function(i) perm_fst_site[.x, .y, i])), TRUE))))


write_csv(pairwise_fst_site, str_c(outDir, '/', str_extract(genpop, '[A-Za-z0-9\\._]+$') %>% str_remove('.haps.genepop'), '_pairwiseFst_site.csv'))

#### Remove Shoals with 1 fish ####
setPop(genotypes) <- ~site/shoal

large_shoals <- seppop(genotypes) %>%
  tibble(pop = names(.),
         genotypes = .) %>%
  rowwise(pop) %>%
  summarise(n = nInd(genotypes),
            .groups = 'drop') %>%
  filter(n >= MIN_SHOAL) %>%
  pull(pop)

message('Using ', length(large_shoals), ' Shoal')
  

#### Pairwise Fst Shoal ####
obs_fst_shoal <- pairwise.fst(genotypes[pop = large_shoals, loc = loci_hwe], FALSE)
perm_fst_shoal <- mc_replicate(BOOTS - 1, pairwise.fst(genotypes[pop = large_shoals, loc = loci_hwe], TRUE), 
                               mc.cores = parallel::detectCores(),
                               varlist = c('pairwise.fst', 'genotypes', 'loci_hwe', 'large_shoals'))

pairwise_fst_shoal <- obs_fst_shoal %>%
  as_tibble(rownames = 'shoal1') %>%
  pivot_longer(cols = -shoal1,
               names_to = 'shoal2',
               values_to = 'Fst') %>%
  mutate(p = map2_dbl(shoal1, shoal2, ~mean(c(obs_fst_shoal[.x,.y] < na.omit(sapply(1:(BOOTS - 1), function(i) perm_fst_shoal[.x, .y, i])), TRUE))))


write_csv(pairwise_fst_shoal, str_c(outDir, '/', str_extract(genpop, '[A-Za-z0-9\\._]+$') %>% str_remove('.haps.genepop'), '_pairwiseFst_shoal.csv'))
