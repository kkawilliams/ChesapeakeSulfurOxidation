library(dada2)
set.seed(100)

# command line arguments
args <- commandArgs(TRUE)

# read in location
base_path <- args[1]
seq_ID <- args[2]
fwd_ID <- args[3]
rev_ID <- args[4]
sample_splitter <- args[5]
n_threads <- strtoi(args[6])
# base_path <- "/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/Processed_data_group"
# seq_ID <- "sprehei1_149186"
# fwd_ID <- "_F_filt"
# rev_ID <- "_R_filt"
# sample_splitter <- "_"
# n_threads <- 1

trim_path = file.path(base_path, seq_ID, "FASTQ", "Trim")

# get file & sample names & paths
pre_fnFs <- sort(list.files(trim_path, pattern=fwd_ID))
pre_fnRs <- sort(list.files(trim_path, pattern=rev_ID))
sample.names <- sapply(strsplit(pre_fnFs, sample_splitter), `[`, 1)
fnFs <- file.path(trim_path, pre_fnFs)
fnRs <- file.path(trim_path, pre_fnRs)
names(fnFs) <- sample.names
names(fnRs) <- sample.names

# make intermediate saved files
error_file_name_F = paste(substr(seq_ID, 1, 15), "errorsWS_F.RData", sep="_")
error_file_path_F = file.path(base_path, seq_ID, error_file_name_F)
error_file_name_R = paste(substr(seq_ID, 1, 15), "errorsWS_R.RData", sep="_")
error_file_path_R = file.path(base_path, seq_ID, error_file_name_R)

# fit error model

if (file.exists(error_file_path_F)){
    write(paste("Loading (F) errors from", error_file_name_F), stdout())
    load(error_file_path_F)
} else {
    write(paste("Start Errors (1)", Sys.time() ), stdout())
    errF <- learnErrors(fnFs, nbases = 1e8, multithread=n_threads, randomize=TRUE)
    save(errF, file=error_file_path_F)
    write(paste("Writing (F) errors to", error_file_name_F), stdout())
}

if (file.exists(error_file_path_R)){
    write(paste("Loading (R) errors from", error_file_name_R), stdout())
    load(error_file_path_R)
} else {
    write(paste("Start Errors (2)", Sys.time() ), stdout())
    errR <- learnErrors(fnRs, nbases = 1e8, multithread=n_threads, randomize=TRUE)
    save(errR, file=error_file_path_R)
    write(paste("Writing (R) errors to", error_file_name_R), stdout())
}

# dereplicate & call OTUs
dds <- vector("list", length(sample.names))
names(dds) <- sample.names

for(sam in sample.names) {
   cat("Processing:", sam, "\n")
   write(paste("Start Derep (1)", Sys.time() ), stdout())
   derepF <- derepFastq(fnFs[[sam]])
   write(paste("Start Derep (2)", Sys.time() ), stdout())
   derepR <- derepFastq(fnRs[[sam]])
   write(paste("Start DADA (1)", Sys.time() ), stdout())
   ddF <- dada(derepF, err=errF, multithread=n_threads, pool="pseudo")
   write(paste("Start DADA (2)", Sys.time() ), stdout())
   ddR <- dada(derepR, err=errR, multithread=n_threads, pool="pseudo")
   merger <- mergePairs(ddF, derepF, ddR, derepR)
   dds[[sam]] <- merger
}

otu_tab_name = paste(substr(seq_ID, 1, 15), "raw_tab.RData", sep="_")
otu_tab_path = file.path(base_path, seq_ID, otu_tab_name)
write(paste("Finished DADA", Sys.time() ), stdout())
write(paste("Writing raw table to", otu_tab_path), stdout())
save(dds, file=otu_tab_path)

seq_tab_name = paste(substr(seq_ID, 1, 15), "seqtab_chim.rds", sep="_")
seq_tab_path = file.path(base_path, seq_ID, seq_tab_name)
seqtab1 <- makeSequenceTable(dds)
write(paste("Finished Seq table", Sys.time() ), stdout())
write(paste("Writing seq table to", seq_tab_path), stdout())
saveRDS(seqtab1, seq_tab_path)

#nochim_tab_name = paste(substr(seq_ID, 1, 15), "seqtab_nochim.rds", sep="_")
#nochim_tab_path = file.path(base_path, seq_ID, nochim_tab_name)
#seqtab2 <- removeBimeraDenovo(seqtab1, method="consensus", multithread=n_threads)
#write(paste("Finished Chimera Checking", Sys.time() ), stdout())
#write(paste("Writing dechimera'd table to", nochim_tab_path), stdout())
#saveRDS(seqtab2, nochim_tab_path)