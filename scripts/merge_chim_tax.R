library(dada2)
set.seed(100)
# Merge multiple runs

# command line arguments
args <- commandArgs(TRUE)

# read in location
taxa_path <- args[1]
# taxa_path <- "../data/Silva_DB"
seq_ID_df <- args[2]
# seq_ID_df <- "../data/TrimOTUsData/seq_tabs.csv"
out_dir <- args[3]
# out_dir <- "../data/TrimOTUsData"

id_df = read.csv(seq_ID_df, sep=",")

otu_tabs <- vector("list")

for (i in 1:dim(id_df)[1]) {
	cat(i, as.character(id_df$LibName[[i]]), "\n")
	otu_tabs[[as.character(id_df$LibName[[i]])]] <- readRDS(as.character(id_df$FilePath[[i]]))
}

st.all <- mergeSequenceTables(otu_tabs[[as.character(id_df$LibName[[1]])]], 
	                          otu_tabs[[as.character(id_df$LibName[[2]])]],
	                          otu_tabs[[as.character(id_df$LibName[[3]])]],
	                          otu_tabs[[as.character(id_df$LibName[[4]])]], 
	                          otu_tabs[[as.character(id_df$LibName[[5]])]], 
	                          otu_tabs[[as.character(id_df$LibName[[6]])]], 
	                          otu_tabs[[as.character(id_df$LibName[[7]])]])
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign genus level taxonomy
silva_gf = "silva_nr_v132_train_set.fa.gz"
silva_gp = file.path(taxa_path, silva_gf)
taxa_g <- assignTaxonomy(seqtab, silva_gp, multithread=TRUE, tryRC=TRUE)

# Add species level taxonomy
silva_sf = "silva_species_assignment_v132.fa.gz"
silva_sp = file.path(taxa_path, silva_sf)
taxa_sp <- addSpecies(taxa_g, silva_sp, tryRC=TRUE)

# Write to disk
saveRDS(seqtab, file.path(out_dir, "seqtab_final.rds"))
saveRDS(taxa_g, file.path(out_dir, "tax_g_final.rds"))
saveRDS(taxa_sp, file.path(out_dir, "tax_sp_final.rds"))
