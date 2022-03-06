library(biomformat)
library(data.table)
suppressMessages(library(Biostrings))

# Patients/samples to remove
rm_spls <- c("CM-01-12-V1",
             "CM-01-12-V2",
             "CM-01-12-V3",
             "CM-01-15-V1",
             "CM-01-15-V2",
             "CM-01-24-V1",
             "CM-02-01-V3",
             "CM-04-05-V3",
             "CM-05-10-V1",
             "CM-11-01-V1",
             "CM-11-01-V2",
             "CM-15-01-V1",
             "CM-15-01-V2")

# load meta info
spl_meta <- fread("data/raw/metadata-cu_cortison-clean.txt")[-1,]
setnames(spl_meta, "sample-id","sample")

# saving asv abundance table
biom <- read_biom("data/raw/unzip/feature-table.biom")
biom <- as.matrix(biom_data(biom))
hist(colSums(biom))

# taxonomy predictions
taxi <- fread("data/raw/unzip/taxonomy.tsv")
all(rownames(biom) %in% taxi$`Feature ID`)

# Remove Cyanobacteria/Chloroplasts, Archaea and Mitrochondria
rm_feats <- taxi[grepl("Cyano", Taxon) | grepl("itoch", Taxon) | !grepl("^d__Bacteria", Taxon), `Feature ID`]
biom <- biom[!(rownames(biom) %in% rm_feats),]
taxi <- taxi[!(`Feature ID` %in% rm_feats)]

# remove some samples (see above)
biom <- biom[, !(colnames(biom) %in% rm_spls)]
spl_meta <- spl_meta[!(sample %in% rm_spls)]

# remove zero-only features
rm_zerofeats <- names(which(apply(biom, 1, function(x) all(x == 0))))
biom <- biom[!(rownames(biom) %in% rm_zerofeats),]
taxi <- taxi[`Feature ID` %in% rownames(biom)]

# save sequence fasta of remaining features
seqs <- readDNAStringSet("data/raw/unzip/dna-sequences.fasta")
seqs <- seqs[rownames(biom)]

#--------------------#
# export clean files #
#--------------------#
write.table(biom,"data/clean/asv_tab.tsv")
writeXStringSet(seqs, "data/clean/asv_seq.fasta")
fwrite(taxi,"data/clean/asv_tax.csv")
fwrite(spl_meta, "data/clean/sample_info.csv")

