
data <- as_tibble(consTables$out)


# List consensus files
fileNames = list.files(
  path = "C:/Users/Stefan/Documents/GitHub/umiSwarm/src",
  pattern =".csv",
  full.names = TRUE
) 

sampleNames = list.files(
  path = "C:/Users/Stefan/Documents/GitHub/umiSwarm/src",
  pattern =".csv",
  full.names = FALSE
) 

# loads all sequence tables in a list with data.frames
consTables <- lapply(fileNames, read.csv)

# name the dataframes within the list
names(consTables) <- gsub(
  pattern = "\\.csv$", 
  replacement = "", 
  x = sampleNames
)

# set colnames suitable for dada2
colnames = c("id", "umi", "sequence", "abundance") 
consTables = lapply(consTables, setNames, colnames)
consTables = lapply(consTables, FUN = function(x){x[, c(3:4)]})

# Generate final sequence table
seqtab_final = dada2::makeSequenceTable(
  samples = consTables
)

# Save sequence table to fasta file
dada2::uniquesToFasta(
  unqs = seqtab_final, 
  fout = file.path(getwd(), "rep-seqs.fna"),
  ids = colnames(seqtab_final)
)
