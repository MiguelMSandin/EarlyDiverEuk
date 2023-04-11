 
library(data.table)
library(dplyr)

setwd("")

data = fread("tipNamesIDs.tsv", header=FALSE)
data$V2 = NULL
colnames(data) = c("name")

# First, we extract the species names
data$species = fifelse(grepl("^PacBio_", data$name),
					   sub(".*_Otu\\d+_\\d+_", "", data$name),
					   sub("_Otu\\d+.*$", "", data$name))
data$species = gsub("_(X)+_", "-X_", data$species)
data$species = sub("^.+((_[^_]+){2})$", "\\1", data$species) %>% sub("^_", "", .)

cat("There are", length(unique(data$species)), "unique 'species' names")
# let's now decide whether a sequence it's 'described' or 'environmental' based on its latin name.
# Briefly, if a sequence has a binomial latin name, we assume that given sequence has been morphologically described,
# OR since PR2 is a CURATED database, at least that given sequence has been phylogenetically annotated as a particular species or genus.
# Otherwise, we assume the sequence comes from the environment and has no known associated morphological name.
# However, there are 27614 unique names, so it will be very time-consuming to do this by hand...
# so we will apply Regular Expression to try to find the best compromise

data$description = data$species
# Let's begin by taking profit of the PR2 database structure:
# If a clade has "_X" (or as in here "-X"), means it's dragging the name from a higher hierarchical level, so it's unknown
data$description[grep("-X", data$description)] = "environmental"
data$description[grep("_X", data$description)] = "environmental"

sort(unique(data[description=="environmental"]$species))

# Now we can proceed with the most abundant specific names
sort(table(data$description), decreasing=TRUE)[1:20]
# unique(grep("Supergroup", data$description, value=TRUE))
data$description[grep("Opisthokonta_Fungi", data$description)] = "environmental"
data$description[grep("mycota", data$description)] = "environmental"
data$description[grep("mycetes_|mycetes$", data$description)] = "environmental"
data$description[grep("Dinophyceae|Dinoflagellata", data$description)] = "environmental"
data$description[grep("Dino-Group", data$description)] = "environmental"
data$description[grep("MAST", data$description)] = "environmental"
data$description[grep("dales.{,3}$", data$description)] = "environmental"
data$description[grep("Eukaryota", data$description)] = "environmental"
data$description[grep("phycea", data$description)] = "environmental"
data$description[grep("Polycystinea|Collodaria", data$description)] = "environmental"
data$description[grep("Rhizaria|Cercozoa", data$description)] = "environmental"
data$description[grep("monadida$", data$description)] = "environmental"
data$description[grep("Obazoa|Opisthokonta", data$description)] = "environmental"
data$description[grep("TSAR|Stramenopiles", data$description)] = "environmental"
data$description[grep("Supergroup", data$description)] = "environmental"

# We consider every species name with at least two consecutive capital letters as environmental
unique(data$name[grep("[A-Z]{2,}", data$description)])
data$description[grep("[A-Z]{2,}", data$description)] = "environmental"
# We check for those with two capital letters
unique(grep("[A-Z].*[A-Z]", data$description, value=TRUE))
sort(table(grep("[A-Z].*[A-Z]", data$description, value=TRUE)), decreasing=TRUE)[1:20]
# We change some families now
unique(grep("[A-Z].*[A-Z].*aceae$|aceae_sp$", data$description, value=TRUE))
data$description[grep("[A-Z].*[A-Z].*idae$|idae_sp$", data$description)] = "environmental"
data$description[grep("[A-Z].*[A-Z].*aceae$|aceae_sp$", data$description)] = "environmental"
data$description[grep("Euglenozoa", data$description)] = "environmental"
data$description[grep("Spirotrichea|Hypotrichia", data$description)] = "environmental"
data$description[grep("Kinetoplastida", data$description)] = "environmental"
data$description[grep("Alveolata|Perkinsea", data$description)] = "environmental"
data$description[grep("Radiolaria|Acantharea", data$description)] = "environmental"
data$description[grep("Endomyxa", data$description)] = "environmental"
data$description[grep("Eugregarinorida", data$description)] = "environmental"
data$description[grep("Cercomonadida", data$description)] = "environmental"
data$description[grep("Filosa", data$description)] = "environmental"
data$description[grep("Alveolata|Ciliophora", data$description)] = "environmental"
data$description[grep("Chelicerata|Arachnida", data$description)] = "environmental"
data$description[grep("Ciliophora|Spirotrichea", data$description)] = "environmental"
data$description[grep("Crustacea|Maxillopoda", data$description)] = "environmental"
# And we repeat the game
sort(table(grep("[A-Z].*[A-Z]", data$description, value=TRUE)), decreasing=TRUE)
data$description[grep("Apicomplexa|Gregarinomorphea", data$description)] = "environmental"
data$description[grep("Dolichomastigaceae", data$description)] = "environmental"
data$description[grep("Choanoflagellatea", data$description)] = "environmental"
data$description[grep("Bicosoecida", data$description)] = "environmental"
data$description[grep("Opisthosporidia", data$description)] = "environmental"
data$description[grep("Nucleariida", data$description)] = "environmental"
data$description[grep("Hexapoda_Insecta", data$description)] = "environmental"
data$description[grep("Mucoromycotina|Mucorales", data$description)] = "environmental"
data$description[grep("Gyrista|Diatomeae", data$description)] = "environmental"
data$description[grep("Ichthyosporea_Ichthyophonida", data$description)] = "environmental"
data$description[grep("Hyphochytriomyceta|Hyphochytriales", data$description)] = "environmental"
data$description[grep("Metazoa|Nematoda", data$description)] = "environmental"
data$description[grep("Group-Te", data$description)] = "environmental"
data$description[grep("Clade_Na13", data$description)] = "environmental"
data$description[grep("[A-Z].*[A-Z].*ales$", data$description)] = "environmental"
data$description[grep("Arcellinida", data$description)] = "environmental"
data$description[grep("Apostomatia", data$description)] = "environmental"
data$description[grep("Group", data$description)] = "environmental"
data$description[grep("lineage", data$description)] = "environmental"
data$description[grep("Centramoebia|Acanthopodida", data$description)] = "environmental"
data$description[grep("Rhynchostomatia", data$description)] = "environmental"
data$description[grep("Monomastigaceae", data$description)] = "environmental"
data$description[grep("Crustomastigaceae", data$description)] = "environmental"
data$description[grep("Gastropoda|Heterobranchia", data$description)] = "environmental"
data$description[grep("Ichthyophonidae_Amoebidianae", data$description)] = "environmental"
data$description[grep("clade", data$description)] = "environmental"
data$description[grep("^Clade_D$", data$description)] = "environmental"

# Check a rough proportion
sum(grepl("environmental", data$description))/nrow(data)
sum(!grepl("environmental", data$description))/nrow(data)

data$species = NULL
data$description[which(data$description != "environmental")] = "described"

table(data$description)

data$colour = fifelse(data$description=="environmental", "#648FFF", "#FFB000")
table(data$colour)
data$description = NULL

write.table(data, "tipNames_envVSdes.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

system(paste0("treeColourBranches.py -t YOUR_TREE_OF_INTEREST -c tipNames_envVSdes.tsv"))
