covid <- read.csv("metadata_filtrato.tsv", sep="\t")
any(is.na(covid$missing_data))
sum(is.na(covid$missing_data))
nrow(covid)

# Cancella righe con NA
covid <- na.omit(covid)

mean(covid$missing_data)
median(covid$missing_data)
sd(covid$missing_data)

sum(covid$missing_data > 900)
sum(covid$missing_data > 1000)
sum(covid$missing_data > 5000)
sum(covid$missing_data > 7000)

# Cancelliamo le righe con missing_data > 900 
# covid <- subset(covid, covid$missing_data <= 900)

sum(covid$host != "Homo sapiens")
sum(covid$pango_lineage == "unclassifiable")

# Togliamo le righe con lineage unclassifiable
covid <- subset(covid, covid$pango_lineage != "unclassifiable")
nrow(covid)

# Deleting useless and problematic rows
covid <- subset(covid, covid$date != "?")
covid <- subset(covid, covid$host != "")
covid <- subset(covid, covid$Nextstrain_clade != "")
covid <- subset(covid, covid$pango_lineage != "")

# They are all humans anyway
covid <- covid[,-6]

# Creating covid_animals dataset
covid_animals <- subset(covid, covid$host != "Homo sapiens")

# See unique parameters

lineages <- unique(covid$pango_lineage)
lineages <- as.data.frame(lineages)
clades <- as.data.frame(unique(covid$Nextstrain_clade))
region_country_division <- as.data.frame(unique(covid[,c('region','country','division')]))

# Togliete l'hashtag del comando sotto per salvare il dataframe in un file .csv 
# write.csv(data.frame(covid), file = "metadata_filtered_humans.csv",row.names=FALSE)
# write.csv(data.frame(covid_animals), file = "metadata_animals.csv",row.names=FALSE)
# write.csv(data.frame(lineages), file = "lineages.csv",row.names=FALSE)
# write.csv(data.frame(clades), file = "clades.csv",row.names=FALSE)
# write.csv(data.frame(region_country_division), file = "region_country_division.csv",row.names=FALSE)

# Playing with aaSequences
# aasequence <- covid_animals[1,10]
# aasequence <- as.data.frame(strsplit(aasequence, ","))

# Takes an aaSequence as a one column dataframe
protein_changes <- function(aasequence) {
  proteins <- as.data.frame(matrix(nrow=0,ncol=2)) # Vector that holds all the proteins and number of changes
  for(i in 1:nrow(aasequence)) {
    tuple <- as.data.frame(strsplit(as.character(aasequence[i,1]), ":"))  
    if(tuple[1,1] %in% proteins[,1]) {
      index <- match(tuple[1,1], proteins[,1])
      proteins[index,2] <- proteins[index,2] + 1 # a change is added
    }
    else {
      proteins <- rbind(proteins, list(tuple[1,1],1))
    }
  }
  colnames(proteins) <- c("protein","numOfVariations")
  return(proteins)
}
