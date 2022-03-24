covid <- read.csv("metadata_filtrato.tsv", sep="\t")
any(is.na(covid$missing_data))
sum(is.na(covid$missing_data))
nrow(covid)
# setwd("C://users//carlo//Documents//AppliedStatistics")
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

b117 <- subset(covid, covid$pango_lineage == "B.1.1.7")

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
  
# Fetch all unique mutations in a certain dataset (to use with lineages or clades)
fetch_mutations <- function(covid_df) {
  # Assumes aaSubstitutions is the 9th column!
  # According to our conventions in the metadata_animals.csv data set
  # the 10th column refers to aaSubstitutions 
  aasequences <- as.data.frame(covid_df[,9])
  unique_mutations <- as.data.frame(matrix(nrow=0,ncol=1))
  for(i in 1:nrow(aasequences)) {
    # Must do an inner loop...
    aamutations <- as.data.frame(strsplit(as.character(aasequences[i,1]), ",")) 
    for(j in 1:nrow(aamutations)) {
      #tuple <- as.data.frame(strsplit(as.character(aachanges[j,1]),":"))
      mutation <- aamutations[j,1]
      if(!(mutation %in% unique_mutations[,1])) {
        unique_mutations <- rbind(unique_mutations, list(mutation))
      }
    }
  }
  return(unique_mutations)
}


# Functions that builds the final matrix to be analyzed
build_matrix(dfcovid, mutations) {
  matrixnames <- cbind("lineage", "location", "date", t(mutations))
  for(i in 1:nrow(dfcovid)) {
    entry <- dfcovid[i,]
    
  }
}






