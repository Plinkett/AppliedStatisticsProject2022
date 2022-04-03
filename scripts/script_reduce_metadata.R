# Come creare la matrice binaria
#dfcovid_b117 <- subset(dfcovid, dfcovid$pango_lineage == "B.1.1.7")

#mutations_b117 <- as.data.frame(matrix(nrow=0,ncol=1))	
#mutations_b117 <- fetch_mutations(dfcovid_b117, mutations_b117)

#frequencies_b117 <- mutations_frequencies(dfcovid_b117,mutations_b117)

#relevant_mutations <- filter5percent(frequencies_b117, nrow(dfcovid_b117))

#binmatrix_b117 <- build_matrix(dfcovid_b117, relevant_mutations)

# Per salvare i dati 
# write.csv(data.frame(frequencies_b117), file = "frequencies_b117.csv",row.names=FALSE)

# binmatrix_b117_0k_100k <- build_matrix(dfcovid_b117[1:100000,], relevant_mutations) 
# write.csv(data.frame(binmatrix_b1117_0k_100k), file = "binmatrix_b117_0k_100k.csv",row.names=FALSE)
# binmatrix_b117_100k_200k <- build_matrix(dfcovid_b117[100001:200000],relevant_mutations)
# binmatrix_b117_0k_200k <- rbind(binmatrix_b117_0k_100k, binmatrix_b117_100k_200k)
# 
#
#

# Aggiunge gli indices
add_indices <- function(covid) {
  for(i in 1:nrow(covid)) {
    covid[i,10] <- i
  }
}

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
b117 <- subset(b117, (b117$length - b117$missing_data) > 29000)

# THIS FUNCTION IS OBSOLETE
# Takes an aaSequence as a one column dataframe
#protein_changes <- function(aasequence) {
#  proteins <- as.data.frame(matrix(nrow=0,ncol=2)) # Vector that holds all the proteins and number of changes
#  for(i in 1:nrow(aasequence)) {
#    tuple <- as.data.frame(strsplit(as.character(aasequence[i,1]), ":"))  
#    if(tuple[1,1] %in% proteins[,1]) {
#      index <- match(tuple[1,1], proteins[,1])
#      proteins[index,2] <- proteins[index,2] + 1 # a change is added
#    }
#    else {
#      proteins <- rbind(proteins, list(tuple[1,1],1))
#    }
#  }
#  colnames(proteins) <- c("protein","numOfVariations")
#  return(proteins)
#}
  
# After having the total number of frequencies, fetch the most meaningful mutations
ratios <- cbind(total_frequencies, total_frequencies[,1] / total_samples)
good_ratios <- subset(ratios, ratios$V2 <= 9.5e-01 & ratios$V2 >= 5e-02)
# Fetch all unique mutations in a certain dataset (to use with lineages or clades)
fetch_mutations <- function(covid_df, unique_mutations) {
  # Assumes aaSubstitutions is the 9th column!
  # According to our conventions in the metadata_animals.csv data set
  # the 10th column refers to aaSubstitutions 
  aasequences <- as.data.frame(covid_df[,9])
  #unique_mutations <- as.data.frame(matrix(nrow=0,ncol=1))
  for(i in 1:nrow(aasequences)) {
    # Must do an inner loop...
    aamutations <- as.data.frame(strsplit(as.character(aasequences[i,1]), ",")) 
    for(j in 1:nrow(aamutations)) {
      #tuple <- as.data.frame(strsplit(as.character(aachanges[j,1]),":"))
      mutation <- aamutations[j,1]
      if(length(mutation) != 0) {
        if(!(mutation %in% unique_mutations[,1])) {
          unique_mutations <- rbind(unique_mutations, list(mutation))
        }
      }
    }
  }
  return(unique_mutations)
}

# Build mutation frequencies, both in dataframe forms
# mutations as a single column dataframe
mutation_frequencies <- function(dfcovid, mutations) {
  freq <-  as.data.frame(matrix(nrow=1,ncol=nrow(mutations)))
  freq[1,] <- integer(nrow(mutations))
  colnames(freq) <- mutations[,1]
  aasequences <- as.data.frame(dfcovid[,9])
  for(i in 1:nrow(aasequences)) {
    aamutations <- as.data.frame(strsplit(as.character(aasequences[i,1]), ","))
    for(j in 1:nrow(aamutations)) {
      mutation_i <- match(aamutations[j,1], mutations[,1])
      freq[1,mutation_i] <- freq[1,mutation_i] + 1      
    }
  }
  return(cbind(mutations,t(freq)))
}

# Function that filters the mutations that appear less than 5% of the time
filter5percent <- function(frequencies, totalsamples) {
  frequencies[,3] <- frequencies[,2] / totalsamples
  frequencies <- subset(frequencies, frequencies[,3] > 0.05)
  frequencies <- frequencies[,-3]
  return(as.data.frame(frequencies))
}

# sapply(1:10, list_mutations) OR
# sapply(1:10, function(mutations) strsplit((aasequence),","))
# Second version with sapply
#list_mutations <- function(aasequence, index) {
#  tokens <- as.data.frame(strsplit((as.character(aasequence)),","))
#  return(t(as.data.frame(sapply(tokens[,1], function(token, id) cbind(id, token), id=index))))
#}

# Second version of split_by_sequence, using sapply
#split_by_sequence_2_0 <- function(dfcovid) {
#  matrices <- mapply(list_mutations, dfcovid[,2], dfcovid[,1])
#}

# First version of the split_mutations_by_sequence
#split_mutations_by_sequence_2_0 <- function(dfcovid) {
#  final_df <- as.data.frame(matrix(nrow=0,ncol=2))
#  for(i in 1:nrow(dfcovid)) {
#    mutations <- list_mutations(dfcovid[i,2], dfcovid[i,1])
#    final_df <- rbind(final_df, mutations)
#  }
# return(final_df)
#``}


split_mutations_by_sequence <- function(dfcovid) {
  final_df <- as.data.frame(matrix(nrow=0,ncol=2))
  for(i in 1:nrow(dfcovid)) {
    mutations <- as.data.frame(strsplit(as.character(dfcovid[i,9]), ","))  
    for(j in 1:nrow(mutations)) {
      final_df <- rbind(final_df, c(dfcovid[i,10],mutations[j,1]))  
    }
  }
  return(final_df)
}


# Function that builds the binary matrix of mutations 
# column names
# [lineage, country, date, mutation1, mutation2, ..., mutationN]
# Input:
#   dfcovid <- dataframe as presented in the original metadata_filtrato.csv dataset 
#              (minus the "host" column, which is assumed to be "Homo sapiens").
#   important_mutations <- the set of mutations deemed relevant for our analysis, that is the mutations
#                          that appear less than 95% and more than 5% of the time.
build_matrix <- function(dfcovid, important_mutations) {
  # Initial declaration of the matrix, empty matrix
  finalmatrix <- as.data.frame(matrix(nrow=0, ncol = (3 + nrow(important_mutations))))
  # Outer for that scans the dataframe one row at a time
  for(i in 1:nrow(dfcovid)) {
    # the "toadd" dataframe is a one row dataframe containing the row to be added to the matrix
    # we fill it with the lineage, country, date + a number of columns equal to that of the 
    # relevant mutations
    toadd <- as.data.frame(matrix(nrow=1, ncol = (3 + nrow(important_mutations))))
    toadd[1,1] <- dfcovid[i,]$pango_lineage;
    toadd[1,2] <- dfcovid[i,]$country
    toadd[1,3] <- dfcovid[i,]$date
    # We initialize the 4th column up to nrow(important_mutations) with zeroes, remember
    # that these are binary variables!
    toadd[1,4:ncol(toadd)] <- t(as.data.frame(integer(nrow(important_mutations))))
    # The aaSubstitutions column is split by "," and converted into a dataframe to be scanned
    # for example:
    #   N:G204R,ORF1a:Q3966R,S:K1191N (one row)
    #   becomes this:
    #     - N:G204R
    #     - ORF1a:Q3966R
    #     - S:K1191N
    aamutations <- as.data.frame(strsplit(as.character(dfcovid[i,9]), ","))
    # Now we scan the mutations
    for(j in 1:nrow(aamutations)) {
      # We search for the index of the current mutation in the important_mutations row
      # If we don't find it then it means that the mutations is not relevant and we don't consider it
      # this also means that match() will return a NA (absence of a value)
      mutation_i <- match(aamutations[j,1], important_mutations[,1])
      if(!is.na(mutation_i)) {
        # If we find the mutation then we mark it with a 1 (the mutation occurred in this sample)
        toadd[1, 3 + mutation_i] <- 1
      }
    }
    # We add the "toadd" row to the matrix (merging the two with rbind())
    finalmatrix <- rbind(finalmatrix, toadd[1,])
  }
  # We rename the column names
  colnames(finalmatrix) <- c("lineage", "country", "date", t(important_mutations[,1]))
  return(finalmatrix)
}






