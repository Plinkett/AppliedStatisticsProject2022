for(i in 790:800) {
    currentLineage <- lineages_frequencies$lineages[i]
    numberOfSamples <- lineages_frequencies$V2[i]
    dfCovidLineage <- subset(dfcovid, dfcovid$pango_lineage == currentLineage)
    mutations <- as.data.frame(matrix(nrow=0,ncol=1))
    mutations <- fetch_mutations(dfCovidLineage, mutations)
    frequencies <- mutation_frequencies(dfCovidLineage, mutations)
    relevant_mutations <- filter05percent(frequencies, numberOfSamples)
    binmatrix <- build_matrix(dfCovidLineage, relevant_mutations)
    filename <- tolower(currentLineage)
    filename <- gsub('\\.','',filename)
    dir.create(filename)
    dirname <- filename
    filename <- paste(filename,".csv",sep="")
    write.csv(data.frame(frequencies), file = paste(dirname,"/frequencies_", filename, sep=""),row.names=FALSE)
    write.csv(data.frame(binmatrix), file = paste(dirname,"/binmatrix_", filename, sep=""),row.names=FALSE)
    remove(dfCovidLineage)
    gc()
}

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

# Function that filters the mutations that appear less than 0.05% of the time

filter05percent <- function(frequencies, totalsamples) {
  frequencies[,3] <- frequencies[,2] / totalsamples
  frequencies <- subset(frequencies, frequencies[,3] > 0.005)
  frequencies <- frequencies[,-3]
  frequencies <- as.data.frame(frequencies)
  frequencies[,1] <- gsub("\\.", ":", frequencies[,1])
  return(frequencies)
}

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
