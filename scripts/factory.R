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
    filename <- paste(filename,".csv",sep="")
    write.csv(data.frame(frequencies), file = paste("frequencies_", filename, sep=""),row.names=FALSE)
    write.csv(data.frame(binmatrix), file = paste("binmatrix_", filename, sep=""),row.names=FALSE)
    remove(dfCovidLineage)
    gc()
}
