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
covid <- subset(covid, covid$missing_data <= 900)

sum(covid$host != "Homo sapiens")
sum(covid$pango_lineage == "unclassifiable")

# Togliamo le righe con lineage unclassifiable
covid <- subset(covid, covid$pango_lineage != "unclassifiable")
nrow(covid)

# Togliete l'hashtag del comando sotto per salvare il dataframe in un file .csv 
# write.csv(data.frame(covid), file = "new_metadata.csv",row.names=TRUE)
