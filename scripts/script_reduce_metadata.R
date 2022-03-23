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


