# Kako narediš filter za NGS?
# V datoteki "input_primers_tags.xlsx" zavihku "tags" se nahaja kombinacija Forward in
# Reverse tagov za 8 plošč.
# Te tage ponovimo kolikor je lokusov in k njim pripnemo imena vzorcev in lokuse iz
# zavihkov "prims" in "plate".
# Izpiše .ngsfilter datoteko, ki bi morala izgledat tako kot se šika.

library(readxl)
library(stringr)

tags <- as.data.frame(read_excel("gatc_trial/input_primers_tags.xlsx",
                    sheet = "tags"))

prims <- as.data.frame(read_excel("gatc_trial/input_primers_tags.xlsx",
                                  sheet = "primers"))

plate <- as.data.frame(read_excel("gatc_trial/input_primers_tags.xlsx",
                                  sheet = "plate"))

plate <- gsub("_", replacement = ".", plate$samplename)
plate <- paste(plate, str_pad(1:96, 3, pad = "0"), sep = "_")
# vsak well pomnoži tolikokrat kolikokrat je lokusov
num.loci <- length(unique(prims$locusname))
wells <- rep(1:nrow(tags), each = num.loci)

ngsfilter <- tags[wells, ]
ngsfilter$samplename <- rep(plate, each = num.loci)
ngsfilter$samplename <- paste(ngsfilter$samplename, "_", ngsfilter$plate, sep = "")
ngsfilter$locusname <- prims$locusname
ngsfilter$primer1 <- prims$primer1
ngsfilter$primer2 <- prims$primer2
ngsfilter$well <- wells

rownames(ngsfilter) <- NULL

ngsfilter[ngsfilter$well == 8*96, ]

# write file

out <- ngsfilter[, c("locusname", "samplename", "tagcombo", "primer1", "primer2")]
head(out)
out$F <- "F"
out$at <- "@"

write.table(out, file = "./gatc_trial/UA_gatc_8_plate.ngsfilter",
            col.names = FALSE, row.names = FALSE, sep = "\t",
            quote = FALSE)
