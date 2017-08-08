# negativna kontrola: vodiš jo že od ekstrakcije, le, da se notri ne doda vzorca (se da pa pufre, segreva...)
# pozitivna kontrola: vzorci, ki so nam znani in so imeli na klasičnem sekvenatorju kvaliteten (neinvaziven) genotip

library(data.table)
library(ggplot2)
library(reshape2)
source("custom_functions.R")
# hiseq run #1 (17.7.2017)
xy1 <- fread("./DAB/data/dab_hiseq1_genotypes.txt", stringsAsFactors = FALSE,
             colClasses = list(character = c(1, 4, 5, 7, 9, 11, 12),
                               numeric = c(2, 3, 6)))

# hiseq run #2 (25.7.2017)
xy2 <- fread("./DAB/data/dab_hiseq2_genotypes.txt", stringsAsFactors = FALSE,
             colClasses = list(character = c(1, 4, 5, 7, 9, 11, 12),
                               numeric = c(2, 3, 6)))
# number of unique samples
length(unique(xy1[, Sample_Name]))

# are there any blanks?
xy1[grepl("blank", Sample_Name), ]

# ali je PCR za plate 3, 4, 5, 6 uspel (za vseh 24 knjižnic)?
# hiseq1
sapply(sort(unique(xy1$Run_Name)), FUN = function(i, x){
  table(x[Run_Name == i & Plate %in% c(3, 4, 5, 6), ][order(Sample_Name, Marker, Plate, length)][, Plate])
}, x = xy1)
sapply(sort(unique(xy1$Run_Name)), FUN = function(i, x){
  table(x[Run_Name == i & Plate %in% c(1, 2, 7, 8), ][order(Sample_Name, Marker, Plate, length)][, Plate])
}, x = xy1)

# hiseq2
sapply(sort(unique(xy2$Run_Name)), FUN = function(i, x){
  table(x[Run_Name == i & Plate %in% c(3, 4, 5, 6), ][order(Sample_Name, Marker, Plate, length)][, Plate])
}, x = xy2)
sapply(sort(unique(xy2$Run_Name)), FUN = function(i, x){
  table(x[Run_Name == i & Plate %in% c(1, 2, 7, 8), ][order(Sample_Name, Marker, Plate, length)][, Plate])
}, x = xy2)


fwrite(xy2[Run_Name %in% c("DAB14") & Plate %in% c(3, 4, 5, 6), ][order(Sample_Name, Marker, -length, Sequence)],
       file = "dab14.txt", sep = "\t")

# poglejmo vzorce iz ponovitev 3, ki imajo več kot 2 alela. če ima več, to pomeni, da je prišlo do
# prilivanja iz drugih vzorcev
fwrite(xy2[Sample_Name %in% xy2[Run_Name == "DAB13" & Plate %in% c(3, 4, 5, 6), ][, .N, by = .(Sample_Name, Marker, Plate)][N == 16, Sample_Name], ],
       file = "16.txt", sep = "\t")

xy2[Sequence == "aacttaccaacaaactaatctatctatctatctatctatctatctatctatctatctatctatctatctatctatctacatatatg", ]
xy2[Sequence == "aacttaccaacaaactaatctatctatctatctatctatctatctatctatctatctatctatctatctatctatctacatatata", ]
xy2[Sequence == "aacttaccaacaaactaatctatctatctatctatctatctatctatctatctatctatctatcgatctatctatctatctatctatctacatatata", ]

xy2[Sequence == "", ]
xy2[Sequence == "", ]

# najdi negativne kontrole
# AK negativne pri tkivnih
# AC negativne pri neinvazivnih
fwrite(xy2[grepl("(^AC\\..*$|^AK\\..*$)", Sample_Name), ][order(Sequence, Run_Name, Sample_Name, Plate, Marker)],
       file = "negkonthiseq2.txt", sep = "\t")
fwrite(xy1[grepl("^NeKo.*", Sample_Name), ][order(Sequence, Run_Name, Sample_Name, Plate, Marker)],
       file = "negkonthiseq1.txt", sep = "\t")

fwrite(xy2[Sequence == "aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctatctatctc", ],
       file = "kontaminacija1.txt", sep = "\t")
# najdi pozitivne kontrole
# to so na H12 (096)
# plate 45, 46 ni blankov in pozitivnih kontrol, na H12 je negativna kontrola
# pri tkivnih
fwrite(xy2[Position == "096", ][order(Marker, Sequence, Sample_Name, Run_Name), ],
       file = "pozitivnekontrole.txt", sep = "\t")

# katere pozitivne kontrole izstopajo?
# ali so plate, na katerih PK izstopajo (npr. majhno število readov), podobne tem PK (ali so "čudne")? na primer,
# ali je cela plata slabo PCR-irana/sekvenirana

# For each library and for each plate, calculate sum (or mean, see code) number of reads per position (sample).
# Result: For each library, 8 plates are printed and number of reads (and sample name) displayed.
#### hiseq1
mc <- xy1
mc[, iscontrol := FALSE]
mc[grepl("(^AC\\..*$|^AK\\..*$|NeKo.*$)", Sample_Name), iscontrol := TRUE]
# find how many reads negative control has compared to the rest
mc[, .(mean.count = sum(Read_Count)), by = .(iscontrol)]

showSumsByLibrary(xy1)

#### hiseq2
mc <- xy2
mc[, iscontrol := FALSE]
mc[grepl("(^AC\\..*$|^AK\\..*$)", Sample_Name), iscontrol := TRUE]
mc[, .(mean.count = sum(Read_Count)), by = .(iscontrol)]

showSumsByLibrary(xy2)

# najdi prazne
xy2[grepl("blank", Sample_Name), ]

