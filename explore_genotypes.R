library(data.table)

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
xy2[Sequence == "", ]

# najdi negativne kontrole
fwrite(xy2[grepl("(^AC\\..*$|^AK\\..*$)", Sample_Name), ],
       file = "negkonthiseq2.txt", sep = "\t")
# najdi pozitivne kontrole
xy2[grepl("POS", Sample_Name), ]
# najdi prazne
xy2[grepl("blank", Sample_Name), ]
