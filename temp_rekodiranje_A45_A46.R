library(plyr)
library(data.table)

xy2 <- fread("./DAB/data/dab_hiseq2_genotypes.txt", stringsAsFactors = FALSE,
             colClasses = list(character = c(1, 4, 5, 7, 9, 11, 12),
                               numeric = c(2, 3, 6)))

ab <- xy2[Run_Name == "DAB23", ]
ap18 <- ab[Plate %in% c(1, 8), ]
ap27 <- ab[Plate %in% c(2, 7)]

rc18 <- fread("./DAB/aliquot_plates/HiSEQ_run2/plate_45_recode.csv", sep = ",")
rc27 <- fread("./DAB/aliquot_plates/HiSEQ_run2/plate_46_recode.csv", sep = ",")

ap18$Sample_Name <- mapvalues(x = ap18$Sample_Name, from = rc18$before,
                               to = rc18$after, warn_missing = TRUE)

ap27$Sample_Name <- mapvalues(x = ap27$Sample_Name, from = rc27$before,
                               to = rc27$after, warn_missing = TRUE)
# The following `from` values were not present in `x`: M1Y3A
# M1Y3A se pretvori iz prave F3 v G8 (AK.05AJ), ki pa je v naši poziciji blank

ap18[Position == "034", -c("TagCombo", "Sequence", "called", "flag", "length")]

ap <- rbindlist(list(ap18, ap27))

ap[Position == "034", ]
fwrite(ap, file = "./DAB/data/dab_hiseq2_A45_A46_new_genotypes.txt",
       sep = "\t")
