# Glossary:
# PP = primer plate, combination of tags and primers. Tags are unique for each well. Eight plates per library.
# AP = aliquot plate, each well holds its own sample.

# Make sure you change locations of files to read according to the plates you're using.

library(readxl)

dir.AP <- "./DAB/aliquot_plates/HiSEQ_run2/" # folder where aliquote plates are located

# 1. Load PP and AP data
primers <- as.data.frame(read_excel("./DAB/aliquot_plates/input_primers_tags.xlsx", sheet = "primers"))
PP <- as.data.frame(read_excel("./DAB/aliquot_plates/UrsusNGSPrimersCrosbreeding.Ljubljana.18Jul2017.xlsx", 
                               sheet = "Tag crossbreeding", skip = 42))
# PP has columns position, slo, PP1, PP2, ... PP8
# specify folder where to look for aliquot plates
AP <- data.frame(location = list.files(dir.AP, pattern = "_DAB_A\\d+\\.xls$", full.names = TRUE))
AP$name <- gsub("^.*(A\\d+)\\.xls$", "\\1", AP$location)

# 2. Find which AP is added to which PP.
pa.loc <- as.data.frame(read_excel("./DAB/aliquot_plates/NGS.Plates.Barcodes.DAB2017.Final.xlsx", 
                                   sheet = "PCR_Plates_4Reps"))

# select which libraries you wish to run through this script
pa.loc <- droplevels(pa.loc[pa.loc$Library_BC %in% sprintf("DAB%02d", 13:24, sep = ""), ])

pa.loc <- split(pa.loc, f = pa.loc$Library_BC)

# For library, split by plate...
out <- sapply(pa.loc, FUN = function(x, PP, AP, primers) {
  x.split <- split(x, f = 1:nrow(x))
  
  # ... and for each plate, construct NGS filter
  out <- sapply(x.split, FUN = function(y, PP, AP, primers) {
    # First find primers for PP in library
    find.pp <- y[, "Primer Plate"]
    PP.x <- na.omit(PP[, names(PP) %in% find.pp, drop = FALSE]) # columns plate and tagcombo
    stopifnot(ncol(PP.x) == 1)
    
    find.ap <- y[, "Aliquot Plate"]
    AP.x <- AP[AP$name %in% find.ap, ] # columns location and name
    AP.x <- read_excel(as.character(AP.x$location))
    
    pos.number <- sprintf("%03d", rep(1:nrow(AP.x), each = nrow(primers)))
    sample.pos.plate <- rep(AP.x$SPositionBC, each = nrow(primers))
    sample.pos.plate <- paste(sample.pos.plate, pos.number, find.pp, sep = "_")
    
    locus <- rep(primers$locusname, times = nrow(AP.x))
    
    tagcombo <- rep(PP.x[, 1], each = nrow(primers))
    
    primer1 <- rep(primers$primer1, times = nrow(PP.x))
    primer2 <- rep(primers$primer2, times = nrow(PP.x))
    
    out <- data.frame(locus, sample.pos.plate, tagcombo, primer1, primer2, ef = "F", at = "@")
  }, PP = PP, AP = AP, primers = primers, simplify = FALSE)
  
  out <- do.call(rbind, out)
  write.table(out, file = sprintf("%s.ngsfilter", unique(x$Library_BC)), row.names = FALSE,
              col.names = FALSE, quote = FALSE, sep = "\t", fileEncoding = "UTF-8")
  out
}, PP = PP, AP = AP, primers = primers, simplify = FALSE)
