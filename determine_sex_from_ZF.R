# This script will extract alleles X and Y and call genetic sex. It produces 2-3 files, depending on your wishes.
# Go through comments to see the logic of the way genetic sex is called.
# Input is the folder where _ZF.uniq.tab files (as produced by the last command in the pipeline - obitab).

library(data.table)

#### INPUT ####
# Point this line to the folder where .uniqe.tab files of ZF are located.
xy.file <- list.files("./DAB/zf_hiseq1/", full.names = TRUE)
ngs.files <- list.files("./DAB/1_ngsfilters_hiseq1/", pattern = ".ngsfilter", full.names = TRUE)
# xy.file <- list.files("./DAB/zf_hiseq2/", full.names = TRUE)
# ngs.files <- list.files("./DAB/1_ngsfilters_hiseq2/", pattern = ".ngsfilter", full.names = TRUE)

# Name of the file(s) into which the results are to be written. There is one more 
# possible intermediate file towards the end.
# Switch between HISEQ1 and HiSEQ2 runs.
file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq1.txt"
file.called <- "./DAB/data/dab_called_genetic_sex_hiseq1.txt"
file.freqs <- "./DAB/data/frequency_of_sequences_by_marker_hiseq1_sex.txt"

# file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq2.txt"
# file.called <- "./DAB/data/dab_called_genetic_sex_hiseq2.txt"
# file.freqs <- "./data/frequency_of_sequences_by_marker_hiseq2_sex.txt"
countTS <- 10 # reads that do not match exact sequence should have this number of reads,
              # otherwise they are discarded
#### INPUT ####

# Exact sequences from X and Y chromosomes
seqUAX <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccgtcagagctcaaaaagcacatgcg"
seqUAY <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccatcagagctcaaaaagcacatgcg"

# For each file/sample, extract alleles along with Position, run, sample name, Run_Name and order properly.
system.time(xy <- sapply(xy.file, FUN = function(fl) {
  out.cols <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", 
                "length", "Position", "called", "flag", "stutter", "Sequence")
  # read in the data
  x <- fread(fl)
  
  # find sample names
  samplecols <- names(x)[grepl("sample:", names(x))]
  # make sure total count is correct sum of all samples
  x[, count := Reduce("+", .SD), .SDcols = samplecols] 
  
  # remove zero alleles and unneeded columns
  x <- x[count != 0, ]
  x[, c("count", "id") := NULL]
  names(x)[grepl("sequence", names(x))] <- "Sequence"
  names(x)[grepl("seq_length", names(x))] <- "length"
  
  # reshape the data from wide to long format
  x <- melt(x, id.vars = c("length", "Sequence"), value.name = "Read_Count",
            variable.name = "Sample_Name")
  
  # add information all extra information
  x[, Position := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\2", Sample_Name)]
  x[, Plate := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\3", Sample_Name)]
  x[, Sample_Name := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\1", Sample_Name)]
  x[, Run_Name := gsub("^.*_(DAB\\d+)_.*$", "\\1", basename(fl))]
  x[, Marker := "ZF"]
  x[, called := TRUE]
  x[, flag := ""]
  x[, stutter := FALSE]
  x[, TagCombo := ""]

  x <- x[, out.cols, with = FALSE]
  
  # put everything in proper order
  x <- x[i = order(-Plate, Read_Count, decreasing = TRUE), 
         j = out.cols, with = FALSE, 
         by = .(Plate)]
  x
}, simplify = FALSE))
# Combine the result into one data.table and exclude zero counts.
xy <- rbindlist(xy)
xy <- xy[Read_Count > 0, ]

# Sequences of X and Y, sent in by Stéphane Lobreaux (Université Grenoble Alpes). Use them to replace with
# human readable format.
uaxy <- c(seqUAX, seqUAY)

# Extract only 100% sequences or with sufficient count.
sexy <- xy[Sequence %in% uaxy | Read_Count > countTS, ]

# add TagCombo
tc <- sapply(ngs.files, FUN = read.table, simplify = FALSE)
tc <- rbindlist(tc)
rownames(tc) <- NULL
tc <- as.data.table(tc[, c(1, 2, 3)])
names(tc) <- c("V1", "V2", "TagCombo")

# create columns by which to merge
sexy[, fn := sprintf("%s_%s_PP%s", Sample_Name, Position, Plate)]
sexy[, fl := sprintf(sprintf("UA_MxRout1_%s", Marker))]

sexy <- merge(sexy, tc, by.x = c("fn", "fl"), by.y = c("V2", "V1"))
sexy[, fn := NULL]
sexy[, fl := NULL]

# TODO:pustim length != 104 notr? tiste k so ekstremno kratki?
# TODO: večina bo imela flag D - a to res hočmo?
sexy[Sample_Name %in% sample(sexy$Sample_Name, 1), -"Sequence"]

freq.of.seq.by.marker <- sexy[, .N, by = .(Marker, Sequence)][order(-Marker, N, decreasing = TRUE), .(Marker, N, Sequence)]
fwrite(freq.of.seq.by.marker, file = file.freqs, sep = "\t")
freq.of.seq.by.marker


# ####### od tu dalje nespremenjeno ##############
# 
# # Add human readable designations for X and Y.
# sexy[, seq := ""]
# sexy[Sequence %in% seqUAX, seq := "X"]
# sexy[Sequence %in% seqUAY, seq := "Y"]
# 
# # write intermediate (raw) result
# fwrite(sexy, file = file.out.seq)
# 
# # Raw sequence is of no importance (from here on), so we remove it. It can always be
# # reconstructed from X and Y
# sexy[, Sequence := NULL]
# 
# # If you are skipping the above steps, you can always read in file.out.seq and go from here.
# # sexy <- fread(file.out.seq)
# 
# # Find samples with more than two alleles which may be causing problems. Uncommend fwrite line to
# # output this to a readable file. Some samples may be salveageable, but only manually.
# tm <- sexy[Sample_Name %in% sexy[, .N, by = .(Sample_Name, Plate, Run_Name)][N > 2, Sample_Name], ]
# tm <- tm[order(Sample_Name, Plate, Run_Name), ]
# # fwrite(tm, file = "./DAB/data/samples_with_more_than_two_alleles.txt", sep = "\t")
# 
# # Exclude samples with potential flaws - handle them manually if critical.
# sexy <- sexy[!(Sample_Name %in% unique(tm$Sample_Name)), ]
# 
# # Recast the data so that X, Y and count go from long to wide format.
# # before:        after:
# # seq count         
# #   X   513       X  Y
# #   Y     1     513  1
# # It is fairly easy to read genotypes in this format.
# sexy <- dcast(sexy, Plate + Sample_Name + length + Position + Run_Name ~ seq, 
#       value.var = c("count"))
# 
# # Since some may not have Y (or X), we turn NAs to 0. This and ordering may somewhat
# # improve readability.
# sexy <- sexy[order(Run_Name, Sample_Name, Plate)]
# sexy[, X := ifelse(is.na(X), 0, X)]
# sexy[, Y := ifelse(is.na(Y), 0, Y)]
# 
# # It would appear males have either balanced X and Y or skewed number of reads for
# # Y. Females have very skewed number of reads toward X. We use this piece of
# # information to call genetic sex. Using coefficient of variation we judge whether
# # X or Y is significantly higher than the other (dispersed far from the mean). If 
# # not, it indicates balanced X and Y and not otherwise.
# 
# # Calculate coefficient of variation - this will tell us if the difference
# # between X and Y is "large" or not. If you plot cv, you will notice that it 
# # clusters around zero and above 1. This can be used to define a cut-off.
# sexy[, mean := apply(sexy[, .(X, Y)], 1, mean)]
# sexy[, cv := apply(sexy[, .(X, Y)], 1, sd) / mean]
# 
# # Small values indicate small differences between X and Y --> calling it male.
# cutoff.cv <- 0.75
# # Define a threshold below which absolute number of reads is considered as 
# # untrustworthy and thus forloring genetic sex calling.
# cutoff.abs <- 5
# 
# # Deem male:
# # If coefficient of variation is very small, this indicates X and Y are balanced.
# sexy[cv < cutoff.cv, sex := "M"] # this also takes care of cases where X == Y
# # If Y is amplified more than X, deem this male.
# sexy[Y > mean & X < mean & cv > cutoff.cv, sex := "M"]
# # Deem female:
# # If X is significantly higher than Y, deem this female.
# sexy[X > mean & Y < mean & cv > cutoff.cv, sex := "F"]
# # ID of samples with really weak, reads are deemed unreliable.
# sexy[X < cutoff.abs & Y < cutoff.abs, sex := NA]
# 
# # This just rounds CV to improve readability.
# sexy[, cv := round(cv, 2)]
# # Order before finally writing the result to disk.
# sexy <- sexy[, .(Sample_Name, Plate, Position, Run_Name, mean, cv, X, Y, sex)]
# fwrite(sexy, file = file.called, sep = "\t")
