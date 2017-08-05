sexy <- fread("./DAB/data/dab_genetic_sex_hiseq1.txt")

### learn of ratio significance
lrn <- fread("./DAB/barkode_spol.txt")
# find samples where sex is known
known.sex <- sexy[sexy$Sample_Name %in% lrn$Sample_Name, ]
# merge with genetic data
known.sex <- merge(known.sex, lrn, by = "Sample_Name")
# known.sex <- known.sex[order(Sample_Name, run, sex)]
# reshape data to better fit our "wide" model
known.sex <- dcast(known.sex, run + Sample_Name + seq_length + position + library + truesex ~ sex, 
                   value.var = c("count"))
# # sort so that X and Y are consistently placed
known.sex <- known.sex[order(Sample_Name, run), ]
# replace missing values with zeros - useful later
known.sex[, X := ifelse(is.na(X), 0, X)]
known.sex[, Y := ifelse(is.na(Y), 0, Y)]
# calculate ratio and round to two digits
known.sex[, ratio := round(X/Y, 2)]

known.sex[, pval := chisq.test(c(X, Y), p = c(0.5, 0.5))$p.value, by = 1:nrow(known.sex)]





# write (intermediate) result to file
fwrite(known.sex, file = "./DAB/genetic_sex_compare_to_true.txt")

ggplotly(ggplot(known.sex, aes(x = truesex, y = ratio, group = Sample_Name)) +
           theme_bw() +
           geom_jitter()
)

chisq.test(c(370, 353), p = c(0.5, 0.5))
chisq.test(c(700, 0), p = c(0.999, 0.001))
chisq.test(c(1200, 0), p = c(0.999, 0.001))
chisq.test(c(199, 155), p = c(0.5, 0.5))$p.value




