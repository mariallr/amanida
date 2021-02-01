
# Exemple 

filename <- "~/OneDrive - URV/metaanalysis_pack/dataset2.xlsx"

datafile <- data.read(filename, coln = c("Compound Name", "P-value", "Fold-change", "N total", "Behaviour",  "References"))

res.met <- metmet(datafile)

write.csv(res.met@down, "~/OneDrive - URV/Review celia/down.csv")

metaplot(res.met, method = "onebyone", cutoff = c(0.05, 4))

voteplot(res.met)

filename <- "~/OneDrive - URV/Review celia/data_prepost.xlsx"

datafile <- data.read(filename, coln = c("Common Name", "p-value", "fold-change", "N", "trend", "Reference"))

res.met <- metmet(datafile)

metaplot(res.met, method = "alltogether", cutoff = c(0.05, 4))

write.csv(res.met@up, "~/OneDrive - URV/Review celia/up_all.csv")

filename <- "~/OneDrive - URV/Review celia/volcano_all.xlsx"

datafile <- data.read(filename, coln = c("Compound Name", "P-value", "Fold-change", "N total", "Behaviour",  "References"))

res.met <- metmet(datafile)

metaplot(res.met, method = "alltogether", cutoff = c(0.05, 4))