biocLite("GenomicFeatures")
library(GenomicFeatures)

fdb_file <- system.file("extdata", "FeatureDb.sqlite",
                        package="GenomicFeatures")
fdb <- loadDb(fdb_file)
fdb
