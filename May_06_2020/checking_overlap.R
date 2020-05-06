# FUNCTIONS ###########

# split the studies included in from van Klink et al. (2020) and 
# from other recent articles that refer to "insect decline"
# this function extracts only the titles to compare overlap
split_references <- function(refs, order = 1) {
  # make it easy to detect where years appear
  for (i in 0:9) {
    refs <- gsub(i, "X", refs)
  }
  refs <- gsub("\\)", "", gsub("\\(", "", refs))
  refs <- refs[grep("XXXX", refs)]
  refs <- gsub("\\.", "", refs)
  refs <- gsub(", ", "", refs)
  split_refs <- unlist(lapply(refs, function(x) {
    strsplit(strsplit(x, "XXXX ")[[1]][order], "\\.")[[1]][]
  }))
  return(trimws(split_refs))
}

# given the list of EntoGEM hits and those mentioned, checks overlap
# adapted from litsearchr v0.3.0
check_recall <- function (true_hits, retrieved) {
  matches <- lapply(true_hits, synthesisr::fuzzdist, b = retrieved)
  similarity_table <-
    cbind(true_hits, retrieved[unlist(lapply(matches,which.min))], 
          1 - unlist(lapply(matches, min, na.rm = TRUE)))
  colnames(similarity_table) <-
    c("Title", "Best_Match", "Similarity")
  return(similarity_table)
}


# LOAD DATA ###########
# load in EntoGEM screening results for 10+-year studies as of May 06, 2020
load("./entogem_May06_2020.rda")

# remove articles that a third party could not verify were 10+ years
entogem <- entogem$Title[entogem$verified == 1]

# load the list of studies listed in the KMB dataset accompanying van Klink et al. (2020)
load("./vanKlink_included_articles.rda")
vanklink <- split_references(vanklink, order = 2)


load("./insect_decline_references.rda")

references <-
  strsplit(paste(insect_decline$references, collapse = "\\."), "\\.")[[1]]
references <- split_references(references, order = 1)

# clean up some of the random noise introduced by the format of exported references
references <-
  references[-unique(c(grep("Export Date", references), grep("; and", references)))]
references <- unique(references)

# merge together the studies included in van Klink et al. with references 
# from other recent (2019-2020) papers mentioning insect decline
all_references <- unique(tolower(append(vanklink, references)))

# CHECK OVERLAP ###########
checked <- check_recall(tolower(entogem), all_references)

load("./preliminary_new_articles.rda")

# only 6 of 117 10+-year studies identified by EntoGEM so far have been referenced
head(checked[order(checked[, 'Similarity'], decreasing = TRUE), ], 10)
