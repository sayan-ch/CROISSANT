library(readr)
library(Matrix)


### DBLP
edge_ac <- readr::read_csv("raw_data/dblp-ac-compact.csv")
edge_ap <- readr::read_csv("raw_data/dblp-ap-compact.csv")
author_label <- read_csv("raw_data/author_label.csv")

A_ac <- Matrix::sparseMatrix(i = edge_ac$V1,
                  j = edge_ac$a,
                  dims = c(4057, 4057))

A_ac <- A_ac | t(A_ac)

A_ap <- sparseMatrix(i = edge_ap$V1,
                     j = edge_ap$a,
                     dims = c(4057, 4057))
A_ap <- A_ap | t(A_ap)

DBLP <- list(A.AC = A_ac,
             A.AP = A_ap,
             label = author_label
             )

save(DBLP, file = "data/DBLP.rda")

####################################
### Twitch
edge <- read_csv("raw_data/twitch_edges.csv")
comm <- read_csv("raw_data/twitch_comm.csv")

A <- sparseMatrix(i = edge$user1, j = edge$user2,
                  dims = c(32407,32407))

Twitch <- list(A = A,
               label = comm)

save(Twitch, file = "data/Twitch.rda")
