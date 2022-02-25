### Adjacency matrices for triads of interest

complete <- matrix(
    c(0, 1, 1,
      1, 0, 1,
      1, 1, 0),
    byrow = TRUE, nrow = 3
)

feedback <- matrix(
    c(0, 1, 0,
      0, 0, 1,
      1, 0, 0),
    byrow = TRUE, nrow = 3
)

feedforward <- matrix(
    c(0, 1, 1,
      0, 0, 1,
      0, 0, 0),
    byrow = TRUE, nrow = 3
)

fromthemiddleout <- matrix(
    c(0, 0, 0,
      1, 0, 1,
      0, 0, 0),
    byrow = TRUE, nrow = 3
)

linemotif <- matrix(
    c(0, 1, 0,
      0, 0, 1,
      0, 0, 0),
    byrow = TRUE, nrow = 3
)

tothemiddle <- matrix(
    c(0, 1, 0,
      0, 0, 0,
      0, 1, 0),
    byrow = TRUE, nrow = 3
)

### A vector of names of empirical data sets stored in this package
empiricals <- c("mine", "jpr", "sn_auth", "pira", "karate", "jazz",
                "weaverbirds", "linux", "taro", "wiring",
                "dolphins", "netsci")
