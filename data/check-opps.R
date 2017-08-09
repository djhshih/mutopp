x <- read.table("opps.tsv", sep="\t", header=TRUE);

# T002.1 is the same as T001.1 but on the reverse stand
# all mutation opportunities should be the same
stopifnot(x[1, -(1:2)] == x[3, -(1:2)]);

# T002.2 is the same as T001.2 but on the reverse stand
stopifnot(x[2, -(1:2)] == x[4, -(1:2)]);

