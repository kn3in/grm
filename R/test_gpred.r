test_me <- function(ground_truth, my_code) {
	gt <- read.table(ground_truth)
	mc <- read.table(my_code)


}

gt <- read.table("../bin/test1_out_gpredictor.txt")
mc <- read.table("../data/test1_predictor_k")
all.equal(gt[ ,1], mc[ ,1])
all.equal(gt[ ,2], mc[ ,2])
all.equal(gt[ ,3], mc[ ,3])
plot(gt[ ,3], mc[ ,3])