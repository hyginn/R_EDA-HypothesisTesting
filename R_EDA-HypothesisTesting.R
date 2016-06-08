# R_EDA-HypothesisTesting.R
#
# Purpose:  Introduction to hypothesis testing for biological data.
#
# Version: 1.0
#
# Date:    2016  06  01
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First code
#
# TODO:
#
#
# == HOW TO WORK WITH THIS FILE ======================================
#
#  Go through this script line by line to read and understand the
#  code. Execute code by typing <cmd><enter>. When nothing is
#  selected, that will execute the current line and move the cursor to
#  the next line. You can also select more than one line, e.g. to
#  execute a block of code, or less than one line, e.g. to execute
#  only the core of a nested expression.
#
#  Edit code, as required, experiment with options, or just play.
#  Especially play.
#
#  DO NOT simply source() this whole file!
#
#  If there are portions you don't understand, use R's help system,
#  Google for an answer, or ask me. Don't continue if you don't
#  understand what's going on. That's not how it works ...
#
#  This is YOUR file. Write your notes in here and add as many
#  comments as you need to understand what's going on here when you
#  come back to it in a year. That's the right way: keep code, data
#  and notes all in one place.
#
# ====================================================================


# ==================================================
# Basic statistical tests
# ==================================================

# Two sample t-test:
# In the two-sample t-test, we take multiple
# observations of our samples and ask: is the
# mean of the observations of the samples different
# from each other.

# t-test a single gene of the GSE26922 dataset for
# differential expression. Remember that these are
# triplicates. Therefore our first sample is characterized
# by three replicates at T0 (columns 1:3), and the
# second sample comprises all other columns lumped
# together. Remember: under the null hypothesis
# they should all be the same, therefore we do not
# need to consider the different time points
# separately.

# Let's look at the top and the bottom differentially
# expressed genes in tT as determined by the
# eBayes() function of the limma package?

# No. 1 in tT
tT$ID[1]
ntT1 <- which(rownames(ex) == tT$ID[1])

# No. 250 in tT
tT$ID[nrow(tT)]
ntT250 <- which(rownames(ex) == tT$ID[nrow(tT)])

# retrieve the actual log(expression) values
# from ex and plot them
gtT1 <- ex[ntT1, ]
gtT250 <- ex[ntT250, ]
plot(gtT1,    ylim=c(8,14), type="b", col="darkseagreen")
lines(gtT250,               type="b", col="limegreen")

# perform a sample-against-sample t-test for the cycle
# values against the T0 values
g1 <- t.test(ex[ntT1, 1:3], ex[ntT1, 4:18])
g1

g250 <- t.test(gtT250[1:3], gtT250[4:18])
g250

# What is in the results object?
typeInfo(g1)
g1$p.value       # this is the p-value of H0
# if it's less than 0.5, we
# generally reject H0

# Let's compute p-values for
# all the 32,321 rows of ex
tAll <- numeric(nrow(ex))
for (i in 1:nrow(ex)) {
    tAll[i] <- t.test(ex[i, 1:3], ex[i, 4:18])$p.value
}

hist(tAll, breaks =40)
abline(v=0.05, col="firebrick2", lwd=4)
range(tAll)

# How many rows have a p less than 0.05?
sum(tAll < 0.05)   # Crafty code! Why does this work?
# Hint try: as.numeric(TRUE)
#           as.numeric(FALSE)

gpMin <- which(tAll == min(tAll))
gpMax <- which(tAll == max(tAll))

plot(ex[gpMin, ],   type="b", col="black", ylim=c(8,14))
lines(ex[gpMax, ],  type="b", col="red")
lines(ex[ntT1, ],   type="b", col="darkseagreen")
lines(ex[ntT250, ], type="b", col="limegreen")

# what are the t-values for these four genes?
tVal    <- t.test(ex[ntT1,   1:3], ex[ntT1,   4:18])$statistic
tVal[2] <- t.test(ex[ntT250, 1:3], ex[ntT250, 4:18])$statistic
tVal[3] <- t.test(ex[gpMin,  1:3], ex[gpMin,  4:18])$statistic
tVal[4] <- t.test(ex[gpMax,  1:3], ex[gpMax,  4:18])$statistic
names(tVal) <- c("ntT1", "ntT250", "gpMin", "gpMax")
tVal


# How do the observed p-values relate to the
# t-statistic?
x <- seq(min(tVal)-2, max(tVal)+2, 0.1)
f <- dt(x, df=16)
plot(x, f, xlab="x", ylab="density", type="l", lwd=2)
abline(h = 0, lwd = 0.3)
segments(tVal, -0.01, tVal, 0.05, col="red", lwd=1)
text(tVal, 0.052, labels=names(tVal), adj=c(0, 0.5), srt=90, col="red")

# Note that the underlying assumptions for t-tests
# are oK for microarrays, but less suitable for
# NGS data. See here:
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0123658
# for alternatives.

# =============================================
# Wilcoxon test

# Non-parametric test principle illustrated.

# generate two random data samples with slightly
# different means
set.seed(53)
n <- 25
M <- matrix(nrow = n+n, ncol=2)
for (i in 1:n) {
    M[i,1]   <- rnorm(1, 10, 1)
    M[i,2]   <- 1
    M[i+n,1] <- rnorm(1, 11, 1)
    M[i+n,2] <- 2
}
boxplot(M[1:25,1], M[26:50,1]) # this is how different the two groups are
plot(M[,1], col=M[,2]) # we can see the trend in the values, plotted by index
plot(sample(1:(2*n)), M[,1], col=M[,2]) # however, with random indices it is
# quite hard to tell that these are
# two different distributions ...

# The Wilcoxon test counts for each observation of one
# group how many observations from the other group are
# ranked below it.

# To illustrate:
o <- order(M[,1])
plot(M[o,1], col=M[o,2])

wilcox.test(M[1:n,1], M[(1:n)+n,1])

# TASK:
#   repeat this test for a single random population with
#   mean 10.5.
#   1 - plot with black and red circles as before, study: does it look differnet
#   2 - are the group differences significant
#   3 - repeat this 1000 times and plot the distribution of p-values
#   4 - add the value we just got as an abline()



# Exercise - try this on expression profiles ...
#   - what profiles would you like to check for similarity?
#   - pull the data out of "ex"
#   - should you use all 18 datapoints, or compare the
#     six means of three replicates?
# What is H0 (your null hypothesis) in this case?



# =============================================
# Multiple testing problem

# When we repeat a test that can have a false
# positive outcome, there is a small chance
# that it will indeed turn out falsely positive. When
# we repeat the test many times, it is virtually
# gueranteed that we will have a number - albeit
# small - of false positives. When we test for
# significant differential expression across the
# genome, we are repeating our test 32,000 times.


set.seed(11358)
n <- 32321

synthVals <- as.matrix(rep(0, n*18))
dim(synthVals) <- c(n, 18)
for (i in 1:n) {
    synthVals[i,1:18] <- sample(ex[i, ])
}

# To do: spike in some "real data"!
# To discuss: is sample a good approach?

myt.test <- function(n){
    p <- t.test(synthVals[n,1:3], synthVals[n,4:18])$p.value
    return(p)
}

pVals <- rep(0, n)
for (i in 1:n) {
    pVals[i] <- myt.test(i)
}

hist(pVals, breaks=50)
abline(v=0.05, col="red")


as.numeric(TRUE)  # 1 ... TRUE is coerced into 1
sum(as.numeric(pVals < 0.05)) # ... so summing
# over all TRUE
# values counts
# them.
min(pVals)
rowMin <- which(pVals == min(pVals))
plot(synthVals[rowMin, ],
     type="b",
     col="darkseagreen")

# "Bonferroni correction"
# The solution is to rescale the p-value cutoff
# for significance, by dividing it by the number
# of observations. Our cutoff of p = 0.05
# becomes ...
0.05 / n  # 1.546982e-06
# ... i.e. a very small number instead.

# What we have done is: replaced false positive
# problems with false negatives.

# But let's consider what that means, applied to
# our real data, and our synthetic data:

# our real data ...
sum(as.numeric(tAll < 0.05))       # as is
sum(as.numeric(tAll < 0.05 / n))   # Bonferroni corrected

# the synthetic data ...
sum(as.numeric(pVals < 0.05))      # as is
sum(as.numeric(pVals < 0.05 / n))  # Bonferroni corrected


# =============================================

# FDR
# A modern new concept is to use FDR
# (False Discovery Rate) instead, to control p-Values

# ...

# Let's illustrte this with sample data
set.seed(100)
N <- 10000
alpha <- 0.05
y1 <- matrix(rnorm(9000*4, 0, 1), 9000, 4)
y2 <- matrix(rnorm(1000*4, 5, 1), 1000, 4)
y <- rbind(y1, y2)

myt.test <- function(y){
    t.test(y, alternative="two.sided")$p.value
}
P <- apply(y, 1, myt.test)
sum(P<alpha)

# ...

Psorted <- sort(P)
plot(Psorted, xlim=c(0, 1000), ylim=c(0, .01))
abline(b=alpha/N, a=0, col=2)

p <- p.adjust(P, method="bonferroni")
sum(p<0.05)
p <- p.adjust(P, method="fdr")
sum(p<0.05)

# Calculate the true FDR
sum(p[1:9000]<0.05)/sum(p<0.05)



# [END]
