#R

context("Test 8 reps data")

test_that("8rep sanity check", {

	x <- msqc1_8rep
	x <- x[grepl("[by]", x$Fragment.Ion) & x$Peptide.Sequence %in% msqc1::peptides$Peptide.Sequence, ]
	x.sum <- aggregate(Area ~ instrument + Isotope.Label.Type + relative.amount + Peptide.Sequence + File.Name.Id, data=x, FUN=sum)

	x.sum.cv <-  aggregate(Area ~ instrument + Isotope.Label.Type + Peptide.Sequence, data = x.sum, FUN = function(x){100 * sd(x) / mean(x)})


	expect_true( table(x$relative.amount) == nrow(x) )
	expect_true( sum(unique(x$Peptide.Sequence) %in% peptides$Peptide.Sequence) == 14)
	expect_true( abs(0.95 * 5 * 2 * 14) <= nrow(x.sum.cv))
})