# this code is to do a natural log transformation of the mtDNA dataset and verify that it did look normally distributed after transformation

# import libraries
library(olsrr)
library(SciViews) # required to use ln()
# import data
md = read.csv(file = "Rep_5_and_6_final_values.csv")

# creating a new column for the ln_mtDNA values
md$ln_mtDNA = ln(md$mtDNA)

# creating a histogram to show the distribution of ln_mtDNA values
jpeg('post_schneider_normality/ln_mtDNA_hist.jpeg')
hist(md$ln_mtDNA)
dev.off()

# creating a qqplot for ln_mtDNA
jpeg('post_schneider_normality/qq_ln_mtDNA.jpeg')
model = lm(ln_mtDNA ~ dev_stage, data = md)
plot.new()
ols_plot_resid_qq(model)
dev.off()

# performing a shapiro wilk test for ln_mtDNA
model = lm(ln_mtDNA ~ dev_stage, data = md)
print(ols_test_normality(model))


