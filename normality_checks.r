# doing a normality check of the ATP and mtDNA data. Need to know if mtDNA data needs to be transformed.

# import libraries
library(olsrr)
# import data
md = read.csv(file = "Rep_5_and_6_final_values.csv")

# checking the normality of the residuals for ATP by generating a qqplot
jpeg('qq_ATP.jpeg')
model = lm(ATP_pg ~ dev_stage, data = md)
plot.new()
ols_plot_resid_qq(model)
dev.off()

# checking the normality of the residuals for mtDNA by generating a qqplot
jpeg('qq_mtDNA.jpeg')
model = lm(mtDNA ~ dev_stage, data = md)
plot.new()
ols_plot_resid_qq(model)
dev.off()

# doing a test for normality that will give us a shapiro wilk statistic and pvalue
sink(file = "resid_normality_check.txt", append = TRUE, type = c("output"), split = TRUE)
print("ATP x dev_stage")
model = lm(ATP_pg ~ dev_stage, data = md)
print(ols_test_normality(model))
print('------------------------------------------------------------------')
print('------------------------------------------------------------------')
print("mtDNA x dev_stage")
model = lm(mtDNA ~ dev_stage, data = md)
print(ols_test_normality(model))
print('------------------------------------------------------------------')
print('------------------------------------------------------------------')
sink()

# creating a histogram to show the distribution of ATP(pg) and mtDNA values
jpeg('ATP_hist.jpeg')
hist(md$ATP_pg)
dev.off()

jpeg('mtDNA_hist.jpeg')
hist(md$mtDNA)
dev.off()

# creating a qq plot for ATP(pg) x mtDNA
jpeg('qq_ATPxmtDNA.jpeg')
model = lm(ATP_pg ~ mtDNA, data = md)
plot.new()
ols_plot_resid_qq(model)
dev.off()

# also doing a shapiro wilk for ATP(pg) x mtDNA:
sink(file = "resid_normality_check.txt", append = TRUE, type = c("output"), split = TRUE)
print("ATP x mtDNA")
model = lm(ATP_pg ~ mtDNA, data = md)
print(ols_test_normality(model))
print('------------------------------------------------------------------')
print('------------------------------------------------------------------')
sink()

# ATP was found to be normally distributed and mtDNA needed to be natural log (ln) transformed

