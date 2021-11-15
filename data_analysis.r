
# import libraries
library(SciViews) # required to use ln()
library(ggplot2) # required to make figures
library(lme4)
library(nlme)

# import data
md = read.csv(file = "Rep_5_and_6_final_values.csv")

# creating a new column for the ln_mtDNA values
md$ln_mtDNA = ln(md$mtDNA)

# removing any values that lie outside of our standard curves (ie greater than 4 million or less than 10 mtDNA copy numbers and greater than .5 or less than .0025 uM of ATP)
md2 = data.frame(matrix(ncol = 5, nrow = 0))
colnames(md2) = colnames(md)
for (i in 1:nrow(md)) {
  if (((md$mtDNA[i] < 4000000) & (md$mtDNA[i] > 10)) & ((md$ATP_uM[i] < 0.5) & (md$ATP_uM[i] > 0.0025))) {
    new_row = md[i,]
    md2 = rbind(md2, new_row)
  } else next
}

# row 20 was removed

# redoing the row numbers
rownames(md2) = 1:nrow(md2)

md = md2

sink(file = "post_schneider_analyses/summary_stats.txt", append = TRUE, type = c("output"), split = TRUE)
print("mtDNA x dev_stage")
print(aggregate(mtDNA~dev_stage,md, summary))
aggregate(mtDNA ~ dev_stage, md, function(x) c(mean = mean(x), sd = sd(x)))
print("---")
print("---")
print("---")
print("ATP_pg x dev_stage")
aggregate(ATP_pg~dev_stage,md, summary)
aggregate(ATP_pg ~ dev_stage, md, function(x) c(mean = mean(x), sd = sd(x)))
print("---")
print("---")
print("---")
print("ATP_uM x dev_stage")
aggregate(ATP_uM~dev_stage,md, summary)
aggregate(ATP_uM ~ dev_stage, md, function(x) c(mean = mean(x), sd = sd(x)))
print("---")
print("---")
print("---")
sink()

# performing an ANOVA for ln_mtDNA and ATP_pg x dev_stage and sinking the output to a .txt file
sink(file = "post_schneider_analyses/ANOVA.txt", append = TRUE, type = c("output"), split = TRUE)
print("ln_mtDNA x dev_stage")
mt = aov(ln_mtDNA ~ dev_stage, data = md)
print(mt)
print(summary(mt))
print("---")
print("---")
print("---")
print("ATP_pg x dev_stage")
atp = aov(ATP_pg ~ dev_stage, data = md)
print(atp)
print(summary(atp))
print("---")
print("---")
print("---")
sink()

# running a correlation analyses to see if there is a relationship between ATP and ln_mtDNA
sink(file = "post_schneider_analyses/lm.txt", append = TRUE, type = c("output"), split = TRUE)
print("ATP x ln_mtDNA")
model = lm(ATP_pg ~ ln_mtDNA, data = md)
print(model)
print(summary(model))
print("---")
print("---")
print("---")
# doing ATP by mtdna with stage as a random effect
print("ATP_pg ~ ln_mtDNA with dev stage as random effect")
random_model = lme(ATP_pg ~ ln_mtDNA, random = (~1|dev_stage), data = md)
print(random_model)
print(summary(random_model))
print("---")
print("---")
print("---")
# atp ~ mtdna with stage as a fixed effect
print("ATP_pg ~ ln_mtDNA with dev stage as fixed effect")
fixed_model = lm(ATP_pg ~ ln_mtDNA + dev_stage, data = md)
print(fixed_model)
print(summary(fixed_model))
sink()


# creating the violin plot for ln_mtDNA vs developmental stage
md$dev_stage = as.factor(md$dev_stage)
text_p = "p = 0.035"
jpeg("post_schneider_analyses/ln_mtDNA_violin.jpeg")
    mtDNA_plot = ggplot(md, aes(x=dev_stage, y=ln_mtDNA)) + 
    geom_violin() + 
    geom_point(size = 3) + 
    xlab("Developmental Stage") + 
    ylab("ln(mtDNA Copy Number)") + 
    annotate(geom = 'text', label = text_p, x = 0.5, y = 15.0, hjust = 0, vjust = 1, cex = 7) + 
    ggtitle("A") + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05),
    axis.title = element_text(size = 20, color = "black", face = "bold")) + 
    stat_summary(fun.y=mean, geom="crossbar",width = 0.2, color = "#0072B2")
    print(mtDNA_plot)
dev.off()

# creating a violin plot for ATP vs developmental stage
text_p = "p < 0.001"
jpeg("post_schneider_analyses/ATP_violin.jpeg")
    ATP_plot = ggplot(md, aes(x=dev_stage, y=ATP_pg)) + 
    geom_violin() + 
    geom_point(size = 3) + 
    xlab("Developmental Stage") + 
    ylab("ATP (pg)") + 
    annotate(geom = 'text', label = text_p, x = 0.5, y = 1600, hjust = 0, vjust = 1, cex = 7) + 
    ggtitle("B") + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05),
    axis.title = element_text(size = 20, color = "black", face = "bold")) + 
    stat_summary(fun.y=mean, geom="crossbar",width = 0.2, color = "#0072B2")
    print(ATP_plot)
dev.off()


# creating a plot of the relationship between ATP and ln_mtDNA
lm_p = "p = 0.23"
jpeg("post_schneider_analyses/atpxln_mtDNA.jpeg")
    atpxmtdna_plot = ggplot(md, aes(x= ln_mtDNA, y=ATP_pg, shape = dev_stage, color = dev_stage)) + 
    ggtitle("C") + 
    xlab("ln(mtDNA Copy Number)") + 
    ylab("ATP (pg)") + 
    annotate(geom = 'text', label = lm_p, x = 9.5, y = 1800, hjust = 0, vjust = 1, cex = 7) + 
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
    geom_point(size = 3) + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05), 
    axis.text = element_text(size=20, color = "black"),
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    legend.title = element_blank())
    print(atpxmtdna_plot)
dev.off()