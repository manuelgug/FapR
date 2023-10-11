
##validation of HFS22_01's DD2 phasing results

phased_haplos <- read.csv("phased_haplos.csv")

#formatting
dd2 <- phased_haplos[grepl("DD2", phased_haplos$SampleID, ignore.case = TRUE), ]
dd2$haplotype <- paste(dd2$dhps_437, dd2$dhps_540, dd2$dhfr_51, dd2$dhfr_59, dd2$dhfr_108, sep = "_")


# sort rows by DD2 concentration in the mix
library(stringr)

numeric_values <- as.numeric(str_extract(dd2$SampleID, "(?<=k13_)[0-9]+(?=_S)"))
index <- order(numeric_values)
sorted_dd2 <- dd2[index, ]


#plot
library(ggplot2)

# factor with the desired order based on SampleID
sorted_dd2$SampleID <- factor(sorted_dd2$SampleID, levels = unique(sorted_dd2$SampleID))

# stacked barplot
a<-ggplot(sorted_dd2, aes(x = SampleID, y = HAPLO_FREQ_RECALC, fill = haplotype)) +
  geom_bar(stat = "identity")+
  labs(title = "Haplotype Frequencies on 3D7-DD2 mixes",
       x = "SampleID",
       y = "Haplo Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  guides(fill = guide_legend(title = "Haplotype"))

ggsave("haplo_profiles_DD2.png", a, width = 12, height = 9) 


#conclusions so far:
# 1) dfhr108 is less sensitive, creating "middle" haplos that doesn't exist. removing dhfr108 leaves the gradient as expected.
# 2) thus, a filtering step may be needed given this limits of detection if all resmarkers should be included in order to return true haplos. 
# 3) dhfr108 may be taken into account when found in high freq tho (between 25 to 60% according to stacked barplots)
# 4) if dhfr108 is monoallelic (TENDS TO BE THAT WAY APPARENTLY! [haplo_profile_multiallelic.png]), all good with the rest.
# 5) rest of markers' limit of detection seem to be around 1-2%
