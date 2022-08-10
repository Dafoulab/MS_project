##############################
##############################
### MS GSM4114296 ANALYSIS ###
##############################
##############################

##############################################################
### WM NAWM AL CA AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
##############################################################

#install.packages("circlize")

library(circlize)

sectors = c("Control WM", "MS NAWM", "MS AL", "MS CA")
fsectors = factor(sectors, levels=c("Control WM", "MS NAWM", "MS AL", "MS CA"))

circos.par("clock.wise" = TRUE, start.degree = 145, gap.after = c("MS NAWM" = 0, "Control WM" = 0, "MS AL" = 0, "MS CA" = 57))

circos.initialize(fsectors, xlim = c(0, 6))

circos.track(fsectors, ylim = c(0, 10), track.height = mm_h(5))

### Type of Editing Sites
circos.track(ylim = c(0,17), track.height = mm_h(22))

circos.barplot(sector.index = "Control WM", value = c("intronic" = 2275/160, "3'5'UTR" = 563/160, "downstream" = 64/160, "upstream" = 6/160, "lncRNA" = 346/160, "exonic" = 4/160), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "MS NAWM", value = c("intronic" = 2152/160, "3'5'UTR" = 644/160, "downstream" = 72/160, "upstream" = 16/160, "lncRNA" = 326/160, "exonic" = 4/160), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "MS AL", value = c("intronic" = 1355/160, "3'5'UTR" = 989/160, "downstream" = 119/160, "upstream" = 10/160, "lncRNA" = 394/160, "exonic" = 3/160), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "MS CA", value = c("intronic" = 1193/160, "3'5'UTR" = 452/160, "downstream" = 78/160, "upstream" = 16/160, "lncRNA" = 306/160, "exonic" = 3/160), 1:6 - 0.5, track.index = 2, col = 1:6)

circos.yaxis(side = "left", at = c(0, 500/160, 1000/160, 1500/160, 2000/160, 2500/160), labels = c("0","500", "1000", "1500","2000","2500"), sector.index = "Control WM", track.index = 2, col = "red", labels.cex = 0.6)

circos.text(x=0.5,y=16.5, sector.index = "Control WM", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)
circos.text(x=0.5,y=15.5, sector.index = "MS NAWM", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)
circos.text(x=0.5,y=13.5, sector.index = "MS AL", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)
circos.text(x=0.5,y=11.5, sector.index = "MS CA", "intronic", facing = 'bending.inside', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=1.5,y=10.5, sector.index = "Control WM", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=1.5,y=10.5, sector.index = "MS NAWM", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=1.5,y=12, sector.index = "MS AL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=1.5,y=10, sector.index = "MS CA", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=2.5,y=9.3, sector.index = "Control WM", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=2.5,y=9.3, sector.index = "MS NAWM", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=2.5,y=9.3, sector.index = "MS AL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=2.5,y=9.3, sector.index = "MS CA", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=3.5,y=8.5, sector.index = "Control WM", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=3.5,y=8.5, sector.index = "MS NAWM", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=3.5,y=8.5, sector.index = "MS AL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=3.5,y=8.5, sector.index = "MS CA", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=4.5,y=9, sector.index = "Control WM", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=4.5,y=9, sector.index = "MS NAWM", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=4.5,y=9, sector.index = "MS AL", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=4.5,y=9, sector.index = "MS CA", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=5.5,y=8, sector.index = "Control WM", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=5.5,y=8, sector.index = "MS NAWM", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=5.5,y=8, sector.index = "MS AL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=5.5,y=8, sector.index = "MS CA", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

### Frequency of Editing Sites
circos.track(ylim = c(0,10), track.height = mm_h(22))

circos.barplot(sector.index = "Control WM", value = c("intronic" = (2275/3258)*10, "3'5'UTR" = (563/3258)*10, "downstream" = (64/3258)*10, "upstream" = (6/3258)*10, "lncRNA" = (346/3258)*10, "exonic" = (4/3258)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "MS NAWM", value = c("intronic" = (2152/3214)*10, "3'5'UTR" = (644/3214)*10, "downstream" = (72/3214)*10, "upstream" = (16/3214)*10, "lncRNA" = (326/3214)*10, "exonic" = (4/3214)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "MS AL", value = c("intronic" = (1355/2874)*10, "3'5'UTR" = (989/2874)*10, "downstream" = (119/2874)*10, "upstream" = (10/2874)*10, "lncRNA" = (394/2874)*10, "exonic" = (3/2874)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "MS CA", value = c("intronic" = (1193/2051)*10, "3'5'UTR" = (452/2051)*10, "downstream" = (78/2051)*10, "upstream" = (16/2051)*10, "lncRNA" = (306/2051)*10, "exonic" = (3/2051)*10), 1:6 - 0.5, track.index = 3, col = 1:6)

circos.text(x=0.5,y=8.5, sector.index = "Control WM", "69.82%", facing = 'bending.inside', track.index = 3,cex = 0.6)
circos.text(x=0.5,y=8.5, sector.index = "MS NAWM", "66.95%", facing = 'bending.inside', track.index = 3,cex = 0.6)
circos.text(x=0.5,y=6.5, sector.index = "MS AL", "47.14%", facing = 'bending.inside', track.index = 3,cex = 0.6)
circos.text(x=0.5,y=7.5, sector.index = "MS CA", "58.16%", facing = 'bending.inside', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=1.5,y=6, sector.index = "Control WM", "17.28%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=1.5,y=6, sector.index = "MS NAWM", "20.03%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=1.5,y=7, sector.index = "MS AL", "34.41%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=1.5,y=6, sector.index = "MS CA", "22.03%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=2.5,y=5, sector.index = "Control WM", "1.96%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=2.5,y=5, sector.index = "MS NAWM", "2.24%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=2.5,y=5, sector.index = "MS AL", "4.14%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=2.5,y=5, sector.index = "MS CA", "3.80%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=3.5,y=5, sector.index = "Control WM", "0.28%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=3.5,y=5, sector.index = "MS NAWM", "0.47%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=3.5,y=5, sector.index = "MS AL", "0.55%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=3.5,y=5, sector.index = "MS CA", "0.49%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=4.5,y=5, sector.index = "Control WM", "3.05%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=4.5,y=5, sector.index = "MS NAWM", "2.94%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=4.5,y=5, sector.index = "MS AL", "3.41%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=4.5,y=5, sector.index = "MS CA", "3.57%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=5.5,y=5, sector.index = "Control WM", "3.21%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=5.5,y=5, sector.index = "MS NAWM", "2.80%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=5.5,y=5, sector.index = "MS AL", "3.57%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=5.5,y=5, sector.index = "MS CA", "3.91%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

### Number of Editing Sites
circos.track(ylim = c(0,17), track.height = mm_h(18))

circos.barplot(sector.index = "Control WM", value = 3258/220 , 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "MS NAWM", value = 3214/220, 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "MS AL", value = 2874/220, 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "MS CA", value = 2051/220, 3, track.index = 4, col = 7, bar_width = 6)

circos.text(x=3,y=10, sector.index = "Control WM", "3258", facing = 'bending.inside', track.index = 4,cex = 0.7)
circos.text(x=3,y=10, sector.index = "MS NAWM", "3214", facing = 'bending.inside', track.index = 4,cex = 0.7)
circos.text(x=3,y=9, sector.index = "MS AL", "2874", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.7)
circos.text(x=3,y=6, sector.index = "MS CA", "2051", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.7)

### Add Header Sector Names
circos.text(x=3,y=5, sector.index = "Control WM", "Control WM", facing = 'bending.inside', track.index = 1)
circos.text(x=3,y=5, sector.index = "MS NAWM", "MS NAWM", facing = 'bending.inside', track.index = 1)
circos.text(x=3,y=5, sector.index = "MS AL", "MS AL", facing = 'bending.inside', track.index = 1)
circos.text(x=3,y=5, sector.index = "MS CA", "MS CA", facing = 'bending.inside', track.index = 1)

circos.clear()


#FOOTNOTE: 
#The distribution, percentage and number of A-to-I RNA editing events across MS sample types.







######################################
### EAE NAIVE PRECLINICAL CLINICAL ###
######################################

#install.packages("circlize")

library(circlize)

sectors = c("EAE NAIVE", "EAE PRECLINICAL", "EAE CLINICAL")
fsectors = factor(sectors, levels=c("EAE NAIVE", "EAE PRECLINICAL", "EAE CLINICAL"))

circos.par("clock.wise" = TRUE, start.degree = 150, gap.after = c("EAE NAIVE" = 0, "EAE PRECLINICAL" = 0, "EAE CLINICAL" = 70))

circos.initialize(fsectors, xlim = c(0, 6))


circos.track(fsectors, ylim = c(0, 10), track.height = mm_h(5))

### Type of Editing Sites
circos.track(ylim = c(0,17), track.height = mm_h(22))

circos.barplot(sector.index = "EAE NAIVE", value = c("intronic" = 408/30, "3'5'UTR" = 185/30, "downstream" = 32/30, "upstream" = 4/30, "intergenic" = 0/30, "exonic" = 146/30), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "EAE PRECLINICAL", value = c("intronic" = 399/30, "3'5'UTR" = 180/30, "downstream" = 17/30, "upstream" = 3/30, "intergenic" = 105/30, "exonic" = 144/30), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "EAE CLINICAL", value = c("intronic" = 177/30, "3'5'UTR" = 222/30, "downstream" = 19/30, "upstream" = 4/30, "intergenic" = 0/30, "exonic" = 124/30), 1:6 - 0.5, track.index = 2, col = 1:6)

circos.yaxis(side = "left", at = c(0, 100/30, 200/30, 300/30, 400/30, 500/30), labels = c("0","100", "200", "300","400","500"), sector.index = "EAE NAIVE", track.index = 2, col = "red", labels.cex = 0.6)

circos.text(x=0.5,y=16, sector.index = "EAE NAIVE", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)
circos.text(x=0.5,y=16, sector.index = "EAE PRECLINICAL", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)
circos.text(x=0.5,y=9.5, sector.index = "EAE CLINICAL", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.6)

circos.text(x=1.5,y=12, sector.index = "EAE NAIVE", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=1.5,y=12, sector.index = "EAE PRECLINICAL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=1.5,y=13, sector.index = "EAE CLINICAL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=2.5,y=9.3, sector.index = "EAE NAIVE", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=2.5,y=9.3, sector.index = "EAE PRECLINICAL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=2.5,y=9.3, sector.index = "EAE CLINICAL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=3.5,y=8, sector.index = "EAE NAIVE", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=3.5,y=8, sector.index = "EAE PRECLINICAL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=3.5,y=8, sector.index = "EAE CLINICAL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=4.5,y=9, sector.index = "EAE NAIVE", "intergenic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=4.5,y=11, sector.index = "EAE PRECLINICAL", "intergenic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=4.5,y=9, sector.index = "EAE CLINICAL", "intergenic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

circos.text(x=5.5,y=10, sector.index = "EAE NAIVE", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=5.5,y=10, sector.index = "EAE PRECLINICAL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)
circos.text(x=5.5,y=10, sector.index = "EAE CLINICAL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.6)

### Frequency of Editing Sites
circos.track(ylim = c(0,10), track.height = mm_h(18))


circos.barplot(sector.index = "EAE NAIVE", value = c("intronic" = (408/775)*10, "3'5'UTR" = (185/775)*10, "downstream" = (32/775)*10, "upstream" = (4/775)*10, "intergenic" = (0/775)*10, "exonic" = (146/775)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "EAE PRECLINICAL", value = c("intronic" = (399/743)*10, "3'5'UTR" = (180/743)*10, "downstream" = (17/743)*10, "upstream" = (3/743)*10, "intergenic" = (105/743)*10, "exonic" = (144/743)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "EAE CLINICAL", value = c("intronic" = (177/546)*10, "3'5'UTR" = (222/546)*10, "downstream" = (19/546)*10, "upstream" = (4/546)*10, "intergenic" = (0/546)*10, "exonic" = (124/546)*10), 1:6 - 0.5, track.index = 3, col = 1:6)

circos.text(x=0.5,y=8, sector.index = "EAE NAIVE", "52.64%", facing = 'bending.inside', track.index = 3,cex = 0.6)
circos.text(x=0.5,y=8, sector.index = "EAE PRECLINICAL", "53.70%", facing = 'bending.inside', track.index = 3,cex = 0.6)
circos.text(x=0.5,y=6, sector.index = "EAE CLINICAL", "32.41%", facing = 'bending.inside', track.index = 3,cex = 0.6)

circos.text(x=1.5,y=6.5, sector.index = "EAE NAIVE", "23.87%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=1.5,y=6.5, sector.index = "EAE PRECLINICAL", "24.22%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=1.5,y=7.5, sector.index = "EAE CLINICAL", "40.65%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=2.5,y=5, sector.index = "EAE NAIVE", "4.12%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=2.5,y=5, sector.index = "EAE PRECLINICAL", "2.28%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=2.5,y=5, sector.index = "EAE CLINICAL", "3.47%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=3.5,y=5, sector.index = "EAE NAIVE", "0.51%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=3.5,y=5, sector.index = "EAE PRECLINICAL", "0.40%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=3.5,y=5, sector.index = "EAE CLINICAL", "0.73%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=4.5,y=5, sector.index = "EAE NAIVE", "0.00%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=4.5,y=6, sector.index = "EAE PRECLINICAL", "14.13%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=4.5,y=5, sector.index = "EAE CLINICAL", "0.00%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

circos.text(x=5.5,y=6, sector.index = "EAE NAIVE", "18.83%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=5.5,y=6, sector.index = "EAE PRECLINICAL", "19.38%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)
circos.text(x=5.5,y=6, sector.index = "EAE CLINICAL", "22.71%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.6)

### Number of Editing Sites
circos.track(ylim = c(0,17), track.height = mm_h(15))

circos.barplot(sector.index = "EAE NAIVE", value = 775/50, 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "EAE PRECLINICAL", value = 294/50 , 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "EAE CLINICAL", value = 177/50, 3, track.index = 4, col = 7, bar_width = 6)

circos.text(x=3,y=9, sector.index = "EAE NAIVE", "775", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)
circos.text(x=3,y=12, sector.index = "EAE PRECLINICAL", "294", facing = 'bending.inside', track.index = 4,cex = 0.8)
circos.text(x=3,y=10, sector.index = "EAE CLINICAL", "177", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)

### Add Header Sector Names
circos.text(x=3,y=5, sector.index = "EAE NAIVE", "EAE NAIVE", facing = 'bending.inside', track.index = 1)
circos.text(x=3,y=5, sector.index = "EAE PRECLINICAL", "EAE PRECLINICAL", facing = 'bending.inside', track.index = 1)
circos.text(x=3,y=5, sector.index = "EAE CLINICAL", "EAE CLINICAL", facing = 'bending.inside', track.index = 1)


circos.clear()


#FOOTNOTE: 
#The distribution, percentage and number of A-to-I RNA editing events across EAE sample types.










######################################################################################
###### MSvsEAE SITES FROM COMMON GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ######
######################################################################################

#install.packages("circlize")

library(circlize)

sectors = c("WMvsNAIVE", "NAWMvsPRECLINICAL", "ALvsCLINICAL", "CAvsCLINICAL")
fsectors = factor(sectors, levels = c("WMvsNAIVE","NAWMvsPRECLINICAL", "ALvsCLINICAL", "CAvsCLINICAL"))

circos.par("clock.wise" = TRUE, start.degree = 145, gap.after = c("NAWMvsPRECLINICAL" = 0, "WMvsNAIVE" = 0, "ALvsCLINICAL" = 0, "CAvsCLINICAL" = 57))

circos.initialize(fsectors, xlim = c(0, 6))


circos.track(fsectors, ylim = c(0, 10), track.height = mm_h(5))

### Type of Editing Sites
circos.track(ylim = c(0,14), track.height = mm_h(22))

circos.barplot(sector.index = "WMvsNAIVE", value = c("intronic" = 259/25, "3'5'UTR" = 48/25, "downstream" = 0/25, "upstream" = 3/25, "lncRNA" = 1/25, "exonic" = 0/25), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "NAWMvsPRECLINICAL", value = c("intronic" = 169/25, "3'5'UTR" = 62/25, "downstream" = 0/25, "upstream" = 0/25, "lncRNA" = 7/25, "exonic" = 0/25), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "ALvsCLINICAL", value = c("intronic" = 32/25, "3'5'UTR" = 54/25, "downstream" = 1/25, "upstream" = 0/25, "lncRNA" = 0/25, "exonic" = 0/25), 1:6 - 0.5, track.index = 2, col = 1:6)
circos.barplot(sector.index = "CAvsCLINICAL", value = c("intronic" = 34/25, "3'5'UTR" = 31/25, "downstream" = 0/25, "upstream" = 0/25, "lncRNA" = 0/25, "exonic" = 0/25), 1:6 - 0.5, track.index = 2, col = 1:6)

circos.yaxis(side = "left", at = c(0, 100/25, 150/25, 200/25, 250/25, 300/25), labels = c("0","100", "150", "200","250","300"), sector.index = "WMvsNAIVE", track.index = 2, col = "red", labels.cex = 0.6)

circos.text(x=0.5,y=13, sector.index = "WMvsNAIVE", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.7)
circos.text(x=0.5,y=10, sector.index = "NAWMvsPRECLINICAL", "intronic", facing = 'bending.inside', track.index = 2,cex = 0.7)
circos.text(x=0.5,y=7, sector.index = "ALvsCLINICAL", "intronic", facing = 'clockwise', niceFacing = FALSE, track.index = 2,cex = 0.7)
circos.text(x=0.5,y=7, sector.index = "CAvsCLINICAL", "intronic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

circos.text(x=1.5,y=8, sector.index = "WMvsNAIVE", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=1.5,y=8, sector.index = "NAWMvsPRECLINICAL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=1.5,y=8, sector.index = "ALvsCLINICAL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=1.5,y=8, sector.index = "CAvsCLINICAL", "3'5'UTR", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

circos.text(x=2.5,y=7.8, sector.index = "WMvsNAIVE", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=2.5,y=7.8, sector.index = "NAWMvsPRECLINICAL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=2.5,y=7.8, sector.index = "ALvsCLINICAL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=2.5,y=7.8, sector.index = "CAvsCLINICAL", "downstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

circos.text(x=3.5,y=8, sector.index = "WMvsNAIVE", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=3.5,y=8, sector.index = "NAWMvsPRECLINICAL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=3.5,y=8, sector.index = "ALvsCLINICAL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=3.5,y=8, sector.index = "CAvsCLINICAL", "upstream", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

circos.text(x=4.5,y=8, sector.index = "WMvsNAIVE", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=4.5,y=8, sector.index = "NAWMvsPRECLINICAL", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=4.5,y=8, sector.index = "ALvsCLINICAL", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=4.5,y=8, sector.index = "CAvsCLINICAL", "lncRNA", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

circos.text(x=5.5,y=8, sector.index = "WMvsNAIVE", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=5.5,y=8, sector.index = "NAWMvsPRECLINICAL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=5.5,y=8, sector.index = "ALvsCLINICAL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)
circos.text(x=5.5,y=8, sector.index = "CAvsCLINICAL", "exonic", facing = 'clockwise', niceFacing = TRUE, track.index = 2,cex = 0.7)

### Frequency of Editing Sites
circos.track(ylim = c(0,11), track.height = mm_h(22))

circos.barplot(sector.index = "WMvsNAIVE", value = c("intronic" = (259/311)*10, "3'5'UTR" = (48/311)*10, "downstream" = (0/311)*10, "upstream" = (3/311)*10, "lncRNA" = (1/311)*10, "exonic" = (0/311)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "NAWMvsPRECLINICAL", value = c("intronic" = (169/238)*10, "3'5'UTR" = (62/238)*10, "downstream" = (0/238)*10, "upstream" = (0/238)*10, "lncRNA" = (7/238)*10, "exonic" = (0/238)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "ALvsCLINICAL", value = c("intronic" = (32/87)*10, "3'5'UTR" = (54/87)*10, "downstream" = (1/87)*10, "upstream" = (0/87)*10, "lncRNA" = (0/87)*10, "exonic" = (0/87)*10), 1:6 - 0.5, track.index = 3, col = 1:6)
circos.barplot(sector.index = "CAvsCLINICAL", value = c("intronic" = (34/65)*10, "3'5'UTR" = (31/65)*10, "downstream" = (0/65)*10, "upstream" = (0/65)*10, "lncRNA" = (0/65)*10, "exonic" = (0/65)*10), 1:6 - 0.5, track.index = 3, col = 1:6)

circos.text(x=0.5,y=10, sector.index = "WMvsNAIVE", "83.27%", facing = 'bending.inside', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=0.5,y=9, sector.index = "NAWMvsPRECLINICAL", "71.00%", facing = 'bending.inside', track.index = 3,cex = 0.7)
circos.text(x=0.5,y=8, sector.index = "ALvsCLINICAL", "36.78%", facing = 'clockwise', track.index = 3,cex = 0.7)
circos.text(x=0.5,y=8.5, sector.index = "CAvsCLINICAL", "52.30%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

circos.text(x=1.5,y=7, sector.index = "WMvsNAIVE", "15.43%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=1.5,y=7, sector.index = "NAWMvsPRECLINICAL", "26.05%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=1.5,y=9, sector.index = "ALvsCLINICAL", "62.06%", facing = 'bending.inside', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=1.5,y=8.2, sector.index = "CAvsCLINICAL", "47.69%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

circos.text(x=2.5,y=6, sector.index = "WMvsNAIVE", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=2.5,y=6, sector.index = "NAWMvsPRECLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=2.5,y=6, sector.index = "ALvsCLINICAL", "1.14%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=2.5,y=6, sector.index = "CAvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

circos.text(x=3.5,y=6, sector.index = "WMvsNAIVE", "0.96%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=3.5,y=6, sector.index = "NAWMvsPRECLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=3.5,y=6, sector.index = "ALvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=3.5,y=6, sector.index = "CAvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

circos.text(x=4.5,y=6, sector.index = "WMvsNAIVE", "0.32%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=4.5,y=6, sector.index = "NAWMvsPRECLINICAL", "2.94%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=4.5,y=6, sector.index = "ALvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=4.5,y=6, sector.index = "CAvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

circos.text(x=5.5,y=6, sector.index = "WMvsNAIVE", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=5.5,y=6, sector.index = "NAWMvsPRECLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=5.5,y=6, sector.index = "ALvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)
circos.text(x=5.5,y=6, sector.index = "CAvsCLINICAL", "0%", facing = 'clockwise', niceFacing = TRUE, track.index = 3,cex = 0.7)

### Number of Editing Sites
circos.track(ylim = c(0,14), track.height = mm_h(18))

circos.barplot(sector.index = "WMvsNAIVE", value = 311/25 , 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "NAWMvsPRECLINICAL", value = 238/25, 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "ALvsCLINICAL", value = 87/25, 3, track.index = 4, col = 7, bar_width = 6)
circos.barplot(sector.index = "CAvsCLINICAL", value = 65/25, 3, track.index = 4, col = 7, bar_width = 6)

circos.text(x=3,y=7, sector.index = "WMvsNAIVE", "311", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)
circos.text(x=3,y=5, sector.index = "NAWMvsPRECLINICAL", "238", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)
circos.text(x=3,y=9, sector.index = "ALvsCLINICAL", "87", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)
circos.text(x=3,y=8, sector.index = "CAvsCLINICAL", "65", facing = 'bending.inside', niceFacing = TRUE, track.index = 4,cex = 0.8)

### Add Header Sector Names
circos.text(x=3,y=28, sector.index = "WMvsNAIVE", "WMvsNAIVE Common Targets", facing = 'bending.inside', adj = c(0.5, 2.6), track.index = 1)
circos.text(x=3,y=28, sector.index = "NAWMvsPRECLINICAL", "NAWMvsPRECLINICAL Common Targets", facing = 'bending.inside', adj = c(0.5, 2.6), track.index = 1)
circos.text(x=3,y=28, sector.index = "ALvsCLINICAL", "ALvsCLINICAL Common Targets", facing = 'bending.inside', adj = c(0.5, 2.6), track.index = 1)
circos.text(x=3,y=28, sector.index = "CAvsCLINICAL", "CAvsCLINICAL Common Targets", facing = 'bending.inside', adj = c(0.5, 2.6), track.index = 1)



circos.clear()


#FOOTNOTE: 
#The distribution, percentage and number of common A-to-I RNA editing events across EAE and MS sample types.






