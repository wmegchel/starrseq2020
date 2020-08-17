####################################################################################
#
# Riverplot of the P53-C1, P53-C2 classes in 2iL and SL
#
# Wout Megchelenbrink
# January 22, 2020
#
####################################################################################
library(riverplot)
source("include/style.R")

# ALL
nodes <- c("A", "B", "C", "D")
edges <- list( A=list(C= 2251, D=249 ), B=list( C= 134, D=1302))
r <- makeRiver( nodes, edges, node_xpos= c( 1,1,2,2 ),
                node_labels= c( A= "C1", B= "C2", C= "C1s", D= "C2s" ),
                node_styles= list( A= list( col=colors.basic[2] ), B= list( col=colors.basic[1] ), C= list( col=colors.basic[2]), D= list( col= colors.basic[1] )))


pdf(file = "IMG/Fig4C_P53_sanky_C1_C2_ALL.pdf", width = 3, height = 4)
riverplot(r, fix.pdf = T)
dev.off()
