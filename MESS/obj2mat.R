library(R.matlab)
library(spdep)

data <- rgdal::readOGR("./", "baltim.csv")
data_mat <- as.matrix(data@data[, c(1:2, 11:21)])
writeMat("data.mat", x = data_mat[, -3], y = data_mat[, 3])

spt_mat <- listw2mat(nb2listw(poly2nb(data)))
writeMat("spt_mat", spt_mat = spt_mat)

library(data.table)
data <- fread("baltim.csv")
Ws <- readMat("./spt.mat")$Ws
listw_obj <- mat2listw(Ws)
moran.mc(data$PRICE, listw_obj, nsim = 10000)


corrplot::corrplot(cor(data[, c(2, 3, 5, 12, 15)]))
impacts <- matrix(c(0.255865, 0.133945,  0.389809,
                    6.091157, 3.188708,  9.279866,
                    5.643839, 2.954538,  8.598377,
                    8.734409, 4.572444, 13.306853,
                    8.544398, 4.472974, 13.017373,
                    7.578163, 3.967152, 11.545315,
                    4.051991, 2.121209,  6.173200,
                    -6.468940,-3.386476, -9.855416,
                    6.085944, 3.185979,  9.271923,
                    0.057127, 0.029906,  0.087033,
                    10.585379, 5.541423, 16.126802,
                    0.436741, 0.228633,  0.665374),nrow = 12, ncol = 3, byrow = TRUE)
impacts <- as.data.table(imapcts)
colnames(impacts) <- c("Direct", "Indirect", "Total")
impacts$Rnm <- c("NROOM", "DWELL", "NBATH", "PATIO", "FIREPL",
                "AC", "BMENT", "NSTOR", "GAR", "AGE", "CITCOU", "SQFT")
impacts <- melt(impacts, id.vars = "Rnm", measure.vars = c("Direct", "Indirect", "Total"), variable.name = "Impact")
library(ggplot2)
ggplot(impacts, aes(x = reorder(Rnm, -value), y = value, fill = Impact))+
  geom_bar(stat = "identity", position = "dodge")+
  xlab("Variables")+
  ylab("change in House Sales Price")+
  labs(title = "Impacts")+
  coord_flip()+
  theme(axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
