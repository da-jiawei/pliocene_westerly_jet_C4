rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)

## read data ----
D47 = read.csv("output/D47_d18sw.csv")[,-1]
D47= D47[order(D47$age),] 
Dp17sw = read.csv("output/Dp17sw.csv")
ob.am = read.csv("data/insolation/ob.am.csv")
ob = read_xlsx("data/insolation/Pliocene_orbital_and_insolation_data.xlsx")
sig.am = read.csv("data/insolation/sig.am.csv")


## time series plot ----
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
site1 = pal[factor(D47$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
site2 = pal[factor(Dp17sw$site, levels = c("Lantian", "Shilou", "Jiaxian"))]

png("figures/FigS2.soil_water_isotopes_timeseries_paired_data.png", 4.9, 5.1, units = "in", res = 400)
par(mar = c(3, 4, 0, 4))
plot(0, 0, ylim = c(0, 5), xlim = c(2.5, 7.5), axes = FALSE,
     xlab = "", ylab = "")
axis(1, at = seq(2.5, 7.5, 0.5), cex = 1, mgp = c(1, 0.7, 0))
mtext("Age (Ma)", 1, line = 2)

yext = range(ob$obliquity)
tix = seq(floor(min(yext*10)), 
          ceiling(max(yext*10+1)), by = 5)/10
ob.rs = cbind(ob$age, 
              4 + (ob$obliquity - min(tix)) / diff(range(tix)))
lines(ob.rs[, 1], ob.rs[, 2], col = pal[1])
ob.am.rs = cbind(ob.am$age,
                 4 + (ob.am$ob.am - min(tix)) / diff(range(tix)))
lines(ob.am.rs[, 1], ob.am.rs[, 2], col = pal[2], lwd = 2)
axis(2, 4 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("obliquity", 2, line = 2, at = 4.5)

yext = range(ob$gradient)
tix = seq(ceiling(max(yext)), 
          floor(min(yext-2.9)), by = -5)
sig.rs = cbind(ob$age,
               4 - (ob$gradient - min(tix)) / diff(range(tix)))
lines(sig.rs[, 1], sig.rs[, 2], col = "grey50")
sig.am.rs = cbind(sig.am$age,
                  4 - (sig.am$sig.am - min(tix)) / diff(range(tix)))
lines(sig.am.rs[, 1], sig.am.rs[, 2], col = "grey20", lwd = 2)
axis(4, 4 - (tix - min(tix)) / diff(range(tix)), tix)
mtext("SIG", 4, line = 2, at = 3.5)

D47$temp.low = D47$temp - D47$temp.sd 
D47$temp.high = D47$temp + D47$temp.sd
yext = range(D47$temp.low, D47$temp.high)
tix = seq(floor(min(yext)+1), 
          ceiling(max(yext)+2), by = 5)
D47.rs = cbind(D47$age,
               2 + (D47$temp - min(tix)) / diff(range(tix)),
               2 + (D47$temp.low - min(tix)) / diff(range(tix)),
               2 + (D47$temp.high - min(tix)) / diff(range(tix)))
arrows(D47.rs[, 1], D47.rs[, 3], D47.rs[, 1], D47.rs[, 4], col = "black",
       angle=90, length=0, code = 0)
points(D47.rs[, 1], D47.rs[, 2], col = "black", bg = site1, pch = 21, cex = 1.5)
data = as.data.frame(D47.rs[,1:2]) 
names(data) = c("age", "value")
fit = loess(value ~ age, data = data, span = 0.3)
age = seq(2.5, 7.5, 0.1)
pr = predict(fit, newdata = data.frame(x = age), se = TRUE)
pred = as.data.frame(cbind(age, "mean" = pr$fit)) |> drop_na()
lines(pred$age, pred$mean, lwd = 2, col = "gray")
axis(2, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste("T"[Delta*"47"]*" (", degree, "C)")), 2, line = 2, at = 2.5)

Dp17sw$d18sw.low = Dp17sw$d18sw - Dp17sw$d18sw.se
Dp17sw$d18sw.high = Dp17sw$d18sw + Dp17sw$d18sw.se
yext = range(Dp17sw$d18sw.low, Dp17sw$d18sw.high)
tix = seq(floor(min(yext)),
          ceiling(max(yext)), by = 2)
d18sw.rs = cbind(Dp17sw$age,
                 1 + (Dp17sw$d18sw - min(tix)) / diff(range(tix)),
                 1 + (Dp17sw$d18sw.low - min(tix)) / diff(range(tix)),
                 1 + (Dp17sw$d18sw.high - min(tix)) / diff(range(tix)))
arrows(d18sw.rs[, 1], d18sw.rs[, 3], d18sw.rs[, 1], d18sw.rs[, 4], col = "black",
       angle=90, length=0, code = 0)
points(d18sw.rs[, 1], d18sw.rs[, 2], col = "black", bg = site2, pch = 22, cex = 1.5)
data = as.data.frame(d18sw.rs[,1:2]) 
names(data) = c("age", "value")
fit = loess(value ~ age, data = data, span = 0.3)
age = seq(2.5, 7.5, 0.1)
pr = predict(fit, newdata = data.frame(x = age), se = TRUE)
pred = as.data.frame(cbind(age, "mean" = pr$fit)) |> drop_na()
lines(pred$age, pred$mean, lwd = 2, col = "gray")
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[sw]*" (\u2030)"), 4, line = 2, at = 1.5)

Dp17sw$Dp17sw.low = Dp17sw$Dp17sw - Dp17sw$Dp17sw.se
Dp17sw$Dp17sw.high = Dp17sw$Dp17sw + Dp17sw$Dp17sw.se
yext = range(Dp17sw$Dp17sw.low, Dp17sw$Dp17sw.high)
tix = seq(-110, 70, 20)
Dp17sw.rs = cbind(Dp17sw$age,
                  1 - (Dp17sw$Dp17sw - min(tix)) / diff(range(tix)),
                  1 - (Dp17sw$Dp17sw.low - min(tix)) / diff(range(tix)),
                  1 - (Dp17sw$Dp17sw.high - min(tix)) / diff(range(tix)))
arrows(Dp17sw.rs[, 1], Dp17sw.rs[, 3], Dp17sw.rs[, 1], Dp17sw.rs[, 4], col = "black",
       angle=90, length=0, code = 0)
points(Dp17sw.rs[, 1], Dp17sw.rs[, 2], col = "black", bg = site2, pch = 21, cex = 1.5)
axis(2, 1 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta^"'17"*"O"[sw]*" (per meg)"), 2, line = 2, at = 0.5)

text(5.8, 1.8, "Lantian", col = pal[1], cex = 0.8)
text(6.5, 1.8, "Shilou", col = pal[2], cex = 0.8)
text(7.2, 1.8, "Jiaxian", col = pal[3], cex = 0.8)
text(3, 4.8, "a", font = 2)
text(3, 3.8, "b", font = 2)
text(3, 3.8, "b", font = 2)
text(3, 2.8, "c", font = 2)
text(3, 1.8, "d", font = 2)
text(2.7, .8, "e", font = 2)

dev.off()
