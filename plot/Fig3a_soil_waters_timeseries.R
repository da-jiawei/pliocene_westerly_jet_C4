rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)

## read data ----
D47 = read.csv("output/D47_d18sw.csv")[,-1]
D47= D47[order(D47$age),] 
d18sw = read.csv("output/D47_d18sw_interpolation.csv")[,-1]
Dp17sw = read.csv("output/Dp17sw.csv")
ob.am = read.csv("data/insolation/ob.am.csv")
ob = read_xlsx("data/insolation/Pliocene_orbital_and_insolation_data.xlsx")
sig.am = read.csv("data/insolation/sig.am.csv")
iron = read_xlsx("data/regional_records/freeiron.xlsx", sheet = "Pianguan")
smi = read_xlsx("data/regional_records/SMI.xlsx") %>%
  mutate(age = age / 1000) %>%
  filter(age > 2 & age < 7.5) %>%
  drop_na(SMI)


## time series plot ----
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
site1 = pal[factor(D47$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
site2 = pal[factor(d18sw$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
site3 = pal[factor(Dp17sw$site, levels = c("Lantian", "Shilou", "Jiaxian"))]

png("figures/Fig3a.soil_water_isotopes_timeseries.png", 9, 4, units = "in", res = 500)
par(mar = c(4, 3, 4, 3))
plot(0, 0, xlim = c(0, 7), ylim = c(2.5, 7.5), axes = FALSE,
     xlab = "", ylab = "")
axis(2, at = seq(2.5, 7.5, 1), cex = 1, mgp = c(1, .7, 0))
mtext("Age (Ma)", 2, line = 2)
rect(xleft = 0, xright = 7, ybottom = 3.3, ytop = 3.8, border = NA, col = "grey90")
rect(xleft = 0, xright = 7, ybottom = 4.7, ytop = 5.3, border = NA, col = "grey90")
rect(xleft = 0, xright = 7, ybottom = 6, ytop = 6.4, border = NA, col = "grey90")
abline(h = 5.32, col = "black", lty = 2)

yext = range(ob$obliquity)
tix = seq(floor(min(yext*10)), 
          ceiling(max(yext*10+1)), by = 5)/10
ob.rs = cbind(0 + (ob$obliquity - min(tix)) / diff(range(tix)),
              ob$age)
lines(ob.rs[, 1], ob.rs[, 2], col = pal[1])
ob.am.rs = cbind(0 + (ob.am$ob.am - min(tix)) / diff(range(tix)),
                 ob.am$age)
lines(ob.am.rs[, 1], ob.am.rs[, 2], col = pal[2], lwd = 2)
axis(3, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("obliquity", 3, line = 2, at = 0.5)

yext = range(ob$gradient)
tix = seq(ceiling(max(yext)), 
          floor(min(yext-2.9)), by = -5)
sig.rs = cbind(2 - (ob$gradient - min(tix)) / diff(range(tix)),
               ob$age)
lines(sig.rs[, 1], sig.rs[, 2], col = "grey50")
sig.am.rs = cbind(2 - (sig.am$sig.am - min(tix)) / diff(range(tix)),
                  sig.am$age)
lines(sig.am.rs[, 1], sig.am.rs[, 2], col = "grey20", lwd = 2)
axis(1, 2 - (tix - min(tix)) / diff(range(tix)), tix)
mtext("SIG", 1, line = 2, at = 1.5)

D47$temp.low = D47$temp - D47$temp.sd 
D47$temp.high = D47$temp + D47$temp.sd
yext = range(D47$temp.low, D47$temp.high)
tix = seq(floor(min(yext)+1), 
          ceiling(max(yext)+2), by = 5)
D47.rs = cbind(D47$age,
               2 + (D47$temp - min(tix)) / diff(range(tix)),
               2 + (D47$temp.low - min(tix)) / diff(range(tix)),
               2 + (D47$temp.high - min(tix)) / diff(range(tix)))
arrows(D47.rs[, 3], D47.rs[, 1], D47.rs[, 4], D47.rs[, 1], col = "black",
       angle=90, length=0, code = 0)
points(D47.rs[, 2], D47.rs[, 1], col = "black", bg = site1, pch = 21, cex = 1.5)
data = as.data.frame(D47.rs[,1:2]) 
names(data) = c("age", "value")
fit = loess(value ~ age, data = data, span = 0.3)
age = seq(2.5, 7.5, 0.1)
pr = predict(fit, newdata = data.frame(x = age), se = TRUE)
pred = as.data.frame(cbind(age, "mean" = pr$fit)) |> drop_na()
lines(pred$mean, pred$age, lwd = 3)
axis(3, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste("T"[Delta*"47"]*" (", degree, "C)")), 3, line = 2, at = 2.5)

d18sw$d18sw.low = d18sw$d18sw - d18sw$d18sw.se
d18sw$d18sw.high = d18sw$d18sw + d18sw$d18sw.se
yext = range(d18sw$d18sw.low, d18sw$d18sw.high)
tix = seq(floor(min(yext)),
          ceiling(max(yext)), by = 2)
d18sw.rs = cbind(d18sw$age,
                3 + (d18sw$d18sw - min(tix)) / diff(range(tix)),
                3 + (d18sw$d18sw.low - min(tix)) / diff(range(tix)),
                3 + (d18sw$d18sw.high - min(tix)) / diff(range(tix)))
arrows(d18sw.rs[, 3], d18sw.rs[, 1], d18sw.rs[, 4], d18sw.rs[, 1], col = "black",
       angle=90, length=0, code = 0)
points(d18sw.rs[, 2], d18sw.rs[, 1], col = "black", bg = site2, pch = 22, cex = 1.2)
axis(1, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[sw]*" (\u2030)"), 1, line = 2, at = 3.5)

Dp17sw$Dp17sw.low = Dp17sw$Dp17sw - Dp17sw$Dp17sw.se
Dp17sw$Dp17sw.high = Dp17sw$Dp17sw + Dp17sw$Dp17sw.se
yext = range(Dp17sw$Dp17sw.low, Dp17sw$Dp17sw.high)
tix = seq(-120, 60, 30)
Dp17sw.rs = cbind(Dp17sw$age,
                5 - (Dp17sw$Dp17sw - min(tix)) / diff(range(tix)),
                5 - (Dp17sw$Dp17sw.low - min(tix)) / diff(range(tix)),
                5 - (Dp17sw$Dp17sw.high - min(tix)) / diff(range(tix)))
arrows(Dp17sw.rs[, 3], Dp17sw.rs[, 1], Dp17sw.rs[, 4], Dp17sw.rs[, 1], col = "black",
       angle=90, length=0, code = 0)
points(Dp17sw.rs[, 2], Dp17sw.rs[, 1], col = "black", bg = site3, pch = 21, cex = 1.5)
abline(v = 5 - (-28 - min(tix)) / diff(range(tix)), col = "black", lty = 2)
axis(3, 5 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta^"'17"*"O"[sw]*" (per meg)"), 3, line = 2, at = 4.5)
yext = range(smi$SMI)
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 0.2)
smi.rs = cbind(5 + (smi$SMI - min(tix)) / diff(range(tix)),
               smi$age)
lines(smi.rs[, 1], smi.rs[, 2], col = pal[1], pch = 21, cex = 1.2)
axis(1, 5 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("SMI"), 1, line = 2, at = 5.5)
yext = range(ob.am$ob.am)
tix = seq(floor(min(yext*10)), 
          ceiling(max(yext*10+3)), by = 5)/10
ob.am.rs = cbind(5 + (ob.am$ob.am - min(tix)) / diff(range(tix)),
                 ob.am$age)
lines(ob.am.rs[, 1], ob.am.rs[, 2], col = "tomato", lwd = 2)

yext = range(iron$Fe)
tix = seq(floor(min(yext * 10)), 
          ceiling(max(yext * 10)), by = 1) /10
iron.rs = cbind(6 + (iron$Fe - min(tix)) / diff(range(tix)),
                iron$age)
lines(iron.rs[, 1], iron.rs[, 2], col = pal[3], pch = 21, cex = 1.2)
axis(3, 6 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("Fe"[2]*"O"[3]*"(f)/Fe"[2]*"O"[3]*"(t)"), 3, line = 2, at = 6.5)
yext = range(ob.am$ob.am)
tix = seq(floor(min(yext*10)), 
          ceiling(max(yext*10+3)), by = 5)/10
ob.am.rs = cbind(6 + (ob.am$ob.am - min(tix)) / diff(range(tix)),
                 ob.am$age)
lines(ob.am.rs[, 1], ob.am.rs[, 2], col = "tomato", lwd = 2)

axis(4, at = seq(2.5, 7.5, 1), cex = 1, mgp = c(1, .7, 0))
mtext("Age (Ma)", 4, line = 2)

text(3, 7.5, "Lantian", col = pal[1])
text(3.5, 7.5, "Shilou", col = pal[2])
text(4, 7.5, "Jiaxian", col = pal[3])
text(0, 5.1, "MPB", col = "black", cex = .8)
text(0, 2.5, "a", cex = 1, col = "black", font = 2)
text(1, 2.5, "b", cex = 1, col = "black", font = 2)
text(2.8, 2.5, "c", cex = 1, col = "black", font = 2)
text(3.2, 2.5, "d", cex = 1, col = "black", font = 2)
text(4.2, 2.5, "e", cex = 1, col = "black", font = 2)
text(5.2, 2.5, "f", cex = 1, col = "black", font = 2)
text(7, 2.5, "g", cex = 1, col = "black", font = 2)

dev.off()
