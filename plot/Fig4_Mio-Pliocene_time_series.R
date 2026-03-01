rm(list = ls())
pacman::p_load(tidyverse, readxl)
pal = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")
pal_rc = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")

# load data ----
jiaxian = read.csv("data/regional_records/d13_org/processed/jiaxian.csv")|> filter(age <= 8 & age >= 2.5)
lingtai = read.csv("data/regional_records/d13_org/processed/lingtai.csv") |> drop_na() |> filter(age <= 8 & age >= 2.5)
northern_china = read.csv("data/regional_records/d13_org/processed/northern_china.csv") |> filter(age <= 8 & age >= 2.5)
japan = read.csv("data/regional_records/d13_org/processed/japan_sea.csv") |> filter(age <= 8 & age >= 2.5)

benthic = read_xlsx("data/global_records/LR04.xlsx") |> mutate(age = age / 1e3) |> filter(age <= 8 & age >= 2.5)
co2 = read.csv("data/global_records/100kyrCO2.csv") |> filter(ages <= 8 & ages >= 2.5) |>
  mutate(mean = exp(X50.),
         low = exp(X2.5.),
         high = exp(X97.5.)) |>
  select(ages, mean, low, high)
DSST = read.csv("output/DSST.csv") |> filter(age <= 8 & age >= 2.5)
d18sw = read.csv("output/D47_d18sw_interpolation.csv")[,-1] |> filter(age <= 8 & age >= 2.5)
D47 = read.csv("output/D47_d18sw.csv") |> filter(age <= 8 & age >= 2.5) |> 
  select(site, age, temp, temp.sd) |>
  mutate(low = temp - temp.sd,
         high = temp + temp.sd)
dust_1208 = read_csv("data/regional_records/d13_org/processed/dust_1208.csv")
dust_885 = read_csv("data/regional_records/d13_org/processed/dust_885.csv")

ob.am = read.csv("data/insolation/ob.am.csv") |> filter(age <= 8 & age >= 2.5)

# plot 
pdf("figures/Fig4.Miocene_Pliocene_time_series.pdf", 11, 5)
par(mar = c(4, 3, 4, 3))
plot(0, 0, ylim = rev(c(2.5, 8)), xlim = c(0, 10), axes = FALSE,
     xlab = "", ylab = "")
axis(2, at = seq(2.5, 8, 0.5), cex = 1, mgp = c(1, .7, 0))
mtext("Age (Ma)", 2, line = 2)

xext = range(benthic$d18O)
tix = seq(min(xext*10-2.3), max(xext*10+0.7), 5) /10
benthic.rs = cbind(1 - (benthic$d18O - min(tix)) / diff(range(tix)), benthic$age)
lines(benthic.rs[, 1], benthic.rs[, 2], col = pal[1])
axis(3, 1 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("Benthic "*delta^"18"*"O (\u2030)"), 3, line = 2, at = 0.5)

xext = range(japan$d13o)
tix = seq(floor(min(xext)), ceiling(max(xext)), 1)
japan.rs = cbind(1 + (japan$d13o - min(tix)) / diff(range(tix)), japan$age)
lines(japan.rs[, 1], japan.rs[, 2], col = pal[2])
points(japan.rs[, 1], japan.rs[, 2], pch = 21, col = pal[2], bg = "white", cex = .8)
axis(1, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[BC]*" (\u2030)"), 1, line = 2.5, at = 1.5)
text(1.3, 2.5, label = "U1430", cex = .9)

nc = northern_china |> filter(age >= 3.5)
xext = range(nc$d13o) 
tix = seq(-26, -22, 1)
g3.rs = cbind(2 + (nc$d13o - min(tix)) / diff(range(tix)), nc$age)
lines(g3.rs[, 1], g3.rs[, 2], col = "black")
points(g3.rs[, 1], g3.rs[, 2], pch = 21, bg = pal[3], col = "black", cex = .8)
axis(3, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[SOM]*" (\u2030)"), 3, line = 2, at = 2.5)
text(2.7, 2.5, label = "G3", cex = .9)

xext = range(jiaxian$d13o)
tix = seq(floor(min(xext)), ceiling(max(xext)), 1)
jiaxian.rs = cbind(3 + (jiaxian$d13o - min(tix)) / diff(range(tix)), jiaxian$age)
lines(jiaxian.rs[, 1], jiaxian.rs[, 2], col = pal[5])
points(jiaxian.rs[, 1], jiaxian.rs[, 2], pch = 21, bg = pal[5], cex = .8)
axis(1, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[OOM]*" (\u2030)"), 1, line = 2.5, at = 3.5)
text(3.5, 7.8, label = "Jiaxian", cex = .9)

xext = range(lingtai$d13)
tix = seq(floor(min(xext)), ceiling(max(xext)), 2)
lingtai.rs = cbind(4 + (lingtai$d13 - min(tix)) / diff(range(tix)), lingtai$age)
lines(lingtai.rs[, 1], lingtai.rs[, 2], col = pal[5])
axis(3, 4 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[BC]*" (\u2030)"), 3, line = 2, at = 4.5)
text(4.5, 7.8, label = "Lingtai", cex = .9)

xext = range(dust_1208$flux)
tix = seq(0, 2.5, 0.5)
dust.rs = cbind(6 - (dust_1208$flux - min(tix)) / diff(range(tix)), dust_1208$age)
lines(dust.rs[, 1], dust.rs[, 2], col = pal[1])
axis(1, 6 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("flux (g/cm"^"2"*"/kyr)"), 1, line = 2.5, at = 5.5)
text(5.5, 4.7, "ODP 1208", cex = .9)

xext = range(dust_885$flux)
tix = seq(0, 0.25, 0.05)
dust.rs = cbind(7 - (dust_885$flux - min(tix)) / diff(range(tix)), dust_885$age)
lines(dust.rs[, 1], dust.rs[, 2], col = pal[1])
axis(3, 7 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("flux (g/cm"^"2"*"/kyr)"), 3, line = 2, at = 6.5)
text(6.5, 4.7, "ODP 885", cex = .9)

site = pal_rc[factor(d18sw$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
xext = range(d18sw$d18sw)
tix = seq(floor(min(xext)), ceiling(max(xext)), by = 2)
soil_water.rs = cbind(7 + (d18sw$d18sw - min(tix)) / diff(range(tix)), d18sw$age)
points(soil_water.rs[, 1], soil_water.rs[, 2], col = "black", bg = site, pch = 21, cex = .8)
axis(1, 7 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[sw]*" (\u2030)"), 1, line = 2.5, at = 7.5)

xext = range(ob.am$ob.am)
tix = seq(23.5, 24.4, .2)
ob.rs = cbind(7 + (ob.am$ob.am - min(tix)) / diff(range(tix)), ob.am$age)
lines(ob.rs[, 1], ob.rs[, 2], col = pal[2], lwd = 2)

site = pal_rc[factor(D47$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
xext = range(D47$low, D47$high)
tix = seq(floor(min(xext)+1), ceiling(max(xext)), by = 5)
D47.rs = cbind(D47$age, 
               8 + (D47$temp - min(tix)) / diff(range(tix)),
               8 + (D47$low - min(tix)) / diff(range(tix)),
               8 + (D47$high - min(tix)) / diff(range(tix)))
arrows(D47.rs[, 3], D47.rs[, 1], D47.rs[, 4], D47.rs[, 1], col = "black",
       angle=90, length=0, code = 0)
points(D47.rs[, 2], D47.rs[, 1], col = "black", bg = site, pch = 22, cex = 1)
axis(3, 8 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste("T"[Delta*"47"]*" (", degree, "C)")), 3, line = 2, at = 8.5)

xext = range(ob.am$ob.am)
tix = seq(23.5, 24.4, .2)
ob.rs = cbind(8 + (ob.am$ob.am - min(tix)) / diff(range(tix)), ob.am$age)
lines(ob.rs[, 1], ob.rs[, 2], col = pal[2], lwd = 2)

DSST_pm = DSST |> drop_na()
xext = range(DSST_pm$DTP_HL)
tix = seq(floor(min(xext)), ceiling(max(xext)), 1)
DSST.rs = cbind(10 - (DSST_pm$DTP_HL - min(tix)) / diff(range(tix)), DSST_pm$age)
lines(DSST.rs[, 1], DSST.rs[, 2], col = pal[6], lwd = 3)
axis(1, 10 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste(Delta*"SST (", degree, "C)")), 1, line = 2.5, at = 9.5)

# xext = range(co2$low, co2$high)
# tix = seq(230, 410, 20)
# co2.rs = cbind(co2$ages,
#                10 + (co2$mean - min(tix)) / diff(range(tix)),
#                10 + (co2$low - min(tix)) / diff(range(tix)),
#                10 + (co2$high - min(tix)) / diff(range(tix)))
# polygon(c(co2.rs[, 3], rev(co2.rs[, 4])), c(co2.rs[, 1], rev(co2.rs[, 1])), col = pal[3], border = NA)
# lines(co2.rs[, 2], co2.rs[, 1], col = pal[1], lwd = 3)
# axis(4, 0 + (tix - min(tix)) / diff(range(tix)), tix)
# mtext(expression("CO"[2]*" (ppmv)"), 4, line = 2.5, at = 0.5)

text(0.3, 7.5, label = "a", font = 2)
text(1.6, 7.5, "b", font = 2)
text(2.5, 7.5, "c", font = 2)
text(3.5, 7.5, "d", font = 2)
text(4.5, 7.5, "e", font = 2)
text(5.5, 7.5, "f", font = 2)
text(6.5, 7.5, "g", font = 2)
text(7.5, 7.5, "h", font = 2)
text(8.5, 7.5, "i", font = 2)
text(9.3, 7.5, "j", font = 2)
axis(4, at = seq(2.5, 8, 0.5), cex = 1, mgp = c(1, .7, 0))
mtext("Age (Ma)", 4, line = 2)

dev.off()
