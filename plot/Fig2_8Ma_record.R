rm(list = ls())
pacman::p_load(tidyverse, readxl)
pal = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")
pal_rc = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")

# load data ----
lingtai = read_csv("data/regional_records/d13_org/processed/lingtai.csv") |> drop_na()
northern_china = read_csv("data/regional_records/d13_org/processed/northern_china.csv")
japan = read_csv("data/regional_records/d13_org/processed/japan_sea.csv") |> filter(age <= 8)

benthic = read_xlsx("data/global_records/LR04.xlsx") |> mutate(age = age / 1e3) |> filter(age <= 8)
co2 = read_csv("data/global_records/100kyrCO2.csv") |> filter(ages <= 8) |>
  mutate(mean = exp(`50%`),
         low = exp(`2.5%`),
         high = exp(`97.5%`)) |>
  select(ages, mean, low, high)

DSST = read_csv("output/DSST.csv") 
DSST_HL = DSST |> drop_na(DTP_HL)
d18sw = read_csv("output/D47_d18sw_interpolation.csv")[,-1]
ob.am = read.csv("data/insolation/ob.am.csv") |> filter(age > 2.8)

# Pliocene - Pleistocene plot ----
# png("figures/Fig5.time_series_8Ma.png", width = 5.3, height = 7.2, units = "in", res = 500)
pdf("figures/Fig2.time_series_8Ma.pdf", width = 5.3, height = 7.2)
par(mar = c(2, 4, 2, 4))
plot(-10, 0, xlim = c(0, 8), ylim = c(0, 7), axes = FALSE,
     xlab = "", ylab = "")
axis(3, at = seq(0, 8, 1), cex = 1, mgp = c(0, -0.3, -1))
mtext("Age (Ma)", 3, line = 1)

yext = range(co2$low, co2$high)
tix = seq(200, 450, 50)
co2.rs = cbind(co2$ages, 
               6 + (co2$mean - min(tix)) / diff(range(tix)),
               6 + (co2$low - min(tix)) / diff(range(tix)),
               6 + (co2$high - min(tix)) / diff(range(tix)))
polygon(c(co2.rs[, 1], rev(co2.rs[, 1])), c(co2.rs[, 3], rev(co2.rs[, 4])), col = pal[3], border = NA)
lines(co2.rs[, 1], co2.rs[, 2], col = pal[1], lwd = 3)
axis(2, 6 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("CO"[2]*" (ppmv)"), 2, line = 2.5, at = 6.5)

yext = range(benthic$d18O)
tix = seq(floor(min(yext*10)), ceiling(max(yext*10)+2), 5) / 10 
benthic.rs = cbind(benthic$age, 6 - (benthic$d18O - min(tix)) / diff(range(tix)))
lines(benthic.rs[, 1], benthic.rs[, 2], col = pal[6], lwd = .7)
axis(4, 6 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("Benthic "*delta^"18"*"O (\u2030)"), 4, line = 2.5, at = 5.5)

tix = seq(0, 1, 0.1)
DSST.rs = cbind(DSST_HL$age, 
                5.2 - (DSST_HL$HL_norm - min(tix)) / diff(range(tix)))
lines(DSST.rs[, 1], DSST.rs[, 2], col = pal[2], lwd = 3)
DSST.rs = cbind(DSST$age, 
                5.2 - (DSST$ML_norm - min(tix)) / diff(range(tix)))
lines(DSST.rs[, 1], DSST.rs[, 2], col = pal[5], lwd = 3)
axis(2, 5.2 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta*"SST"[norm]), 2, line = 2.5, at = 4.7)

yext = range(japan$d13o)
tix = seq(floor(min(yext)), ceiling(max(yext)), 2)
japan.rs = cbind(japan$age, 3 + (japan$d13o - min(tix)) / diff(range(tix)))
lines(japan.rs[, 1], japan.rs[, 2], col = pal[2])
points(japan.rs[, 1], japan.rs[, 2], pch = 21, col = pal[2], bg = "white", cex = .8)
axis(4, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[BC]*" (\u2030)"), 4, line = 2.5, at = 3.5)
text(7.5, 3.9, label = "U1430", cex = 1)

yext = range(northern_china$d13o)
tix = seq(floor(min(yext)), ceiling(max(yext)), 2)
g3.rs = cbind(northern_china$age, 2.3 + (northern_china$d13o - min(tix)) / diff(range(tix)))
lines(g3.rs[, 1], g3.rs[, 2], col = "black")
points(g3.rs[, 1], g3.rs[, 2], pch = 21, bg = pal[3], col = "black", cex = .8)
axis(2, 2.3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[SOM]*" (\u2030)"), 2, line = 2.5, at = 2.8)
text(7.5, 2.8, label = "G3", cex = 1)

yext = range(lingtai$d13)
tix = seq(floor(min(yext)), ceiling(max(yext)), 2)
lingtai.rs = cbind(lingtai$age, 1.3 + (lingtai$d13 - min(tix)) / diff(range(tix)))
lines(lingtai.rs[, 1], lingtai.rs[, 2], col = pal[5])
axis(4, 1.3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[BC]*" (\u2030)"), 4, line = 2.5, at = 1.8)
text(7.5, 2.1, label = "Lingtai", cex = 1)

site = pal_rc[factor(d18sw$site, levels = c("Lantian", "Shilou", "Jiaxian"))]
yext = range(d18sw$d18sw)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 2)
soil_water.rs = cbind(d18sw$age, 0.3 + (d18sw$d18sw - min(tix)) / diff(range(tix)))
points(soil_water.rs[, 1], soil_water.rs[, 2], col = "black", bg = site, pch = 21, cex = .8)
axis(2, 0.3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[sw]*" (\u2030)"), 2, line = 2.5, at = 0.8)

yext = range(ob.am$ob.am)
tix = seq(23.5, 24.4, .2)
ob.rs = cbind(ob.am$age, 0.3 + (ob.am$ob.am - min(tix)) / diff(range(tix)))
lines(ob.rs[, 1], ob.rs[, 2], col = pal[2], lwd = 2)

axis(1, at = seq(0, 8, 1), cex = 1, mgp = c(0, -0.3, -1))
mtext("Age (Ma)", 1, line = 1)

text(.5, 6.8, "a", font = 2)
text(.5, 5.9, "b", font = 2)
text(7.5, 4.5, "c", font = 2)
text(.5, 4, "d", font = 2)
text(.5, 3, "e", font = 2)
text(.5, 1.3, "f", font = 2)
text(.5, 0.3, "g", font = 2)
# text(.5, .5, "h", font = 2)

dev.off()
