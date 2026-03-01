# This script calculates the d18O and Dp17O of soil waters, with errors propagated using Monte Carlo random sampling. 
rm(list = ls())
pacman::p_load(tidyverse, patchwork, readxl)
nsyth = 1e4
set.seed(42)

# load data ----
d18c = read_xlsx("data/stable_isotope.xlsx")[1:289,]
names(d18c) = c("site", "age", "d18c", "d18c.se")
d18c = d18c |>
  mutate(d18c.se = as.numeric(d18c.se)) |>
  mutate(d18c.se = ifelse(is.nan(d18c.se), 0.25, d18c.se))

D47c = read_xlsx("data/clumped_isotope.xlsx")
colnames(D47c)[3] = "age"
Dp17c = read_xlsx("data/triple_oxygen_isotope.xlsx")
colnames(Dp17c)[3] = "age"
Dp17c = Dp17c |>
  mutate(Dp17c.se = as.numeric(Dp17c.se)) |>
  mutate(Dp17c.se = ifelse(is.nan(Dp17c.se), 18, Dp17c.se))

# calculate soil water isotope using paired data ----
# soil water d18O
for (i in 1:nrow(D47c)) {
  d18c.vpdb = rnorm(nsyth, D47c$d18[i], D47c$d18.se[i])
  st = rnorm(nsyth, D47c$temp[i], D47c$temp.sd[i])
  d18c.vsmow = (1.03091 * d18c.vpdb) + 30.91
  alpha = exp((16.1*1000/(273.15 + st)-24.6)/1000) # Tremaine 2011
  d18sw = (1000 + d18c.vsmow) / alpha - 1000
  D47c$d18sw[i] = mean(d18sw)
  D47c$d18sw.se[i] = sd(d18sw)
}

ggplot(D47c, aes(x = age, y = d18sw, fill = site)) +
  geom_point(size = 3, shape = 21)
write.csv(D47c, "output/D47_d18sw.csv")

# soil water Dp17O
Dp17c = left_join(Dp17c, D47c[, c("age", "site", "D47", "D47.se", "temp", "temp.sd")], by = c("age", "site"))
for (i in 1:nrow(Dp17c)) {
  d18c.vpdb = rnorm(nsyth, Dp17c$d18c[i], Dp17c$d18c.se[i])
  Dp17c.smow = rnorm(nsyth, Dp17c$Dp17c[i], Dp17c$Dp17c.se[i])
  st = rnorm(nsyth, Dp17c$temp[i], Dp17c$temp.sd[i])
  
  # R18c = (d18c.vpdb / 1e3 + 1) * 0.0020672
  # dp17c.vpdb = Dp17c.vpdb / 1e3 + 0.528 * log(d18c.vpdb / 1e3 + 1) * 1e3
  # R17c = exp(dp17c.vpdb / 1e3) * 0.0003860
  
  d18c.smow = (d18c.vpdb + 29.98) / 0.97002
  dp18c.smow = 1e3 * log(d18c.smow/1e3 + 1)
  dp17c.smow = Dp17c.smow/1e3 + 0.528 * dp18c.smow
  
  alpha18_c_w_eq = exp((1.61e4 / (st + 273.15) - 24.6) / 1e3) # Wostbrock (2020) 
  theta_c_w = 0.5305 - 1.39 / (st + 273.15)
  alpha17_c_w_eq = alpha18_c_w_eq ^ theta_c_w 
  
  d18sw = ((d18c.smow/1e3 + 1) / alpha18_c_w_eq - 1) * 1e3 
  dp18sw = 1e3 * log(d18sw/1e3 + 1)
  dp17sw = log(exp(dp17c.smow/1e3) / alpha17_c_w_eq) * 1e3
  Dp17sw = (dp17sw - 0.528 * dp18sw) * 1e3
  
  Dp17c$d18sw[i] = mean(d18sw) 
  Dp17c$d18sw.se[i] = sd(d18sw)
  Dp17c$Dp17sw[i] = mean(Dp17sw) 
  Dp17c$Dp17sw.se[i] = sd(Dp17sw)
}

ggplot(Dp17c, aes(x = age, y = Dp17sw, fill = site)) +
  geom_point(size = 3, shape = 21)
write.csv(Dp17c, "output/Dp17sw.csv")

# calculate soil water isotope using interpolated data ----
d18c$temp = approx(D47c$age, D47c$temp, xout = d18c$age)$y
d18c$temp.sd = approx(D47c$age, D47c$temp.sd, xout = d18c$age)$y
d18c = d18c |> drop_na()

for (i in 1:nrow(d18c)) {
  d18c.vpdb = rnorm(nsyth, d18c$d18c[i], d18c$d18c.se[i])
  st = rnorm(nsyth, d18c$temp[i], d18c$temp.sd[i])
  d18c.vsmow = (1.03091 * d18c.vpdb) + 30.91
  alpha = exp((16.1*1000/(273.15 + st)-24.6)/1000) # Tremaine 2011
  d18sw = (1000 + d18c.vsmow) / alpha - 1000
  d18c$d18sw[i] = mean(d18sw)
  d18c$d18sw.se[i] = sd(d18sw)
}

ggplot(d18c, aes(x = age, y = d18sw, fill = site)) +
  geom_point(size = 3, shape = 21)
write.csv(d18c, "output/D47_d18sw_interpolation.csv")

