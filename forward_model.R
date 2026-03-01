rm(list = ls())

dat = read.csv("output/Dp17sw.csv")
dat$D47c = 0.0391e6 / (dat$temp + 273.15) ^ 2 + 0.154

R18smow = 0.0020052
R17smow = 0.0003799
R18vpdb = 0.0020672
R17vpdb = 0.0003860
alpha18_c_w_eq = exp((1.61e4 / (dat$temp + 273.15) - 24.6) / 1e3) # Wostbrock (2020)
theta_c_w = 0.5305 - 1.39 / (dat$temp + 273.15)
alpha17_c_w_eq = alpha18_c_w_eq ^ theta_c_w
R18s = (dat$d18sw / 1e3 + 1) * R18smow
dp17sw = dat$Dp17sw / 1e3 + 0.528 * log(dat$d18sw/1e3 + 1) * 1e3
R17s = exp(dp17sw/1e3) * R17smow
R18c = R18s * alpha18_c_w_eq
R17c = R17s * alpha17_c_w_eq
d18c = (R18c / R18vpdb - 1) * 1e3
dp18c = log(R18c / R18smow) * 1e3
dp17c = log(R17c / R17smow) * 1e3
Dp17c = (dp17c - 0.528 * dp18c) * 1e3 # per meg
D47c = 0.0391e6 / Tsoil.K ^ 2 + 0.154 # Andersen (2021)