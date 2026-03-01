rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

# load data ----
rc = read.csv("output/Dp17sw.csv")[,-1]
rc$site = factor(rc$site, levels = c("Lantian", "Shilou", "Jiaxian"))
rc = rc |>
  mutate(dp18sw = 1e3 * log(d18sw/1e3 + 1),
         dp18sw.low = 1e3 * log((d18sw - d18sw.se)/1e3 + 1),
         dp18sw.high = 1e3 * log((d18sw + d18sw.se)/1e3 + 1))
mw = read_xlsx("data/isotope_records/ModernWater_China_tian2019.xlsx")
sw = read_xlsx("data/isotope_records/Kelson2023.xlsx", sheet = "data")
bayes_data = read_csv("output/BASS_bayes_vary_evap_ver2.csv")
cat("\014")

# all modern waters ----
p1 = ggplot(rc, aes(x = dp18sw, y = Dp17sw)) +
  geom_point(data = mw, aes(x = dp18O, y = Dp17O), color = "grey", size = 3, shape = 21) +
  geom_errorbar(aes(ymin = Dp17sw - Dp17sw.se, ymax = Dp17sw + Dp17sw.se), linewidth = 0.2) +
  geom_errorbar(aes(xmin = dp18sw.low, xmax = dp18sw.high), linewidth = 0.2) +
  geom_point(data = sw, aes(x = dp18sw, y = Dp17sw), color = "tomato", shape = 3, size = 2) +
  geom_point(aes(fill = temp, shape = site), size = 4) +
  scale_fill_distiller(palette = "RdBu") +
  scale_shape_manual(values = c(21,22,23)) +
  theme_bw() + theme +
  annotate("text", x = -18, y = -110, label = "h", fontface = "bold", size = 5) +
  labs(x = expression(delta^"'18"*"O (\u2030, VSMOW)"),
       y = expression(Delta^"'17"*"O (per meg, VSMOW)"),
       fill = expression(paste("T"[Delta][47]*" (", degree, "C)")), shape = "")
p1

# joint proxy inversion of steady state model ----
bayes_data$site = factor(bayes_data$site, levels = c("Lantian", "Shilou", "Jiaxian"))
rc2 = rc |> select(age, temp)
colnames(rc2)[2] = "T47"
rc2 = rc2 |> arrange(age)
bayes_data = cbind(bayes_data, T47 = rc2$T47)

p2 = ggplot(bayes_data, aes(x = 1e3 * (exp(d18sw / 1e3) - 1), y = post_d18p)) +
  geom_errorbar(aes(ymin = post_d18p - post_d18p_sd, ymax = post_d18p + post_d18p_sd), 
                linewidth = 0.2) +
  geom_point(aes(fill = T47, shape = site), size = 4) +
  scale_fill_distiller(palette = "RdBu") +
  scale_shape_manual(values = c(21,22,23)) +
  theme_bw() + theme +
  annotate("text", x = -6, y = -19, label = "i", fontface = "bold", size = 5) +
  labs(x = expression(delta^"18"*"O"[sw]*" (\u2030, VSMOW)"),
       y = expression(delta^"'18"*"O"[p]*" (\u2030, VSMOW)"),
       fill = expression(paste("T"[Delta][47]*" (", degree, "C)")), shape = "")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv",
          common.legend = TRUE, legend = "right")
ggsave("figures/Fig3b.dp18sw_Dp17sw.jpg", width = 9, height = 4, dpi = 500)

