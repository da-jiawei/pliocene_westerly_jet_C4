rm(list = ls()) 
pacman::p_load(tidyverse, readxl, IsoplotR)

dat = read_csv("output/Dp17sw.csv")

dat$site = factor(dat$site, levels = c("Lantian", "Shilou", "Jiaxian"))
m1 = lm(temp ~ d18sw, data = dat)
summary(m1)

dat2 = na.omit(dat[, c("d18sw","d18sw.se","temp","temp.sd")])
m2 = york(dat2)
m2$p.value
ggplot(dat, aes(x = d18sw, y = temp)) +
  geom_errorbar(aes(xmin = d18sw - d18sw.se, xmax = d18sw + d18sw.se),
                width = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = temp - temp.sd, ymax = temp + temp.sd),
                width = 0, linewidth = 0.2) +
  # geom_smooth(method = "lm", fill = "grey50", color = "black", linewidth = 1, linetype = "dashed") +
  geom_abline(slope = m2$b, intercept = m2$a, linewidth = 1, linetype = "dashed") +
  geom_point(aes(fill = site), shape = 21, size = 3) +
  scale_fill_brewer(palette = "Paired") +
  # annotate("text", x = -12, y = 30, label = expression("R"^"2"*"=0.74")) +
  annotate("text", x = -11.5, y = 28, label = expression(italic(p)*"<0.001")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(2, "mm"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2)) +
  labs(x = expression(delta^"18"*"O"[sw]*" (\u2030, VSMOW)"),
       y = expression(paste("T"[Delta][47]*" (", degree, "C)")),
       fill = "")
ggsave("figures/TD47_vs_d18sw.png", width = 3.5, height = 3.7, dpi = 300)
