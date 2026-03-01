rm(list = ls())
pacman::p_load(tidyverse, patchwork)
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

post_data = read.csv("output/BASS_bayes_vary_evap_ver2.csv")

m1 = lm(post_d18p ~ temp, data = post_data)
summary(m1)
p1 = ggplot(post_data, aes(x = temp, y = post_d18p)) +
  geom_errorbar(aes(xmin = temp - temp.sd, xmax = temp + temp.sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = post_d18p - post_d18p_sd, ymax = post_d18p + post_d18p_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = site, alpha = d18p_eff),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = expression(delta^"18"*"O"[p]*" (\u2030, VSMOW)"),
       fill = "")

p2 = ggplot(post_data, aes(x = temp, y = post_RH * 1e2)) +
  geom_errorbar(aes(xmin = temp - temp.sd, xmax = temp + temp.sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = (post_RH - post_RH_sd) * 1e2, ymax = (post_RH + post_RH_sd) * 1e2),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = site, alpha = RH_eff),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = "RH (%)",
       fill = "")

p3 = ggplot(post_data, aes(x = temp, y = post_MAST)) +
  geom_errorbar(aes(xmin = temp - temp.sd, xmax = temp + temp.sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = post_MAST - post_MAST_sd, ymax = post_MAST + post_MAST_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = site, alpha = d18p_eff),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = expression(paste("MAST (",degree, "C)")),
       fill = "")

p1 + p2 + p3 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0.3, 0.9),
    legend.position = "bottom"
  )

ggsave("figures/FigS3_Bayes_BASS.jpg", width = 8, height = 3.5, dpi = 300)
