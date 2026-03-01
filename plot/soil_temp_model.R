rm(list = ls())
library(tidyverse, patchwork)
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())
bs = read_csv("output/BASS_bayes_vary_evap_ver2.csv")

# Hillel 1980
soil_temp = function(T_avg, T_amp, month, z) {
  t = (month - 3.5) / 12 
  w = 2 * pi
  d = 1.53
  temp = T_avg + T_amp * sin(w*t - z/d) / exp(z/d)
}

# time series
month = seq(1, 12, 1)
temp = soil_temp(12, 33, month, 0.5)
plot(month, temp)

# depth profile
z = seq(0.1, 1, 0.01)

# min T_avg and T_amp
for (i in seq_along(month)) {
  temp = soil_temp(12, 33, month[i], z)
  result = data.frame(depth = z, temp = temp)
  result$month = month[i]
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}

p1 = ggplot(sims, aes(x = temp, y = depth * 1e2, group = month, color = month)) +
  geom_path() +
  scale_color_gradientn(colors = c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
                                   "#F4A582", "#D6604D", "#B2182B", "#D6604D", "#F4A582", "#F7F7F7",
                                   "#92C5DE", "#4393C3", "#2166AC"), limits = c(1, 12)) +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(2, "mm")) +
  labs(x = expression(paste("T (", degree, "C)")),
       y = "Depth (cm)")

p2 = ggplot(bs, aes(x = temp, y = post_MAST)) +
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

p1 + p2 +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag.position = c(0.7, 0.2)
  )
ggsave("figures/soil_temp_model.png", width = 8.5, height = 3.6, dpi = 300)

for (i in seq_along(month)) {
  temp_s = sims |> 
    filter(depth == 0.4 & month == i)
  temp_d = sims |>
    filter(depth == 0.7 & month == i)
  temp_diff = temp_s$temp - temp_d$temp
  result = data.frame(month = i, DT = temp_diff)
  if (i == 1) {
    sims2 = result
  } else {
    sims2 = rbind(sims2, result)
  }
}
