rm(list = ls())
library(tidyverse)
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

# Hillel 1980
soil_temp = function(T_avg, T_amp, month, z) {
  t = (month - 3.5) / 12 
  w = 2 * pi
  d = 1.53
  temp = T_avg + T_amp * sin(w*t - z/d) / exp(z/d)
}

# time series
month = seq(1, 12, 1)
temp = soil_temp(12, 17, month, 0.5)
plot(month, temp)

# depth profile
z = seq(0.1, 1, 0.01)

# min T_avg and T_amp
for (i in seq_along(month)) {
  temp = soil_temp(12, 17, month[i], z)
  result = data.frame(depth = z, temp = temp)
  result$month = month[i]
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}

ggplot(sims, aes(x = temp, y = depth * 1e2, group = month, color = month)) +
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
ggsave("figures/soil_temp_model.png", width = 4, height = 3.6, dpi = 300)

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
