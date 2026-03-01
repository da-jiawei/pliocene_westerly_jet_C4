rm(list = ls())
pacman::p_load(R2jags, tidyverse, readxl, patchwork)
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

# load functions
source('JPI/helpers.R')
cat("\014")

## Read and groom data ----
clp = read_csv("output/Dp17sw.csv")[,-1]
clp = clp[order(clp$age),]

## Parse data into series
d18c = clp[c("d18c", "d18c.se")]
Dp17c = clp[c("Dp17c", "Dp17c.se")]
D47c = clp[c("D47", "D47.se")]
ages = clp$age

## MCMC ----
# priors
d = list(ages = ages, 
         d18c.obs = d18c, 
         Dp17c.obs = Dp17c,
         D47c.obs = D47c)

parms = c("d18p", "z", "RH", "evap", "MAST", "tsc", "MAP")

system.time({post.clp = jags.parallel(d, NULL, parms, "JPI/BASS_bayes.R",
                                      n.iter = 4e5, n.chains = 3, n.burnin = 2e5)})

View(post.clp$BUGSoutput$summary)

for (i in 1:length(parms)) {
  param = parms[i]
  plot.jpi(ages, post.clp$BUGSoutput$sims.list[[param]], n = 100, ylab = param)
}

post_data = clp |>
  select(age, site, temp, temp.sd, d18sw, d18sw.se)
for (i in 1:length(parms)) {
  post_data$Rhat = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "Rhat"]
  post_data$n.eff = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "n.eff"]
  post_data[, paste0("post_", parms[i])] = post.clp$BUGSoutput$mean[parms[i]]
  post_data[, paste0("post_", parms[i], "_sd")] = post.clp$BUGSoutput$sd[parms[i]]
  for (p in 1:nrow(post_data)) {
    if (post_data$Rhat[p] <= 1.05 & post_data$n.eff[p] >= 200){
      post_data[p, paste0(parms[i], "_eff")] = "positive"
    } else {
      post_data[p, paste0(parms[i], "_eff")] = "negative"
    }
  }
}
# write.csv(post_data, file = "output/BASS_bayes_vary_evap_ver2.csv")

## plot ----
# post_data = read.csv("output/BASS_bayes_vary_evap_ver2.csv")
m1 = lm(post_d18p ~ d18sw, data = post_data)
summary(m1)
ggplot(post_data, aes(x = d18sw, y = post_d18p)) +
  geom_errorbar(aes(xmin = d18sw - d18sw.se, xmax = d18sw + d18sw.se),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = post_d18p - post_d18p_sd, ymax = post_d18p + post_d18p_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = site, alpha = d18p_eff),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none")


ggplot(post_data, aes(x = temp, y = post_evap)) +
  geom_errorbar(aes(xmin = temp - temp.sd, xmax = temp + temp.sd),
                linewidth = .2, width = 0, color = "grey80") +
  # geom_errorbar(aes(ymin = post_evap - post_evap_sd, ymax = post_evap + post_evap_sd),
  #               linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = site),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = expression(delta^"'18"*"O"[sw]*" (\u2030, VSMOW)"),
       fill = "")

