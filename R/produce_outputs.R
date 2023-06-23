# load libraries
library(tidyverse)    # v2.0.0
library(tidybayes)    # v3.0.4
library(bayesplot)    # v1.10.0
library(brms)         # v2.19.0
library(broom.mixed)  # v0.2.9.4
library(performance)  # v0.10.4
library(patchwork)    # v1.1.2
library(tidygraph)    # v1.2.3
library(igraph)       # v1.5.0
library(ggraph)       # v2.1.0
library(ggpubr)       # v0.6.0
library(ggforce)      # v0.4.1
library(geomtextpath) # v0.1.1
library(vegan)        # v2.6-4
library(RColorBrewer) # v1.1-3
library(ggdist)       # v3.3.0
library(ggtext)       # v0.1.2
library(gtsummary)    # v1.7.1
library(extrafont)    # v0.19
loadfonts()

# load fcts
source("R/functions.R")

# choose model
modelname <- "main"
# ... and its modeled outcome
outc <- "PAR 30"
outc2 <- str_remove_all(str_to_lower(outc), " ")

# import model and underlying data
modelobj  <- readRDS(file = paste0("modelfiles/m_", modelname, ".rds"))
df        <- read_csv(paste0("modelfiles/data_model_unscaled_", modelname, ".csv"))
df_ld     <- read_csv(paste0("modelfiles/df_ld_", modelname, ".csv"))
ld        <- read_csv(paste0("modelfiles/ld_", modelname, ".csv"))
df_scaled <- modelobj$data
stopifnot(nrow(df) == nrow(df_scaled))

# set required number of cluster colors
nr.cols  <- length(unique(df$ling_cluster))
mycolors <- rev(colorRampPalette(brewer.pal(9, "Set1"), space = "Lab")(nr.cols))

# size of subsets of posterior draws for visualizations (saving time)
few  <- 100
some <- 500

# make posterior predictions
pp_few <- posterior_predict(modelobj, ndraws = few, seed = 12345)
pp     <- posterior_predict(modelobj, seed = 12345)

# outcome vector
y  <- as.numeric(modelobj$data[,outc2])

# relevant range of cooperative norms for plotting (untransformed scale)
range_coop <- c(6,10)

# get variable names
varnames <- read_csv2("utils/varnames.csv")
# ... and create labelling string
gt_labelling <- varnames %>% 
  filter(var %in% names(df_scaled) & !var %in% c("country.name", "ling_cluster")) %>% 
  summarise(expr = paste0("list(",
                          paste(
                            paste0(var, "~'", label, "'"), 
                            collapse = ","),
                          ")")) %>% pull(expr)

#############################
#### dot chart countries ####
#############################
p1 <- df %>%
  filter(!duplicated(country.name)) %>%
  mutate(country.name = fct_reorder(factor(country.name), cooperation)) %>%
  ggplot(aes(x = country.name,
             y = cooperation)) +
  geom_point(size = 2.5) +
  geom_segment(aes(yend = 1,
                   xend = country.name)) +
  geom_text(
    aes(label = format(round(cooperation, 2), nsmall = 2)),
    hjust = -1 / 2,
    family = fontfam,
    size = 3
  ) +
  theme_classic(14) +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    text = element_text(family = fontfam)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(1, 10, 1),
    limits = c(1, 11)
  ) +
  labs(y = "Cooperative norms",
       x = "")

p2 <- df %>%
  filter(!duplicated(country.name)) %>%
  mutate(
    country.name = fct_reorder(factor(country.name), cooperation),
    ling_cluster = str_remove(ling_cluster, "cluster_")
  ) %>%
  ggplot() +
  geom_label(
    aes(
      y = country.name,
      x = -1,
      label = country.name,
      fill = ling_cluster
    ),
    size = 3.5,
    hjust = "left",
    label.padding = unit(0.175, "lines"),
    label.r = unit(0, "lines"),
    label.size = 0,
    alpha = alpha_hgh,
    family = fontfam
  ) +
  theme_classic(14) +
  scale_x_continuous(limits = c(-1, 1)) +
  theme(
    axis.text = element_blank(),
    legend.position = "left",
    legend.key.size = unit(0.5, "cm"),
    text = element_text(family = fontfam),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 12, hjust = 0),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(hjust = 0),
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
    
  ) +
  scale_fill_manual(values = mycolors) +
  labs(fill = "Cultural\ncluster",
       x = "Country",) +
  guides(fill = guide_legend(override.aes = list(size = 0)))

p3 <- df %>%
  filter(!duplicated(country.name)) %>%
  mutate(country.name = fct_reorder(factor(country.name), cooperation)) %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = country.name, x = 0, label = coop_srvy_years),
    size = 3,
    hjust = 0,
    family = fontfam
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(family = fontfam),
    axis.title.x = element_text(hjust = 4 / 5),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0)) +
  labs(x = "Survey years")

p4 <- df %>%
  mutate(country.name = fct_reorder(factor(country.name), cooperation)) %>%
  group_by(country.name) %>%
  summarise(n = n()) %>%
  ggplot() +
  theme_classic(14) +
  geom_histogram(
    aes(x = country.name, y = n),
    stat = "identity",
    fill = greycol_lgt,
    col = greycol_lgt,
    width = 0.75
  ) +
  geom_text(
    aes(x = country.name, y = n, label = n),
    hjust = -1 / 5,
    size = 3,
    family = fontfam
  ) +
  labs(x = "", y = "Observations") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(family = fontfam),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 200, 50),
    limits = c(0, 240)
  ) +
  coord_flip()

png(paste0("outputs/fig_dotchart_", modelname, ".png"), width = 3800, height = 4100, res = 390)
p2 + p1 + p3 + p4 + plot_layout(design = "1#2#34
                                          ",
                                widths = c(1,-0.15,1.5,-0.5,1.5,0.7))
dev.off()

#######################
#### density plots ####
#######################
rawdf <- df %>% 
  pivot_longer(cols = c(par30, share_grp, cooperation), 
               names_to = "var", values_to = "val") %>% 
  mutate(
    var = fct_reorg(factor(var),
                    "Cooperative norms" = "cooperation",
                    "Group lending share" = "share_grp",
                    "PAR 30" = "par30"))

p5 <- rawdf %>%
  filter(var == "Cooperative norms") %>%
  ggplot() +
  geom_density(
    aes(x = val),
    col = "black",
    fill = greycol,
    alpha = alpha_med,
    linewidth = 1,
    bounds = c(1, 10)
  ) +
  theme_classic(14) +
  theme(
    text = element_text(family = fontfam),
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  ) +
  facet_wrap(~ var , scales = "free") +
  scale_fill_manual(values = twocolors) +
  scale_color_manual(values = twocolors) +
  labs(x = "",
       y = "Density",
       fill = "",
       col = "")

p6 <- rawdf %>%
  filter(var != "Cooperative norms") %>%
  ggplot() +
  geom_density(
    aes(x = val),
    col = "black",
    fill = greycol,
    alpha = alpha_med,
    linewidth = 1,
    bounds = c(0, 1)
  ) +
  theme_classic(14) +
  theme(
    text = element_text(family = fontfam),
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  ) +
  facet_wrap(~ var , scales = "free") +
  scale_fill_manual(values = twocolors) +
  scale_color_manual(values = twocolors) +
  labs(x = "",
       y = "",
       fill = "",
       col = "")

p7 <- df %>%
  ggplot() +
  geom_point(
    aes(x = cooperation, y = par30,  fill = share_grp),
    col = "white",
    size = 2.75,
    pch = 21
  ) +
  theme_classic(14) +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = fontfam),
    plot.caption = element_text(
      family = fontfam,
      size = 11,
      hjust = 0
    ),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = c(0.15, 0.75)
  ) +
  scale_x_continuous(breaks = seq(6, 10, 1), limits = c(5.5, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Cooperative norms", y = "PAR 30", fill = "Group lending share") +
  scale_fill_gradient(low = twocolors[2], high = twocolors[1]) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(caption = paste0("N = ", format(nrow(df), big.mark = ",")))

png(paste0("outputs/fig_distr_", modelname, ".png"), width = 3200, height = 2900, res = 480)
(p5 + p6 + plot_layout(widths = c(1,2))) / p7 + 
  plot_annotation(tag_levels = c("A")) + plot_layout(heights = c(1,1.5))
dev.off()

#################################
#### linguistic network plot ####
#################################
nodes <- df_ld %>% mutate(ling_cluster = str_remove(ling_cluster, "cluster_"))
ld_prox <- 1-ld
edges <- ld_prox %>% 
  mutate(e1 = names(.)) %>% 
  pivot_longer(-e1, names_to = "e2", values_to= "prox") %>% 
  mutate(prox = ifelse(prox < 0.01, NA, prox))

net_clust <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

for(i in seq(1, length(E(net_clust)), 1)){
  E(net_clust)$weight[i] <- ifelse(df_ld$ling_cluster[get.edgelist(net_clust)[i,1]] == 
                                     df_ld$ling_cluster[get.edgelist(net_clust)[i,2]], 
                                   180, 1)}
leg <- get_legend(p2)

set.seed(12345)
p8 <- ggraph(net_clust, layout = "fr") + 
  scale_edge_width(range = c(0.1,3)) +
  geom_edge_link(aes(width = prox), alpha = 0.08, color = "black", show.legend = FALSE) + 
  geom_node_point(aes(color = as.factor(ling_cluster)), size = 4, alpha = alpha_hgh) + 
  geom_node_text(aes(label = country.name), 
                 repel = TRUE, 
                 size = 4,
                 family = fontfam) +
  scale_color_manual(values = mycolors) +
  labs(color = "Cultural\ncluster") +
  theme_graph(background = "white") +
  theme(legend.position = "none",
        text = element_text(family = fontfam),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        plot.caption = element_text(size = 10, face = "plain", 
                                    family = fontfam, hjust = 0),
        plot.margin = unit(c(0,0,0.25,0), "cm")) + 
  guides(color = guide_legend(ncol = 1))

png(paste0("outputs/fig_network_", modelname, ".png"), width = 2100, height = 2000, res = 350)
as_ggplot(leg)+p8+plot_layout(widths = c(1,10))
dev.off()

#####################################
#### compute mantel correlations ####
#####################################
df_cultmap <- read_csv2("data/Cultural_Map.csv") %>% 
  mutate(country.name = countrycode::countrycode(Country, origin = "country.name", destination= "country.name")) %>% 
  select(country.name, Zone)

df_ld <- df_ld %>% 
  left_join(df %>% select(country.name, region) %>% distinct) %>% 
  left_join(df_cultmap)

clusterma <- matrix(NA, nrow = nrow(df_ld), ncol = nrow(df_ld))
rownames(clusterma) <- df_ld$country.name
colnames(clusterma) <- df_ld$country.name

regionma <- matrix(NA, nrow = nrow(df_ld), ncol = nrow(df_ld))
rownames(regionma) <- df_ld$country.name
colnames(regionma) <- df_ld$country.name

zonema <- matrix(NA, nrow = nrow(df_ld), ncol = nrow(df_ld))
rownames(zonema) <- df_ld$country.name
colnames(zonema) <- df_ld$country.name

rando <- matrix(NA, nrow = nrow(df_ld), ncol = nrow(df_ld))

for(i in 1:nrow(df_ld)){
  for(j in 1:nrow(df_ld)){
    clusterma[i,j] <- ifelse(df_ld$ling_cluster[i]==df_ld$ling_cluster[j], 1, 0)
    regionma[i,j]  <- ifelse(df_ld$region[i]==df_ld$region[j], 1, 0)
    zonema[i,j]    <- ifelse(df_ld$Zone[i]==df_ld$Zone[j], 1, 0)
    rando[i,j]     <- sample(c(0,1), size = 1, prob = c(0.9,0.1))
  }
}

mcor_ling <- round(mantel(1-ld, clusterma, permutations = 2000)$statistic, 2)
mcor_reg  <- round(mantel(regionma, clusterma, permutations = 2000)$statistic, 2)
mcor_zon  <- round(mantel(zonema, clusterma, permutations = 2000)$statistic, 2)

saveRDS(mcor_ling, file = "outputs/mcor_ling.RDS")
saveRDS(mcor_reg, file =  "outputs/mcor_reg.RDS")
saveRDS(mcor_zon, file =  "outputs/mcor_zon.RDS")

#################################################
# posterior summaries of parameters of interest #
#################################################
ests_all <- as_draws_df(modelobj, regex = TRUE, variable = c("share", "cooperation", "sd_")) %>% 
  select(-contains("phi"),
         -contains("agricul")) %>% 
  mutate(b_cooperation_for_grp = b_cooperation + `b_share_grp:cooperation`,
         b_cooperation_for_halfgrp = b_cooperation + 0.5 * `b_share_grp:cooperation`) %>% 
  select(contains("b_") | contains ("sd_")) %>% 
  mutate(across(.cols = everything(), exp, .names = "{.col}_or")) %>% 
  pivot_longer(cols = everything(), names_to = "var", values_to = "val") %>% 
  group_by(var) %>% 
  summarise(mean = mean(val),
            q05 = quantile(val,0.05),
            q95 = quantile(val,0.95)) %>% 
  mutate(across(.cols = where(is.numeric), ~format(round(.x, 2), nsmall = 2, trim = TRUE))) %>% 
  mutate(verbose = paste0(mean, " (90% CI: [", q05, ", ", q95, "])")) %>% 
  mutate(var = janitor::make_clean_names(var))

saveRDS(split(ests_all, ests_all$var), file = paste0("outputs/ests_all_", modelname, ".rds"))

#####################
# median odds ratio #
#####################
mor <- as_draws_df(modelobj, regex = TRUE, variable = c("sd_ling_cluster:country.name__Intercept")) %>% 
  mutate(mor = exp(`sd_ling_cluster:country.name__Intercept` * qnorm(3/4)),
         mor_inv = 1/mor) %>% 
  select(contains("mor")) %>% 
  summarise(across(.cols = everything(), ~format(round(mean(.x), 2), nsmall = 2, trim = TRUE)))

saveRDS(mor, file = paste0("outputs/mor_", modelname, ".rds"))

############
# bayes R2 #
############
r2 <- as_tibble(r2_bayes(modelobj, ci = 0.9, robust = FALSE))
saveRDS(r2, file = paste0("outputs/r2_", modelname, ".rds"))

################################################
# posterior probability of negative b1 + b3 ####
################################################
prob_neg <- as_draws_df(modelobj, regex = TRUE, variable = c("cooperation")) %>% 
  as.matrix %>% {.[,c(1,2)]} %>% rowSums %>% {mean(.<0)}
saveRDS(prob_neg, file = paste0("outputs/prob_neg_", modelname, ".rds"))

#################################################################
# variance components and country-level effect visualization ####
#################################################################
plot_at_group_share <- 1

m_dr <- modelobj %>%
  spread_draws(
    `b_Intercept`,
    `r_ling_cluster`[cluster, term],
    `r_ling_cluster:country.name`[country, term],
    `b_share_grp:cooperation`,
    `b_cooperation`,
    `b_share_grp`
  ) %>%
  mutate(country.name = str_remove(country, "cluster_[0-9][0-9]_")) %>%
  mutate(country.name = str_replace_all(country.name, "\\.", " ")) %>%
  left_join(df_scaled %>% select(country.name, cooperation) %>%
              distinct(),
            by = "country.name") %>%
  filter(cluster == str_sub(country, 1, 10))

global_est <- tibble(
  intercept = tidy(modelobj) %>% filter(term %in% c("(Intercept)", "share_grp")) %>%
    summarise(sum(
      estimate[1], estimate[2] * plot_at_group_share
    )) %>% pull,
  slope = tidy(modelobj) %>% filter(term %in% c("cooperation", "share_grp:cooperation")) %>%
    summarise(sum(
      estimate[1], estimate[2] * plot_at_group_share
    )) %>% pull
)

m_dr_sum <- m_dr %>%
  mutate(
    alpha_k =
      `r_ling_cluster:country.name` +
      b_Intercept +
      b_share_grp * plot_at_group_share +
      (
        b_cooperation + `b_share_grp:cooperation` * plot_at_group_share
      ) * cooperation,
    alpha_kc =
      `r_ling_cluster` +
      `r_ling_cluster:country.name` +
      b_Intercept +
      b_share_grp * plot_at_group_share +
      (
        b_cooperation + `b_share_grp:cooperation` * plot_at_group_share
      ) * cooperation
  ) %>%
  group_by(country.name) %>%
  summarise(
    est_k = mean(alpha_k),
    lwr_k = est_k - sd(alpha_k),
    upr_k = est_k + sd(alpha_k),
    est_kc = mean(alpha_kc),
    lwr_kc = est_kc - sd(alpha_kc),
    upr_kc = est_kc + sd(alpha_kc),
    ling_cluster = unique(cluster),
    cooperation = unique(cooperation)
  )

p8 <- ggplot() + 
  geom_mark_ellipse(data = m_dr_sum,
                    aes(x = cooperation, y = est_kc,
                        fill = str_sub(ling_cluster, 9, 10)),
                    alpha = alpha_med,
                    show.legend = FALSE,
                    col = NA,
                    expand = unit(1.25, "mm")) +
  geom_abline(data = global_est,
              aes(intercept = intercept, slope = slope)) +
  geom_errorbar(data = m_dr_sum, aes(x = cooperation, ymin = lwr_kc, ymax = upr_kc), 
                alpha = alpha_low) + 
  geom_point(data = m_dr_sum, aes(x = cooperation, y = est_kc, 
                                  fill = str_sub(ling_cluster, 9, 10)),
             key_glyph = "rect",
             pch = 21,
             col = "black",
             stroke = 0.25,
             size = 1.75) + 
  theme_classic(14) + 
  labs(x = "Cooperative norms (std.)",
       col = "Cultural\ncluster",
       fill = "Cultural\ncluster",
       y = paste0("Combined country-level effect", 
                  "<br><span style = 'font-size:10pt'> 
                                             on linear predictor for a group lender</span>")) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        panel.grid = element_blank(),
        axis.title.y = element_markdown(),
        text = element_text(family = fontfam)) + 
  scale_color_manual(values = mycolors) + 
  scale_fill_manual(values = mycolors)

v <- as_draws_df(modelobj, regex = TRUE, variable = c("sd_"))  %>% pivot_longer(
  cols = contains("sd_"),
  names_to = "var",
  values_to = "val") %>%  mutate(var = case_when(
    
    grepl("ling_cluster_", var)              ~ "Cultural\ncluster",
    grepl("ling_cluster:country.name_", var) ~ "Country"
  )) %>% 
  mutate(var = factor(var, ordered = TRUE,
                      levels = c("Country", "Cultural\ncluster")))

p9 <- ggplot(v)  +
  stat_halfeye(
    aes(x = val, y = var),
    fill = greycol_lgt,
    density = density_bounded(bounds = c(0, NA)),
    .width = iv_lims,
    point_interval = iv_type
  ) +
  
  geom_vline(xintercept = 0,
             lty = 2,
             col = "black") +
  theme_classic(14) +
  scale_thickness_shared() +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  labs(y = "Posterior density", x = "Standard deviation") +
  theme(
    text = element_text(family = fontfam),
    panel.grid.major.y = element_line(),
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  ) +
  scale_fill_manual(values = c(greycol, greycol_lgt))

png(paste0("outputs/fig_varcomp_", modelname, ".png"), width = 3700, height = 1800, res = 460)
p9 + p8 + plot_layout(widths = c(1, 1.75)) +
  plot_annotation(tag_levels = c("A"))
dev.off()

###########################################################
#### plot posterior predictions for selected countries ####
###########################################################
ctries <- c("Philippines", "Mongolia", "Pakistan")
times  <- list(c(2008:2014), c(2008:2014), c(2008:2014))

newdata <- map2_dfr(ctries, times, 
                    function(ctry, time){
                      df_scaled %>% 
                        mutate(country.name = as.character(country.name)) %>% 
                        filter(country.name == ctry,
                               year %in% time
                        ) %>% 
                        mutate(country.name = fct_reorder(factor(country.name), cooperation),
                               cooperation_obs = cooperation) %>% 
                        select(-cooperation, -share_grp) %>% 
                        expand_grid(cooperation = seq(range_coop[1], range_coop[2], 0.1),
                                    share_grp = c(0,1)) %>% 
                        mutate(cooperation = (cooperation-mean(df$cooperation))/sd(df$cooperation))
                    })

preds <- epred_draws(object = modelobj, newdata = newdata, seed = 12345, ndraws = few)

plotdata <- preds %>%
  mutate(
    cooperation_bt = cooperation * sd(df$cooperation) + mean(df$cooperation),
    cooperation_obs_bt = cooperation_obs * sd(df$cooperation) + mean(df$cooperation)
  ) %>%
  group_by(country.name, cooperation_bt, .draw , share_grp) %>%
  summarise(
    cooperation_obs_bt = unique(cooperation_obs_bt),
    marginalized_pred = mean(.epred),
    .groups = "drop"
  ) %>%
  mutate(country.name = fct_reorder(country.name, cooperation_obs_bt))

p10 <- plotdata %>%
  ggplot() +
  geom_vline(aes(xintercept = cooperation_obs_bt, linetype = "Observed value"), col = greycol) +
  geom_line(aes(x = cooperation_bt, 
                y = marginalized_pred, 
                group = paste0(.draw, share_grp),
                col = share_grp), alpha = alpha_low) + 
  labs(x = "Cooperative norms", 
       y = paste0("Expected ", outc),
       col = "Group lending share",
       linetype = "") + 
  facet_wrap(~country.name, nrow = 1) +
  theme_classic(14) + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 11),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        text = element_text(family = fontfam)) + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.25,0.05)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(6,10,1)) +
  scale_color_gradient(low = twocolors[2], high = twocolors[1], breaks = c(0,1)) + 
  scale_linetype_manual(values = "dashed") +
  coord_cartesian(ylim = c(0,0.25)) +
  guides(col = guide_legend(override.aes = list(lwd = 1, alpha = 1),
                            title.theme = element_text(family = fontfam, size = 11)))

png(paste0("outputs/fig_predcoop_", modelname, ".png"), width = 3100, height = 1400, res = 400)
p10
dev.off()

###############################################
#### multi-panel figure with model results ####
###############################################
# densities
d <- as_draws_df(modelobj, variable = c("b_cooperation", "b_share_grp", "b_share_grp:cooperation")) 
d <- d %>% pivot_longer(
  cols = contains("b_"),
  names_to = "var",
  values_to = "val") %>%
  mutate(var = case_when(
    grepl("share_grp:cooperation", var) ~ "Cooperative norms x\nGroup lending share",
    grepl("share_grp", var)             ~ "Group lending share",
    grepl("cooperation", var)           ~ "Cooperative norms"
  )) %>% 
  mutate(var = factor(var, ordered = TRUE,
                      levels = c("Cooperative norms x\nGroup lending share",
                                 "Group lending share",
                                 "Cooperative norms" )))

dd <- d
levels(dd$var) <- c("COOP \u00B7\nGROUP", "GROUP", "COOP")

p11 <- ggplot(dd) + 
  stat_halfeye(aes(x = val, y = var),
               fill = greycol_lgt,
               .width = iv_lims,
               point_interval = iv_type) +
  
  geom_vline(xintercept = 0, lty = 2, col = "black") +
  theme_classic(14) + 
  scale_thickness_shared() +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  labs(y = "Posterior density", x = "Coefficient") + 
  theme(text = element_text(family = fontfam),
        panel.grid.major.y = element_line(),
  ) + 
  scale_fill_manual(values = c(greycol, greycol_lgt))

ddd <- d %>% 
  pivot_wider(values_from = val, names_from = var) %>% 
  expand_grid(x = seq(0,1,0.01)) %>% 
  mutate(link_pred = `Cooperative norms` * 1 +
           `Cooperative norms x\nGroup lending share` * 1 * x) 

# linear predictor effect modification
p12 <- ggplot(ddd) +
  stat_lineribbon(
    aes(x = x, y = link_pred),
    .width = iv_lims,
    point_interval = iv_type,
    show.legend = FALSE
  ) +
  theme_classic(14) +
  theme(text = element_text(family = fontfam),
        axis.title.y = element_markdown()) +
  geom_hline(yintercept = 0,
             lty = 2,
             col = "black") +
  labs(
    x = "Group lending share",
    y = paste0(
      "Marginal effect",
      "<br><span style = 'font-size:10pt'>
                                             of one SD increase in cooperative norms<br>on linear predictor </span>"
    )
  ) +
  coord_cartesian(ylim = c(-0.35, 0.1)) +
  scale_fill_manual(values = c(greycol_lgt, greycol))

# average marginal coop effect by decile
newdata <- modelobj$data %>% mutate(share_grp_group = cut(share_grp,
                                                          seq(0,1,0.1),
                                                          include.lowest = TRUE,
                                                          dig.lab = 1),
                                    cooperation = cooperation + 1,
                                    cooperation_at = "plus_1_sd",
                                    rowid = 1:nrow(.)) %>% 
  rbind(modelobj$data %>% mutate(share_grp_group = cut(share_grp, 
                                                       seq(0,1,0.1),
                                                       include.lowest = TRUE,
                                                       dig.lab = 1),
                                 cooperation_at = "observed",
                                 rowid = 1:nrow(.))
  ) %>% 
  group_by(share_grp_group) %>% 
  mutate(n_obs = n()/2) %>% 
  ungroup

dddd <- epred_draws(object = modelobj,
                    newdata = newdata,
                    seed = 12345,
                    ndraws = some) %>% 
  as_tibble %>% 
  select(rowid, .draw, .epred, cooperation_at, share_grp, share_grp_group, n_obs) %>%
  pivot_wider(names_from = "cooperation_at", values_from = ".epred") %>% 
  mutate(eff = 
           (
             (plus_1_sd / observed) - 1
           ) * 100,
  ) %>% 
  group_by(share_grp_group) %>% 
  mutate(n_obs = length(unique(rowid)))

p13 <- ggplot(dddd) +
  stat_pointinterval(aes(y = eff, x = factor(share_grp_group)),
                     point_interval = iv_type,
                     .width = iv_lims) +
  
  geom_hline(yintercept = 0,
             lty = 2,
             col = "black") +
  geom_text(
    aes(
      x = factor(share_grp_group),
      y = 17.5,
      label = paste0("N =\n", n_obs)
    ),
    check_overlap = TRUE,
    size = 3,
    family = fontfam
  ) +
  theme_classic(14) +
  labs(
    y = paste0(
      "Average marginal effect<br>
                   <span style = 'font-size:10pt'>
                    of one SD increase in cooperative norms<br>as relative difference in expected ",
      outc,
      "</span>"
    ),
    x = "Group lending share",
    col = "Credible interval"
  ) +
  theme(
    text = element_text(family = fontfam),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(-30, 30, 10),
    labels = function(x)
      paste0(x, "%")
  ) +
  coord_cartesian(ylim = c(-30, 20)) +
  scale_color_manual(values = c(greycol_lgt, greycol)) +
  guides(col = guide_legend(nrow = 1, reverse = TRUE))

# contour plot
newdata_contour <- df_scaled %>% 
  
  select(-share_grp, -cooperation) %>% 
  expand_grid(cooperation = seq(range_coop[1], range_coop[2], 0.2)) %>%
  mutate(cooperation = (cooperation-mean(df$cooperation))/sd(df$cooperation)) %>%
  expand_grid(share_grp = seq(0,1,0.1))

epreds_contour <- epred_draws(object = modelobj, newdata = newdata_contour, 
                              ndraws = some,
                              seed = 12345) %>% 
  mutate(cooperation_bt = cooperation*sd(df$cooperation) + mean(df$cooperation)) %>% 
  group_by(cooperation_bt, share_grp) %>% 
  summarise(mean_epred = mean(.epred), .groups = "drop")

avg_coop <- df %>%
  filter(!duplicated(country)) %>%
  summarise(avg_coop = mean(cooperation),
            n=n()) %>%
  pull(avg_coop)

p14 <- ggplot(epreds_contour) +
  geom_contour_filled(
    aes(x = share_grp,
        y = cooperation_bt,
        z = mean_epred),
    alpha = alpha_hgh,
    binwidth = 0.005
  ) +
  geom_texthline(
    data = tibble(avg_coop = avg_coop),
    aes(yintercept = avg_coop, label = "Average country"),
    hjust = 0.15,
    col = "grey20",
    size = 2.75,
    lty = 3,
    lwd = 1 / 3,
    family = "Segoe UI"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = 6:10,
                     limits = c(range_coop[1], range_coop[2])) +
  theme_classic(14) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    text = element_text(family = fontfam),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.margin = margin(l = 0, r = 0, unit = "pt")
  ) +
  labs(y = "Cooperative norms",
       x = "Group lending share",
       fill = paste0("Expected\n", outc)) +
  scale_fill_manual(values = rev(palcolors))

png(paste0("outputs/fig_params_", modelname, ".png"), width = 4100, height = 3100, res = 450)
(p11 + p12) / (p13 + p14) + 
  plot_annotation(tag_levels = c("A"))
dev.off()

##############################################
#### marginal posterior predictions check ####
##############################################
pred_comp <- as_tibble(t(pp_few)) %>% 
  pivot_longer(cols = everything(),
               names_to = "var",
               values_to = "val") %>% 
  mutate(origin = "Model") %>% 
  bind_rows(tibble(
    var = "V1",
    val = y,
    origin = "Observed"))

p15 <- ggplot(pred_comp) +
  geom_density(
    aes(
      x = val,
      group = paste0(origin, var),
      col = origin,
      linewidth = origin
    ),
    key_glyph = draw_key_path,
    bounds = c(0, 1)
  ) +
  theme_classic(14) +
  labs(y = "Density",
       x = outc,
       col = "") +
  theme(legend.position = c(0.75, 0.95),
        text = element_text(family = fontfam)) +
  scale_color_manual(values = c(greycol_lgt, "black")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linewidth_discrete(range = c(0.25, 1)) +
  guides(linewidth = "none",
         col = guide_legend(reverse = TRUE))

png(paste0("outputs/fig_overlay_", modelname, ".png"), width = 1600, height = 1500, res = 400)
p15
dev.off()

####################################################
#### country median posterior predictions check ####
####################################################
country.name.n <- df %>%
  group_by(country.name) %>%
  mutate(country.name.n = paste0(country.name, " (N = ", n(), ")"),
         .groups = "drop") %>% pull(country.name.n)

color_scheme_set("darkgray")

p16 <- ppc_stat_grouped(
  y = y,
  yrep = pp,
  group = country.name.n,
  facet_args = list(scales = "free_y", ncol = 5),
  binwidth = 0.0025,
  stat = "median"
) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  theme_classic(14) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    text = element_text(family = fontfam, size = 14),
    strip.text = element_text(size = 6.5, face = "bold"),
    strip.background = element_blank()
  ) +
  labs(x = outc) +
  coord_cartesian(xlim = c(0, 0.25))

png(paste0("outputs/fig_postcheckctry_", modelname, ".png"), width = 3100, height = 4800, res = 400)
p16
dev.off()

###########################
#### descriptive table ####
###########################
tab1 <- df %>% 
  mutate(wdi_gdp = wdi_gdp / 1000) %>%
  select(any_of(names(df_scaled))) %>% 
  select(-country.name, -ling_cluster) %>% 
  select(-contains("_sq")) %>% 
  relocate(par30, share_grp) %>% 
  tbl_summary(digits    = list(all_continuous() ~ 2,
                               all_categorical() ~ c(0, 1)),
              statistic = list(all_continuous() ~ "{mean},{sd},{median},{p25},{p75}"),
              label     = eval(parse(text = gt_labelling))) %>% 
  modify_header(
    label = "Variable",
    stat_0 = "val"
  ) %>% 
  as_tibble() %>% 
  separate(col = val, sep = ",", into = c("Mean (or N)", "SD", "Median", "P25", "P75")) %>% 
  mutate(across(.cols = everything(), function(x)ifelse(is.na(x), "", x)))

saveRDS(tab1, file = paste0("outputs/tab_descriptives_", modelname, ".rds"))

##################################################
#### main interest parameter posteriors table ####
##################################################
shape_par <- as_tibble(summary(modelobj, 
                               robust = FALSE,
                               prob = 0.9)$spec_pars) %>% 
  select(Mean = Estimate,
         SD = Est.Error,
         `0.05 quantile` = `l-90% CI`,
         `0.95 quantile` = `u-90% CI`) %>% 
  mutate(Parameter = c(
    "phi",
    "cut1",
    "cut2"
  ))

est_pars <- tidy(modelobj,
                 robust = FALSE,
                 conf.level = 0.9,
                 conf.method = "quantile") %>%
  filter(term %in% c("(Intercept)", 
                     "cooperation", 
                     "share_grp", 
                     "share_grp:cooperation",
                     "sd__(Intercept)")) %>% 
  mutate(term = 
           case_when(group == "ling_cluster" ~ "SD Cultural cluster",
                     group == "ling_cluster:country.name" ~ "SD Country",
                     TRUE ~ term)) %>% 
  select(Parameter = term,
         Mean = estimate,
         SD = std.error,
         `0.05 quantile` = conf.low,
         `0.95 quantile` = conf.high) %>% 
  mutate(Parameter = factor(Parameter, ordered = TRUE,
                            levels = c(
                              "(Intercept)", 
                              "cooperation", 
                              "share_grp", 
                              "share_grp:cooperation",
                              "SD Cultural cluster",
                              "SD Country"
                            ))) %>% 
  arrange(Parameter) %>% 
  mutate(Parameter = as.character(Parameter))

tab2 <- rbind(est_pars , shape_par) %>% 
  mutate(Description =
           case_when(Parameter == "(Intercept)"           ~ "Global regression intercept",
                     Parameter == "cooperation"           ~ "Regression coefficient for $\\text{COOP}$",
                     Parameter == "share_grp"             ~ "Regression coefficient for $\\text{GROUP}$",
                     Parameter == "share_grp:cooperation" ~ "Regression coefficient for $\\text{COOP} \\cdot \\text{GROUP}$",
                     Parameter == "SD Cultural cluster"   ~ Parameter,
                     Parameter == "SD Country"            ~ Parameter,
                     Parameter == "phi"                   ~ "Precision parameter of Beta distribution",
                     Parameter == "cut1"                  ~ "First cutpoint parameter",
                     Parameter == "cut2"                  ~ "Second cutpoint parameter"
           ), .after = "Parameter",
         Parameter =
           case_when(Parameter == "(Intercept)"           ~ "$\\beta_0$",
                     Parameter == "cooperation"           ~ "$\\beta_1$",
                     Parameter == "share_grp"             ~ "$\\beta_2$",
                     Parameter == "share_grp:cooperation" ~ "$\\beta_3$",
                     Parameter == "SD Cultural cluster"   ~ "$\\sigma_{\\psi}$",
                     Parameter == "SD Country"            ~ "$\\sigma_{\\zeta}$",
                     Parameter == "phi"                   ~ "$\\phi$",
                     Parameter == "cut1"                  ~ "$c_1$",
                     Parameter == "cut2"                  ~ "$c_2$"
           )
  ) 

saveRDS(tab2, file = paste0("outputs/tab_model_", modelname, ".rds"))

#############################################
#### varying intercepts posteriors table ####
#############################################
tab3 <- tidy(modelobj, effect = "ran_vals",
             robust = FALSE,
             conf.level = 0.9,
             conf.method = "quantile") %>% 
  mutate(level = case_when(
    group == "ling_cluster" ~ str_remove(level, "cluster_"),
    group == "ling_cluster:country.name" ~ str_sub(str_replace_all(level, "\\.", " "), 12, 100)
  )) %>% 
  select(Parameter = level,
         Mean = estimate,
         SD = std.error,
         `0.05 quantile` = conf.low,
         `0.95 quantile` = conf.high) %>% 
  arrange(Parameter)

saveRDS(tab3, file = paste0("outputs/tab_varyicepts_", modelname, ".rds"))

#####################################################
#### all regression coefficient posteriors table ####
#####################################################
tab4 <- tidy(modelobj, effect = "fixed",
             robust = FALSE, # mean/SD instead of median/MAD of posterior
             conf.level = 0.9,
             conf.method = "quantile") %>% 
  filter(!term %in% c(
    "sd__(Intercept)"
  )) %>% 
  left_join(varnames, by = c("term" = "var")) %>% 
  mutate(label = case_when(
    term == "wdi_gdp_sq" ~ paste0(label[term=="wdi_gdp"], " squared"),
    term == "(Intercept)" ~ "Intercept",
    grepl("year", term) ~ as.character(extract_numeric(term)),
    grepl("inst_type", term) ~ str_sub(term, 10,50),
    TRUE ~ label)) %>%
  mutate(
    label = ifelse(grepl("CreditUnion", label), "Credit Union/Cooperative", label),
    label = ifelse(grepl("RuralBank", label), "Rural Bank", label),
  ) %>% 
  mutate(label = ifelse(term == "share_grp:wdi_agricul",
                        paste0(
                          varnames$label[varnames$var=="share_grp"],
                          " $\\cdot$ ",
                          varnames$label[varnames$var=="wdi_agricul"]
                        ),
                        label)) %>% 
  mutate(label = ifelse(term == "share_grp:cooperation",
                        paste0(
                          varnames$label[varnames$var=="share_grp"],
                          " $\\cdot$ ",
                          varnames$label[varnames$var=="cooperation"]
                        ),label)) %>% 
  select(Parameter = label,
         Mean = estimate,
         SD = std.error,
         `0.05 quantile` = conf.low,
         `0.95 quantile` = conf.high) %>% 
  mutate(Parameter = str_replace_all(Parameter, "\\%", "\\\\%"))

saveRDS(tab4, file = paste0("outputs/tab_adjcoeff_", modelname, ".rds"))

############################################
#### predictor correlation matrix table ####
############################################
corm <- df_scaled %>% 
  select(-contains("_sq"), 
         -contains("country.name"),
         -ling_cluster, 
         -par30,
         -year,
         -inst_type) %>% 
  cor %>% 
  {format(round(., 2), nsmall=2)}

corm[lower.tri(corm)] <- "-"

tab5 <- as_tibble(corm, rownames = "var") %>% 
  left_join(varnames) %>% 
  relocate(Variable = label) %>% 
  mutate(ID = LETTERS[1:nrow(.)],
         Variable = paste0(ID, ": ", Variable)) %>% 
  select(-var, -ID)

names(tab5) <- c("", LETTERS[1:(ncol(tab5)-1)])

saveRDS(tab5, file = paste0("outputs/tab_corrm_", modelname, ".rds"))

######################################
#### sensitivity analyses outputs ####
######################################
# collect model objects
file_names <- list.files("modelfiles", ".rds", full.names = TRUE)
file_names <- file_names[!grepl("prior", file_names)]
modlist <- list(NA)
for(i in 1:length(file_names)){
  modlist[[i]] <- readRDS(file_names[i])
}
modlist_names <- str_remove(str_remove(file_names, "modelfiles/m_"), ".rds")
names(modlist) <- modlist_names

drws <- map_dfr(modlist, .id = "model", function(modelobj) {
  spread_draws(modelobj,
               `b_share_grp:cooperation`,
               `b_cooperation`) %>%
    mutate(b_cooperation_for_grp = `b_cooperation` + 1 * `b_share_grp:cooperation`)
}) %>%
  mutate(
    model = fct_reorg(
      model,
      "S9"   = "noitcpts",
      "S8"   = "phi",
      "S7"   = "ranslopes",
      "S6"   = "mfi",
      "S5"   = "par90",
      "S4"   = "allintact",
      "S3"   = "otheradj",
      "S2"   = "nohaiti",
      "S1"   = "weaker",
      "Main" = "main"
    )
  ) %>%
  pivot_longer(
    cols = contains("b_cooperation"),
    names_to = "share_grp",
    values_to = "b"
  ) %>%
  mutate(share_grp = ifelse(share_grp == "b_cooperation_for_grp", 1, 0))

p19 <- ggplot(drws) +
  stat_pointinterval(
    aes(
      y = b,
      x = factor(model),
      col = factor(share_grp)
    ),
    position = "dodge",
    point_interval = iv_type,
    .width = iv_lims
  ) +
  geom_hline(yintercept = 0,
             lty = 2,
             col = "black") +
  theme_classic(14) +
  labs(
    y = paste0(
      "Marginal effect",
      "<br><span style = 'font-size:10pt'>
                  of one SD increase in cooperative norms on linear predictor</span>"
    ),
    x = "Model specification",
    col = "Group lending share"
  ) +
  theme(
    text = element_text(family = fontfam),
    axis.title.x = element_markdown(),
    legend.position = "bottom"
  ) +
  coord_flip(ylim = c(NA, 0.1)) +
  scale_color_manual(values = rev(twocolors)) +
  guides(colour = guide_legend(
    override.aes = list(linetype = 0, size = 3),
    size = "none",
    linetype = "none"
  ))

nobs <- map_dfr(names(modlist), function(modelname) {
  dat <-
    read_csv(paste0(
      "modelfiles/data_model_unscaled_",
      modelname,
      ".csv"
    ))
  tibble(
    model = modelname,
    n_obs = format(nrow(dat), big.mark = ","),
    n_countries = length(unique(dat$country.name))
  )
}) %>%
  mutate(
    model = fct_reorg(
      model,
      "S9"   = "noitcpts",
      "S8"   = "phi",
      "S7"   = "ranslopes",
      "S6"   = "mfi",
      "S5"   = "par90",
      "S4"   = "allintact",
      "S3"   = "otheradj",
      "S2"   = "nohaiti",
      "S1"   = "weaker",
      "Main" = "main"
    )
  ) %>% arrange(model)

p20 <- nobs  %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = model, x = 0, label = n_obs),
    size = 3,
    hjust = 1 / 2,
    family = fontfam
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 1 / 2),
    text = element_text(family = fontfam),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0), position = "top") +
  labs(x = "Observations")

p21 <- nobs  %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = model, x = 0, label = n_countries),
    size = 3,
    hjust = 1 / 2,
    family = fontfam
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 1 / 2),
    text = element_text(family = fontfam),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0), position = "top") +
  labs(x = "Countries")


prbs_r2 <- map_dfr(modlist, .id = "model", function(modelobj) {
  pr_neg_fullgrp = as_draws_df(modelobj,
                               regex = TRUE,
                               variable = c("cooperation")) %>%
    as.matrix %>% {
      .[, c(1, 2)]
    } %>% rowSums %>% {
      mean(. < 0)
    }
  pr_neg_idv = as_draws_df(modelobj, variable = c("b_cooperation")) %>%
    as.matrix %>% {
      .[, c(1)]
    } %>% {
      mean(. < 0)
    }
  
  
  r2 = as_tibble(performance:::r2_bayes(modelobj, ci = 0.9, robust = FALSE))
  if (nrow(r2) > 1)
    r2 <- r2 %>% filter(Component == "conditional")
  
  tibble(
    pr_neg_fullgrp = pr_neg_fullgrp,
    pr_neg_idv     = pr_neg_idv,
    r2 = r2 %>% pull(R2)
  )
}) %>%
  mutate(
    model = fct_reorg(
      model,
      "S9"   = "noitcpts",
      "S8"   = "phi",
      "S7"   = "ranslopes",
      "S6"   = "mfi",
      "S5"   = "par90",
      "S4"   = "allintact",
      "S3"   = "otheradj",
      "S2"   = "nohaiti",
      "S1"   = "weaker",
      "Main" = "main"
    ),
    across(contains("pr_neg"), ~ format(round(.x, 3), nsmall = 3)),
    r2 = format(round(r2, 3), nsmall = 3)
  ) %>% arrange(model)

p22 <- prbs_r2 %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = model, x = 0, label = pr_neg_fullgrp),
    col = twocolors[1],
    size = 3,
    hjust = 1 / 2,
    family = fontfam,
    fontface = "bold"
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 1 / 2),
    text = element_text(family = fontfam),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0), position = "top") +
  labs(x = "Pr negative")

p23 <- prbs_r2 %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = model, x = 0, label = pr_neg_idv),
    col = twocolors[2],
    size = 3,
    hjust = 1 / 2,
    family = fontfam,
    fontface = "bold"
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 1 / 2),
    text = element_text(family = fontfam),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0), position = "top") +
  labs(x = "")

p24 <- prbs_r2 %>%
  ggplot() +
  theme_classic(14) +
  geom_text(
    aes(y = model, x = 0, label = r2),
    size = 3,
    hjust = 1 / 2,
    family = fontfam
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 1 / 2),
    text = element_text(family = fontfam),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 0), position = "top") +
  labs(x = expression(Bayesian ~ R ^ 2))

png("outputs/fig_sensitivity.png", width = 3900, height = 2000, res = 450)
p19 + p22 + p23 + p24 + p20 + p21 + plot_layout(widths = c(1,0.25,0.25, 0.25, 0.25,0.25))
dev.off()


######################################
#### visualize prior simulations #####
######################################
modelobj  <- readRDS(file = paste0("modelfiles/prior_m_", modelname, ".rds"))

d <- as_draws_df(modelobj, variable = c("b_cooperation", "b_share_grp", "b_share_grp:cooperation")) 
d <- d %>% pivot_longer(
  cols = contains("b_"),
  names_to = "var",
  values_to = "val"
) %>% 
  mutate(var = case_when(
    grepl("share_grp:cooperation", var) ~ "Cooperative norms x\nGroup lending share",
    grepl("share_grp", var)             ~ "Group lending share",
    grepl("cooperation", var)           ~ "Cooperative norms"
  )) %>% 
  mutate(var = factor(var, ordered = TRUE,
                      levels = c("Cooperative norms x\nGroup lending share",
                                 "Group lending share",
                                 "Cooperative norms" ))
  )
dd <- d
levels(dd$var) <- c("COOP \u00B7\nGROUP", "GROUP", "COOP")

ddd <- d %>% 
  pivot_wider(values_from = val, names_from = var) %>% 
  expand_grid(x = seq(0,1,0.01)) %>% 
  mutate(link_pred = `Cooperative norms` * 1 +
           `Cooperative norms x\nGroup lending share` * 1 * x) 

p17 <- ggplot(ddd) +
  stat_lineribbon(
    aes(x = x, y = link_pred),
    .width = iv_lims,
    point_interval = iv_type,
    show.legend = FALSE
  ) +
  theme_classic(14) +
  theme(text = element_text(family = fontfam),
        axis.title.y = element_markdown()) +
  geom_hline(yintercept = 0,
             lty = 2,
             col = "black") +
  labs(
    x = "Group lending share",
    y = paste0(
      "Marginal effect",
      "<br><span style = 'font-size:10pt'>
      of one SD increase in cooperative norms<br>on linear predictor </span>"
    )
  ) +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_fill_manual(values = c(greycol_lgt, greycol))

m_dr <- modelobj %>%
  spread_draws(
    `b_Intercept`,
    `r_ling_cluster`[cluster, term],
    `r_ling_cluster:country.name`[country, term],
    `b_share_grp:cooperation`,
    `b_cooperation`,
    `b_share_grp`
  ) %>%
  mutate(country.name = str_remove(country, "cluster_[0-9][0-9]_")) %>%
  mutate(country.name = str_replace_all(country.name, "\\.", " ")) %>%
  left_join(df_scaled %>% select(country.name, cooperation) %>%
              distinct(),
            by = "country.name") %>%
  filter(cluster == str_sub(country, 1, 10))

global_est <- tibble(
  intercept = tidy(modelobj) %>% filter(term %in% c("(Intercept)", "share_grp")) %>%
    summarise(sum(
      estimate[1], estimate[2] * plot_at_group_share
    )) %>% pull,
  slope = tidy(modelobj) %>% filter(term %in% c("cooperation", "share_grp:cooperation")) %>%
    summarise(sum(
      estimate[1], estimate[2] * plot_at_group_share
    )) %>% pull
)

m_dr_sum <- m_dr %>%
  mutate(
    alpha_k =
      `r_ling_cluster:country.name` +
      b_Intercept +
      b_share_grp * plot_at_group_share +
      (
        b_cooperation + `b_share_grp:cooperation` * plot_at_group_share
      ) * cooperation,
    alpha_kc =
      `r_ling_cluster` +
      `r_ling_cluster:country.name` +
      b_Intercept +
      b_share_grp * plot_at_group_share +
      (
        b_cooperation + `b_share_grp:cooperation` * plot_at_group_share
      ) * cooperation
  ) %>%
  group_by(country.name) %>%
  summarise(
    est_k = mean(alpha_k),
    lwr_k = est_k - sd(alpha_k),
    upr_k = est_k + sd(alpha_k),
    est_kc = mean(alpha_kc),
    lwr_kc = est_kc - sd(alpha_kc),
    upr_kc = est_kc + sd(alpha_kc),
    ling_cluster = unique(cluster),
    cooperation = unique(cooperation)
  )

p18 <- ggplot() +
  geom_mark_ellipse(
    data = m_dr_sum,
    aes(
      x = cooperation,
      y = est_kc,
      fill = str_sub(ling_cluster, 9, 10)
    ),
    alpha = alpha_med,
    show.legend = FALSE,
    col = NA,
    expand = unit(1.25, "mm")
  ) +
  geom_abline(data = global_est,
              aes(intercept = intercept, slope = slope)) +
  geom_errorbar(data = m_dr_sum,
                aes(x = cooperation, ymin = lwr_kc, ymax = upr_kc),
                alpha = alpha_low) +
  geom_point(
    data = m_dr_sum,
    aes(
      x = cooperation,
      y = est_kc,
      fill = str_sub(ling_cluster, 9, 10)
    ),
    key_glyph = "rect",
    pch = 21,
    col = "black",
    stroke = 0.25,
    size = 1.75
  ) +
  theme_classic(14) +
  labs(
    x = "Cooperative norms (std.)",
    col = "Cultural\ncluster",
    fill = "Cultural\ncluster",
    y = paste0(
      "Combined country-level effect",
      "<br><span style = 'font-size:10pt'>
      on linear predictor for a group lender</span>"
    )
  ) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    axis.title.y = element_markdown(),
    text = element_text(family = fontfam)
  ) +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)


png(paste0("outputs/fig_priorsims_", modelname, ".png"), width = 4000, height = 1700, res = 460)
p17 + p18 + plot_annotation(tag_levels = c("A"))
dev.off()

