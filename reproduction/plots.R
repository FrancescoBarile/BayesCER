library(tidyverse)
library(ggh4x)

n_vec = c( 10, 100, 1000)
N_vec = c( 10, 100, 1000, 2000)
alpha_vec = c(0.1, 0.2, 0.3)

n_labels <- setNames(
  paste0("n = ", n_vec),
  as.character(n_vec)
)

N_labels <- setNames(
  paste0("N = ", N_vec),
  as.character(N_vec)
)

alpha_labels <- setNames(
  paste0("alpha == ", alpha_vec),
  as.character(alpha_vec)
)

library(BayesCER)
data("sim_data_result", package = "BayesCER")

output_df_long = sim_data_result
output_df_long$method <- factor(output_df_long$method,levels=c("EB","Gibbs"))

paper_colors <- c(EB='#F8766D',Gibbs='#00BFC4')
paper_outlines <- setNames(colorspace::darken(unname(paper_colors),.2),names(paper_colors))


## PLOT MAIN
# ECC
pp <- ggplot(
  output_df_long %>% filter(alpha==0.2),
  aes(
    x = as.factor(n),
    y = ecc_hat,
    fill = method, colour = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    . ~ N,
    scales = "free_y",
    independent = "y",
    axes = "y",
    labeller = labeller(
      N = as_labeller(N_labels)
    )
  ) +
  geom_hline(
    data=distinct(output_df_long,alpha,N,ecc_true) %>% filter(alpha==0.2),
    aes(yintercept = ecc_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "", y = "", fill = "Method") +
  theme_bw() +
  theme(
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14)
  ) +
  theme(plot.margin = margin(2, 2, 0, 2),
        panel.spacing.y = unit(0, "mm")
  )
pp <- pp + BayesCER:::paper_style(main=TRUE) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
pp

# EDGE DENSITY
pp1 <- ggplot(
  output_df_long %>% filter(alpha==0.2),
  aes(
    x = as.factor(n),
    y = density,
    fill = method, colour = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    . ~ N,
    scales = "free_y",
    independent = "y",
    axes = "y",
    labeller = labeller(
      N = function(x) rep("", length(x))
    )
  ) +
  geom_hline(
    data=distinct(output_df_long,alpha,N,density_true) %>% filter(alpha==0.2),
    aes(yintercept = density_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14)
  ) +
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(-20, 2, 2, 2),
    panel.spacing.y = unit(0, "mm")
  )
pp1 <- pp1 + BayesCER:::paper_style(main=TRUE) +
  theme(strip.text.x=element_blank(),strip.background.x=element_blank(), legend.position="none")
pp1



library(grid)
library(gridExtra)
library(cowplot)

leg <- get_legend(
  pp + theme(
    legend.position = "right"
  )
)

g1 <- ggplotGrob(pp + theme(legend.position = "none"))
g2 <- ggplotGrob(pp1 + theme(legend.position = "none"))

panel1 <- unique(g1$layout$t[grepl("^panel", g1$layout$name)])
panel2 <- unique(g2$layout$t[grepl("^panel", g2$layout$name)])

strip2 <- unique(g2$layout$t[grepl("^strip-t", g2$layout$name)])

if (length(strip2) > 0) {
  g2$heights[strip2] <- unit(0, "mm")
}

panel_height <- unit(40, "mm")
g1$heights[panel1] <- panel_height
g2$heights[panel2] <- panel_height

w <- unit.pmax(g1$widths, g2$widths)
g1$widths <- w
g2$widths <- w

g1$heights[length(g1$heights)] <- unit(0, "mm")
g2$heights[1] <- unit(0, "mm")

plots <- arrangeGrob(
  g1,
  g2,
  ncol = 1,
  heights = unit.c(
    sum(g1$heights),
    sum(g2$heights)
  )
)

p <- arrangeGrob(
  plots,
  leg,
  ncol = 2,
  widths = unit.c(
    unit(1, "null"),
    grobWidth(leg) + unit(2, "mm")
  )
)

grid.newpage()
grid.draw(p)



## PLOTS SUPPLEMENT

# DISCREPANCY OF ESTIMATE OF C
pp <- ggplot(
  output_df_long,
  aes(
    x = as.factor(n),
    y = discrepancy_C_hat,
    fill = method, colour = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    alpha ~ N,
    scales = "free_y",
    remove = "y",
    axes = "margins",
    labeller = labeller(
      N = as_labeller(N_labels),
      alpha = as_labeller(alpha_labels, default = label_parsed)
    )
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp <- pp + BayesCER:::paper_style(error=TRUE)
pp

# DISCREPANCY OF ESTIMATE OF ALPHA
pp <- ggplot(
  output_df_long,
  aes(
    x = as.factor(n),
    y = discrepancy_alpha_hat,
    fill = method, colour = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    alpha ~ N,
    scales = "free_y",
    remove = "y",
    axes = "margins",
    labeller = labeller(
      N = as_labeller(N_labels),
      alpha = as_labeller(alpha_labels, default = label_parsed)
    )
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp <- pp + BayesCER:::paper_style(error=TRUE)
pp

# EXPECTED CLUSTERING COEFFICIENT
pp <- ggplot(
  output_df_long %>% filter(alpha!=0.2),
  aes(
    x = as.factor(n),
    y = ecc_hat,
    fill = method, colour = method
  )) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    alpha ~ N,
    scales = "free_y",
    remove = "y",
    axes = "margins",
    labeller = labeller(
      N = as_labeller(N_labels),
      alpha = as_labeller(alpha_labels, default = label_parsed)
    )
  ) +
  geom_hline(
    data=distinct(output_df_long,alpha,N,ecc_true) %>% filter(alpha!=0.2),
    aes(yintercept = ecc_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp <- pp + BayesCER:::paper_style(error=FALSE)
pp

# EXPECTED EDGE DENSITY
pp <- ggplot(
  output_df_long %>% filter(alpha!=0.2),
  aes(
    x = as.factor(n),
    y = density,
    fill = method, colour = method
  )) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.62,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(position=position_jitterdodge(jitter.width=.13,jitter.height=0,dodge.width=.8,seed=110),
             size=.55,alpha=.22,show.legend=FALSE) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.2,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_grid2(
    alpha ~ N,
    scales = "free_y",
    remove = "y",
    axes = "margins",
    labeller = labeller(
      N = as_labeller(N_labels),
      alpha = as_labeller(alpha_labels, default = label_parsed)
    )
  ) +
  geom_hline(
    data=distinct(output_df_long,alpha,N,density_true) %>% filter(alpha!=0.2),
    aes(yintercept = density_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(1,1,1,1, unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp <- pp + BayesCER:::paper_style(error=FALSE)
pp


