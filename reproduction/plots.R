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

output_df_long = readRDS( paste0(getwd(), "/data/output_df_long.rds"))


## PLOT MAIN
# ECC
pp <- ggplot(
  output_df_long %>% filter(alpha==0.2),
  aes(
    x = as.factor(n),
    y = ecc_hat,
    fill = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = FALSE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
    aes(yintercept = ecc_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "", y = "", fill = "Method") +
  theme_bw() +
  theme(
    plot.margin = margin(rep(1,4), unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14)
  ) +
  theme(plot.margin = margin(2, 2, 0, 2),
        panel.spacing.y = unit(0, "mm")
  ) +
  scale_y_continuous(
    labels = function(x) sprintf("%.3f", x)
  )
pp

# EDGE DENSITY
pp1 <- ggplot(
  output_df_long %>% filter(alpha==0.2),
  aes(
    x = as.factor(n),
    y = total_degree_norm,
    fill = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = FALSE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
    aes(yintercept = total_degree_true_norm),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    plot.margin = margin(rep(1,4), unit="mm"),
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
  ) +
  scale_y_continuous(
    labels = function(x) sprintf("%.3f", x)
  )
pp1


library(grid)
library(gridExtra)

g1 <- ggplotGrob(pp)
g2 <- ggplotGrob(pp1)

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

p <- arrangeGrob(
  g1,
  g2,
  ncol = 1,
  heights = unit.c(
    sum(g1$heights),
    sum(g2$heights)
  )
)

grid.newpage()
grid.draw(p)

ggsave( paste0(out_path, "/ecc_ecd.pdf"),
        plot = p,
        width = 200,
        height = 145,
        units = "mm",
        device = cairo_pdf
)

## PLOTS SUPPLEMENT

# DISCREPANCY OF ESTIMATE OF C
pp <- ggplot(
  output_df_long,
  aes(
    x = as.factor(n),
    y = discrepancy_C_hat_norm,
    fill = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
  facetted_pos_scales(
    y = list(
      alpha == 0.1 ~ scale_y_continuous(limits = c(0, 0.003)),
      TRUE ~ scale_y_continuous()
    )
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(rep(1,4), unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp

ggsave( paste0(out_path, "/discrepancy_C_hat.pdf"),
        plot = pp,
        width = 200,
        height = 150,
        units = "mm",
        device = cairo_pdf
)


# DISCREPANCY OF ESTIMATE OF ALPHA
pp <- ggplot(
  output_df_long,
  aes(
    x = as.factor(n),
    y = discrepancy_alpha_hat,
    fill = method
  )
) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
    plot.margin = margin(rep(1,4), unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp

ggsave( paste0(out_path, "/discrepancy_alpha_hat.pdf"),
        plot = pp,
        width = 200,
        height = 150,
        units = "mm",
        device = cairo_pdf
)


# EXPECTED CLUSTERING COEFFICIENT
pp <- ggplot(
  output_df_long %>% filter(alpha!=0.2),
  aes(
    x = as.factor(n),
    y = ecc_hat,
    fill = method
  )) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
    aes(yintercept = ecc_true),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(rep(1,4), unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp

ggsave( paste0(out_path, "/ECC.pdf"),
        plot = pp,
        width = 200,
        height = 100,
        units = "mm",
        device = cairo_pdf
)


# EXPECTED EDGE DENSITY
pp <- ggplot(
  output_df_long %>% filter(alpha!=0.2),
  aes(
    x = as.factor(n),
    y = total_degree_norm,
    fill = method
  )) +
  geom_boxplot(
    aes(colour = after_scale(colorspace::darken(fill, 0.2))),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.3,
    show.legend = c(colour = FALSE, fill = TRUE)
  ) +
  stat_boxplot(
    aes(colour = method),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    show.legend = FALSE
  ) +
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
    aes(yintercept = total_degree_true_norm),
    colour = "black",
    linetype = "dashed"
  ) +
  labs(x = "Sample size n", y = "", fill = "Method") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "mm"),
    plot.margin = margin(rep(1,4), unit="mm"),
    strip.text = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text.y.right  = element_text(size = 14),
    strip.text.x.top    = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
pp

ggsave(  paste0(out_path, "/degree.pdf"),
         plot = pp,
         width = 200,
         height = 100,
         units = "mm",
         device = cairo_pdf
)


