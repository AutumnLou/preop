# DCA

## Load Data
roc_dat <- read.csv("data/roc.csv", header = TRUE, row.names = 1)

## Load Libraries
library(dcurves)
library(gridExtra)

# Set Colour Palette and Theme
Palette <- c("#CC79A7", "#E69F00", "gold", "red", "blue", "black")
ggplot2::theme_set(theme_bw(16))

# DCA at 3 years
dca3 <- dca(Surv(surv, time) ~ Cox + LASSO + Boosted + RSF,
            data = roc_dat,
            time = 36,
            thresholds = 1:50 / 100)
# Plot DCA at 3 years
g1 <- as_tibble(dca3) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label), show.legend = FALSE) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", 
              span = 0.2) + ggplot2::scale_colour_manual(values=Palette) + 
  coord_cartesian(ylim = c(-0.10606219431217232, 0.206219431217232
  )) + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "DCA at 3 years", tag = "A", x = "Threshold Probability", y = "Net Benefit", color = "") +
  ggplot2::theme(title = ggplot2::element_text(face = "bold", size = 20),
                 axis.text = ggplot2::element_text(face = "bold", size = 20))

# DCA at 5 years
dca5 <- dca(Surv(surv, time) ~ Cox + LASSO + Boosted + RSF,
            data = roc_dat,
            time = 60,
            thresholds = 1:50 / 100)
# Plot DCA at 5 years
g2 <- as_tibble(dca5) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", 
              span = 0.2) + ggplot2::scale_colour_manual(values=Palette) + 
  coord_cartesian(ylim = c(-0.10606219431217232, 0.206219431217232
  )) + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "DCA at 5 years", tag = "B", x = "Threshold Probability", y = "", color = "") +
  ggplot2::theme(title = ggplot2::element_text(face = "bold", size = 20),
                 axis.text = ggplot2::element_text(face = "bold", size = 20),
                 legend.text = ggplot2::element_text(size =20))

# Plot DCA
grid.arrange(g1 + ggplot2::theme(legend.position = "none"), g2, ncol = 2, widths = c(1, 1.35))
