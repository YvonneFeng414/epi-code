################################################################################
# Controlled I0 experiment: I0 = 2000 against I0 = 500.
#
# Everything except I0 is held fixed - same truth per arm, set.seed(123),
# scoring over days 14-253, the +/-60% grid, the gamma_bin assumed-SI family,
# 100 replicates, and filter_start = 8 for EpiFilter in BOTH arms.
#
# WHY THE QUESTION. I0 = 2000 is the seed study's sens_I0[1] and is what the
# aligned runs use. It has one known cost: a single 2000-case seeded day makes
# day 1 a spike and day 2 a trough, and at assumed mean +60% the resulting
# total infectiousness is near zero on days 3-6 while incidence is already in
# the hundreds. The implied Rt reaches 8794 against an EpiFilter grid capped at
# 10, the Poisson likelihood underflows on every grid point, and the recursion
# returns NaN. Dropping I0 removes the cliff - but it also shrinks the whole
# simulated epidemic, and the repo's own notes warn that a series far below the
# observed scale stops describing Queens at all.
#
# So the two arms answer: how much does the epidemic scale change, and how much
# do the metrics move with it?
#
# Colour is the estimator, linetype and point fill are I0, so a reader can see
# both dimensions at once rather than flipping between two figures.
#
# Run from queens_rt/:  Rscript combined/experiments/make_i0_panel.R
################################################################################

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })

base    <- "combined/experiments"
out_dir <- base
fam     <- "gamma_bin"
i0_vals <- c(2000, 500)

arms <- list(
  list(label = "EpiEstim",     path = "epiestim_arm_metrics.csv",            arm = NULL),
  list(label = "EpiLPS (MAP)", path = "epilps_map/si_misspec_metrics.csv",   arm = NULL),
  list(label = "EpiFilter",    path = "epifilter/si_misspec_metrics.csv",    arm = "smoother")
)

take <- function(i0, a) {
  f <- file.path(base, paste0("I0-", i0), a$path)
  if (!file.exists(f)) stop("Missing ", f, " - run run_i0_experiment.sh ", i0, " first.")
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[grepl(paste0("^", fam, "_"), d$key), , drop = FALSE]
  if (!is.null(a$arm)) d <- d[d$arm == a$arm, , drop = FALSE]
  data.frame(I0 = i0, method = a$label, scenario = d$scenario, rel_error = d$rel_error,
             Bias = d$Bias, RMSE = d$RMSE, Coverage = d$Coverage95,
             CIWidth = d$MeanCIWidth, stringsAsFactors = FALSE)
}

combined <- do.call(rbind, lapply(i0_vals, function(i0)
  do.call(rbind, lapply(arms, function(a) take(i0, a)))))

key_order <- c("EpiEstim", "EpiFilter", "EpiLPS (MAP)")
combined$method   <- factor(combined$method, levels = key_order)
combined$I0       <- factor(combined$I0, levels = i0_vals,
                            labels = paste0("I0 = ", i0_vals))
combined$scenario <- factor(combined$scenario, levels = c("vary_mean", "vary_sd"),
                            labels = c("SI mean misspecified (SD fixed)",
                                       "SI SD misspecified (mean fixed)"))
combined$rel_pct  <- 100 * combined$rel_error

write.csv(combined, file.path(out_dir, "i0_comparison_metrics.csv"), row.names = FALSE)

# Any cell that failed to produce a number, named rather than left as a gap.
full <- expand.grid(I0 = levels(combined$I0), method = key_order,
                    scenario = levels(combined$scenario),
                    rel_pct = sort(unique(combined$rel_pct)), stringsAsFactors = FALSE)
have <- paste(combined$I0, combined$method, combined$scenario, combined$rel_pct)
gaps <- full[!paste(full$I0, full$method, full$scenario, full$rel_pct) %in% have, ]
if (nrow(gaps) > 0L) {
  cat("\nMISSING CELLS (gaps in the line):\n")
  for (i in seq_len(nrow(gaps)))
    cat(sprintf("  %-12s %-14s %-32s %+.0f%%\n", gaps$I0[i], gaps$method[i],
                gaps$scenario[i], gaps$rel_pct[i]))
} else cat("\nNo missing cells - every setting produced a number in both arms.\n")

method_colours <- c("EpiEstim" = "#1f77b4", "EpiFilter" = "#ff7f0e",
                    "EpiLPS (MAP)" = "#2ca02c")
method_shapes  <- c("EpiEstim" = 16, "EpiFilter" = 15, "EpiLPS (MAP)" = 17)
i0_lty         <- c("solid", "22")

row_fill <- c("SI mean misspecified (SD fixed)" = "#FBE4D5",
              "SI SD misspecified (mean fixed)" = "#DEEAF6")

metrics <- list(
  list(col = "Bias",     short = "Bias",            y = expression("Bias of " * R[t]), hline = 0),
  list(col = "RMSE",     short = "RMSE",            y = expression("RMSE of " * R[t]), hline = NA),
  list(col = "Coverage", short = "95% CI coverage", y = "95% CI coverage",             hline = 0.95)
)
lims <- lapply(metrics, function(m) {
  r <- range(combined[[m$col]], na.rm = TRUE)
  if (!is.na(m$hline)) r <- range(c(r, m$hline))
  r + c(-1, 1) * 0.08 * diff(r)
})
x_breaks <- sort(unique(combined$rel_pct))

panel_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

make_panel <- function(m, sc, letter, lim) {
  d <- combined[combined$scenario == sc, ]
  d$strip <- sprintf("%s. %s (%s)", letter, m$short,
                     sub(" \\(.*", "", sub("SI ", "SI ", sc)))
  p <- ggplot(d, aes(x = rel_pct, y = .data[[m$col]], colour = method,
                     shape = method, linetype = I0, group = interaction(method, I0)))
  if (!is.na(m$hline))
    p <- p + geom_hline(yintercept = m$hline, linetype = "dashed",
                        colour = "grey45", linewidth = 0.35)
  p +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45", linewidth = 0.35) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
    facet_wrap(~ strip) +
    scale_colour_manual(values = method_colours, limits = key_order, drop = FALSE) +
    scale_shape_manual(values = method_shapes, limits = key_order, drop = FALSE) +
    scale_linetype_manual(values = i0_lty, limits = levels(combined$I0), drop = FALSE) +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(ylim = lim) +
    labs(x = "Relative error in the assumed SI parameter (%)", y = m$y) +
    panel_theme +
    theme(strip.background = element_rect(fill = row_fill[[sc]], colour = "grey30",
                                          linewidth = 0.4),
          strip.text = element_text(size = 8.5, face = "bold", hjust = 0,
                                    margin = margin(3, 3, 3, 5))) +
    guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid")),
           shape  = guide_legend(order = 1),
           linetype = guide_legend(order = 2, override.aes = list(colour = "grey20", shape = NA)))
}

scs  <- levels(combined$scenario)
lets <- matrix(LETTERS[1:6], nrow = 2, byrow = TRUE)
rows <- lapply(seq_along(scs), function(i)
  wrap_plots(lapply(seq_along(metrics), function(j)
    make_panel(metrics[[j]], scs[i], lets[i, j], lims[[j]])), nrow = 1))

page <- wrap_plots(rows, ncol = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), legend.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 2, 0))

pdf(file.path(out_dir, "i0_comparison_panel.pdf"), width = 12, height = 7.2)
print(page); invisible(dev.off())
png(file.path(out_dir, "i0_comparison_panel.png"), width = 12, height = 7.2,
    units = "in", res = 200)
print(page); invisible(dev.off())

cat("\n===== Control point (correct SI) by I0 =====\n")
c0 <- combined[combined$rel_error == 0 &
               combined$scenario == "SI SD misspecified (mean fixed)", ]
c0 <- c0[order(c0$method, c0$I0), ]
print(data.frame(method = c0$method, I0 = c0$I0, Coverage = round(c0$Coverage, 4),
                 Bias = round(c0$Bias, 4), RMSE = round(c0$RMSE, 4),
                 Width = round(c0$CIWidth, 4)), row.names = FALSE)

cat("\n===== What moved when I0 went 2000 -> 500 =====\n")
w <- combined %>%
  tidyr::pivot_wider(id_cols = c(method, scenario, rel_error),
                     names_from = I0, values_from = c(Coverage, RMSE, CIWidth))
names(w) <- make.names(names(w))
w %>% group_by(method) %>%
  summarise(dCoverage = median(`Coverage_I0...500` - `Coverage_I0...2000`, na.rm = TRUE),
            RMSE_ratio = median(`RMSE_I0...500` / `RMSE_I0...2000`, na.rm = TRUE),
            Width_ratio = median(`CIWidth_I0...500` / `CIWidth_I0...2000`, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.frame() %>%
  (function(x) print(data.frame(lapply(x, function(v)
    if (is.numeric(v)) round(v, 4) else v)), row.names = FALSE))

cat(sprintf("\nWrote %s/i0_comparison_panel.pdf, .png and i0_comparison_metrics.csv\n", out_dir))
