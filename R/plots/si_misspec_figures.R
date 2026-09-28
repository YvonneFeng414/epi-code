# ==============================================================================
# si_misspec_figures.R
# Base-R figures for the SI-misspecification studies. Each function writes one
# file, PNG or PDF by the extension of `file`. The panels are parameterised so
# the studies draw the same figures with their own overlays: the MALA study
# adds the HPD interval as a dashed line, the EpiFilter study adds its two arms
# as line types.
#
# `metrics` is the expanded table from expand_si_metrics(); the moment-sweep
# figures restrict it to the generating family (gamma_discr) so they show
# exactly what they showed before the family dimension existed.
# ==============================================================================

# x is the relative error in the swept moment, not its absolute value, so the
# two sweeps share an axis and the correctly specified point is at 0 in both.
rel_error_xlab <- "relative error in the assumed SI moment (%)"

sweep_scenarios <- c("vary_sd", "vary_mean")
sweep_legend <- c(vary_sd = "assumed sd wrong", vary_mean = "assumed mean wrong")
sweep_short <- c(vary_sd = "sd wrong", vary_mean = "mean wrong")

# Panel specifications: which metric column, its axis label, panel title, the
# reference line (NA for none) and where the legend goes.
misspec_panels_bias <- list(
  list(col = "Coverage95",  lab = "95% CI coverage", main = "coverage", href = 0.95),
  list(col = "Bias",        lab = "mean bias in Rt", main = "bias",     href = 0),
  list(col = "MeanCIWidth", lab = "mean CI width",   main = "CI width", href = NA)
)

misspec_panels_rmse <- list(
  list(col = "Coverage95",  lab = "95% CI coverage", main = "Coverage", href = 0.95, legend = "bottomright"),
  list(col = "RMSE",        lab = "RMSE of Rt",      main = "RMSE",     href = NA,   legend = "topright"),
  list(col = "MeanCIWidth", lab = "mean CI width",   main = "CI width", href = NA,   legend = "topleft")
)

# The HPD companions of the equal-tailed columns, for hpd = TRUE.
hpd_companion <- c(Coverage95 = "CoverageHPD95", MeanCIWidth = "MeanHPDWidth")

# Open the device the file name asks for. The figures were laid out in pixels
# for png(); a .pdf extension gets the same layout in inches at the same
# resolution, so the two renderings of one figure are the same picture. The
# caller closes the device with dev.off() as before.
open_figure_device <- function(file, width, height, res) {
  if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
    pdf(file, width = width / res, height = height / res)
  } else {
    png(file, width = width, height = height, res = res)
  }
}

# ------------------------------------------------------------------------------
# Each sweep on its own axes, one row per sweep, one panel per metric
# ------------------------------------------------------------------------------
# Coverage carries +/- 2 MCSE bars so the sd sweep can be read against its own
# noise. hpd = TRUE overlays the HPD interval's coverage and width as a dashed
# open-symbol line. `arms` draws one line per arm, distinguished by arm_lty /
# arm_pch, all in the sweep's colour.
plot_misspec_curves <- function(metrics, file, panels, caption, scenario_colors,
                                hpd = FALSE, arms = NULL, arm_lty = NULL,
                                arm_pch = NULL, family = "gamma_discr",
                                width = 1500, height = 1000, res = 130,
                                caption_cex = 0.8) {
  base <- metrics[metrics$family == family, ]
  arm_set <- if (is.null(arms)) list(NULL) else as.list(arms)

  open_figure_device(file, width, height, res)
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3, 1), oma = c(0, 0, 2, 0))

  for (sc in sweep_scenarios) {
    d_all <- base[base$scenario == sc, ]

    for (p in panels) {
      ys <- d_all[[p$col]]
      if (p$col == "Coverage95") {
        ys <- c(d_all$Coverage95 - 2 * d_all$MCSE_Coverage,
                d_all$Coverage95 + 2 * d_all$MCSE_Coverage, 0.95)
      }
      if (hpd && p$col %in% names(hpd_companion)) {
        ys <- c(ys, d_all[[hpd_companion[[p$col]]]])
      }

      plot(NA, xlim = range(d_all$rel_error) * 100, ylim = range(ys, na.rm = TRUE),
           xlab = rel_error_xlab, ylab = p$lab, main = paste0(sc, ": ", p$main))
      if (!is.na(p$href)) abline(h = p$href, lty = 2, col = "grey40")
      abline(v = 0, lty = 3, col = "grey25")

      for (arm in arm_set) {
        d <- if (is.null(arm)) d_all else d_all[d_all$arm == arm, ]
        d <- d[order(d$rel_error), ]
        x <- d$rel_error * 100
        pch <- if (is.null(arm)) 19 else arm_pch[[arm]]
        lty <- if (is.null(arm)) 1 else arm_lty[[arm]]

        if (p$col == "Coverage95") {
          arrows(x, d$Coverage95 - 2 * d$MCSE_Coverage,
                 x, d$Coverage95 + 2 * d$MCSE_Coverage,
                 angle = 90, code = 3, length = 0.03, col = scenario_colors[[sc]])
        }
        lines(x, d[[p$col]], type = "b", pch = pch, lty = lty,
              col = scenario_colors[[sc]])
        if (hpd && p$col %in% names(hpd_companion)) {
          lines(x, d[[hpd_companion[[p$col]]]], type = "b", pch = 1, lty = 2,
                col = scenario_colors[[sc]])
        }
      }

      if (p$col == "Coverage95") {
        if (hpd) {
          legend("bottomright", legend = c("equal-tailed", "HPD"),
                 pch = c(19, 1), lty = c(1, 2), col = scenario_colors[[sc]],
                 bty = "o", bg = "white", box.col = NA, cex = 0.8)
        }
        if (!is.null(arms)) {
          legend("bottom", bty = "n", cex = 0.8, horiz = TRUE, legend = arms,
                 lty = unlist(arm_lty[arms]), pch = unlist(arm_pch[arms]),
                 col = scenario_colors[[sc]])
        }
      }
    }
  }

  mtext(caption, outer = TRUE, cex = caption_cex)
  par(op)
  invisible(dev.off())
}

# ------------------------------------------------------------------------------
# Both sweeps on one axis - the direct "which moment matters more" comparison
# ------------------------------------------------------------------------------
# Bias is left to the per-sweep figure: it stays tiny across the whole grid, so
# it explains almost none of the coverage collapse. Coverage tracks CI width
# instead, which is why the width panel is here.
plot_misspec_relative_error <- function(metrics, file, panels, relative_error,
                                        scenario_colors, arms = NULL,
                                        arm_lty = NULL, arm_pch = NULL,
                                        family = "gamma_discr") {
  base <- metrics[metrics$family == family, ]
  arm_set <- if (is.null(arms)) list(NULL) else as.list(arms)

  open_figure_device(file, 1500, 500, 130)
  op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 3, 1))

  for (p in panels) {
    plot(NA, xlim = range(relative_error) * 100,
         ylim = range(base[[p$col]], p$href, na.rm = TRUE),
         xlab = rel_error_xlab, ylab = p$lab, main = p$main)
    if (!is.na(p$href)) abline(h = p$href, lty = 2, col = "grey40")
    abline(v = 0, lty = 3, col = "grey25")

    for (arm in arm_set) {
      for (sc in sweep_scenarios) {
        d <- base[base$scenario == sc, ]
        if (!is.null(arm)) d <- d[d$arm == arm, ]
        d <- d[order(d$rel_error), ]
        lines(d$rel_error * 100, d[[p$col]], type = "b",
              pch = if (is.null(arm)) 19 else arm_pch[[arm]],
              lty = if (is.null(arm)) 1 else arm_lty[[arm]],
              col = scenario_colors[[sc]])
      }
    }

    if (is.null(arms)) {
      legend(p$legend, legend = unname(sweep_legend[sweep_scenarios]),
             col = scenario_colors[sweep_scenarios], lty = 1, pch = 19, bty = "n")
    }
  }

  if (!is.null(arms)) {
    combos <- expand.grid(sc = sweep_scenarios, arm = arms, stringsAsFactors = FALSE)
    legend("topleft", bty = "n", cex = 0.85,
           legend = paste0(sweep_short[combos$sc], ", ", combos$arm),
           col = scenario_colors[combos$sc],
           lty = unlist(arm_lty[combos$arm]), pch = unlist(arm_pch[combos$arm]))
  }

  par(op)
  invisible(dev.off())
}

# ------------------------------------------------------------------------------
# The kernels themselves, before any result is read
# ------------------------------------------------------------------------------
# No fitting involved. If the families turn out to score alike, this figure is
# what says whether that is because they ARE alike or in spite of not being.
# The short-lag detail is on a log scale: all five agree on mean and sd, so this
# is where they actually differ - and lags 1-3 are the weight the most recent
# days carry into the infectiousness.
plot_si_family_kernels <- function(grid, family_colors, file) {
  fit_families <- grid$fit_families
  true_keys <- vapply(fit_families,
                      function(f) sprintf("%s_%.2f_%.2f", f, grid$mean_si, grid$sd_si),
                      character(1))
  true_sis <- grid$si_list[true_keys]

  open_figure_device(file, 1400, 600, 130)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3, 1), oma = c(0, 0, 2, 0))

  plot(NA, xlim = c(0, 20), ylim = c(0, max(vapply(true_sis, max, numeric(1)))),
       xlab = "serial interval (days)", ylab = "P(SI = lag)",
       main = "Assumed SI at the true moments")
  for (i in seq_along(true_sis)) {
    si <- true_sis[[i]]
    lines(seq.int(0L, length(si) - 1L), si, type = "b", pch = 19, cex = 0.6,
          col = family_colors[[fit_families[i]]])
  }
  legend("topright", legend = fit_families, col = family_colors[fit_families],
         lty = 1, pch = 19, bty = "n", cex = 0.8)

  short <- 1:6
  ymin <- max(1e-5, min(vapply(true_sis, function(s) min(s[short + 1L][s[short + 1L] > 0]),
                               numeric(1))))
  plot(NA, xlim = range(short),
       ylim = c(ymin, max(vapply(true_sis, function(s) max(s[short + 1L]), numeric(1)))),
       log = "y", xlab = "serial interval (days)", ylab = "P(SI = lag), log scale",
       main = "Short lags: where the families disagree")
  for (i in seq_along(true_sis)) {
    y <- true_sis[[i]][short + 1L]
    y[y <= 0] <- NA  # uniform has an exact zero at lag 1; a log axis cannot show it
    lines(short, y, type = "b", pch = 19, cex = 0.7,
          col = family_colors[[fit_families[i]]])
  }
  legend("bottomright", legend = fit_families, col = family_colors[fit_families],
         lty = 1, pch = 19, bty = "n", cex = 0.8)

  mtext(sprintf("Same mean (%.1f) and sd (%.1f), four different shapes - uniform has NO lag-1 mass",
                grid$mean_si, grid$sd_si),
        outer = TRUE, cex = 0.8)
  par(op)
  invisible(dev.off())
}

# ------------------------------------------------------------------------------
# Coverage / RMSE / width by family across both sweeps
# ------------------------------------------------------------------------------
# Open symbols mark settings whose realised moments drifted past drift_tol; the
# x-axis is the REQUESTED moment, and for those points it overstates the error.
# `arm` restricts a multi-arm table to one arm (the EpiFilter study draws the
# smoother alone: five families x two arms would be ten lines a panel).
plot_si_family_curves <- function(metrics, file, panels, grid, family_colors,
                                  caption, arm = NULL) {
  fit_families <- grid$fit_families
  if (!is.null(arm)) metrics <- metrics[metrics$arm == arm, ]

  open_figure_device(file, 1500, 1000, 130)
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3, 1), oma = c(0, 0, 2, 0))

  for (sc in sweep_scenarios) {
    d_all <- metrics[metrics$scenario == sc, ]

    for (p in panels) {
      plot(NA, xlim = range(d_all$rel_error) * 100,
           ylim = range(c(d_all[[p$col]], p$href), na.rm = TRUE),
           xlab = rel_error_xlab, ylab = p$lab, main = paste0(sc, ": ", p$lab))
      if (!is.na(p$href)) abline(h = p$href, lty = 2, col = "grey40")
      abline(v = 0, lty = 3, col = "grey25")

      for (fam in fit_families) {
        d <- d_all[d_all$family == fam, ]
        if (nrow(d) == 0L) next
        d <- d[order(d$rel_error), ]
        lines(d$rel_error * 100, d[[p$col]], type = "l", col = family_colors[[fam]])
        points(d$rel_error * 100, d[[p$col]], pch = ifelse(d$drift_ok, 19, 1),
               cex = 0.8, col = family_colors[[fam]])
      }
      if (identical(p$col, "Coverage95")) {
        legend("bottomright", legend = fit_families, col = family_colors[fit_families],
               lty = 1, pch = 19, bty = "n", cex = 0.7)
      }
    }
  }

  mtext(caption, outer = TRUE, cex = 0.8)
  par(op)
  invisible(dev.off())
}
