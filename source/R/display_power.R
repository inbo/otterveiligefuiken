# grafische weergave van power voor 3 modellen
display_power <- function(sims, power = 0.9, predicted, soortnaam) {
  stopifnot(
    require("dplyr"), require("ggplot2"), require("purrr"), require("scales"), require("patchwork")
  )
  if (missing(predicted)) {
    p <- ggplot(
      sims,
      aes(x = .data$n_datum, ymin = .data$lcl, ymax = .data$ucl)
    ) +
      geom_hline(yintercept = power, linetype = 2) +
      geom_errorbar(aes(colour = .data$sims)) +
      geom_point(aes(y = .data$simpower, colour = .data$sims)) +
      scale_x_continuous("datum", limits = c(0, NA)) +
      scale_y_continuous("power", limits = c(0, 1), labels = percent) +
      scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000)) +
      facet_wrap(~modeltype)

    print(p)
    return(invisible(NULL))
  }

  pred_high <- min(predicted[["lengte"]]$n_datum[power <= predicted[["lengte"]]$lcl])
  pred_low <- min(predicted[["lengte"]]$n_datum[power <= predicted[["lengte"]]$ucl])
  pred_fit <- min(predicted[["lengte"]]$n_datum[power <= predicted[["lengte"]]$fit])

  pl <- ggplot(
    sims |> filter(modeltype == "1_lengte"),
    aes(x = .data$n_datum, ymin = .data$lcl, ymax = .data$ucl)
  ) +
    geom_hline(yintercept = power, linetype = 2) +
    geom_errorbar(aes(colour = .data$sims)) +
    geom_point(aes(y = .data$simpower, colour = .data$sims)) +
    scale_x_continuous("datum", limits = c(0, NA)) +
    scale_y_continuous("power", limits = c(0, 1), labels = percent) +
    scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000), guide = "none") +
    geom_rect(
      xmin = pred_low,
      xmax = pred_high,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.05
    ) +
    geom_vline(xintercept = pred_fit, linetype = 3) +
    geom_ribbon(data = predicted[["lengte"]], alpha = 0.2, fill = "darkgreen") +
    geom_line(data = predicted[["lengte"]], aes(y = .data$fit), colour = "darkgreen") +
    ggtitle(
      sprintf(
        #"kleinste aantal datums voor aantonen verschil in lengte bij %s: %.2f%% (%.2f%%; %.2f%%)",
        "lengte bij %s:\n %.0f (%.0f; %.0f)",
        soortnaam,
        ceiling(pred_fit),
        floor(pred_low),
        ceiling(pred_high)
      )
    )

  if (!is.null(predicted[["lengte_aantal"]])) {
    pred_high <- min(predicted[["lengte_aantal"]]$n_datum[power <= predicted[["lengte_aantal"]]$lcl])
    pred_low <- min(predicted[["lengte_aantal"]]$n_datum[power <= predicted[["lengte_aantal"]]$ucl])
    pred_fit <- min(predicted[["lengte_aantal"]]$n_datum[power <= predicted[["lengte_aantal"]]$fit])

    pln <- ggplot(
      sims |> filter(modeltype == "2_lengte_aantal"),
      aes(x = .data$n_datum, ymin = .data$lcl, ymax = .data$ucl)
    ) +
      geom_hline(yintercept = power, linetype = 2) +
      geom_errorbar(aes(colour = .data$sims)) +
      geom_point(aes(y = .data$simpower, colour = .data$sims)) +
      scale_x_continuous("datum", limits = c(0, NA)) +
      scale_y_continuous("power", limits = c(0, 1), labels = percent) +
      scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000), guide = "none") +
      geom_rect(
        xmin = pred_low,
        xmax = pred_high,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.05
      ) +
      geom_vline(xintercept = pred_fit, linetype = 3) +
      geom_ribbon(data = predicted[["lengte_aantal"]], alpha = 0.2, fill = "darkgreen") +
      geom_line(data = predicted[["lengte_aantal"]], aes(y = .data$fit), colour = "darkgreen") +
      ggtitle(
        sprintf(
          #"kleinste aantal datums voor aantonen verschil in lengte en aantal bij %s: %.2f%% (%.2f%%; %.2f%%)",
          "lengte en aantal bij \n%s:\n %.0f (%.0f; %.0f)",
          soortnaam,
          ceiling(pred_fit),
          floor(pred_low),
          ceiling(pred_high)
        )
      )
  }

  if (!is.null(predicted[["aantal"]])) {
    pred_high <- min(predicted[["aantal"]]$n_datum[power <= predicted[["aantal"]]$lcl])
    pred_low <- min(predicted[["aantal"]]$n_datum[power <= predicted[["aantal"]]$ucl])
    pred_fit <- min(predicted[["aantal"]]$n_datum[power <= predicted[["aantal"]]$fit])

    pn <- ggplot(
      sims |> filter(modeltype == "3_aantal"),
      aes(x = .data$n_datum, ymin = .data$lcl, ymax = .data$ucl)
    ) +
      geom_hline(yintercept = power, linetype = 2) +
      geom_errorbar(aes(colour = .data$sims)) +
      geom_point(aes(y = .data$simpower, colour = .data$sims)) +
      scale_x_continuous("datum", limits = c(0, NA)) +
      scale_y_continuous("power", limits = c(0, 1), labels = percent) +
      scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000)) +
      geom_rect(
        xmin = pred_low,
        xmax = pred_high,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.05
      ) +
      geom_vline(xintercept = pred_fit, linetype = 3) +
      geom_ribbon(data = predicted[["aantal"]], alpha = 0.2, fill = "darkgreen") +
      geom_line(data = predicted[["aantal"]], aes(y = .data$fit), colour = "darkgreen") +
      ggtitle(
        sprintf(
          #"kleinste aantal datums voor aantonen verschil in aantal bij %s: %.2f%% (%.2f%%; %.2f%%)",
          "aantal bij %s:\n %.0f (%.0f; %.0f)",
          soortnaam,
          ceiling(pred_fit),
          floor(pred_low),
          ceiling(pred_high)
        )
      )
  }

  if (exists("pln") & exists("pn")) {
    print(pl + pln + pn)
  } else {
    print(pl)
  }

  return(invisible(NULL))
}
