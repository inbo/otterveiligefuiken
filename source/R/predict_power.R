predict_power <- function(sims, step_size = 1) {
  stopifnot(require("dplyr"))
  sims |>
    mutate(not = .data$sims - .data$significant) |>
    glm(formula = cbind(significant, not) ~ n_datum, family = binomial) -> model

  # stukje code toegevoegd om ook negatieve trends in rekening te brengen
  # TO
  #direction <- if (max(sims$trend) < 0) -1 else 1   # alle trends negatief? → -1
  direction <- ifelse(max(sims$n_datum) < 0, -1, 1)
  # trend_seq <- seq(
  #   direction * step_size,
  #   direction * 2 * max(abs(sims$trend[sims$simpower < 1])),
  #   by = direction * step_size
  # )
  trend_seq <- #direction *
    seq(step_size, 2 * max(abs(sims$n_datum[sims$simpower < 1])), by = step_size)
  # vermijd trendwaarden die kleiner zijn dan -1 omdat die onmogelijk zijn
  trend_seq <- trend_seq[trend_seq > -1]
  new_data <- data.frame(n_datum = trend_seq)

  new_data |>
    predict(object = model, se.fit = TRUE) -> preds
  new_data |>
    mutate(
      fit = plogis(preds$fit),
      lcl = qnorm(0.025, preds$fit, preds$se.fit) |>
        plogis(),
      ucl = qnorm(0.975, preds$fit, preds$se.fit) |>
        plogis()
    )
}
