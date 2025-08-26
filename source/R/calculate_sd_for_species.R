# Deze functie berekent de nodige sd voor een soort die als input kunnen dienen voor een simulatie

calculate_sd_for_species <- function(
    soortnaam,
    dataset_aantal
) {
  # model voor aantal vissen per waterlichaam, vispunt en datum
  n_fish <- empirical_aantal |>
    filter(.data$soort == soortnaam) |>
    count(.data$waterlichaam, .data$vispunt_id, .data$datum_id) |>
    complete(
      nesting(waterlichaam, datum_id, vispunt_id),
      fill = list(n = 0)
    )
  model_count <- glmmTMB(
    n ~ (1 | waterlichaam) + (1 | datum_id) + (1 | vispunt_id),
    data = n_fish,
    family = poisson(link = "log")
  )
  print(model_count)
  count_intercept <- fixef(model_count)$cond
  count_sd_waterlichaam <- VarCorr(model_count)$cond$waterlichaam
  count_sd_datum <- VarCorr(model_count)$cond$datum_id
  count_sd_vispunt <- VarCorr(model_count)$cond$vispunt_id

  return(
    c(
      count_intercept = count_intercept,
      count_sd_waterlichaam = count_sd_waterlichaam,
      count_sd_datum = count_sd_datum,
      count_sd_vispunt = count_sd_vispunt
    )
  )
}
