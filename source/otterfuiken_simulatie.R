library(tidyverse)
library(glmmTMB)
library(duckdb)
library(showtext)
library(scales)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

#### FUNCTIES VOOR SIMULATIES ####

sim_linear <- function(
    # meetnetparameters
    fuiktype = c("SF", "SFO"),     # 2 fuiktypes, waarbij SFO otterproof
    n_waterlichaam = 3,   # aantal waterlichamen waarin wordt gevangen
    n_sample = 5,     # aantal fuiken per waterlichaam

    # populatieparameters
    initialdensity = exp(8.864),   # model‑intercept
    trend = -0.174,                 # fuiktypeSFO

    # ruisparameters
    #overdispersion = 4.87,
    sd_waterlichaam = sqrt(5.375),

    # simulatieparameters
    n_sim = 100
) {
  stopifnot(require(glmmTMB), require(tidyverse))

  # design
  expand.grid(
    fuiktype = fuiktype,
    waterlichaam_id = seq_len(n_waterlichaam)
  ) -> design

  # simulatie
  replicate(n_sim, {
    design |>
      mutate(
        eta = log(initialdensity) + log(1 + trend) * .data$fuiktype +     # eta: verwachte waarde in log-schaal
          rnorm(n_waterlichaam, 0, sd_waterlichaam)[.data$waterlichaam_id],         # effect van waterlichaam
        count = rnbinom(n(), mu = exp(.data$eta), size = overdispersion)
      ) -> dataset

    # fit glmm op gesimuleerde data
    m <- try(
      glmmTMB(
        aantal_ind ~ fuiktype + (1 | waterlichaam_id),   # simpel model
        data = dataset, family = nbinom2()
      )
    )

    if (inherits(m, "try-error") || any(is.na(m$sdr$cov.fixed))) {
      return(numeric(0))
    }
    if (!inherits(m, "glmmTMB")) {
      stop(class(m))
    }

    # haal p-waarde voor jaar-effect uit summary van model
    coef(summary(m))$cond["year", "Pr(>|z|)"]
  }) |>
    unlist() |>
    na.omit()
}


# wat doet bovenstaande functie? op basis van opgegeven parameters wordt data gesimuleerd, een model wordt gefit op de data en er wordt gekeken of het jaar-effect gedetecteerd wordt via de p-waarde. Resultaat van de functie is een lijst met n p-waarden (n is het aantal replicaties). Deze lijst met p-waarden kan vervolgens gebruikt worden om een power te berekenen als het aandeel significante p-waarden op het totaal aantal simulaties. vb. data werd gesimuleerd met een gekende trend van -5%. We simuleren de data 100 keer en zien dat het model slechts 30 keer een significant jaar-effect detecteert. De power van de proefopzet om de trend van -5% te meten is dus gelijk aan 30%

display_power <- function(sims, power = 0.9, predicted) {
  stopifnot(
    require("dplyr"), require("ggplot2"), require("purrr"), require("scales")
  )
  p <- ggplot(
    sims,
    aes(x = .data$trend, ymin = .data$lcl, ymax = .data$ucl)
  ) +
    geom_hline(yintercept = power, linetype = 2) +
    geom_errorbar(aes(colour = .data$sims)) +
    geom_point(aes(y = .data$simpower, colour = .data$sims)) +
    #scale_x_continuous("trend", limits = c(0, NA), labels = percent) +
    scale_y_continuous("power", limits = c(0, 1), labels = percent) +
    scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000))
  if (missing(predicted)) {
    print(p)
    return(invisible(NULL))
  }

  if (max(predicted$trend) < 0) {
    pred_high <- max(predicted$trend[power <= predicted$ucl])
    pred_low <- max(predicted$trend[power <= predicted$lcl])
    pred_fit <- max(predicted$trend[power <= predicted$fit])
  } else {
    pred_high <- min(predicted$trend[power <= predicted$lcl])
    pred_low <- min(predicted$trend[power <= predicted$ucl])
    pred_fit <- min(predicted$trend[power <= predicted$fit])
  }

  p <- p +
    geom_rect(
      xmin = pred_low,
      xmax = pred_high,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.05
    ) +
    geom_vline(xintercept = pred_fit, linetype = 3) +
    geom_ribbon(data = predicted, alpha = 0.2, fill = "darkgreen") +
    geom_line(data = predicted, aes(y = .data$fit), colour = "darkgreen") +
    ggtitle(
      sprintf(
        "smallest detectable trend: %.2f%% (%.2f%%; %.2f%%)",
        100 * pred_fit,
        100 * pred_low,
        100 * pred_high
      )
    )
  print(p)
  return(invisible(NULL))
}

predict_power <- function(sims, step_size = 0.001) {
  stopifnot(require("dplyr"))
  sims |>
    mutate(not = .data$sims - .data$significant) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model

  # stukje code toegevoegd om ook negatieve trends in rekening te brengen
  # TO
  #direction <- if (max(sims$trend) < 0) -1 else 1   # alle trends negatief? → -1
  direction <- ifelse(max(sims$trend) < 0, -1, 1)
  # trend_seq <- seq(
  #   direction * step_size,
  #   direction * 2 * max(abs(sims$trend[sims$simpower < 1])),
  #   by = direction * step_size
  # )
  trend_seq <- direction *
    seq(step_size, 2 * max(abs(sims$trend[sims$simpower < 1])), by = step_size)
  # vermijd trendwaarden die kleiner zijn dan -1 omdat die onmogelijk zijn
  trend_seq <- trend_seq[trend_seq > -1]
  new_data <- data.frame(trend = trend_seq)

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


# Wat voor trend kunnen we wel zien?
## Gegeven een bepaalde meetopzet (aantal jaren, transecten, meetfout, enz.), wat is de minimale trend die met een opgegeven power (hier gebruiken we een power van 90%) kan worden gedetecteerd?

### VOOR POSITIEVE TRENDS
estimate_trend_pos <- function(
    # meetnetparameters
    fuiktype = c("SF", "SFO"),     # 2 fuiktypes, waarbij SFO otterproof
    n_waterlichaam = 3,   # aantal waterlichamen waarin wordt gevangen
    n_sample = 5,     # aantal fuiken per waterlichaam

    # populatieparameters
    initialdensity = exp(8.864),   # model‑intercept

    # ruisparameters
    #overdispersion = 4.87,
    sd_waterlichaam = sqrt(5.375),

    # powerparameters
    alpha = 0.10,   # significantieniveau (α)
    power = 0.90,   # gewenste statistische power
    step_size = 0.001,

    connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = FALSE
  )
) {
  stopifnot(
    require("assertthat"), require("DBI"), require("dplyr"), require("duckdb"),
    require("purrr"), require("tidyr")
  )

  # valideer inputwaarden
  assert_that(

    # discrete variabelen
    is.character(fuiktype),
    all(fuiktype %in% c("SF", "SFO")),
    is.count(n_sample),
    is.count(n_waterlichaam),

    # continue variabelen
    is.number(initialdensity),
    is.number(sd_waterlichaam), noNA(sd_waterlichaam), sd_waterlichaam > 0,
    #is.number(overdispersion), noNA(overdispersion), overdispersion > 0,

    # alfa en power begrensd tussen 0 en 1
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )

  # tabel "trend" in DuckDB initialiseren indien deze nog niet bestaat
  if (!"trend" %in% dbListTables(conn = connection)) {
    #   dbRemoveTable(conn = connection, name = "trend")
    data.frame(
      trend = numeric(0),
      fuiktype = character(0),
      n_waterlichaam = integer(0),
      n_sample = integer(0),
      initialdensity = numeric(0),
      overdispersion = numeric(0),
      sd_waterlichaam = numeric(0),
      #sd_day = numeric(0),
      p = numeric(0)
    ) |>
      dbCreateTable(conn = connection, name = "trend")
  }

  # qery: reeds bestaande simulaties op basis van gespecifieerde simulatieparameters ophalen
  # per trend wordt het totale aantal simulaties (sims) en het aantal significante simulaties p < alpha gegeven
  # omdat sd_transect, sd_day, initialdensity, enz. komma-getallen zijn (float), gebruik je een tolerantie — anders afrondingsfouten. Zoek "ongeveer gelijk aan", binnen een marge van 0.0001.
  query <- sprintf(
    "SELECT trend,
          SUM(p < %f) AS significant,
          COUNT(trend) AS sims
   FROM trend
   WHERE ABS(sd_waterlichaam    - %f) < 1e-4
     --AND ABS(sd_day         - %f) < 1e-4
     AND n_waterlichaam         = %d
     AND n_sample           = %d
     AND fuiktype           = %s
     AND ABS(initialdensity - %f) < 1e-4
     --AND ABS(overdispersion - %f) < 1e-4
   GROUP BY trend
   ORDER BY trend;",
    alpha,
    sd_waterlichaam,
    #sd_day,
    n_waterlichaam,
    n_sample,
    fuiktype,
    initialdensity #,
    #overdispersion
  )

  dbGetQuery(conn = connection, statement = query) |>
    mutate(simpower = .data$significant / .data$sims) -> sims

  # start simulatie indien query geen resultaat gaf
  if (nrow(sims) == 0) {
    message("initializing to 0.05")
    data.frame(
      trend = 0.05,
      fuiktype = fuiktype,
      n_waterlichaam = n_waterlichaam,
      n_sample = n_sample,
      initialdensity = initialdensity,
      #overdispersion = overdispersion,
      sd_waterlichaam = sd_waterlichaam,
      #sd_day = sd_day,
      p = sim_linear(
        trend = 0.05,
        fuiktype = fuiktype,
        n_waterlichaam = n_waterlichaam,
        n_sample = n_sample,
        initialdensity = initialdensity,
        #overdispersion = overdispersion,
        sd_waterlichaam = sd_waterlichaam
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  # OPSCHALEN: zolang de hoogste power < 90%, verdubbel maximale trend
  while (max(sims$simpower) <= power) {
    extra <- 2 * max(sims$trend)
    message("increasing to ", extra)
    data.frame(
      trend = extra,
      fuiktype = fuiktype,
      n_waterlichaam = n_waterlichaam,
      n_sample = n_sample,
      #overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_waterlichaam = sd_waterlichaam,
      p = sim_linear(
        trend = 0.05,
        fuiktype = fuiktype,
        n_waterlichaam = n_waterlichaam,
        n_sample = n_sample,
        initialdensity = initialdensity,
        #overdispersion = overdispersion,
        sd_waterlichaam = sd_waterlichaam
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  # AFSCHALEN - als de kleinste, gesimuleerde trend een relatief hoge power (power / 2) heeft, wordt deze gehalveerd.
  while (
    min(sims$trend) >= 2 * step_size && min(sims$simpower) >= (power / 2)
  ) {
    extra <- round(min(sims$trend) / 2, -floor(log10(step_size)))
    stopifnot(extra > 0)
    message("decreasing to ", extra)
    data.frame(
      trend = extra,
      fuiktype = fuiktype,
      n_waterlichaam = n_waterlichaam,
      n_sample = n_sample,
      #overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_waterlichaam = sd_waterlichaam,
      p = sim_linear(
        trend = 0.05,
        fuiktype = fuiktype,
        n_waterlichaam = n_waterlichaam,
        n_sample = n_sample,
        initialdensity = initialdensity,
        #overdispersion = overdispersion,
        sd_waterlichaam = sd_waterlichaam
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  stopifnot(nrow(sims) > 1)
  predicted <- predict_power(sims = sims, step_size = step_size)

  # Voor welke trendwaarden overlapt het CI nog met 90%; hier extra simulaties uitvoeren om de schatting te verfijnen
  predicted |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      predicted |>
        filter(.data$lcl > power) |>   # selecteer trendwaarden waarvan het 95%-interval de drempelwaarde (power) overlapt
        slice_min(.data$trend, n = 1),   # voeg de eerstvolgende trend toe waarbij het hele interval boven de powerdrempel ligt
      sims |>
        filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
        select("trend")
    ) |>
    left_join(
      sims |>
        select(-"lcl", -"ucl"),
      by = "trend"
    ) |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>   # bereken hoeveel simulaties er nog bij moeten om op 5 uit te komen
    filter(.data$to_do > 0) -> candidate

  while (nrow(candidate) > 0) {
    display_power(sims = sims, power = power, predicted = predicted)
    if (nrow(candidate) == 1) {
      extra <- candidate$trend
    } else {
      extra <- sample(candidate$trend, size = 1, prob = candidate$to_do)
    }
    message("interpolate ", extra)
    data.frame(
      trend = extra,
      fuiktype = fuiktype,
      n_waterlichaam = n_waterlichaam,
      n_sample = n_sample,
      #overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_waterlichaam = sd_waterlichaam,
      p = sim_linear(
        trend = 0.05,
        fuiktype = fuiktype,
        n_waterlichaam = n_waterlichaam,
        n_sample = n_sample,
        initialdensity = initialdensity,
        #overdispersion = overdispersion,
        sd_waterlichaam = sd_waterlichaam
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims

    predicted <- predict_power(sims = sims, step_size = step_size)
    predicted |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        predicted |>
          filter(.data$lcl > power) |>
          slice_min(.data$trend, n = 1),
        sims |>
          filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
          select("trend")
      ) |>
      left_join(
        sims |>
          select(-"lcl", -"ucl"),
        by = "trend"
      ) |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }

  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  predicted <- predict_power(sims = sims, step_size = step_size)

  display_power(sims = sims, power = power, predicted = predicted)
  return(
    c(
      estimate = min(predicted$trend[predicted$fit >= power]),
      lcl = min(predicted$trend[predicted$ucl >= power]),
      ucl = min(predicted$trend[predicted$lcl >= power])
    )
  )
}

### VOOR NEGATIEVE TRENDS
estimate_trend_neg <- function(
    duration = 7,
    n_transect = 5,
    n_sample = 3,
    initialdensity = exp(0.466),
    trend = -0.05,
    overdispersion = 4.87,
    sd_transect = sqrt(0.24520),
    sd_day = sqrt(0.03633),
    alpha = 0.10,
    power = 0.90,
    step_size = 0.001,

    connection = duckdb::dbConnect(
      duckdb::duckdb(),
      #TO: nooit absolute paden naar databases gebruiken, anders werkt het niet op andere computers
      #dbdir = "power_cc_neg.duckdb",
      dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc_neg.duckdb",
      read_only = FALSE # ik sla deze resultaten op in een aparte database
    ),
    #TO
    max_sim = 1000
    ) {

  stopifnot(
    require("assertthat"), require("DBI"), require("dplyr"), require("duckdb"),
    require("purrr"), require("tidyr")
    )

  assert_that(
    is.count(duration),
    is.count(n_sample),
    is.count(n_transect),
    is.number(initialdensity),
    is.number(sd_transect), noNA(sd_transect), sd_transect > 0,
    is.number(sd_day), noNA(sd_day), sd_day > 0,
    is.number(overdispersion), noNA(overdispersion), overdispersion > 0,
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )

  if (!"trend" %in% dbListTables(conn = connection)) {
    data.frame(
      trend = numeric(0),
      duration = integer(0),
      n_transect = integer(0),
      n_sample = integer(0),
      initialdensity = numeric(0),
      overdispersion = numeric(0),
      sd_transect = numeric(0),
      sd_day = numeric(0),
      p = numeric(0)
    ) |>
      dbCreateTable(conn = connection, name = "trend")
  }

  query <- sprintf(
    "SELECT trend,
          SUM(p < %f) AS significant,
          COUNT(trend) AS sims
   FROM trend
   WHERE ABS(sd_transect    - %f) < 1e-4
     AND ABS(sd_day         - %f) < 1e-4
     AND n_transect         = %d
     AND n_sample           = %d
     AND duration           = %d
     AND ABS(initialdensity - %f) < 1e-4
     AND ABS(overdispersion - %f) < 1e-4
   GROUP BY trend
   ORDER BY trend;",
    alpha,
    sd_transect,
    sd_day,
    n_transect,
    n_sample,
    duration,
    initialdensity,
    overdispersion
  )

  dbGetQuery(conn = connection, statement = query) |>
    mutate(simpower = .data$significant / .data$sims) -> sims

  if (nrow(sims) == 0) {
    message("initializing to -0.05")
    data.frame(
      trend = -0.05,
      duration = duration,
      n_transect = n_transect,
      n_sample = n_sample,
      initialdensity = initialdensity,
      overdispersion = overdispersion,
      sd_transect = sd_transect,
      sd_day = sd_day,
      p = sim_linear(
        trend = -0.05,
        duration = duration,
        n_transect = n_transect,
        n_sample = n_sample,
        initialdensity = initialdensity,
        overdispersion = overdispersion,
        sd_transect = sd_transect,
        sd_day = sd_day
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  # AFSCHALEN - zolang de hoogste power < 90 %, (negatieve) kleinste trend * 2
  while (max(sims$simpower) <= power) {
    extra <- 2 * min(sims$trend)          # min(sims$trend) is de kleinste negatieve trend
    message("decreasing to ", extra)
    data.frame(
      trend = extra,
      duration = duration,
      n_transect = n_transect,
      n_sample = n_sample,
      overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_transect = sd_transect,
      sd_day = sd_day,
      p = sim_linear(
        trend = extra,
        duration = duration,
        n_transect = n_transect,
        n_sample = n_sample,
        initialdensity = initialdensity,
        overdispersion = overdispersion,
        sd_transect = sd_transect,
        sd_day = sd_day
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  # OPSCHALEN - als de kleinste, gesimuleerde trend een relatief hoge power heeft (power / 2), wordt deze gehalveerd
  while (
    max(sims$trend) <= -2 * step_size &&
    min(sims$simpower) >= (power / 2)
  ) {
    extra <- round(max(sims$trend) / 2, -floor(log10(step_size)))
    stopifnot(extra < 0)
    message("increasing to ", extra)
    data.frame(
      trend = extra,
      duration = duration,
      n_transect = n_transect,
      n_sample = n_sample,
      overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_transect = sd_transect,
      sd_day = sd_day,
      p = sim_linear(
        trend = extra,
        duration = duration,
        n_transect = n_transect,
        n_sample = n_sample,
        overdispersion = overdispersion,
        initialdensity = initialdensity,
        sd_transect = sd_transect,
        sd_day = sd_day
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  stopifnot(nrow(sims) > 1)
  predicted <- predict_power(sims = sims, step_size = step_size)

  predicted |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      predicted |>
        filter(.data$lcl > power) |>
        slice_min(.data$trend, n = 1),
      sims |>
        filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
        select("trend")
    ) |>
    left_join(
      sims |>
        select(-"lcl", -"ucl"),
      by = "trend"
    ) |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
    filter(.data$to_do > 0) -> candidate

  while (nrow(candidate) > 0) {
    if (nrow(candidate) == 1) {
      extra <- candidate$trend
    } else {
      extra <- sample(candidate$trend, size = 1, prob = candidate$to_do)
    }
    message("interpolate ", extra)
    data.frame(
      trend = extra,
      duration = duration,
      n_transect = n_transect,
      n_sample = n_sample,
      overdispersion = overdispersion,
      initialdensity = initialdensity,
      sd_transect = sd_transect,
      sd_day = sd_day,
      p = sim_linear(
        trend = extra,
        duration = duration,
        n_transect = n_transect,
        n_sample = n_sample,
        overdispersion = overdispersion,
        initialdensity = initialdensity,
        sd_transect = sd_transect,
        sd_day = sd_day
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims

    predicted <- predict_power(sims = sims, step_size = step_size)
    predicted |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        predicted |>
          filter(.data$lcl > power) |>
          slice_min(.data$trend, n = 1),
        sims |>
          filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
          select("trend")
      ) |>
      left_join(
        sims |>
          select(-"lcl", -"ucl"),
        by = "trend"
      ) |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }

  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  predicted <- predict_power(sims = sims, step_size = step_size)

  display_power(sims = sims, power = power, predicted = predicted)
  return(
    c(
      estimate = min(predicted$trend[predicted$fit >= power]),
      lcl = min(predicted$trend[predicted$ucl >= power]),
      ucl = min(predicted$trend[predicted$lcl >= power])
    )
  )
}


#### SIMULATIE SCENARIO'S ####

# BASIS SCENARIO obv PARTRIDGE tellingen in een gemiddelde site
estimate_trend_pos(
  # meetnetparameters
  fuiktype = c("SF", "SFO"),
  n_sample = 5,   # 5 herhalingen per plas
  n_waterlichaam = 3,

  # populatieparameters
  initialdensity = exp(8.864),   # gemiddelde densiteit volgens intercept in model "blankvoorn"

  # ruisparameters
  #overdispersion = 4.87,
  sd_waterlichaam = sqrt(5.375)
)

# estimate_trend_neg(
#   # meetnetparameters
#   duration = 7,   # 2017-2023
#   n_sample = 3,   # 3 tellingen per jaar
#   n_transect = 5,   # 5 x 1km transect per site ~ 1km / 100 ha -> 500 ha
#
#   # populatieparameters
#   initialdensity = exp(0.466),   # gemiddelde densiteit volgens intercept in model "spring_count_data.Rmd"
#
#   # ruisparameters
#   overdispersion = 4.87,
#   sd_transect = sqrt(0.24520),
#   sd_day = sqrt(0.03633)
# )

# simuleer hypothetisch maximum scenario (maximale duur, maximaal gebiedsopp, max aantal tellingen per jaar)
estimate_trend_pos(
  initialdensity = 4.5,
  duration = 10,
  n_sample = 6,
  n_transect = 100)

estimate_trend_neg(
  initialdensity = 4.5,
  duration = 10,
  n_sample = 6,
  n_transect = 100)


#### EVALUEREN VAN TELMETHODE

##### LAGE DENSITEIT
# effect van duration
for (i in c(4,5,6,7,8,9,10)) {
  estimate_trend_pos(
    initialdensity = 0.8,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 0.5,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )
}

# effect van transecten
for (i in c(5,10,20,40,80,100)) {
  estimate_trend_pos(
    initialdensity = 0.8,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )

  estimate_trend_neg(
    initialdensity = 0.5,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )
}

# effect aantal tellingen
for (i in c(2,3,4,5,6)) {
  estimate_trend_pos(
    initialdensity = 0.8,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 0.5,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )
}

##### MEDIUM DENSITEIT
# effect van duration
for (i in c(4,5,6,7,8,9,10)) {
  estimate_trend_pos(
    initialdensity = 1.5,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 1.5,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )
}

# effect van transecten
for (i in c(5,10,20,40,80,100)) {
  estimate_trend_pos(
    initialdensity = 1.5,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )

  estimate_trend_neg(
    initialdensity = 1.5,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )
}

# effect aantal tellingen
for (i in c(2,3,4,5,6)) {
  estimate_trend_pos(
    initialdensity = 1.5,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 1.5,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )
}

##### HOGE DENSITEIT
# effect van duration
for (i in c(4,5,6,7,8,9,10)) {
  estimate_trend_pos(
    initialdensity = 4.5,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 4.5,
    duration = i,
    n_sample = 6,
    n_transect = 100
  )
}

# effect van transecten
for (i in c(5,10,20,40,80,100)) {
  estimate_trend_pos(
    initialdensity = 4.5,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )

  estimate_trend_neg(
    initialdensity = 4.5,
    n_transect = i,
    duration = 10,
    n_sample = 6
  )
}

# effect aantal tellingen
for (i in c(2,3,4,5,6)) {
  estimate_trend_pos(
    initialdensity = 4.5,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )

  estimate_trend_neg(
    initialdensity = 4.5,
    n_sample = i,
    duration = 10,
    n_transect = 100
  )
}

### default density
# effect van duration
for (i in c(4,5,6,7,8,9,10)) {
  estimate_trend_pos(
    duration = i,
    n_sample = 6,
    n_transect = 100
  )

  estimate_trend_neg(
    duration = i,
    n_sample = 6,
    n_transect = 100
  )
}

# effect van transecten
for (i in c(5,10,20,40,80,100)) {
  estimate_trend_pos(
    n_transect = i,
    duration = 10,
    n_sample = 6
  )

  estimate_trend_neg(
    n_transect = i,
    duration = 10,
    n_sample = 6
  )
}

# effect aantal tellingen
for (i in c(2,3,4,5,6)) {
  estimate_trend_pos(
    n_sample = i,
    duration = 10,
    n_transect = 100
  )

  estimate_trend_neg(
    n_sample = i,
    duration = 10,
    n_transect = 100
  )
}

##### RESULTATEN
# Max scenario (default waarden)
max_n_transect <- 100
max_n_sample   <- 6
max_duration   <- 10
overdispersion <- 4.87
sd_transect <- sqrt(0.24520)
sd_day <- sqrt(0.03633)
alpha <- 0.10
power <- 0.90
step_size <- 0.001


con <- dbConnect(duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = TRUE)

# Resultaten uit de database
trend_results <- dbGetQuery(con, "
  SELECT trend, duration, n_transect, n_sample, initialdensity, overdispersion, sd_transect, sd_day, p
  FROM trend
")

summary_results <- trend_results %>%
  group_by(duration, n_transect, n_sample, initialdensity, overdispersion, sd_transect, sd_day, trend) %>%
  filter(initialdensity %in% c(0.8,1.5,4.5)) %>%
  summarise(
    sims = n(),
    significant = sum(p < 0.1),
    simpower = significant / sims,
    .groups = "drop"
  )

estimate_from_sims <- function(dat, power = 0.9, step_size = 0.001) {
  # Zorg dat trend positief/negatief wordt opgeschaald zoals in jouw script
  direction <- ifelse(max(dat$trend) < 0, -1, 1)
  # Binomiaal GLM voor power ~ trend
  m <- glm(cbind(significant, sims-significant) ~ trend, family = binomial, data = dat)
  trend_seq <- direction * seq(step_size, 2 * max(abs(dat$trend[dat$simpower < 1])), by = step_size)
  trend_seq <- trend_seq[trend_seq > -1]
  new_data <- data.frame(trend = trend_seq)
  preds <- predict(m, newdata = new_data, se.fit = TRUE)
  fit <- plogis(preds$fit)
  lcl <- plogis(preds$fit - 1.96 * preds$se.fit)
  ucl <- plogis(preds$fit + 1.96 * preds$se.fit)

  out <- data.frame(trend = trend_seq, fit = fit, lcl = lcl, ucl = ucl)
  estimate = min(out$trend[out$fit >= power])
  lcl_val = min(out$trend[out$ucl >= power])
  ucl_val = min(out$trend[out$lcl >= power])
  c(estimate = estimate, lcl = lcl_val, ucl = ucl_val)
}

# Selecteer relevante scenario's
scenarios <- bind_rows(
  summary_results %>% filter(n_transect == 100, n_sample == 6)   %>% mutate(scenario = "duration",   value = duration),
  summary_results %>% filter(duration == 10, n_sample == 6)      %>% mutate(scenario = "n_transect", value = n_transect),
  summary_results %>% filter(duration == 10, n_transect == 100)  %>% mutate(scenario = "n_sample",   value = n_sample)
)

# Loop per scenario/density/value
detect_trends <- scenarios %>%
  group_by(scenario, value, initialdensity, overdispersion, sd_transect, sd_day) %>%
  group_modify(~{
    res <- estimate_from_sims(.x, power = 0.9, step_size = 0.001)
    tibble(estimate = res[1], lcl = res[2], ucl = res[3])
  }) %>%
  ungroup() %>%
  rename(density = initialdensity)

detect_trends <- detect_trends %>%
  mutate(
    density = factor(density, levels = c(0.8, 1.5, 4.5),
                     labels = c("Low (0.8)", "Medium (1.5)", "High (4.5)")),
    scenario = factor(scenario, levels = c("duration", "n_sample", "n_transect"),
                      labels = c("Number of years", "Counts per year", "Number of transects"))
  )

duration_plot_data   <- detect_trends %>% filter(scenario == "Number of years")
n_sample_plot_data   <- detect_trends %>% filter(scenario == "Counts per year")
n_transect_plot_data <- detect_trends %>% filter(scenario == "Number of transects")

all_y <- c(
  100*detect_trends$lcl,
  100*detect_trends$ucl,
  100*detect_trends$estimate
)
ymin <- floor(min(all_y, na.rm = TRUE))
ymax <- ceiling(max(all_y, na.rm = TRUE))

p_duration <- ggplot(duration_plot_data,
                     aes(x = value, y = 100*estimate, ymin = 100*lcl, ymax = 100*ucl,
                         color = density, fill = density, linetype = density, group = density, shape = density)) +
  geom_ribbon(aes(ymin = 100 * lcl, ymax = 100 * ucl), alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_point(size = 2.7) +
  scale_x_continuous(breaks = 4:10, labels = 4:10) +
  labs(x = "Number of years", y = "Minimum detectable trend (%)") +
  scale_color_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  scale_fill_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(2, 13, by = 2))

p_nsample <- ggplot(n_sample_plot_data,
                    aes(x = value, y = 100*estimate, ymin = 100*lcl, ymax = 100*ucl,
                        color = density, fill = density, linetype = density, group = density, shape = density)) +
  geom_ribbon(aes(ymin = 100 * lcl, ymax = 100 * ucl), alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_point(size = 2.7) +
  scale_x_continuous(breaks = 2:6, labels = 2:6) +
  labs(x = "Counts per year", y = NULL) +
  scale_color_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  scale_fill_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(2, 13, by = 2))

p_ntransect <- ggplot(n_transect_plot_data,
                      aes(x = value, y = 100*estimate, ymin = 100*lcl, ymax = 100*ucl,
                          color = density, fill = density, linetype = density, group = density, shape = density)) +
  geom_ribbon(aes(ymin = 100 * lcl, ymax = 100 * ucl), alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_point(size = 2.7) +
  scale_x_continuous(breaks = c(5, 10, 20, 40, 80, 100), labels = c(5, 10, 20, 40, 80, 100)) +
  labs(x = "Number of transects", y = NULL) +
  scale_color_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  scale_fill_manual(values = c("#4C6A92", "#8B98B5", "#D7D8E9")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(2, 13, by = 2), name = NULL)

final_plot <- p_duration + p_nsample + p_ntransect +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
final_plot

write.csv2(detect_trends, "./output/detect_trends.csv")

# Optimal sampling design (based on simulation outcome)
for (i in c(4.5)) {
  estimate_trend_pos(
    initialdensity = i,
    duration = 7,
    n_sample = 4,
    n_transect = 20
  )
}

# Simulaties obv optimal sampling design
## Lage densiteit
for (i in c(4,6,8,10,15,20)) {
  estimate_trend_pos(
    initialdensity = 0.8,
    duration = i,
    n_sample = 4,
    n_transect = 20
  )
}

## Medium densiteit
for (i in c(4,6,8,10,15,20)) {
  estimate_trend_pos(
    initialdensity = 1.5,
    duration = i,
    n_sample = 4,
    n_transect = 20
  )
}

## Hoge densiteit
for (i in c(4,6,8,10,15,20)) {
  estimate_trend_pos(
    initialdensity = 4.5,
    duration = i,
    n_sample = 4,
    n_transect = 20
  )
}

