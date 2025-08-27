# hetzelfde als estimate_trend_pos, maar dan voor simulaties voor 1 waterlichaam

# Op basis van opgegeven parameters wordt data gesimuleerd, 3 modellen worden gefit op de data en de p-waardes worden teruggegeven (zie functie sim_3p). Resultaat van de functie is een lijst met n p-waarden (n is het aantal replicaties). Deze lijst met p-waarden kan vervolgens gebruikt worden om een power te berekenen als het aandeel significante p-waarden op het totaal aantal simulaties. vb. data werd gesimuleerd met een gekende trend van -5%. We simuleren de data 100 keer en zien dat het model slechts 30 keer een significant jaar-effect detecteert. De power van de proefopzet om de trend van -5% te meten is dus gelijk aan 30%

# Wat voor trend kunnen we wel zien?
## Gegeven een bepaalde meetopzet (aantal jaren, transecten, meetfout, enz.), wat is de minimale trend die met een opgegeven power (hier gebruiken we een power van 90%) kan worden gedetecteerd?

### VOOR POSITIEVE TRENDS
estimate_trend_pos_waterlichaam <- function(
    soortnaam = "snoekbaars",
    lengte_soort = empirical_length |> filter(soort == soortnaam),
    intercept_waterlichaam = 1,
    n_datum = 2,
    n_vispunt = 5,
    drempel_lengte = 20,
    voorkeur_standaard = 0.5,
    count_intercept = 0,
    count_sd_datum = 5e-12,
    count_sd_vispunt = 1e-10,
    n_sim = 2,

    # powerparameters
    alpha = 0.10,   # significantieniveau (α)
    power = 0.90,   # gewenste statistische power
    step_size = 1,   # kleinste afstand tussen stappen in design (eerst ruw, dan fijner")

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
    is.character(soortnaam),
    is.count(n_datum),
    is.count(n_vispunt),

    # continue variabelen
    is.number(drempel_lengte), noNA(drempel_lengte), drempel_lengte > 0,
    is.number(voorkeur_standaard), noNA(voorkeur_standaard),
    voorkeur_standaard >= 0, voorkeur_standaard <= 1,
    is.number(count_intercept), noNA(count_intercept), #count_intercept >= 0,
    is.number(intercept_waterlichaam), noNA(intercept_waterlichaam),
    #intercept_waterlichaam > 0,
    is.number(count_sd_datum), noNA(count_sd_datum), count_sd_datum >= 0,
    is.number(count_sd_vispunt), noNA(count_sd_vispunt), count_sd_vispunt >= 0,

    # alfa en power begrensd tussen 0 en 1
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )

  # tabel "trend" in DuckDB initialiseren indien deze nog niet bestaat
  if (!"trend" %in% dbListTables(conn = connection)) {
    #   dbRemoveTable(conn = connection, name = "trend")
    data.frame(
      methode = character(0),
      soortnaam = character(0),
      n_waterlichaam = integer(0),
      n_datum = integer(0),
      n_vispunt = integer(0),
      drempel_lengte = numeric(0),
      voorkeur_standaard = numeric(0),
      intercept = numeric(0),
      intercept_waterlichaam = numeric(0),
      sd_datum = numeric(0),
      sd_vispunt = numeric(0),
      p1v1 = numeric(0),
      p1v2 = numeric(0),
      p2v1 = numeric(0),
      p2v2 = numeric(0),
      p3 = numeric(0)
    ) |>
      dbCreateTable(conn = connection, name = "trend")
  }

  # qery: reeds bestaande simulaties op basis van gespecifieerde simulatieparameters ophalen
  # per trend wordt het totale aantal simulaties (sims) en het aantal significante simulaties p < alpha gegeven
  # omdat sd_transect, sd_day, initialdensity, enz. komma-getallen zijn (float), gebruik je een tolerantie — anders afrondingsfouten. Zoek "ongeveer gelijk aan", binnen een marge van 0.0001.
  query <- sprintf(
    "SELECT n_datum,
          SUM(p1v1 < %f) AS significant_1_lengte,
          SUM(p2v1 < %f) AS significant_2_lengte_aantal,
          SUM(p3 < %f) AS significant_3_aantal,
          COUNT(n_datum) AS sims
   FROM trend
   WHERE soortnaam = '%s'
     AND methode = '1 waterlichaam'
     AND n_waterlichaam         = 1
     AND n_vispunt           = %d
     AND ABS(drempel_lengte - %f) < 1e-4
     AND ABS(voorkeur_standaard - %f) < 1e-4
     AND ABS(intercept - %f) < 1e-4
     AND ABS(intercept_waterlichaam - %f) < 1e-4
     AND ABS(sd_datum - %f) < 1e-4
     AND ABS(sd_vispunt - %f) < 1e-4
   GROUP BY n_datum
   ORDER BY n_datum;",
    alpha, alpha, alpha,
    soortnaam,
    n_vispunt,
    drempel_lengte,
    voorkeur_standaard,
    count_intercept,
    intercept_waterlichaam,
    count_sd_datum,
    count_sd_vispunt
  )

  dbGetQuery(conn = connection, statement = query) |>
    pivot_longer(
      cols = starts_with("significant"),
      names_to = "modeltype",
      names_pattern = "significant_(\\d_\\w*)",
      values_to = "significant"
    ) |>
    filter(!is.na(significant)) |>
    mutate(
      simpower = .data$significant / .data$sims
    ) -> sims

  # start simulatie indien query geen resultaat gaf
  if (nrow(sims) == 0) {
    message("initializing simulation")
    p_sim <- sim_3p_waterlichaam(
      soortnaam = soortnaam,
      lengte_soort = lengte_soort,
      n_datum = n_datum,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      count_intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      count_sd_datum = count_sd_datum,
      count_sd_vispunt = count_sd_vispunt,
      n_sim = n_sim
    )
    data.frame(
      methode = "1 waterlichaam",
      soortnaam = soortnaam,
      n_waterlichaam = 1,
      n_datum = n_datum,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      sd_datum = count_sd_datum,
      sd_vispunt = count_sd_vispunt,
      p1v1 = p_sim["p1v1.version 1:", ],
      p1v2 = p_sim["p1v2.version 2:", ],
      p2v1 = p_sim["p2v1.version 1:", ],
      p2v2 = p_sim["p2v2.version 2:", ],
      p3 = p_sim["p3", ]
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      pivot_longer(
        cols = starts_with("significant"),
        names_to = "modeltype",
        names_pattern = "significant_(\\d_\\w*)",
        values_to = "significant"
      ) |>
      filter(!is.na(significant)) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }

  # OPSCHALEN: zolang de hoogste power < 90%, verdubbel aantal datums
  while (max(sims$simpower, na.rm = TRUE) <= power) {
    extra <- 2 * max(sims$n_datum)
    message("increasing to ", extra)
    p_sim <- sim_3p_waterlichaam(
      soortnaam = soortnaam,
      lengte_soort = lengte_soort,
      n_datum = extra,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      count_intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      count_sd_datum = count_sd_datum,
      count_sd_vispunt = count_sd_vispunt,
      n_sim = n_sim
    )
    data.frame(
      methode = "1 waterlichaam",
      soortnaam = soortnaam,
      n_waterlichaam = 1,
      n_datum = extra,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      sd_datum = count_sd_datum,
      sd_vispunt = count_sd_vispunt,
      p1v1 = p_sim["p1v1.version 1:", ],
      p1v2 = p_sim["p1v2.version 2:", ],
      p2v1 = p_sim["p2v1.version 1:", ],
      p2v2 = p_sim["p2v2.version 2:", ],
      p3 = p_sim["p3", ]
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      pivot_longer(
        cols = starts_with("significant"),
        names_to = "modeltype",
        names_pattern = "significant_(\\d_\\w*)",
        values_to = "significant"
      ) |>
      filter(!is.na(significant)) |>
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
    pivot_longer(
      cols = starts_with("significant"),
      names_to = "modeltype",
      names_pattern = "significant_(\\d_\\w*)",
      values_to = "significant"
    ) |>
    filter(!is.na(significant)) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  stopifnot(nrow(sims) > 1)
  predicted_lengte <- predict_power(
    sims = sims |>
      filter(modeltype == "1_lengte"),
    step_size = step_size
  )

  predicted_lengte_aantal <- predict_power(
    sims = sims |>
      filter(modeltype == "2_lengte_aantal"),
    step_size = step_size
  )

  predicted_aantal <- predict_power(
    sims = sims |>
      filter(modeltype == "3_aantal"),
    step_size = step_size
  )

  # Voor welke trendwaarden overlapt het CI nog met 90%; hier extra simulaties uitvoeren om de schatting te verfijnen
  #we doen dit voorlopig enkel om een verschil in de lengte te kunnen waarnemen
  predicted_lengte |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      predicted_lengte |>
        filter(.data$lcl > power) |>   # selecteer trendwaarden waarvan het 95%-interval de drempelwaarde (power) overlapt
        slice_min(.data$n_datum, n = 1),   # voeg de eerstvolgende trend toe waarbij het hele interval boven de powerdrempel ligt
      sims |>
        filter(modeltype == "1_lengte") |>
        filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
        select("n_datum")
    ) |>
    left_join(
      sims |>
        filter(modeltype == "1_lengte") |>
        select(-"lcl", -"ucl"),
      by = "n_datum"
    ) |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>   # bereken hoeveel simulaties er nog bij moeten om op 5 uit te komen
    filter(.data$to_do > 0) -> candidate

  while (nrow(candidate) > 0) {
    display_power(sims = sims, power = power, predicted = list(lengte = predicted_lengte), soortnaam = soortnaam)
    if (nrow(candidate) == 1) {
      extra <- as.integer(round(candidate$n_datum))
    } else {
      extra <- as.integer(round(sample(candidate$n_datum, size = 1, prob = candidate$to_do)))
    }
    message("interpolate ", extra)
    p_sim <- sim_3p_waterlichaam(
      soortnaam = soortnaam,
      lengte_soort = lengte_soort,
      n_datum = extra,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      count_intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      count_sd_datum = count_sd_datum,
      count_sd_vispunt = count_sd_vispunt,
      n_sim = n_sim
    )
    data.frame(
      methode = "1 waterlichaam",
      soortnaam = soortnaam,
      n_waterlichaam = 1,
      n_datum = extra,
      n_vispunt = n_vispunt,
      drempel_lengte = drempel_lengte,
      voorkeur_standaard = voorkeur_standaard,
      intercept = count_intercept,
      intercept_waterlichaam = intercept_waterlichaam,
      sd_datum = count_sd_datum,
      sd_vispunt = count_sd_vispunt,
      p1v1 = p_sim["p1v1.version 1:", ],
      p1v2 = p_sim["p1v2.version 2:", ],
      p2v1 = p_sim["p2v1.version 1:", ],
      p2v2 = p_sim["p2v2.version 2:", ],
      p3 = p_sim["p3", ]
    ) |>
      dbAppendTable(conn = connection, name = "trend")

    dbGetQuery(conn = connection, statement = query) |>
      pivot_longer(
        cols = starts_with("significant"),
        names_to = "modeltype",
        names_pattern = "significant_(\\d_\\w*)",
        values_to = "significant"
      ) |>
      filter(!is.na(significant)) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims

    predicted_lengte <- predict_power(sims = sims |> filter(modeltype == "1_lengte"), step_size = step_size)
    predicted_lengte |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        predicted_lengte |>
          filter(.data$lcl > power) |>
          slice_min(.data$n_datum, n = 1),
        sims |>
          filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
          select("n_datum")
      ) |>
      left_join(
        sims |>
          filter(modeltype == "1_lengte") |>
          select(-"lcl", -"ucl"),
        by = "n_datum"
      ) |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }

  dbGetQuery(conn = connection, statement = query) |>
    pivot_longer(
      cols = starts_with("significant"),
      names_to = "modeltype",
      names_pattern = "significant_(\\d_\\w*)",
      values_to = "significant"
    ) |>
    filter(!is.na(significant)) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims

  predicted_lengte <- predict_power(
    sims = sims |>
      filter(modeltype == "1_lengte"),
    step_size = step_size
  )

  predicted_lengte_aantal <- predict_power(
    sims = sims |>
      filter(modeltype == "2_lengte_aantal"),
    step_size = step_size
  )

  predicted_aantal <- predict_power(
    sims = sims |>
      filter(modeltype == "3_aantal"),
    step_size = step_size
  )

  display_power(
    sims = sims,
    power = power,
    predicted = list(
      lengte = predicted_lengte,
      lengte_aantal = predicted_lengte_aantal,
      aantal = predicted_aantal
    ),
    soortnaam = soortnaam
  )
  return(
    c(
      estimate_lengte = min(predicted_lengte$n_datum[predicted_lengte$fit >= power]),
      lcl_lengte = min(predicted_lengte$n_datum[predicted_lengte$ucl >= power]),
      ucl_lengte = min(predicted_lengte$n_datum[predicted_lengte$lcl >= power]),
      estimate_lengte_aantal =
        min(predicted_lengte_aantal$n_datum[predicted_lengte_aantal$fit >= power]),
      lcl_lengte_aantal =
        min(predicted_lengte_aantal$n_datum[predicted_lengte_aantal$ucl >= power]),
      ucl_lengte_aantal =
        min(predicted_lengte_aantal$n_datum[predicted_lengte_aantal$lcl >= power]),
      estimate_aantal = min(predicted_aantal$n_datum[predicted_aantal$fit >= power]),
      lcl_aantal = min(predicted_aantal$n_datum[predicted_aantal$ucl >= power]),
      ucl_aantal = min(predicted_aantal$n_datum[predicted_aantal$lcl >= power])
    )
  )
}


