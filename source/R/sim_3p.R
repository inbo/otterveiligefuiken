# functie die op basis van simulatie 3 modellen fit en per model 3 p-waarden teruggeeft
# (eigenlijk 5 omwille van 2 versies voor eerste twee p-waarden)
sim_3p <- function(
    soortnaam = "snoekbaars",
    lengte_soort = empirical_length |> filter(soort == soortnaam),
    n_waterlichaam = 1,
    n_datum = 2,
    n_vispunt = 5,
    drempel_lengte = 20,
    voorkeur_standaard = 0.5, # hebben ook kleinere vissen een voorkeur voor een type fuik? (0.5 = nee)
    count_intercept = 0,
    count_sd_waterlichaam = 7e-13,
    count_sd_datum = 5e-12,
    count_sd_vispunt = 1e-10,
    n_sim = 2
) {
  stopifnot(require(glmmTMB), require(tidyverse))

  replicate(n_sim, {
    # simulatie maken met aantal vissen per waterlichaam, datum en vispunt
    # (met per datum enkel fuiken in eenzelfde waterlichaam),
    # en dan voor elk van die vissen een lengte invullen door willekeurige trekking uit dataset  met teruglegging.
    # per vispunt staan 2 fuiken (met en zonder ottergrid) en vissen worden verdeeld over de 2 fuiken:
    # - vissen groter dan drempel_lengte in fuik zonder ottergrid
    # - andere worden willekeurig aan fuik toegewezen volgens verdeling 'voorkeur_standaard'
    expand.grid(
      waterlichaam = seq_len(n_waterlichaam),
      datum = seq_len(n_datum),
      vispunt = seq_len(n_vispunt)
    ) |>
      mutate(
        datum = (.data$waterlichaam - 1) * n_datum + .data$datum,
        vispunt = (.data$datum - 1) * n_vispunt + .data$vispunt,
        eta = count_intercept +
          rnorm(n_waterlichaam, sd = count_sd_waterlichaam)[.data$waterlichaam] +
          rnorm(n_waterlichaam * n_datum, sd = count_sd_datum)[.data$datum] +
          rnorm(n_waterlichaam * n_datum * n_vispunt, sd = count_sd_vispunt)[
            .data$vispunt
          ],
        aantal = rpois(n(), lambda = exp(.data$eta)),
        lengte = map(
          .data$aantal,
          ~ sample(na.omit(lengte_soort$lengte), size = .x, replace = TRUE)
        )
      ) |>
      unnest(lengte) |>
      mutate(
        type = ifelse(
          .data$lengte > drempel_lengte,
          1,
          rbinom(n(), size = 1, prob = voorkeur_standaard)
        ) |>
          factor(levels = c(1, 0), labels = c("standaard", "aangepast"))
      ) -> sim_data
    # ggplot(sim_data, aes(x = lengte, colour = type)) +
    #   stat_ecdf()  # = cumulatieve distributiefunctie

    # Is cumulatieve distributie verschillend?
    # (info: zie https://fisheries.org/docs/books/55049C/9.pdf)
    # hier is een versie 1 en 2 van, in documentatie bekijken wat het beste is
    # version 1: continuous populations (https://stackoverflow.com/questions/31133870/computing-anderson-darling-test-statistics-for-continuous-distributions-in-r)
    # version 2: discrete populations
    # best version 1 gebruiken, maar misschien best ook version 2 opslaan?  (kost amper extra rekenkracht, en vermijdt dat het opnieuw gedaan moet worden)
    m1 <- ad.test(
      sim_data$lengte[sim_data$type == "aangepast"],
      sim_data$lengte[sim_data$type == "standaard"]
    )
    m1
    p1 <- m1$ad[, 3]

    # sim_data |>
    #   distinct(.data$waterlichaam, .data$datum, .data$vispunt, .data$type) |>
    #   count(.data$waterlichaam, .data$datum, .data$vispunt) |>
    #   arrange(n)
    #
    # sim_data |>
    #   distinct(.data$waterlichaam, .data$datum, .data$type) |>
    #   count(.data$waterlichaam, .data$datum) |>
    #   arrange(n)

    # ad.test.combined combineert vergelijkingen tussen afvissingen
    # (er is geaggregeerd per afvissing om problemen te vermijden als een van de fuiken een nulwaarde heeft)
    # de functie heeft lijsten van vectoren nodig
    # deze houdt ook rekening met verschillen in het aantal dieren, ad.test niet.
    m2 <- try(
      sim_data |>
        group_by(.data$waterlichaam, .data$datum, .data$type) |>
        summarise(ecdf = list(.data$lengte), .groups = "drop_last") |>
        summarise(ecdf = list(.data$ecdf), .groups = "drop") |>
        pull(.data$ecdf) |>
        c(method = "exact") |>
        do.call(what = ad.test.combined)
    )
    if (inherits(m2, "try-error")) {
      p2 <- c("version 1:" = NA_real_, "version 2:" = NA_real_)
    } else {
      p2 <- m2$ad.c[, 3]
    }

    # en deze vergelijkt enkel het aantal individuen (en houdt geen rekening met lengtes)
    sim_data |>
      count(.data$waterlichaam, .data$datum, .data$vispunt, .data$type) |>
      glmmTMB(
        formula = n ~ type + (1 | waterlichaam) + (1 | datum) + (1 | vispunt),
        family = poisson(link = "log")
      ) -> model_count2
    # summary(model_count2)

    p3 <- summary(model_count2)$coefficients$cond[2, 4]

    return(
      c(
        "p1v1" = p1[1], "p1v2" = p1[2],
        "p2v1" = p2[1], "p2v2" = p2[2],
        "p3" = p3
      )
    )
  })
}
