library(tidyverse)
library(glmmTMB)
library(duckdb)
library(showtext)
library(scales)
library(patchwork)
library(kSamples)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

source("source/R/calculate_sd_for_species.R")
source("source/R/sim_3p.R")
source("source/R/display_power.R")
source("source/R/predict_power.R")
source("source/R/estimate_trend_pos.R")

source("source/R/sim_3p_waterlichaam.R")
source("source/R/estimate_trend_pos_waterlichaam.R")



afvissingsdata <- read_csv("data/afvissingsdata.csv") |>
  transmute(
    vispunt_id = factor(.data$VispuntID),
    waterlichaam = factor(.data$WaterlichaamNaam),
    datum_id = factor(.data$WaarnemingDatumWID),
    soort = factor(.data$NaamNL),
    aantal = .data$MetingTaxonAantal,
    lengte = .data$MetingTaxonLengteTotaal
  )

# om het aantal vissen per fuik te modelleren,
# gebruiken we alle data waaruit aantallen af te leiden zijn,
# dus ook de fuiken met groepsmetingen en individuen zonder lengtemetingen.
# Om modelfitproblemen door te weinig vissen per fuik te vermijden,
# laten we waterlichamen met minder dan 100 vissen voor een soort weg
empirical_aantal <- afvissingsdata |>
  group_by(.data$vispunt_id, .data$waterlichaam, .data$datum_id, .data$soort) |>
  summarise(
    aantal = sum(.data$aantal)
  ) |>
  ungroup() |>
  complete(
    nesting(waterlichaam, vispunt_id, datum_id), soort,
    fill = list(aantal = 0)
  )

empirical_aantal |>
  group_by(.data$waterlichaam, .data$soort) |>
  summarise(aantal = sum(.data$aantal)) |>
  arrange(soort, desc(aantal))

veel_vissen_in_plas <- empirical_aantal |>
  group_by(.data$waterlichaam, .data$soort) |>
  summarise(aantal = sum(.data$aantal)) |>
  ungroup() |>
  filter(aantal > 100) |>
  select("soort", "waterlichaam")

empirical_aantal <- empirical_aantal |>
  inner_join(veel_vissen_in_plas, by = c("waterlichaam", "soort"))

# voor de simulatie van de lengtematen proberen we een zo representatief
# mogelijke dataset samen te stellen:
# - enkel waterlichamen waarin de soort veel gezien is
# - enkel fuiken waarin alle individuen van de soort opgemeten zijn (geen groepsmetingen)
# - records met NA-waarden voor lengte laten we weg
niet_volledig_gemeten <- afvissingsdata |>
  filter(aantal > 1) |>
  distinct(.data$vispunt_id, .data$soort)

empirical_length <- afvissingsdata |>
  anti_join(niet_volledig_gemeten, by = c("vispunt_id", "soort")) |>
  filter(!is.na(.data$lengte))

empirical_length |>
  count(.data$soort, .data$waterlichaam) |>
  arrange(soort, desc(n))

# met > 100 per waterlichaam hebben we voor veel soorten minstens 1 plas
veel_in_waterlichaam <- empirical_length |>
  count(.data$soort, .data$waterlichaam) |>
  filter(.data$n > 100)

empirical_length <- empirical_length |>
  inner_join(veel_in_waterlichaam, by = c("soort", "waterlichaam"))


# SNOEKBAARS

calculate_sd_for_species(soortnaam = "snoekbaars", dataset_aantal = empirical_aantal)

estimate_trend_pos(
  soortnaam = "snoekbaars",
  lengte_soort = empirical_length |> filter(soort == "snoekbaars"),
  n_waterlichaam = 1,
  n_datum = 2,
  n_vispunt = 10,
  drempel_lengte = 20,
  voorkeur_standaard = 0.5,
  count_intercept = 1.2,
  count_sd_waterlichaam = 2e-8,
  count_sd_datum = 2.2,
  count_sd_vispunt = 1.2,
  n_sim = 100,
  alpha = 0.10,   # significantieniveau (α)
  power = 0.90,   # gewenste statistische power
  step_size = 1,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = FALSE
  )
)

# effect van waterlichaam uitschakelen:
# - bij 1 waterlichaam specifiek kiezen voor welbepaald waterlichaam
# - random intercept uit formule van waterlichaam weghalen, en het intercept van dat specifieke waterlichaam toevoegen -> gekozen waterlichaam ook als parameter in databank opslaan
estimate_trend_pos_waterlichaam(
  soortnaam = "snoekbaars",
  lengte_soort = empirical_length |> filter(soort == "snoekbaars"),
  n_datum = 2,
  n_vispunt = 5,
  drempel_lengte = 20,
  voorkeur_standaard = 0.5,
  count_intercept = 1.2,
  intercept_waterlichaam = 0,
  count_sd_datum = 2.2,
  count_sd_vispunt = 1.2,
  n_sim = 100,
  alpha = 0.10,   # significantieniveau (α)
  power = 0.90,   # gewenste statistische power
  step_size = 1,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = FALSE
  )
)




# deze wordt geoptimaliseerd om het verschil in lengte te kunnen detecteren

# Nadat we hebben kunnen achterhalen hoe groot de steekproef moet zijn om een verschil in lengte te kunnen detecteren, nagaan welk verschil in aantal vissen (voorkeur voor een bepaalde val) we met deze steekproefgrootte zouden kunnen berekenen


# BLANKVOORN

calculate_sd_for_species(soortnaam = "blankvoorn", dataset_aantal = empirical_aantal)

estimate_trend_pos(
  soortnaam = "blankvoorn",
  lengte_soort = empirical_length |> filter(soort == "blankvoorn"),
  n_waterlichaam = 1,
  n_datum = 2,
  n_vispunt = 10,
  drempel_lengte = 20,
  voorkeur_standaard = 0.5,
  count_intercept = 0,
  count_sd_waterlichaam = 3e-13,
  count_sd_datum = 1e-11,
  count_sd_vispunt = 2e-10,
  n_sim = 100,
  alpha = 0.10,   # significantieniveau (α)
  power = 0.90,   # gewenste statistische power
  step_size = 1,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = FALSE
  )
)




# BAARS

calculate_sd_for_species(soortnaam = "baars", dataset_aantal = empirical_aantal)

estimate_trend_pos(
  soortnaam = "baars",
  lengte_soort = empirical_length |> filter(soort == "baars"),
  n_waterlichaam = 1,
  n_datum = 2,
  n_vispunt = 10,
  drempel_lengte = 20,
  voorkeur_standaard = 0.5,
  count_intercept = 0,
  count_sd_waterlichaam = 2e-11,
  count_sd_datum = 5e-11,
  count_sd_vispunt = 2e-10,
  n_sim = 100,
  alpha = 0.10,   # significantieniveau (α)
  power = 0.90,   # gewenste statistische power
  step_size = 1,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "C:/Users/els_lommelen/Documents/data/otterveiligefuiken/power_cc.duckdb", read_only = FALSE
  )
)



