# overall Anderson-Darling

library(tidyverse)
library(glmmTMB)
library(kSamples)
read_csv("data/afvissingsdata.csv") |>
  transmute(
    vispunt_id = factor(.data$VispuntID),
    waterlichaam = factor(.data$WaterlichaamNaam),
    datum_id = factor(.data$WaarnemingDatumWID),
    soort = factor(.data$NaamNL),
    lengte = .data$MetingTaxonLengteTotaal
  ) -> empirical

# model voor aantal vissen per waterlichaam, vispunt en datum
empirical |>
  filter(.data$soort == "snoekbaars") |>
  count(
    waterlichaam = factor(.data$waterlichaam),
    vispunt_id = factor(.data$vispunt_id),
    datum_id = factor(.data$datum_id)
  ) |>
  complete(
    nesting(waterlichaam, datum_id, vispunt_id),
    fill = list(n = 0)
  ) -> n_fish
glmmTMB(
  n ~ (1 | waterlichaam) + (1 | datum_id) + (1 | vispunt_id),
  data = n_fish,
  family = poisson(link = "log")
) -> model_count
count_intercept <- fixef(model_count)$cond
count_sd_waterlichaam <- VarCorr(model_count)$cond$waterlichaam
count_sd_datum <- VarCorr(model_count)$cond$datum_id
count_sd_vispunt <- VarCorr(model_count)$cond$vispunt_id

# simulatie maken met aantal vissen per waterlichaam, datum en vispunt
# (met per datum enkel fuiken in eenzelfde waterlichaam),
# en dan voor elk van die vissen een lengte invullen door willekeurige trekking uit dataset met teruglegging.
# per vispunt staan 2 fuiken (met en zonder ottergrid) en vissen worden verdeeld over de 2 fuiken:
# - vissen groter dan drempel_lengte in fuik zonder ottergrid
# - andere worden willekeurig aan fuik toegewezen volgens verdeling 'voorkeur_standaard'
n_waterlichaam <- 10
n_datum <- 10
n_vispunt <- 10
drempel_lengte <- 20
voorkeur_standaard <- 0.7
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
      ~ sample(na.omit(empirical$lengte), size = .x, replace = TRUE)
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
ggplot(sim_data, aes(x = lengte, colour = type)) +
  stat_ecdf()  # = cumulatieve distributiefunctie

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

sim_data |>
  distinct(.data$waterlichaam, .data$datum, .data$vispunt, .data$type) |>
  count(.data$waterlichaam, .data$datum, .data$vispunt) |>
  arrange(n)

sim_data |>
  distinct(.data$waterlichaam, .data$datum, .data$type) |>
  count(.data$waterlichaam, .data$datum) |>
  arrange(n)

# ad.test.combined combineert vergelijkingen tussen afvissingen
# (er is geaggregeerd per afvissing om problemen te vermijden als een van de fuiken een nulwaarde heeft)
# de functie heeft lijsten van vectoren nodig
# deze houdt ook rekening met verschillen in het aantal dieren, ad.test niet.
m2 <- sim_data |>
  group_by(.data$waterlichaam, .data$datum, .data$type) |>
  summarise(ecdf = list(.data$lengte), .groups = "drop_last") |>
  summarise(ecdf = list(.data$ecdf), .groups = "drop") |>
  pull(.data$ecdf) |>
  c(method = "exact") |>
  do.call(what = ad.test.combined)
m2
p2 <- m2$ad.c[, 3]

# en deze vergelijkt enkel het aantal individuen (en houdt geen rekening met lengtes)
sim_data |>
  count(.data$waterlichaam, .data$datum, .data$vispunt, .data$type) |>
  glmmTMB(
    formula = n ~ type + (1 | waterlichaam) + (1 | datum) + (1 | vispunt),
    family = poisson(link = "log")
  ) -> model_count2
summary(model_count2)

p3 <- summary(model_count2)$coefficients$cond[2, 4]

# best power berekenen voor alledrie?
# hoeveel vangdagen zijn nodig zodat de aantallen in 90 % significant zijn?
# en als er een effect is, moet dat in 90 % van de gevallen terecht zijn
# dus zowel type I als type II fout 10 %
