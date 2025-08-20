# ontwerpen van simulatie op basis van gemeten lengtes per fuik_id voor baars

# afleiden van gemiddelde en stdev uit echte gegevens
library(glmmTMB)
selectie_soorten |>
  filter(!is.na(.data$MetingTaxonLengteTotaal), NaamNL == "baars") |>
  transmute(lengte =.data$MetingTaxonLengteTotaal, fuik_id = .data$WaarnemingWID) -> ds
summary(ds)

# we vertrekken hierbij vanuit een eenvoudig model
m <- glmmTMB(lengte ~ (1 | fuik_id), data = ds, family = Gamma(link = "log"))
summary(m)

# dan gegevens simuleren op basis van deze variabelen,
# we kiezen voor een gamma-distributie met deze vorm:

#rate = shape / exp(intercept)
shape <- summary(m)$sigma ^ -2

n_fuik <- 9
n_vis <- 1e4
gewoon <- 2.27493
sd <- 0.2279
max_vislengte <- 10

# omdat de simulatie geen rekening houdt met de vislengte,
# kappen we bij de aanpassing de vislengte af op de maximumlengte van vis die door het grid kan
expand.grid(
  vis_id = seq_len(n_vis),
  fuik_id = seq_len(n_fuik),
  aanpassing = c(0, 1)
) |>
  mutate(
    rf_fuik = rnorm(n_fuik, mean = 0, sd = sd)[fuik_id],
    eta = gewoon + rf_fuik,
    lengte = rgamma(n = n_vis*n_fuik * 2, shape = shape, rate = shape / exp(eta))
  ) |>
  filter(aanpassing == 0 | lengte < max_vislengte) -> z
ggplot(z, aes(x = lengte, colour = factor(aanpassing))) + geom_density() + facet_wrap(~fuik_id)

# dit is een goed begin, maar die afkapping is niet ideaal
# daarom testen we enkele nieuwe modellen vertrekkend vanuit het model van de vorige simulatie (opgeslagen in z)

m_sim <- glmmTMB(lengte ~ factor(aanpassing) + (1 | fuik_id), data = z, family = Gamma(link = "log"))
m_sim2 <- glmmTMB(lengte ~ factor(aanpassing) + (1 + aanpassing| fuik_id), data = z, family = Gamma(link = "log"))
summary(m_sim)
summary(m_sim2)

# het laatste model lijkt zeer geschikt,
# omdat er een negatieve correlatie is tussen fuik_id en de aanpassing
# (m.a.w. de grootte van de aanpassing hangt af van de ligging het intercept van de fuik (gemiddelde vislengte))

# we maken een nieuwe simulatie startend vanuit dit model en gebruik makend van de random factors

library(mvtnorm)

random_factors <- rmvnorm(n_fuik, sigma = VarCorr(m_sim2)$cond$fuik_id)

expand.grid(
  vis_id = seq_len(n_vis),
  fuik_id = seq_len(n_fuik),
  aanpassing = c(0, 1)
) |>
  mutate(
    rf_fuik_gewoon = random_factors[fuik_id, 1],
    rf_fuik_aangepast = random_factors[fuik_id, 2],
    eta = 2.23 - 0.15245 * aanpassing + rf_fuik_gewoon * (1 - aanpassing) + rf_fuik_aangepast * aanpassing,
    lengte = rgamma(n = n_vis * n_fuik * 2, shape = shape, rate = shape / exp(eta))
  ) |>
  ggplot(aes(x = lengte, colour = factor(aanpassing))) + geom_density() + facet_wrap(~fuik_id)


# Vraag/test: waarom een simulatie van een simulatie?
# Is een gelijkaardig resultaat mogelijk vertrekkend vanuit basisgegevens waarbij die vislengte afgekapt wordt?

ds2 <- ds |>
  mutate(aanpassing = 0) |>
  bind_rows(
    ds |>
      filter(lengte < max_vislengte) |>
      mutate(aanpassing = 1)
  )

m_sim3 <- glmmTMB(lengte ~ factor(aanpassing) + (1 + aanpassing| fuik_id), data = ds2, family = Gamma(link = "log"))
summary(m_sim3)

random_factors3 <- rmvnorm(n_fuik, sigma = VarCorr(m_sim3)$cond$fuik_id)

expand.grid(
  vis_id = seq_len(n_vis),
  fuik_id = seq_len(n_fuik),
  aanpassing = c(0, 1)
) |>
  mutate(
    rf_fuik_gewoon = random_factors3[fuik_id, 1],
    rf_fuik_aangepast = random_factors3[fuik_id, 2],
    eta = 2.27517 - 0.18903 * aanpassing + rf_fuik_gewoon * (1 - aanpassing) + rf_fuik_aangepast * aanpassing,
    lengte = rgamma(n = n_vis * n_fuik * 2, shape = shape, rate = shape / exp(eta))
  ) |>
  ggplot(aes(x = lengte, colour = factor(aanpassing))) + geom_density() + facet_wrap(~fuik_id)
