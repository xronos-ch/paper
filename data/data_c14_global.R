# data_c14_global.R
# Summarise coverage of 'global' radiocarbon databases (c14bazAAR, IntChron,
# p3k14c, XRONOS)
library(c14)       # https://github.com/joeroe/c14
library(c14bazAAR) # https://github.com/rOpenSci/c14bazAAR
library(countrycode)
library(dplyr, warn.conflicts = FALSE)
library(future)
library(here)
library(purrr)
library(readr)
library(rintchron) # https://github.com/joeroe/rintchron
library(tidyr)
library(xronos)    # https://github.com/xronos-ch/xronos.R

plan(multisession)

# c14bazAAR (excl. p3k14c and XRONOS)
c14bazAAR <- get_c14data("all") |>
  filter(sourcedb != "p3k14c")

# IntChron
egyptdb <- intchron("egyptdb")
nrcf <- intchron("nrcf", tabulate = FALSE)
nrcf <- intchron_tabulate(nrcf[-c(778, 781, 783, 791, 792)]) # bad records
oxa <- intchron("oxa")
sadb <- intchron("sadb")
intchron <- bind_rows(egyptdb, nrcf, oxa, sadb)

# p3k14c
p3k14c <- get_c14data("p3k14c")

# XRONOS
xronos <- chron_data(.everything = TRUE)

# Merge basic information from each database
c14_global <- bind_rows(
  transmute(
    c14bazAAR,  source = "c14bazAAR", 
    lab_id = labnr, c14_age = c14age, c14_error = c14std, country, 
    latitude = lat, longitude = lon
  ),
  transmute(
    intchron, source = "IntChron",
    lab_id = labcode, c14_age = r_date, c14_error = r_date_sigma, 
    country = record_country, latitude = record_latitude, 
    longitude = record_longitude
  ),
  transmute(
    p3k14c, source = "p3k14c",
    lab_id = labnr, c14_age = c14age, c14_error = c14std, country, 
    latitude = lat, longitude = lon
  ),
  transmute(
    xronos, source = "XRONOS",
    lab_id = labnr, c14_age = bp, c14_error = std, country, 
    latitude = as.numeric(lat), longitude = as.numeric(lng)
  )
)

# Clean and deduplicate
c14_global <- c14_global |>
  mutate(
    lab_id = c14_control_lab_id(lab_id, quiet = TRUE, warn_unmatched = FALSE),
    # TODO: something not right here, why is IntChron losing a tonne?
    country = countrycode(country, "country.name", "iso2c", warn = FALSE, nomatch = NULL),
    country = countrycode(country, "iso3c", "iso2c", warn = FALSE, nomatch = NULL),
    country = countrycode(country, "country.name.fr", "iso2c", warn = FALSE, nomatch = NULL),
    country = case_match(
      country,
      "GEQ" ~ "GQ",
      "Channel Isles" ~ "GB",
      "Corsica" ~ "FR",
      "Crete" ~ "GR",
      "England/Wales" ~ "GB",
      "Kosovo" ~ "XK",
      "Mallorca" ~ "ES",
      "Scotland" ~ "GB",
      "SFRY" ~ NA_character_,
      "Sicily" ~ "IT",
      "Yugoslavia" ~ NA_character_,
      "Mauretania" ~ "MR",
      "Meri" ~ NA_character_,
      "Bonaire" ~ "NL",
      "St. Martin" ~ "GB",
      "Alemania" ~ "DE",
      "España" ~ "ES",
      "Francia" ~ "FR",
      "Italia" ~ "IT",
      "Suiza" ~ "CH",
      "Moravia" ~ "CZ",
      "Spaiin" ~ "ES",
      "Danmark" ~ "DK",
      "Switserland" ~ "CH",
      "Denamrk" ~ "DK",
      "Lithuannia" ~ "LT",
      "Tsjech" ~ "CZ",
      "Wales" ~ "GB",
      "Andodrra" ~ "AD",
      "Lituania" ~ "LT",
      "Ris" ~ NA_character_,
      "Unkr" ~ NA_character_,
      "North Sea" ~ NA_character_,
      "Nederland" ~ "NL",
      "IL/PS" ~ NA_character_,
      "CAR" ~ "CF",
      "Israel, Palestine" ~ NA_character_,
      "Israel/West Bank/Gaza Strip" ~ NA_character_,
      .default = country
    ),
    country = na_if(country, "")
  ) |>
  group_by(source) |>
  distinct(lab_id, c14_age, c14_error, .keep_all = TRUE)

# Write basic info
write_tsv(c14_global, here("data", "c14_global.tsv"))

# Sum calibrate by source
c14_global_sum <- c14_global |>
  ungroup() |>
  drop_na(c14_age, c14_error) |>
  filter(c14_age <= 80000, c14_age > 0) |>
  mutate(cal = c14_calibrate(c14_age, c14_error, threshold = 0.001)) |>
  summarise(cal = cal_sum(cal, range = seq(55000, 0, by = -100)), .by = source)

write_tsv(unnest(c14_global_sum, cal), here("data", "c14_global_sum.tsv"))
