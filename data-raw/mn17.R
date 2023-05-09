## code to prepare `mn17` dataset goes here

library(readr)

ct <- cols(
  subject = col_double(),
  rt = col_double(),
  correct = col_double(),
  jitter = col_double(),
  noise = col_double(),
  n200trialerp_c1r = col_double(),
  n200trialerplat_c1r = col_double(),
  p200trialerp_c1n = col_double(),
  p200trialerplat_c1n = col_double(),
  goodtrials = col_double()
)

mn17 <- readr::read_csv("data-raw/nunez_et_al_2017_eeg_decisionmaking.csv",
                        col_types = ct)

usethis::use_data(mn17, overwrite = TRUE)
