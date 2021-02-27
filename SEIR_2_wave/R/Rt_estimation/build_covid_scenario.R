covid_data<-read.csv(('data/cases_mm.csv'), header = TRUE, sep=";")
covid_si_dist<-read.csv('data/dist.csv',header = TRUE, sep=";")
covid_data$Datum <- as.Date(covid_data$Datum)

date_from = "2020-10-31"
idx<-which(covid_data$Datum==date_from)
covid2020r <- list("incidence"=covid_data$Dennych.PCR.prirastkov[-seq(idx)],"si_dist"=covid_si_dist$y)
covid2020r$date <- as.Date(covid_data$Datum[-seq(idx)])
#covid2020r <- list("incidence"=covid_data$Dennych.PCR.prirastkov,"si_dist"=covid_si_dist$y)
#covid2020r$date <- as.Date(covid_data$Datum)

res_covid2020r <- estimate_R(covid2020r$incidence,
                                method="non_parametric_si",
                                config=make_config(list(si_distr=covid2020r$si_dist)))

plot(res_covid2020r,legend=FALSE)

####################################

# si_mean = 6.5
# si_std = 0.62
# res_covid2020rp <- estimate_R(covid2020r$incidence,
#                              method="parametric_si",
#                              config=make_config(list(
#                                mean_si=si_mean*si_std*si_std,
#                                std_si=1/(si_std*si_std))))
# plot(res_covid2020rp,legend=FALSE)

