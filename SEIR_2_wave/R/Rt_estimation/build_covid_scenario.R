covid_data<-read.csv(('data/cases_smooth.csv'), header = TRUE, sep=";")
covid_si_dist<-read.csv('data/dist_shift.csv',header = TRUE, sep=";")
covid_data$Datum <- as.Date(covid_data$Datum)

date_from = "2020-7-1"
idx<-which(covid_data$Datum==date_from)
covid2020r <- list("incidence"=covid_data$Dennych_PCR_prirastkov[-seq(idx)],"si_dist"=covid_si_dist$y)
covid2020r$date <- as.Date(covid_data$Datum[-seq(idx)])

T <- length(covid_data$Dennych_PCR_prirastkov[-seq(idx)])
ilen = 6
t_start0 <- seq(2,T-ilen+1)
t_end0 <- t_start0+ilen-1
res_covid2020r <- estimate_R(covid2020r$incidence,
                             method="non_parametric_si",
                             config=make_config(list(t_start=t_start0,
                                                     t_end=t_end0,
                                                     si_distr=covid2020r$si_dist)))
# res_covid2020r <- estimate_R(covid2020r$incidence,
#                              method="parametric_si",
#                              config=make_config(list(mean_si = 7.0*0.62*0.62,
#                                                     std_si = 1/(0.62*0.62))))

# m = 6.5; m0 = 5; m1 = 7.5;
# s = 0.62; s0 = 0.75; s1 = 0.55;
# config <- make_config(list(mean_si = m*s*s, std_mean_si = 0.2,
#                            min_mean_si = m0*s*s, max_mean_si = m1*s*s,
#                            std_si = 1/(s*s), std_std_si = 0.2,
#                            min_std_si = 1/(s0*s0), max_std_si = 1/(s1*s1)))
# res_covid2020r <- estimate_R(covid2020r$incidence,
#                                method = "uncertain_si",
#                                config = config)

plot(res_covid2020r)

res_covid2020r$date <- as.Date(covid_data$Datum[-seq(idx)])
res_covid2020r$R$Date <-res_covid2020r$date[-seq(ilen)]
res_covid2020r$R$Cases<-res_covid2020r$I[-seq(ilen)]

write.csv(res_covid2020r$R,file=paste("results/output_R_reported.csv",sep=";"),)
