covid_data<-read.csv(('data/cases_mm.csv'), header = TRUE, sep=";")
covid_si_dist<-read.csv('data/dist.csv',header = TRUE, sep=";")
covid_data$Datum <- as.Date(covid_data$Datum)

date_from = "2020-8-15"
idx<-which(covid_data$Datum==date_from)
covid2020r <- list("incidence"=covid_data$Dennych_PCR_prirastkov[-seq(idx)],"si_dist"=covid_si_dist$y)
covid2020r$date <- as.Date(covid_data$Datum[-seq(idx)])

T <- length(covid_data$Dennych_PCR_prirastkov[-seq(idx)])
ilen = 1
t_start0 <- seq(2,T-ilen+1)
t_end0 <- t_start0+ilen-1
res_covid2020r <- estimate_R(covid2020r$incidence,
                             method="non_parametric_si",
                             config=make_config(list(t_start=t_start0,
                                                     t_end=t_end0,
                                                     si_distr=covid2020r$si_dist)))

plot(res_covid2020r,legend=FALSE)

res_covid2020r$date <- as.Date(covid_data$Datum[-seq(idx)])
res_covid2020r$R$Date <-res_covid2020r$date[-seq(ilen)]
res_covid2020r$R$Cases<-res_covid2020r$I[-seq(ilen)]

write.csv(res_covid2020r$R,file=paste("results/output_R_reported.csv",sep=";"),)
