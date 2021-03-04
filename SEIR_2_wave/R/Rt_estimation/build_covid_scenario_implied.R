covid_data_imp<-read.csv(('data/cases_implied.csv'), header = TRUE, sep=";")
covid_si_dist<-read.csv('data/dist.csv',header = TRUE, sep=";")
covid_data_imp$Datum <- as.Date(covid_data_imp$Datum)

date_from = "2020-9-30"
idx<-which(covid_data_imp$Datum==date_from)
covid2020imp <- list("incidence"=covid_data_imp$Dennych_PCR_prirastkov[-seq(idx)],"si_dist"=covid_si_dist$y)
covid2020imp$date <- as.Date(covid_data_imp$Datum[-seq(idx)])

T <- nrow(covid_data_imp$Dennych_PCR_prirastkov[-seq(idx)])
ilen = 0
t_start0 <- seq(2,T-ilen+1)
t_end0 <- t_start0+ilen-1
res_covid2020_imp <- estimate_R(covid2020imp$incidence,
                             method="non_parametric_si",
                             config=make_config(list(t_start=t_start0,
                                                     t_end=t_end0,
                                                     si_distr=covid2020imp$si_dist)))

plot(res_covid2020_imp,legend=FALSE)

res_covid2020_imp$date <- as.Date(covid_data_imp$Datum[-seq(idx)])
res_covid2020_imp$R$Date <-res_covid2020_imp$date[-seq(ilen)]
res_covid2020_imp$R$Cases<-res_covid2020_imp$I[-seq(ilen)]


write.csv(res_covid2020_imp$R,file=paste("results/output_R_implied.csv",sep=""),)
