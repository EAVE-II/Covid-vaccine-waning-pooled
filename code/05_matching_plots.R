##### Final plots
# Colours 
eave_green <- rgb(54, 176, 136, maxColorValue = 255)
eave_blue <- rgb(71,93,167, maxColorValue = 255)
eave_blue2 <- rgb(0,192,209, maxColorValue = 255)
eave_gold <- rgb(255,192,0, maxColorValue = 255)
eave_orange <- rgb(244,143,32, maxColorValue = 255)

###### Cumulative risk plots using ggplot ########

library(survival)
library(survminer)


#### Overall
par(mfrow=c(1,2))
z <- survfit(Surv(time_to_hosp, event) ~ vacc, data=df_cc)

ggsurvplot(z, data=df_cc, fun="event", conf.int = T,
           ggtheme = theme_light(),
           xlab= "Follow-up (days)",
           palette = c(eave_blue, eave_green),
           legend.labs = c("Unvaccinated", "Vaccinated")) 

title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc$vacc), lty=1, col=c(1,2), cex=0.8)

# From 14 days
z <- survfit(Surv(time_to_event14, event) ~ vacc, data=df_cc)
plot(z, fun="event", col=c(1,2), conf.int=T)
title(main=z_title, 
      xlab="14 days from vaccination",ylab="cumulative risk")
legend("topleft",legend=levels(df_cc$vacc), lty=1, col=c(1,2), cex=0.8)

