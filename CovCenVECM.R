library(gtrendsR)
library(tidyverse)
library(aTSA)
library(tsDyn) ## imports vars and urca

p04_data <- read_csv("GTrends.csv")
hosp01 <- read_csv("Hospital.csv")
HBot <- read_csv("HBot.csv")

master_data <- full_join(p04_data, HBot)
master_data <- full_join(master_data, hosp01)
master_data <- master_data %>% replace_na(list(Census = 0))
master_data$Time <- 1:nrow(master_data)

## Spearman verification
step01_I_list <- list()
for(j in 2:3){
info <- matrix(NA, nrow = 15, ncol = 2)
for(i in 0:14){
info[i + 1, 1] <- i
info[i + 1, 2] <- cor(master_data$Census, lag(master_data[[j]], i), use = "pairwise.complete.obs", method = "spearman")
result <- as_tibble(info) %>% filter(abs(V2) == max(abs(V2)))
}
step01_I_list[[j - 1]] <- result
big_result01_I <- bind_rows(step01_I_list)
}
big_result01_I

## Linearization of Xs
master_data$HB <- log(master_data$HB) 
master_data$Testing <- sqrt(master_data$Testing)

## Iterate through Testing and HB
blast02_I <- matrix(NA, nrow = 15, ncol = 3)
blast02_I[,1] <- 0:-14
for(i in 0:14){blast02_I[i + 1, 2] <- cor(master_data$Census, lag(master_data[[2]], i), use = "pairwise.complete.obs")}
for(i in 0:14){blast02_I[i + 1, 3] <- cor(master_data$Census, lag(master_data[[3]], i), use = "pairwise.complete.obs")}
blast02_I <- as_tibble(blast02_I)

## Figure 2 in paper
p21 <- ggplot(blast02_I, aes(x = V1)) 
p21 <- p21 + geom_line(aes(y = V2, colour = "Num1"))
p21 <- p21 + geom_point(aes(y = V2, colour = "Num1"), size = 1.5)
p21 <- p21 + geom_line(aes(y = V3, colour = "Num2"))
p21 <- p21 + geom_point(aes(y = V3, colour = "Num2"), size = 1.5)
p21 <- p21 + ylim(0.50, 1.00)
p21 <- p21 + labs(y = "Correlation", 
                  title = "Longitudinal Cross-Correlation Profiles for Census",
                  subtitle = "2020 Time Period: Google Trends, Feb 21 to Aug 1; Health Bot, Apr 16 to Aug 1")
p21 <- p21 + scale_x_continuous("Time Lag (Days) For Leading Internet Variable", 
                                breaks = 0:-14)
cols <- c("Num1" = "black", "Num2" = "red")
p21 <- p21 + scale_colour_manual(name = "Internet Variable", values = cols, labels = c("Testing", "Health Bot"))
p21 <- p21 + theme(legend.position = c(0.15, .20), legend.key.width = unit(3, "cm"),
                   panel.grid.minor.x = element_blank())
p21

p22_data <- master_data %>% mutate(Scaled_C = 100*(Census/max(Census)), 
                                   Scaled_HB = 100*(HB - min(HB, na.rm = TRUE))/(max(HB, na.rm = TRUE) - min(HB, na.rm = TRUE)),
                                   Scaled_Test = 100*(Testing - min(Testing, na.rm = TRUE))/(max(Testing, na.rm = TRUE) - min(Testing, na.rm = TRUE)))
p22_data$date <- as.Date(p22_data$date, "%m/%d/%y") 

## Figure 3 in paper
p23b <- ggplot(p22_data, aes(x = date)) 
p23b <- p23b + geom_line(aes(y = Scaled_C, colour = "Num1"))
p23b <- p23b + geom_line(aes(y = Scaled_Test, colour = "Num2"))
p23b <- p23b + geom_line(aes(y = Scaled_HB, colour = "Num3"))
p23b <- p23b + labs(y = "Standardized Scale", x = "Date",
                  title = "Scaled Time Series for Census, Testing, and Health Bot",
                  subtitle = "2020 Time Period: Google Trends, Feb 21 to Aug 1; Health Bot, Apr 16 to Aug 1")
p23b <- p23b + scale_x_date(date_breaks = "2 weeks",
                          date_labels = "%b %d",  
                          limits = as.Date(c("2020-2-21", "2020-8-01")))
cols <- c("Num1" = "darkblue", "Num2" = "red", "Num3" = "darkgreen")
p23b <- p23b + scale_colour_manual(name = "Scaled Time Series", values = cols, labels = c("Census", "Testing", "Health Bot"))
p23b <- p23b + theme(legend.position = c(0.15, .80), legend.key.width = unit(2, "cm"),
                   panel.grid.minor.x = element_blank())
p23b <- p23b + guides(linetype = guide_legend(override.aes = list(size = 0.75)))
p23b

## VECM

dat <- cbind(master_data$Testing, master_data$HB, master_data$Census)[56:163,]
colnames(dat) <- c("Testing", "HB", "Census")
VARselect(dat, lag.max = 14, type = "const") ## Inconclusive; accept default K = 2 for ca.jo()
H1 <- ca.jo(dat, type = "trace", ecdet = "const", spec = "transitory")
summary(H1)

stb <- vec2var(H1, r = 2)
Mod(eigen(rbind(cbind(stb$A$A1, stb$A$A2), cbind(diag(3),0,0,0)))$values) ## Looks good

lttest(H1, r = 2) ## Do not include linear trend

## Look at Pi long-run elements to see if they are stationary
## Testing (Health Bot not shown for brevity)
beta_01 <- 0.319786070 + dat[,1] - 0.77046287*dat[,2] - 0.01559538*dat[,3] 
plot(beta_01, type = "l")
trunc(12*(length(beta_01)/100)^(1/4)) ## 12 chosen as upper bound
tryit01 <- ur.df(beta_01, type = "trend", lags = 12, selectlags = "BIC")  
summary(tryit01) ## 2 lags using Holm correction method
tryit02 <- ur.df(beta_01, type = "trend", lags = 2)
summary(tryit02) ## tau_3 and Phi_3 marginally significant
tryit03 <- ur.df(beta_01, type = "drift", lags = 2)
summary(tryit03) ## tau_2 significant 

## Put constant inside the co-integration component
vecm2 <- VECM(dat, lag = 1, estim = "ML", r = 2) ## Include = "const" is default
vecm <- VECM(dat, lag = 1, estim = "ML", r = 2, LRinclude = "const")
AIC(vecm); AIC(vecm2) 

cajorls(H1, r = 2)$beta ## Restricted VECM; Table 4
summary(cajorls(H1, r = 2)$rlm) ## Table 3 

## Check assumptions
serial.test(stb)
normality.test(stb, multivariate.only = FALSE)
acf(resid(stb))

## MAPE (Aug 2 - 15)
Obs <- c(211, 232, 249, 232, 227, 216, 229, 221, 215, 215, 207, 200, 192, 198)
first_try <- predict(stb, n.ahead = 14, ci = 0.95)
Fcst <- first_try$fcst$Census[,1] 
100*mean(abs(Obs - Fcst)/Obs) ## 6.4% 

## Figure 4 in paper 
upper <- c(NA, NA, fitted(vecm, level = "original")[,3] + 2*sqrt(summary(vecm)$sigma[3,3]))
lower <- c(NA, NA, fitted(vecm, level = "original")[,3] - 2*sqrt(summary(vecm)$sigma[3,3]))
my_dates <- seq(as.Date("2020/04/16"), as.Date("2020/08/01"), "days")
poly01 <- tibble(Date = c(my_dates, rev(my_dates)), Bounds = c(upper, rev(lower)))
my_dates02 <- seq(as.Date("2020/08/02"), as.Date("2020/08/15"), "days")
poly02 <- tibble(Date = c(my_dates02, rev(my_dates02)), Bounds = c(first_try$fcst$Census[,3], rev(first_try$fcst$Census[,2])))
p24_graph <- tibble(Testing = c(dat[, "Testing"], rep(NA, 14)), 
                    HB = c(dat[, "HB"], rep(NA, 14)), 
                    Census = c(dat[, "Census"], Obs),  
                    Date = c(my_dates, my_dates02),
                    Fit = c(NA, NA, fitted(vecm, level = "original")[,3], Fcst))

p24 <- ggplot(p24_graph, aes(x = Date))
p24 <- p24 + geom_polygon(data = poly01, aes(x = Date, y = Bounds), fill = "#00BFC4", alpha = 0.4)
p24 <- p24 + geom_polygon(data = poly02, aes(x = Date, y = Bounds), fill = "#F8766D", alpha = 0.4)
p24 <- p24 + geom_line(aes(y = Fit), colour = "red")
p24 <- p24 + geom_point(aes(y = Census), size = 0.5, color = "black")
p24 <- p24 + ylim(0, NA) + xlim(as.Date("2020-04-16"), as.Date("2020-08-15"))
p24 <- p24 + labs(y = "Hospital Covid-19+ Census", 
                  title = "Hospital Covid-19+ Census Internet Model, 14-Day Forecast
Greater Charlotte, NC Region as of August 1, 2020",
                  subtitle = "Fitted Values (red), Observations (black),
95% Prediction Interval for In-Sample One-Step-Ahead Forecasts (blue), 
95% Forecast Interval for Out-of-Sample Forecasts (pink)", 
                  caption = "MAPE (Aug 2 - Aug 15): 6.4%")
p24 <- p24 + theme(plot.caption = element_text(hjust = 0), legend.position = "none")
p24