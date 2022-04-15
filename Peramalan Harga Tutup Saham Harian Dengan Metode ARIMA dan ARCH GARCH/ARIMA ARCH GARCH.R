#Tugas Akhir UTS ADW - Saham ICBP
library(forecast)
library(FitAR)
library(lmtest)
library(tseries)
library(psych)
library(ggplot2)
library(readr)
library(tsoutliers)
library(readr)
library(timeDate) 
library(MASS) 
library(imputeTS) 
library(rugarch)

data <- read.csv(file.choose(), header=F, sep=';')
# data : saham icbp (missing value)

head(data)
dim(data)
data <- data[,2]

#Missing value 
require(imputeTS)
ggplot_na_distribution(data)
data1 <- na_ma(data,k=3,weighting="simple")
round(data1,0)

#Data di time series
series <- ts(data1, start=2020, frequency=365) ; head(series)

#Plot data
ts.plot(series,ylab="Value",main="Plot Indeks Harga Tutup Saham ICBP (sebelum differencing)",col="turquoise",lwd=2)
legend("bottomright",c("Value"),cex=0.8,lty=5,text.font=2,col=c("turquoise"))

#Identifikasi asumsi homogenitas varians (stasioner varians)
lambda <- BoxCox.lambda(series) ;lambda

#Identifikasi asumsi stasioner rata-rata
adf.test(series)

#Plot ACF dan PACF
ggAcf(series, lag=25)+ ggtitle("Plot ACF sebelum differencing")
ggPacf(series, lag=25)+ ggtitle("Plot PACF sebelum differencing")

#Differencing karena tidak stasioner rata-rata
library(urca)
summary(ur.kpss(series))
ndiffs(series)
df.series <- diff(series,1)

#Cek kembali stasioner rata-rata
adf.test(df.series)

#Plot setelah differencing
ts.plot(df.series,ylab="Value",main="Plot Indeks Harga Tutup Saham ICBP (setelah differencing)",col="turquoise",lwd=2)
legend("bottomright",c("Value"),cex=0.8,lty=5,text.font=2,col=c("turquoise"))

#Identifikasi model (Plot ACF dan PACF setelah differencing)
ggAcf(df.series, lag=50)+ ggtitle("Plot ACF setelah differencing")
ggPacf(df.series, lag=50) + ggtitle("Plot PACF setelah differencing")

#Menentukan model dengan lihat cut off & menulis tiga pilihan model
#Q(0,1,2,3) P(0,1,2) d(1)
model1 <- arima(series,order=c(0,1,0)) ;model1 #IMA(q) integrated
model2 <- arima(series,order=c(0,1,1)) ;model2 #ARI(p) integrated
model3 <- arima(series,order=c(0,1,2)) ;model3
model4 <- arima(series,order=c(0,1,3)) ;model4
model5 <- arima(series,order=c(1,1,0)) ;model5
model6 <- arima(series,order=c(1,1,1)) ;model6
model7 <- arima(series,order=c(1,1,2)) ;model7
model8 <- arima(series,order=c(1,1,3)) ;model8
model9 <- arima(series,order=c(2,1,0)) ;model9
model10 <- arima(series,order=c(2,1,1)) ;model10
model11 <- arima(series,order=c(2,1,2)) ;model11
model12 <- arima(series,order=c(2,1,3)) ;model12

#Uji signifikansi parameter
coeftest(model1) #0 parameter
coeftest(model2) #sig
coeftest(model3)
coeftest(model4)
coeftest(model5) #sig
coeftest(model6) 
coeftest(model7) #sig
coeftest(model8)
coeftest(model9) #sig
coeftest(model10)
coeftest(model11)
coeftest(model12)

#Asumsi model 2
r.2 <- residuals(model2)
#1. Normalitas residual
n.2 = length(r.2)
mean.2 = mean(r.2)
sd.2 = sd(r.2)
res.2 = rnorm(n.2,mean.2,sd.2)
normalitas.2 <- ks.test(r.2,res.2) ;normalitas.2
#2. Residual white noise (tidak ada autokorelasi)
WN.2 <- Box.test(r.2, type=c("Ljung-Box")) ;WN.2
#3. Homokedastisitas
h.2 = r.2^2
heteros.2 <- Box.test(h.2, type=c("Ljung-Box")) ;heteros.2

#Asumsi model 5
r.5 <- residuals(model5)
#1. Normalitas residual
n.5 = length(r.5)
mean.5 = mean(r.5)
sd.5 = sd(r.5)
res.5 = rnorm(n.5,mean.5,sd.5)
normalitas.5 <- ks.test(r.5,res.5) ;normalitas.5
#2. Residual white noise (tidak ada autokorelasi)
WN.5 <- Box.test(r.5, type=c("Ljung-Box")) ;WN.5
#3. Homokedastisitas
h.5 = r.5^2
heteros.5 <- Box.test(h.5, type=c("Ljung-Box")) ;heteros.5

#Asumsi model 7
r.7 <- residuals(model7)
#1. Normalitas residual
n.7 = length(r.7)
mean.7 = mean(r.7)
sd.7 = sd(r.7)
res.7 = rnorm(n.7,mean.7,sd.7)
normalitas.7 <- ks.test(r.7,res.7) ;normalitas.7
#2. Residual white noise (tidak ada autokorelasi)
WN.7 <- Box.test(r.7, type=c("Ljung-Box")) ;WN.7
#3. Homokedastisitas
h.7 = r.7^2
heteros.7 <- Box.test(h.7, type=c("Ljung-Box")) ;heteros.7

#Asumsi model 9
r.9 <- residuals(model9)
#1. Normalitas residual
n.9 = length(r.9)
mean.9 = mean(r.9)
sd.9 = sd(r.9)
res.9 = rnorm(n.9,mean.9,sd.9)
normalitas.9 <- ks.test(r.9,res.9) ;normalitas.9
#2. Residual white noise (tidak ada autokorelasi)
WN.9 <- Box.test(r.9, type=c("Ljung-Box")) ;WN.9
#3. Homokedastisitas
h.9 = r.9^2
heteros.9 <- Box.test(h.9, type=c("Ljung-Box")) ;heteros.9

#Pemilihan model terbaik menggunakan AIC
AIC(model2)
AIC(model5)
AIC(model7)
AIC(model9)
summary(model2)
summary(model5)
summary(model7)
summary(model9) #Model 9 arima terbaik

#ARCH GARCH
#Pengujian efek heteroskedastisitas
resk <- (r.9)^2
Box.test(resk,lag=1,type="Ljung-Box")
Box.test(resk,lag=2,type="Ljung-Box")
Box.test(resk,lag=3,type="Ljung-Box")
library(FinTS)
return.archtest <- ArchTest(series, lags=1, demean=TRUE) ;return.archtest

#Identifikasi Model ARCH GARCH
plot(acf(resk), main="Squared Residual ACF")
plot(pacf(resk), main="Squared Residual PACF")
ggAcf(resk, lag=25) + ggtitle("Plot ACF Squared Residuals Model ARIMA(2,1,0)") #Q(0,1,2,3)
ggPacf(resk, lag=25) + ggtitle("Plot PACF Squared Residuals Model ARIMA(2,1,0)") #P(0,1,2,3,4)

#Pemilihan model
library(rugarch)
spec <- ugarchspec()
spec

#Model ARCH pake syntax teh novika 
#ARCH(1)
arch1 <- ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=TRUE),
                   variance.model=list(model="sGARCH",garchOrder=c(1,0)),distribution.model="norm")
fit1 <- ugarchfit(data=c(series), spec=arch1, out.sample=8,solver='hybrid')
fit1 #sig

#ARCH(2)
arch2 <- ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=TRUE),
                    variance.model=list(model="sGARCH",garchOrder=c(2,0)),distribution.model="norm")
fit2 <- ugarchfit(data=c(series), spec=arch2, out.sample=8,solver='hybrid')
fit2 #non sig

#ARCH(3)
arch3 <- ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=TRUE),
                    variance.model=list(model="sGARCH",garchOrder=c(3,0)),distribution.model="norm")
fit3 <- ugarchfit(data=c(series), spec=arch3, out.sample=8,solver='hybrid')
fit3 #non sig

#ARCH(4)
arch4 <- ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=TRUE),
                    variance.model=list(model="sGARCH",garchOrder=c(4,0)),distribution.model="norm")
fit4 <- ugarchfit(data=c(series), spec=arch4, out.sample=8,solver='hybrid')
fit4 #non sig

#GARCH(1,1)
arch5 <- ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=TRUE),
                    variance.model=list(model="sGARCH",garchOrder=c(1,1)),distribution.model="norm")
fit5 <- ugarchfit(data=c(series), spec=arch5, out.sample=8,solver='hybrid')
fit5 #sig

#Forecasting model ARCH(1)
ugarchforecast(fit1, data=NULL, n.ahead=15,n.roll=0, out.sample=8) #Model fit1 terbaik
