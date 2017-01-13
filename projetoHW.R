#--------------------------------------------------------------------------
# Séries Temporais - Suavização exponencial pelo método de Holt-Winters
#--------------------------------------------------------------------------

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Acessando a série que será analisada
data("AirPassengers")

ap <- 3 # Anos a serem Previstos
s <- 12 # frequência Sazonal
h <- s*ap # quant. meses a serem previstos
Zcplt <- AirPassengers # Série completa
Zmdl <- Zcplt[1:(length(Zcplt)-h)] # Séria utilizada na modelagem

#--------------------------------------------------------------------------
# Ajuste da modelagem - Suavização exponencial

# Equações de recorrência - parâmetros iniciais
# Nível incial z
z <- c(rep(NA, s-1), mean(Zmdl[1:s])) # OU c(rep(NA, s-1), 120.4284)

# Tendência inicial t
t <- c(rep(NA, s-1), 0) # OU c(rep(NA, s-1),1.9035) #c(rep(NA, s-1), 0)

# Sazonalidade inicial f
f <- Zmdl[1:s]-z[s] # OU c(22.1981, -43.8653, -17.7495, 12.8625, 43.4763,
#                          50.1649,26.7827, -4.5286, -2.9447, 4.7007,
#                          -26.4387, -20.2622)

# Constantes de suavização
A <- 0.25 # 0.9999
C <- 0.02 # 0.0004
D <- 0.9  # 0.0004

# Equações de estimativa - nível z(t), tendência te(t) e sazonalidade sa(t) 
for (i in (s+1):length(Zmdl)) {
      z[i] <- A*(Zmdl[i]-f[i-s])+(1-A)*(z[i-1]+t[i-1])
      t[i] <- C*(z[i]-z[i-1])+(1-C)*t[i-1]
      f[i] <- D*(Zmdl[i]-z[i])+(1-D)*f[i-s]
}

# Série de ajuste da modelagem Zfit(t)
Zfit <- c(z + t + f, rep(NA, h))

#--------------------------------------------------------------------------
# Previsão sem atualização a cada amostra (outdated)

# Zprev1 se refere a previsão desatualizada; Zprev2, a previsão atualizada
Zprev1 <- Zprev2 <- rep(NA, length(Zmdl))
n <- 1 
for (p in 1:h) {
      if (p > n*s) n <- n + 1
      # zp estima a série p passos a frente do último valor de Zmdl
      zp <- z[length(Zmdl)] + p*t[length(Zmdl)] + f[length(Zmdl)+p-n*s]
      Zprev1 <- c(Zprev1, zp)
}

# Previsão com atualização a cada amostra (updated)

m <- (length(Zmdl)+1):length(Zcplt) # m indica os meses a serem previstos
for (i in m) {
      Znew <- z[i-1] + te[i-1] + sa[i-s]
      z[i] <- A*(Znew-f[i-s])+(1-A)*(z[i-1]+t[i-1])
      t[i] <- C*(z[i]-z[i-1])+(1-C)*t[i-1]
      f[i] <- D*(Znew-z[i])+(1-D)*f[i-s]
}

Zprev2 <- c(Zprev2, (z + t + f)[m])

#--------------------------------------------------------------------------
# Gráficos utilizando o pacote ggplot2

library(ggplot2)
dataAP <- data.frame(
      Date = rep(1:length(Zcplt), 4),
      Data = as.numeric(c(Zcplt, Zfit, Zprev1, Zprev2)),
      Type = gl(4, length(Zcplt),
                labels = c("Dados", "Ajuste", "Prev OTD", "Prev UPD")))

dataAP <- dataAP[complete.cases(dataAP),]

dataAPprevOTD <- subset(dataAP, Type!="Prev UPD")
dataAPprevUPD <- subset(dataAP, Type!="Prev OTD")
dataAPprevBOTH <- subset(dataAP, Type!="Ajuste")

# Gráfico 1 - Previsão sem atualização (outdated)
p1 <- ggplot(dataAPprevOTD, aes(x=Date, y=Data, color=Type)) +
      geom_line(size = 1.1) +
      geom_vline(xintercept = length(Zmdl), lty = 2) +
      scale_color_manual(
            values = c("cyan2", "darkorange", "darkslategray4")) +
      scale_x_continuous(
            breaks = seq(1, length(Zcplt), by = 12),
            labels = paste0(rep("Jan/", end(Zcplt)[1]-start(Zcplt)[1]+1),
                            start(Zcplt)[1]:end(Zcplt)[1])) +
      labs(title="Método de Holt-Winters \n Previsão Desatualizada (OTD)",
           x = "Tempo (mês/ano)", y = "Nº de Passageiros (1000s)")

# Gráfico 2 - Previsão com atualização (updated) 
p2 <- ggplot(dataAPprevUPD, aes(x=Date, y=Data, color=Type)) +
      geom_line(size = 1.1) +
      geom_vline(xintercept = length(Zmdl), lty = 2) +
      scale_color_manual(
            values = c("cyan2", "darkorange", "darkslategray4")) +
      scale_x_continuous(
            breaks = seq(1, length(Zcplt), by = 12),
            labels = paste0(rep("Jan/", end(Zcplt)[1]-start(Zcplt)[1]+1),
                            start(Zcplt)[1]:end(Zcplt)[1])) +
      labs(title="Método de Holt-Winters \n Previsão Atualizada (UPD)",
           x = "Tempo (mês/ano)", y = "Nº de Passageiros (1000s)")

# Gráfico 3 - Comparação entre previsões
p3 <- ggplot(dataAPprevBOTH, aes(x=Date, y=Data, color=Type)) +
      geom_line(size = 1.1) +
      scale_color_manual(
            values = c("cyan2", "darkorange", "darkslategray4")) +
      scale_x_continuous(
            breaks = seq(1, length(Zcplt), by = 4),
            labels = paste(rep(c("Jan","May","Sep"),
                               end(Zcplt)[1]-start(Zcplt)[1]+1),
                           sort(rep(start(Zcplt)[1]:end(Zcplt)[1], 3))),
            limits = c(length(Zmdl), length(Zcplt))) +
      labs(title="Método de Holt-Winters \n Comparação entre previsões",
           x = "Tempo (mês/ano)", y = "Nº de Passageiros (1000s)")

#__________________________________________________________________________
#--------------------------------------------------------------------------
# Comparação com os dados da função ets(ts, model="AAA")
library(forecast)
fit <- ets(ts(Zmdl, start = c(1949,1), frequency = 12), model = "AAA")
pred <- forecast(fit, h)

# Gráfico 4 - Comparação com os resultados do pacote forecast
forecastPrev <- data.frame(Date = m, Data = pred$mean, Type = "Prev FCST")
dataAPforecast <- rbind(dataAPprevBOTH, forecastPrev)
p4 <- ggplot(dataAPforecast, aes(x=Date, y=Data, color=Type)) +
      geom_line(size = 1.1) +
      scale_color_manual(
            values = c("cyan2", "darkorange", "darkslategray4", "red3")) +
      scale_x_continuous(
            breaks = seq(1, length(Zcplt), by = 4),
            labels = paste(rep(c("Jan","May","Sep"),
                               end(Zcplt)[1]-start(Zcplt)[1]+1),
                           sort(rep(start(Zcplt)[1]:end(Zcplt)[1], 3))),
            limits = c(length(Zmdl), length(Zcplt))) +
      labs(title = "Suavização Exponencial - Método de Holt-Winters",
           x = "Tempo (anos)", y = "Nº de Passageiros (1000s)")




