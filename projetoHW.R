#--------------------------------------------------------------------------
# Séries Temporais - Suavização exponencial pelo método de Holt-Winters
#--------------------------------------------------------------------------

# Definindo a biblioteca e os pacotes não padrões a serem utilizados.
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(ggplot2)

data("AirPassengers")

anos.prev <- 3
s <- 12
h <- s*anos.prev
Zcplt <- AirPassengers
Zmdl <- Zcplt[1:(length(Zcplt)-h)]

#--------------------------------------------------------------------------
# Equações de recorrência - parâmetros iniciais

# Nível z(t)
z <- c(rep(NA, s-1), mean(Zdml[1:s]))

# Tendência te(t)
te <- c(rep(NA, s-1), 0)

# Sazonalidade sa(t)
sa <- Zmdl[1:s]/z[s]

#--------------------------------------------------------------------------
# Ajuste da modelagem - Suavização
a = 0.25
c = 0.02
d = 0.95

for (i in (s+1):length(Zmdl)) {
      z[i] <- a*(Zmdl[i]-sa[i-s])+(1-a)*(z[i-1]+te[i-1])
      te[i] <- c*(z[i]-z[i-1])+(1-c)*te[i-1]
      sa[i] <- d*(Zmdl[i]-z[i])+(1-d)*sa[i-s]
}

Zfit <- c(rep(NA, s-1), z + te + sa)

#--------------------------------------------------------------------------
# Previsão desatualizada

Zprev1 <- rep(NA, length(Zmdl))
n <- 1
for (p in 1:h) {
      if (p > n*s) n <- n + 1
      zp <- z[length(Zmdl)] + p*te[length(Zmdl)] + sa[length(Zmdl)+p-n*s]
      Zprev1 <- c(Zprev1, zp)
}

# Previsão atualizada

for (i in (length(Zmdl)+1):length(Zcplt)) {
      Znew <- z[i-1] + te[i-1] + sa[i-s]
      z[i] <- a*(Znew-sa[i-s])+(1-a)*(z[i-1]+te[i-1])
      te[i] <- c*(z[i]-z[i-1])+(1-c)*te[i-1]
      sa[i] <- d*(Znew-z[i])+(1-d)*sa[i-s]
}

Zprev2 <- c(rep(NA,length(Zmdl)),
            (z + te + sa)[(length(Zmdl)+1):length(Zcplt)])


rest <- (length(Zmdl)+1):length(Zcplt)
#--------------------------------------------------------------------------
# Gráficos utilizando ggplot2
dataAP1 <- data.frame(Time = 1:length(Zmdl), Data = Zmdl, isin = "Dados")
dataAP2 <- data.frame(Time = (s+1):length(Zmdl), Dados = Zfit,
                      isin = "Ajuste")
dataAP3 <- data.frame(Time = rest, Dados = Zcplt[rest],
                      isin = "Test Prev")
dataAP4 <- data.frame(Time = rest, Dados = Zprev1,
                      isin = "Prev OTD")
dataAP5 <- data.frame(Time = rest, Dados = Zprev2,
                      isin = "Prev UPD")

dataAPprevOTD <- rbind(dataAP1, dataAP2, dataAP3, dataAP4)
dataAPprevUPD <- rbind(dataAP1, dataAP2, dataAP3, dataAP5)
dataAPprevs <- rbind(dataAP3, dataAP4, dataAP5)

p1 <- ggplot(dataAPprevOTD, aes(x = Time, y = Data)) +
      geom_line(color = "cyan2", size = 1.1) +
      geom_line(mapping = aes(y = Fitted), color = "darkorange",
                size = 1.1) +
      geom_vline(xintercept = length(Zmdl), lty = 2) +
      
      

      
      geom_line(mapping = aes(y = Prev_otd), color = "darkslategray4",
                size = 1.1)
      






