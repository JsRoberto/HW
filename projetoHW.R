#--------------------------------------------------------------------------
# Séries Temporais - Suavização exponencial pelo método de Holt-Winters
#--------------------------------------------------------------------------

# Definindo a biblioteca e os pacotes não padrões a serem utilizados.
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(ggplot2)

data("AirPassengers")

ap <- 3
s <- 12
h <- s*ap
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
dataAP1 <- data.frame(Time = 1:length(Zcplt), Data = Zcplt, isin = "Dados")
dataAP2 <- data.frame(Time = (s+1):length(Zmdl), Data = Zfit,
                      isin = "Ajuste")
dataAP3 <- data.frame(Time = rest, Data = Zprev1,
                      isin = "Prev OTD")
dataAP4 <- data.frame(Time = rest, Data = Zprev2,
                      isin = "Prev UPD")

dataAPprevOTD <- rbind(dataAP1, dataAP2, dataAP3)
dataAPprevUPD <- rbind(dataAP1, dataAP2, dataAP4)
dataAPprevs <- rbind(dataAP1[rest,], dataAP3, dataAP4)

p1 <- ggplot(dataAPprevOTD, aes(x=Time, y=Data, color=isin))
      #geom_smooth(aes(x=Time, y=Data, ymax=Data+2*sd(Data), ymin=Data-2*sd(Data)), 
      #            data = df5, colour = 'cyan2', stat = 'identity') +
      #scale_x_continuous(breaks = seq(1, length(Zcplt), by = 6),
      #                   labels = 
      #                         paste0(rep(c("Jan","Jul"),
      #                                    end(Zcplt)[1]-start(Zcplt)[1]+1),
      #                                seq(start(Zcplt)[1], end(Z)[1]-ap,
      #                                    by = 2))) +

p2 <- ggplot(dataAPprevUPD, aes(x=Time, y=Data, color=isin))

p3 <- ggplot(dataAPprevs, aes(x=Time, y=Data, color=isin))

conf <- scale_color_manual(values = c("cyan2", "darkorange",
                                        "darkslategray4")) +
        geom_line(size = 1.1) +
        geom_vline(xintercept = length(Zmdl), lty = 2) +
        labs(title = "Suavização Exponencial - Método de Holt-Winters",
             x = "Tempo (anos)", y = "Nº de Passageiros (1000s)")      

p1 + conf
p2 + conf
p3 + conf

