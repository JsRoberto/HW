#--------------------------------------------------------------------------
# Séries Temporais - Suavização exponencial pelo método de Holt-Winters
#--------------------------------------------------------------------------

# Definindo a biblioteca e os pacotes não padrões a serem utilizados.
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(ggplot2)

data("AirPassengers")

ap <- 3 # Anos a serem Previstos
s <- 12 # frequência Sazonal
h <- s*ap # quant. meses a serem previstos
Zcplt <- AirPassengers # Série completa
Zmdl <- Zcplt[1:(length(Zcplt)-h)] # Séria utilizada na modelagem

#--------------------------------------------------------------------------
# Equações de recorrência - parâmetros iniciais

# Nível z(t)
z <- c(rep(NA, s-1), mean(Zmdl[1:s]))

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

Zfit <- z + te + sa

#--------------------------------------------------------------------------
# Previsão desatualizada

Zprev1 <- rep(NA, length(Zmdl))
n <- 1
for (p in 1:h) {
      if (p > n*s) n <- n + 1
      zp <- z[length(Zmdl)] + p*te[length(Zmdl)] + sa[length(Zmdl)+p-n*s]
      Zprev1 <- c(Zprev1, zp)
}

#--------------------------------------------------------------------------
# Previsão atualizada

for (i in (length(Zmdl)+1):length(Zcplt)) {
      Znew <- z[i-1] + te[i-1] + sa[i-s]
      z[i] <- a*(Znew-sa[i-s])+(1-a)*(z[i-1]+te[i-1])
      te[i] <- c*(z[i]-z[i-1])+(1-c)*te[i-1]
      sa[i] <- d*(Znew-z[i])+(1-d)*sa[i-s]
}

Zprev2 <- c(rep(NA,length(Zmdl)),
            (z + te + sa)[(length(Zmdl)+1):length(Zcplt)])


tempo <- 1:length(Zcplt)

#--------------------------------------------------------------------------
# Gráficos utilizando ggplot2
dataAP1 <- data.frame(Time = tempo, Data = as.numeric(Zcplt),
                      isin = "Dados")
dataAP2 <- data.frame(Time = 1:length(Zmdl), Data = as.numeric(Zfit),
                      isin = "Ajuste")
dataAP3 <- data.frame(Time = tempo, Data = as.numeric(Zprev1),
                      isin = "Prev OTD")
dataAP4 <- data.frame(Time = tempo, Data = as.numeric(Zprev2),
                      isin = "Prev UPD")

dataAPprevOTD <- rbind(dataAP1, dataAP2[complete.cases(dataAP2),],
                       dataAP3[complete.cases(dataAP3),])
dataAPprevUPD <- rbind(dataAP1, dataAP2[complete.cases(dataAP2),],
                       dataAP4[complete.cases(dataAP4),])
dataAPprevs <- rbind(dataAP1[complete.cases(dataAP3),], dataAP3, dataAP4)

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

# Gráfico 1 - Previsão desatualizada
p1 + scale_color_manual(values = c("cyan2", "darkorange",
                                   "darkslategray4")) +
      geom_line(size = 1.1) +
      geom_vline(xintercept = length(Zmdl), lty = 2) +
      labs(title = "Suavização Exponencial - Método de Holt-Winters",
           x = "Tempo (anos)", y = "Nº de Passageiros (1000s)")

# Gráfico 2 - Previsão atualizada 
p2 + scale_color_manual(values = c("cyan2", "darkorange",
                                   "darkslategray4")) +
      geom_line(size = 1.1) +
      geom_vline(xintercept = length(Zmdl), lty = 2) +
      labs(title = "Suavização Exponencial - Método de Holt-Winters",
           x = "Tempo (anos)", y = "Nº de Passageiros (1000s)")

# Gráfico 3 - Comparação entre previsões
p3 + scale_color_manual(values = c("cyan2", "darkorange",
                                   "darkslategray4")) +
      geom_line(size = 1.1) +
      #geom_vline(xintercept = length(Zmdl), lty = 2) +
      labs(title = "Suavização Exponencial - Método de Holt-Winters",
           x = "Tempo (anos)", y = "Nº de Passageiros (1000s)")

# Observações: 
# 
# 1) O que falta para o código? Falta deixar mais enxuto o código, melhorar
# os comentários e alguns aspectos dos gráficos.
# 
# 2) Fazer comparação dos resultados obtidos pelo nossa implementação do
# algoritmo de Holt-Winters com duas implementações consolidadas: com a 
# função HoltWinters() da instalação básica e a função ets() do pacote 
# "forecast" (escreverei alguns comentários sobre essas implementações para
# o relatório ). Ref.: R in action. pp. 352-9 

