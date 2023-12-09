library(psych)
library(psychTools)
library(dplyr)
library(EFAtools)
library(GPArotation)
library(QuantPsyc)
library(parallel)

# Lendo os dados e fazendo algumas manipulações ----------------------
dados_poms <- read.csv("dados_poms.csv", sep = ";") |>
  filter(group == 2) |>
  dplyr::select(!c(1:3))

## Corrigindo os nomes das colunas
colnames(dados_poms) <- paste0("item", 1:30)

colnames(dados_poms) <- case_when(
  names(dados_poms) %in% paste0("item", c(9, 14, 23, 27, 30)) ~ paste0(names(dados_poms), "_tensao"),
  names(dados_poms) %in% paste0("item", c(5, 17, 20, 24, 28)) ~ paste0(names(dados_poms), "_depressao"),
  names(dados_poms) %in% paste0("item", c(6, 13, 18, 22, 29)) ~ paste0(names(dados_poms), "_raiva"),
  names(dados_poms) %in% paste0("item", c(1, 4, 10, 15, 19)) ~ paste0(names(dados_poms), "_vigor"),
  names(dados_poms) %in% paste0("item", c(3, 7, 11, 16, 25)) ~ paste0(names(dados_poms), "_fatiga"),
  names(dados_poms) %in% paste0("item", c(2, 8, 12, 21, 26)) ~ paste0(names(dados_poms), "_simpatia")
)

## Mudando a ordem das colunas para facilitar a visualização
dados_poms <- dados_poms |>
  dplyr::select(
    ends_with("raiva"), 
    ends_with("vigor"),
    ends_with("fatiga"),
    ends_with("simpatia"), 
    ends_with("depressao"), 
    ends_with("tensao")
  )

## Estimando a matriz de correlações policóricas
cor_poli <- polychoric(dados_poms)
corrplot::corrplot(
  cor_poli$rho,
  method = "color",
  type = "lower",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  diag = TRUE,
  number.cex = 0.6,
  tl.cex = 0.8,
  cl.pos = "r"
)


# Verificando se a análise fatorial é aplicável ---------------------------
## Verificando a plausibilidade do uso da AF pelo critério de Kaiser-Meyer-Olkin
KMO(cor_poli$rho)

## Aplicando o teste de esfericidade Bartlett 
BARTLETT(cor_poli$rho, N = 350)


# Determinando o número de fatores a serem retidos ------------------------
## Pelo critério de Kaiser (autovalores > 1)
autovalores <- eigen(cor_poli$rho)$values
round(autovalores, 2)  # Sugere que 6 fatores devem ser retidos

## Pelo critério do diagrama de declividade
plot(autovalores, type = "b", ylab = "Autovalores", pch = 19, xlab = "Número do fator")
abline(0.81, 0, lty = 2)  # Sugere que 6 fatores devem ser retidos

##Pelo critério da proporção da variância acumulada explicada pelos fatores
prop_explicada <- autovalores/sum(autovalores)
round(cumsum(prop_explicada), 2)  # Sugere que 4 fatores devem ser retidos



# Estimando as cargas fatoriais e fazendo a rotação dos fatores -----------
## Ajustando um modelo fatorial de 6 fatores
fa <- fa(cor_poli$rho, 6, rotate = "oblimin", fm = "minres")

## Loadings estimados
print(loadings(fa), cutoff = 0.3)


#  Verificando se cada grupo de itens é, realmente, unidimensional --------
df_tension <- dados_poms |>
  select_at(vars(matches("tensao")))

df_depresion <- dados_poms |>
  select_at(vars(matches("depressao")))

df_colera <- dados_poms |>
  select_at(vars(matches("raiva")))

df_vigor <- dados_poms |>
  select_at(vars(matches("vigor")))

df_fatiga <- dados_poms |>
  select_at(vars(matches("fatiga")))

df_amistad <- dados_poms |>
  select_at(vars(matches("simpatia")))

par(mfrow = c(3, 2))

## Raiva
plot(
  eigen(polychoric(df_colera)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Raiva"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_colera)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)

## Fatiga
plot(
  eigen(polychoric(df_fatiga)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Fatiga"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_fatiga)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)

## Vigor
plot(
  eigen(polychoric(df_vigor)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Vigor"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_vigor)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)

## Simpatia
plot(
  eigen(polychoric(df_amistad)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Simpatia"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_amistad)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)


## Tensão
plot(
  eigen(polychoric(df_tension)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Tensão"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_tension)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)


## Depressão
plot(
  eigen(polychoric(df_depresion)$rho)$values,
  type = "b", 
  main = 'Para os itens do domínio "Depressão"',
  ylab = "Autovalor", 
  xlab = "Número do fator",
  pch = 19
)

cor_poli <- polychoric(df_depresion)$rho

autovalores <- eigen(cor_poli)$values
round(autovalores, 2)

(prop_explicada <- round(autovalores/sum(autovalores) * 100, 2))
round(cumsum(prop_explicada), 2)


# Calculando as estimativas e ICs dos coeficientes de confiabilide --------
## Criando as funções que geram intervalos BCa
gera_est_boot <- function(x_trans, n, n_iter_boot) {
  cl <- makeCluster(getOption("cl.cores", 2))
  clusterEvalQ(cl, { library(psych) })
  
  boot <- parLapply(cl, 1:n_iter_boot, function(i, x_trans, n, n_iter_boot) {
    #Criando a reamostra bootstrap
    amostra_ruim <- TRUE
    while (amostra_ruim) {
      amostra_ruim <- tryCatch({
        data_boot <- x_trans[sample.int(n, replace = TRUE), ]
        cor_poli <- polychoric(data_boot)$rho
        cor_pearson <- cor(data_boot, use = "pairwise")
        fit_fa <- fa(cor_poli, 1, fm = "minres")
        fit_fa_pearson <- fa(data_boot, 1, fm = "minres")
        amostra_ruim <- FALSE
      },
      warning = function(cond) return(TRUE)
      )
    }
    
    #Calculando o alfa de Cronbach e o alfa ordinal
    n.itens <- dim(cor_pearson)[2]
    alfa <- (1 - tr(cor_pearson)/sum(cor_pearson)) * (n.itens/(n.itens - 1))
    alfa_ordinal <- (1 - tr(cor_poli)/sum(cor_poli)) * (n.itens/(n.itens - 1))
    
    #Calculando o omega e o omega ordinal
    omega <- (sum(fit_fa_pearson$loadings[,1]))^2 / ((sum(fit_fa_pearson$loadings[,1]))^2 + sum(fit_fa_pearson$uniquenesses))
    omega_ordinal <- (sum(fit_fa$loadings[,1]))^2 / ((sum(fit_fa$loadings[,1]))^2 + sum(fit_fa$uniquenesses))
    
    return(
      list(
        alfa = alfa,
        alfa_ordinal = alfa_ordinal,
        omega = omega,
        omega_ordinal = omega_ordinal
      )
    )
  }, x_trans, n, n_iter_boot) 
  
  stopCluster(cl)
  
  estimativas_boot <- unlist(boot)
  
  alfa_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "alfa"])
  alfa_ordinal_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "alfa_ordinal"])
  omega_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "omega"])
  omega_ordinal_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "omega_ordinal"])
  
  return(
    list(
      alfa_boot = alfa_boot,
      alfa_ordinal_boot = alfa_ordinal_boot,
      omega_boot = omega_boot,
      omega_ordinal_boot = omega_ordinal_boot
    )
  )
}

gera_boot_delete <- function(x_trans, n) {
  ### calculating the a_0_hat (acceleration) 
  cl <- makeCluster(getOption("cl.cores", 2))
  clusterEvalQ(cl, { library(psych) })
  boot <- parLapply(cl, 1:n, function(i, x_trans, n_iter_boot, poly) { 
    cor_poli <- polychoric(x_trans[-i, ])$rho
    cor_pearson <- cor(x_trans[-i, ], use = "pairwise")
    fit_fa <- fa(cor_poli, 1, fm = "minres")
    fit_fa_pearson <- fa(x_trans[-i, ], 1, fm = "minres")
    
    #Calculando o alfa de Cronbach e o alfa ordinal
    n.itens <- dim(cor_pearson)[2]
    alfa <- (1 - tr(cor_pearson)/sum(cor_pearson)) * (n.itens/(n.itens - 1))
    alfa_ordinal <- (1 - tr(cor_poli)/sum(cor_poli)) * (n.itens/(n.itens - 1))
    
    #Calculando o omega e o omega ordinal
    omega <- (sum(fit_fa_pearson$loadings[,1]))^2 / ((sum(fit_fa_pearson$loadings[,1]))^2 + sum(fit_fa_pearson$uniquenesses))
    omega_ordinal <- (sum(fit_fa$loadings[,1]))^2 / ((sum(fit_fa$loadings[,1]))^2 + sum(fit_fa$uniquenesses))
    
    return(
      list(
        alfa = alfa,
        alfa_ordinal = alfa_ordinal,
        omega = omega,
        omega_ordinal = omega_ordinal
      )
    )
  }, x_trans, n) 
  
  stopCluster(cl)
  
  estimativas_boot_delete <- unlist(boot)
  
  alfa_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "alfa"])
  alfa_ordinal_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "alfa_ordinal"])
  omega_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "omega"])
  omega_ordinal_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "omega_ordinal"])
  
  return(
    list(
      alfa_boot = alfa_boot_delete,
      alfa_ordinal = alfa_ordinal_boot_delete,
      omega_boot = omega_boot_delete,
      omega_ordinal_boot = omega_ordinal_boot_delete
    )
  )
  
}

gera_bca <- function(estimativa_original, estimativas_boot, estimativas_boot_delete, n_iter_boot) {
  ### calculating z0_hat (bias-correction)
  p1 <- which(estimativas_boot < estimativa_original)
  z0 <- qnorm(length(p1)/n_iter_boot)
  
  a0 <- sum((mean(estimativas_boot_delete)-estimativas_boot_delete)^3)/(6*(sum((mean(estimativas_boot_delete)-estimativas_boot_delete)^2))^(3/2))
  
  #### Obtaing the bootstrap interval of coverage
  alpha1 <- pnorm(z0 + (z0 + qnorm(0.05/2))/(1-a0*(z0 + qnorm(0.05/2))))
  alpha2 <- pnorm(z0 + (z0 + qnorm(1-(0.05/2)))/(1-a0*(z0 + qnorm(1-(0.05/2)))))
  
  #### Final statistics
  as.numeric(c(quantile(estimativas_boot, prob = alpha1), quantile(estimativas_boot, prob = alpha2)))
}

## Calculando as estimativas dos coeficientes de confiabilidade para cada grupo de itens
grupos <- c("colera", "fatiga", "vigor", "amistad", "depresion", "tension")
df_completo <- data.frame()

for (grupo in grupos) {
  set.seed(412)
  
  lista_estimativas_boot <- gera_est_boot(x_trans = get(paste0("df_", grupo)), n = 350, n_iter_boot = 1000)
  lista_estimativas_boot_delete <- gera_boot_delete(x_trans = get(paste0("df_", grupo)), n = 350)
  
  cor_poli <- polychoric(get(paste0("df_", grupo)))$rho
  fit_fa <- fa(cor_poli, 1, rotate = "oblimin", fm = "minres")
  fit_fa_pearson <- fa(get(paste0("df_", grupo)), 1, rotate = "oblimin", fm = "minres")
  
  #Calculando a estimativa e o IC do alfa 
  alfa_est <- as.numeric(alpha(get(paste0("df_", grupo)))$total[2])
  alfa_ci <- gera_bca(
    estimativa_original = alfa_est,
    estimativas_boot = lista_estimativas_boot$alfa_boot,
    estimativas_boot_delete = lista_estimativas_boot_delete$alfa_boot,
    n_iter_boot = 1000
  )
  alfa_lower <- alfa_ci[1]
  alfa_upper <- alfa_ci[2]
  
  
  #Calculando a estimativa e o IC do alfa ordinal 
  alfa_ordinal_est <- as.numeric(alpha(cor_poli)$total[1])
  alfa_ordinal_ci <- gera_bca(
    estimativa_original = alfa_ordinal_est,
    estimativas_boot = lista_estimativas_boot$alfa_ordinal_boot,
    estimativas_boot_delete = lista_estimativas_boot_delete$alfa_ordinal,
    n_iter_boot = 1000
  )
  alfa_ordinal_lower <- alfa_ordinal_ci[1]
  alfa_ordinal_upper <- alfa_ordinal_ci[2]
  
  
  #Calculando a estimativa e o IC do omega
  omega_est <- (sum(fit_fa_pearson$loadings[,1]))^2 / ((sum(fit_fa_pearson$loadings[,1]))^2 + sum(fit_fa_pearson$uniquenesses))
  omega_ci <- gera_bca(
    estimativa_original = omega_est,
    estimativas_boot = lista_estimativas_boot$omega_boot,
    estimativas_boot_delete = lista_estimativas_boot_delete$omega_boot,
    n_iter_boot = 1000
  )
  omega_lower <- omega_ci[1]
  omega_upper <- omega_ci[2]
  
  #Calculando a estimativa e o IC do omega ordinal
  omega_ordinal_est <- (sum(fit_fa$loadings[,1]))^2 / ((sum(fit_fa$loadings[,1]))^2 + sum(fit_fa$uniquenesses))
  omega_ordinal_ci <- gera_bca(
    estimativa_original = omega_ordinal_est,
    estimativas_boot = lista_estimativas_boot$omega_ordinal_boot,
    estimativas_boot_delete = lista_estimativas_boot_delete$omega_ordinal_boot,
    n_iter_boot = 1000
  )
  
  df_grupo <- data.frame(
    grupo = grupo,
    alfa = paste(round(alfa_est, 3), glue::glue("({round(alfa_ci[1], 3)}; {round(alfa_ci[2], 3)})")),
    alfa_ordinal = paste(round(alfa_ordinal_est, 3), glue::glue("({round(alfa_ordinal_ci[1], 3)}; {round(alfa_ordinal_ci[2], 3)})")),
    omega = paste(round(omega_est, 3), glue::glue("({round(omega_ci[1], 3)}; {round(omega_ci[2], 3)})")),
    omega_ordinal = paste(round(omega_ordinal_est, 3), glue::glue("({round(omega_ordinal_ci[1], 3)}; {round(omega_ordinal_ci[2], 3)})"))
  )
  
  df_completo <- rbind(df_completo, df_grupo)
}

write.csv(df_completo, "resultados_aplicacao_questionario.csv", row.names = FALSE)
