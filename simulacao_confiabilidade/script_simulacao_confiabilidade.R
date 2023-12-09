library(MASS)
library(psych)
library(dplyr)
library(MBESS)
library(parallel)

simulacao_confiabilidade <- function(M, inicio, n_iter_boot, vec_n, vec_num_itens, vec_num_pontos, vec_confiabilidade, var_erros, paralelo, vec_nivel_variacao_loadings) {
  #Para uma dada confiabilidade no vetor vec_confiabilidade
  do.call(rbind, lapply(vec_confiabilidade, function(confiabilidade) {
    
    #Para um dado "nível de variação" dos loadings no vetor vec_nivel_variacao_loadings
    do.call(rbind, lapply(vec_nivel_variacao_loadings, function(nivel_variacao_loadings) {
      
      #Para um dado número de itens no vetor num_itens
      do.call(rbind, lapply(vec_num_itens, function(num_itens) {
        
        #Para um dado tamanho de amostra no vetor vec_n
        do.call(rbind, lapply(vec_n, function(n) {
          
          #Para um dado número de pontos na escala no vetor vec_num_pontos
          df_iteracao <- do.call(rbind, lapply(vec_num_pontos, function(num_pontos) {
            
            #Criando a matriz que irá receber as estimativas dos coeficientes
            matriz_coeficientes <- matrix(0, nrow = M, ncol = 5, dimnames = list(1:M, c("alfa", "alfa_std", "alfa_ordinal", "omega", "omega_ordinal")))
            
            #Criando a matriz que irá receber as informações sobre a cobertura e a amplitude dos ICs
            matriz_ics <- matrix(
              0, nrow = M, ncol = 10,
              dimnames = list(
                1:M,
                c("alfa_cobre", "alfa_amplitude", "alfa_std_cobre", "alfa_std_amplitude", "alfa_ordinal_cobre", "alfa_ordinal_amplitude",
                  "omega_cobre", "omega_amplitude", "omega_ordinal_cobre", "omega_ordinal_amplitude")
              )
            )
            
            gera_dados <- function() {
              #Gerando os valores do fator em cada indivíduo
              fator <- rnorm(n)
              
              #Gerando os erros de cada indivíduo e em cada item
              erros <- mvrnorm(n, mu = rep(0, num_itens), Sigma = diag(num_itens)*var_erros)  
              
              if (paralelo == TRUE) {  #Se os itens forem paralelos
                #Calculando o valor único dos loadings
                loadings <- sqrt((confiabilidade*num_itens*var_erros)/(num_itens^2 - num_itens^2*confiabilidade))
                
                #Calculando o valor do escore observado de cada indivíduo em cada item
                x <- matrix(loadings*fator + erros, nrow = n, ncol = num_itens)
                
                #Calculando o desvio padrão da distribuição teórica dos escores observados (será o mesmo em todos os itens)
                sd_x <- sqrt(loadings^2 + var_erros)
                
                #Criando a matriz que irá receber os escores transformados
                x_trans <- matrix(1, nrow(x), ncol(x))
                
                #Populando a matriz de escores transformados a partir da distribuição teórica dos escores observados 
                #(todos os itens tem distribuição N(0, sd_x))
                for (i in 1:nrow(x)){
                  for (j in 1:ncol(x)){
                    quantil_base <- 1/num_pontos
                    for (ponto in 1:(num_pontos - 1)) {
                      if ((x[i, j] > qnorm(quantil_base*ponto, sd = sd_x)) & (x[i, j] <= qnorm(quantil_base*(ponto + 1), sd = sd_x))) {
                        x_trans[i, j] <- ponto + 1
                      }
                    }
                  }
                }
              } else {  #Se os itens não forem paralelos
                #Calculando o valor da soma dos loadings
                soma_loadings <- sqrt((confiabilidade*num_itens*var_erros)/(1 - confiabilidade))
                
                #Definindo o valor do loading de cada item de forma semelhante ao que vimos naquele site
                if (confiabilidade == 0.9) {
                  fator_loadings <- 0.005
                } else {
                  fator_loadings <- 0.025
                }
                
                loadings <- seq(
                  soma_loadings/5 + nivel_variacao_loadings * fator_loadings,
                  soma_loadings/5 - nivel_variacao_loadings * fator_loadings, 
                  length = num_itens
                )
                
                #Calculando o valor do escore observado de cada indivíduo em cada item
                x <- t(matrix(loadings, nrow = num_itens) %*% t(matrix(fator))) + erros
                
                #Calculando os desvios padrões das distribuições teóricas de cada item (cada item terá um desvio padrão diferente)
                sd_x <- sqrt(loadings^2 + var_erros)
                
                #Criando a matriz que irá receber os escores transformados
                x_trans <- matrix(1, nrow(x), ncol(x))
                
                #Populando a matriz de escores transformados a partir da distribuição teórica dos escores observados 
                #cada item tem distribuição N(0, sd_x[j]))
                for (i in 1:nrow(x)){
                  for (j in 1:ncol(x)){
                    quantil_base <- 1/num_pontos
                    for (ponto in 1:(num_pontos - 1)) {
                      if ((x[i, j] > qnorm(quantil_base*ponto, sd = sd_x[j])) & (x[i, j] <= qnorm(quantil_base*(ponto + 1), sd = sd_x[j]))) {
                        x_trans[i, j] <- ponto + 1
                      }
                    }
                  }
                }
              } 
              return(x_trans)
            }
            
            for (k in inicio:M) {
              set.seed(432 + k)
              
              amostra_ruim <- TRUE
              while (amostra_ruim) {
                amostra_ruim <- tryCatch({
                  x_trans <- gera_dados()
                  cor_poli <- polychoric(x_trans)$rho
                  fit_fa <- fa(x_trans, 1, fm = "minres")
                  fit_fa_poli <- fa(cor_poli, 1, fm = "minres")
                  amostra_ruim <- FALSE
                },
                warning = function(cond) return(TRUE)
                )
              }
              
              #Criando a função que será utilizada na obtenção do IC das correlações por bootstrap
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
                      cov_pearson <- cov(data_boot, use = "pairwise")
                      cor_pearson <- cor(data_boot, use = "pairwise")
                      fit_fa_poli <- fa(cor_poli, 1, fm = "minres")
                      fit_fa_pearson <- fa(data_boot, 1, fm = "minres")
                      amostra_ruim <- FALSE
                    },
                    warning = function(cond) return(TRUE)
                    )
                  }
                  
                  #Calculando os alfas
                  n.itens <- dim(cov_pearson)[2]
                  alfa <- (1 - tr(cov_pearson)/sum(cov_pearson)) * (n.itens/(n.itens - 1))
                  alfa_std <- (1 - tr(cor_pearson)/sum(cor_pearson)) * (n.itens/(n.itens - 1))
                  alfa_ordinal <- (1 - tr(cor_poli)/sum(cor_poli)) * (n.itens/(n.itens - 1))
                  
                  #Calculando o omega e o omega ordinal
                  omega <- (sum(fit_fa_pearson$loadings[,1]))^2 / ((sum(fit_fa_pearson$loadings[,1]))^2 + sum(fit_fa_pearson$uniquenesses))
                  omega_ordinal <- (sum(fit_fa_poli$loadings[,1]))^2 / ((sum(fit_fa_poli$loadings[,1]))^2 + sum(fit_fa_poli$uniquenesses))
                  
                  return(
                    list(
                      alfa = alfa,
                      alfa_std = alfa_std,
                      alfa_ordinal = alfa_ordinal,
                      omega = omega,
                      omega_ordinal = omega_ordinal
                    )
                  )
                }, x_trans, n, n_iter_boot) 
                
                stopCluster(cl)
                
                estimativas_boot <- unlist(boot)
                
                alfa_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "alfa"])
                alfa_std_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "alfa_std"])
                alfa_ordinal_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "alfa_ordinal"])
                omega_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "omega"])
                omega_ordinal_boot <- as.vector(estimativas_boot[names(estimativas_boot) == "omega_ordinal"])
                
                return(
                  list(
                    alfa_boot = alfa_boot,
                    alfa_std_boot = alfa_std_boot,
                    alfa_ordinal_boot = alfa_ordinal_boot,
                    omega_boot = omega_boot,
                    omega_ordinal_boot = omega_ordinal_boot
                  )
                )
              }
              
              gera_boot_delete <- function(x_trans, n) {
                cl <- makeCluster(getOption("cl.cores", 2))
                clusterEvalQ(cl, { library(psych) })
                
                boot <- parLapply(cl, 1:n, function(i, x_trans, n_iter_boot, poly) { 
                  cor_poli <- polychoric(x_trans[-i, ])$rho
                  cov_pearson <- cov(x_trans[-i, ], use = "pairwise")
                  cor_pearson <- cor(x_trans[-i, ], use = "pairwise")
                  fit_fa_poli <- fa(cor_poli, 1, fm = "minres")
                  fit_fa_pearson <- fa(x_trans[-i, ], 1, fm = "minres")
                  
                  #Calculando o alfa de Cronbach e o alfa ordinal
                  n.itens <- dim(cov_pearson)[2]
                  alfa <- (1 - tr(cov_pearson)/sum(cov_pearson)) * (n.itens/(n.itens - 1))
                  alfa_std <- (1 - tr(cor_pearson)/sum(cor_pearson)) * (n.itens/(n.itens - 1))
                  alfa_ordinal <- (1 - tr(cor_poli)/sum(cor_poli)) * (n.itens/(n.itens - 1))
                  
                  #Calculando o omega e o omega ordinal
                  omega <- (sum(fit_fa_pearson$loadings[,1]))^2 / ((sum(fit_fa_pearson$loadings[,1]))^2 + sum(fit_fa_pearson$uniquenesses))
                  omega_ordinal <- (sum(fit_fa_poli$loadings[,1]))^2 / ((sum(fit_fa_poli$loadings[,1]))^2 + sum(fit_fa_poli$uniquenesses))
                  
                  return(
                    list(
                      alfa = alfa,
                      alfa_std = alfa_std,
                      alfa_ordinal = alfa_ordinal,
                      omega = omega,
                      omega_ordinal = omega_ordinal
                    )
                  )
                }, x_trans, n) 
                
                stopCluster(cl)
                
                estimativas_boot_delete <- unlist(boot)
                
                alfa_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "alfa"])
                alfa_std_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "alfa_std"])
                alfa_ordinal_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "alfa_ordinal"])
                omega_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "omega"])
                omega_ordinal_boot_delete <- as.vector(estimativas_boot_delete[names(estimativas_boot_delete) == "omega_ordinal"])
                
                return(
                  list(
                    alfa_boot = alfa_boot_delete,
                    alfa_std_boot = alfa_std_boot_delete,
                    alfa_ordinal = alfa_ordinal_boot_delete,
                    omega_boot = omega_boot_delete,
                    omega_ordinal_boot = omega_ordinal_boot_delete
                  )
                )
                
              }
              
              gera_bca <- function(estimativa_original, estimativas_boot, estimativas_boot_delete) {
                #Calculando o z0_hat (bias-correction)
                p1 <- which(estimativas_boot < estimativa_original)
                z0 <- qnorm(length(p1)/n_iter_boot)
                
                #Calculando o a_hat (acceleration)
                a0 <- sum((mean(estimativas_boot_delete)-estimativas_boot_delete)^3)/(6*(sum((mean(estimativas_boot_delete)-estimativas_boot_delete)^2))^(3/2))
                
                #Obtendo os quantis que serão utilizados para a construção do intervalo BCa
                alpha1 <- pnorm(z0 + (z0 + qnorm(0.05/2))/(1-a0*(z0 + qnorm(0.05/2))))
                alpha2 <- pnorm(z0 + (z0 + qnorm(1-(0.05/2)))/(1-a0*(z0 + qnorm(1-(0.05/2)))))
                
                #Obtendo os limites do intervalo
                as.numeric(c(quantile(estimativas_boot, prob = alpha1), quantile(estimativas_boot, prob = alpha2)))
              }
              
              #Gerando as estimativas bootstrap
              lista_estimativas_boot <- gera_est_boot(x_trans = x_trans, n = n, n_iter_boot = n_iter_boot)
              
              #Gerando as estimativas dos coeficientes sem a i-ésima observação
              lista_estimativas_boot_delete <- gera_boot_delete(x_trans = x_trans, n = n)
              
              #Calculando a estimativa e o IC do alfa 
              alfa_est <- as.numeric(alpha(x_trans)$total[1])
              alfa_ci <- gera_bca(
                estimativa_original = alfa_est,
                estimativas_boot = lista_estimativas_boot$alfa_boot,
                estimativas_boot_delete = lista_estimativas_boot_delete$alfa_boot
              )
              alfa_lower <- alfa_ci[1]
              alfa_upper <- alfa_ci[2]
              
              
              #Calculando a estimativa e o IC do alfa std 
              alfa_std_est <- as.numeric(alpha(x_trans)$total[2])
              alfa_std_ci <- gera_bca(
                estimativa_original = alfa_std_est,
                estimativas_boot = lista_estimativas_boot$alfa_std_boot,
                estimativas_boot_delete = lista_estimativas_boot_delete$alfa_std_boot
              )
              alfa_std_lower <- alfa_std_ci[1]
              alfa_std_upper <- alfa_std_ci[2]
              
              
              #Calculando a estimativa e o IC do alfa ordinal 
              alfa_ordinal_est <- as.numeric(alpha(cor_poli)$total[1])
              alfa_ordinal_ci <- gera_bca(
                estimativa_original = alfa_ordinal_est,
                estimativas_boot = lista_estimativas_boot$alfa_ordinal_boot,
                estimativas_boot_delete = lista_estimativas_boot_delete$alfa_ordinal
              )
              alfa_ordinal_lower <- alfa_ordinal_ci[1]
              alfa_ordinal_upper <- alfa_ordinal_ci[2]
              
              
              #Calculando a estimativa e o IC do omega
              omega_est <- (sum(fit_fa$loadings[,1]))^2 / ((sum(fit_fa$loadings[,1]))^2 + sum(fit_fa$uniquenesses))
              omega_ci <- gera_bca(
                estimativa_original = omega_est,
                estimativas_boot = lista_estimativas_boot$omega_boot,
                estimativas_boot_delete = lista_estimativas_boot_delete$omega_boot
              )
              omega_lower <- omega_ci[1]
              omega_upper <- omega_ci[2]
              
              #Calculando a estimativa e o IC do omega ordinal
              omega_ordinal_est <- (sum(fit_fa_poli$loadings[,1]))^2 / ((sum(fit_fa_poli$loadings[,1]))^2 + sum(fit_fa_poli$uniquenesses))
              omega_ordinal_ci <- gera_bca(
                estimativa_original = omega_ordinal_est,
                estimativas_boot = lista_estimativas_boot$omega_ordinal_boot,
                estimativas_boot_delete = lista_estimativas_boot_delete$omega_ordinal_boot
              )
              omega_ordinal_lower <- omega_ordinal_ci[1]
              omega_ordinal_upper <- omega_ordinal_ci[2]
              
              
              #Atribuindo os valores das estimativas à matriz de coeficientes
              matriz_coeficientes[k, "alfa"] <- alfa_est
              matriz_coeficientes[k, "alfa_std"] <- alfa_std_est
              matriz_coeficientes[k, "alfa_ordinal"] <- alfa_ordinal_est
              matriz_coeficientes[k, "omega"] <- omega_est
              matriz_coeficientes[k, "omega_ordinal"] <- omega_ordinal_est
              
              
              #Atribuindo as informações de cobertura e amplitude dos ICs à matriz de ICs
              matriz_ics[k, "alfa_cobre"] <- ifelse(confiabilidade > alfa_lower & confiabilidade < alfa_upper, 1, 0)
              matriz_ics[k, "alfa_amplitude"] <- round(alfa_upper - alfa_lower, 3)
              
              matriz_ics[k, "alfa_std_cobre"] <- ifelse(confiabilidade > alfa_std_lower & confiabilidade < alfa_std_upper, 1, 0)
              matriz_ics[k, "alfa_std_amplitude"] <- round(alfa_std_upper - alfa_std_lower, 3)
              
              matriz_ics[k, "alfa_ordinal_cobre"] <- ifelse(confiabilidade > alfa_ordinal_lower & confiabilidade < alfa_ordinal_upper, 1, 0)
              matriz_ics[k, "alfa_ordinal_amplitude"] <- round(alfa_ordinal_upper - alfa_ordinal_lower, 3)
              
              matriz_ics[k, "omega_cobre"] <- ifelse(confiabilidade > omega_lower & confiabilidade < omega_upper, 1, 0)
              matriz_ics[k, "omega_amplitude"] <- round(omega_upper - omega_lower, 3)
              
              matriz_ics[k, "omega_ordinal_cobre"] <- ifelse(confiabilidade > omega_ordinal_lower & confiabilidade < omega_ordinal_upper, 1, 0)
              matriz_ics[k, "omega_ordinal_amplitude"] <- round(omega_ordinal_upper - omega_ordinal_lower, 3)
              
              if (!file.exists(glue::glue("matrizes_estimativas/matriz_coeficientes_{confiabilidade}_{num_pontos}_{n}_{nivel_variacao_loadings}.csv"))) {
                write.table(cbind(as.data.frame(t(matriz_coeficientes[k, ])), as.data.frame(t(matriz_ics[k, ]))), file = glue::glue("matrizes_estimativas/matriz_coeficientes_{confiabilidade}_{num_pontos}_{n}_{nivel_variacao_loadings}.csv"), row.names = FALSE, sep = ",")
              } else {
                write.table(cbind(as.data.frame(t(matriz_coeficientes[k, ])), as.data.frame(t(matriz_ics[k, ]))), file = glue::glue("matrizes_estimativas/matriz_coeficientes_{confiabilidade}_{num_pontos}_{n}_{nivel_variacao_loadings}.csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
              }
              
              print(k)
            }
            
          }))
        }))
      }))
    }))
  }))
}

system.time(
  resultados_paralelo <- simulacao_confiabilidade(
    M = 1000,
    inicio = 1,
    n_iter_boot = 1000,
    vec_n = c(100, 200, 350, 500),
    vec_num_itens = 5,
    vec_num_pontos = c(3, 5, 7), 
    vec_confiabilidade = c(0.6, 0.8, 0.9),
    var_erros = 0.5,
    paralelo = TRUE,
    vec_nivel_variacao_loadings = 1
  )
)

system.time(
  resultados_congenerico <- simulacao_confiabilidade(
    M = 1000,
    inicio = 1,
    n_iter_boot = 1000,
    vec_n = c(100, 200, 350, 500),
    vec_num_itens = 5,
    vec_num_pontos = c(3, 5, 7), 
    vec_confiabilidade = c(0.6, 0.8, 0.9),
    var_erros = 0.5,
    paralelo = FALSE,
    vec_nivel_variacao_loadings = 10
  )
)


#Gerando os data.frames com as estatísticas de interesse
df_paralelo_completo <- df_congenerico_completo <- data.frame()

for (confiabilidade in c(0.6, 0.8, 0.9)) {
  for (n in c(100, 200, 350, 500)) {
    for (num_pontos in c(3, 5, 7)) {
      matriz_coeficientes_paralelo <- read.csv(glue::glue("matrizes_estimativas/matriz_coeficientes_{confiabilidade}_{num_pontos}_{n}_1.csv"))
      matriz_coeficientes_congenerico <- read.csv(glue::glue("matrizes_estimativas/matriz_coeficientes_{confiabilidade}_{num_pontos}_{n}_10.csv"))
      
      M_paralelo <- nrow(matriz_coeficientes_paralelo)
      M_congenerico <- nrow(matriz_coeficientes_congenerico)
      
      df_estatisticas_paralelo <- data.frame(
        M = rep(M_paralelo, 5),
        confiabilidade = rep(confiabilidade, 5),
        n = rep(n, 5),
        num_pontos = rep(num_pontos, 5),
        coeficiente = c("alfa", "alfa_std", "alfa ordinal", "omega", "omega ordinal"),
        vicio = as.vector(round(colMeans(matriz_coeficientes_paralelo[, c(1:5)], na.rm = TRUE) - confiabilidade, 3)),
        eqm = as.vector(round(apply(matriz_coeficientes_paralelo[, c(1:5)], 2, function(coluna) sum((coluna - confiabilidade)^2, na.rm = TRUE)/M_paralelo), 3)),
        ic_taxa_de_cobertura = as.vector(round(colMeans(matriz_coeficientes_paralelo[, c(6, 8, 10, 12, 14)], na.rm = TRUE), 3)),
        ic_amplitude = as.vector(round(colMeans(matriz_coeficientes_paralelo[, c(7, 9, 11, 13, 15)], na.rm = TRUE), 3))
      )
      
      df_estatisticas_congenerico <- data.frame(
        M = rep(M_congenerico, 5),
        confiabilidade = rep(confiabilidade, 5),
        n = rep(n, 5),
        num_pontos = rep(num_pontos, 5),
        coeficiente = c("alfa", "alfa_std", "alfa ordinal", "omega", "omega ordinal"),
        vicio = as.vector(round(colMeans(matriz_coeficientes_congenerico[, c(1:5)], na.rm = TRUE) - confiabilidade, 3)),
        eqm = as.vector(round(apply(matriz_coeficientes_congenerico[, c(1:5)], 2, function(coluna) sum((coluna - confiabilidade)^2, na.rm = TRUE)/M_congenerico), 3)),
        ic_taxa_de_cobertura = as.vector(round(colMeans(matriz_coeficientes_congenerico[, c(6, 8, 10, 12, 14)], na.rm = TRUE), 3)),
        ic_amplitude = as.vector(round(colMeans(matriz_coeficientes_congenerico[, c(7, 9, 11, 13, 15)], na.rm = TRUE), 3))
      )
      
      df_paralelo_completo <- rbind(df_paralelo_completo, df_estatisticas_paralelo)
      df_congenerico_completo <- rbind(df_congenerico_completo, df_estatisticas_congenerico)
      
    }
  }
}

df_paralelo_completo <- df_paralelo_completo |>
  dplyr::filter(coeficiente != "alfa") |>
  dplyr::mutate(
    pc_amplitude = glue::glue("{ic_taxa_de_cobertura} ({ic_amplitude})"),
    .keep = "unused"
  )

df_congenerico_completo <- df_congenerico_completo |>
  dplyr::filter(coeficiente != "alfa") |>
  dplyr::mutate(
    pc_amplitude = glue::glue("{ic_taxa_de_cobertura} ({ic_amplitude})"),
    .keep = "unused"
  )

write.csv(df_paralelo_completo, "resultados_simulacao_confiabilidade_paralelo.csv", row.names = FALSE)
write.csv(df_congenerico_completo, "resultados_simulacao_confiabilidade_congenerico.csv", row.names = FALSE)
