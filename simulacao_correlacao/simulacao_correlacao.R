library(psych)
library(fMultivar)
library(parallel)

vec_n <- c(100, 350, 500)
vec_rho <- c(0, 0.2, 0.5, 0.8)
vec_num_pontos <- c(3, 5, 7)
inicio <- 1
M <- 1000
n.iter_boot <- 1000

df_estatisticas <- do.call(rbind, lapply(vec_rho, function(rho) {
  do.call(rbind, lapply(vec_n, function(n) {
    df_iteracao <- do.call(rbind, lapply(vec_num_pontos, function(num_pontos) {
      Cor_poli <- Cor_pearson <- data.frame(
        "estimativa" = rep(0, M), 
        "ic_normal_cobre" = rep(0, M),
        "amplitude_ic_normal" = rep(0, M),
        "ic_boot_cobre" = rep(0, M),
        "amplitude_ic_boot" = rep(0, M)
      )
      
      for (k in inicio:M) {
        set.seed(429 + k)
        
        X <- rnorm2d(n, rho = rho)
        X_trans <- matrix(1, nrow(X), ncol(X))
        
        for (i in 1:nrow(X)){
          for (j in 1:ncol(X)){
            quantil_base <- 1/num_pontos
            for (ponto in 1:(num_pontos - 1)) {
              if ((X[i, j] > qnorm(quantil_base*ponto)) & (X[i, j] <= qnorm(quantil_base*(ponto + 1)))) {
                X_trans[i, j] <- ponto + 1
              }
            }
          }
        }
        
        #Criando a função que será utilizada na obtenção do IC das correlações por bootstrap
        gera_est_boot <- function(X_trans, n, n.iter_boot) {
          cl <- makeCluster(getOption("cl.cores", 2))
          clusterEvalQ(cl, { library(psych) })
          
          boot <- parLapply(cl, 1:n.iter_boot, function(i, X_trans, n, n.iter_boot) {
            #Criando a reamostra bootstrap
            data_boot <- X_trans[sample.int(n, replace = TRUE), ]
            
            #Calculando as correlações policórica e de Pearson
            cor_poli <- polychoric(data_boot)$rho[1, 2]
            cor_pearson <- cor(data_boot)[1, 2]
            
            return(list(cor_poli = cor_poli, cor_pearson = cor_pearson))
          }, X_trans, n, n.iter_boot) 
          
          stopCluster(cl)
          
          cor_boot <- unlist(boot)
          
          cor_pearson_boot <- as.vector(cor_boot[names(cor_boot) == "cor_pearson"])
          cor_poli_boot <- as.vector(cor_boot[names(cor_boot) == "cor_poli"])
          
          return(list(cor_pearson_boot = cor_pearson_boot, cor_poli_boot = cor_poli_boot))
        }
      
        gera_bca <- function(X_trans, estimativa_original, estimativas_boot, n.iter_boot, poly) {
          #Calculando o z0_hat (bias-correction)
          p1 <- which(estimativas_boot < estimativa_original)
          z0 <- qnorm(length(p1)/n.iter_boot)
          
          #Calculando o a_hat (acceleration)
          cl <- makeCluster(getOption("cl.cores", 2))
          clusterEvalQ(cl, { library(psych) })
          cor_boot_delete <- parLapply(cl, 1:n, function(i, X_trans, n.iter_boot, poly) { 
            if (poly == TRUE) {
              cor <- polychoric(X_trans[-i, ])$rho[1, 2]
            } else {
              cor <- cor(X_trans[-i, ])[1, 2]
            }
            cor
          }, X_trans, n.iter_boot, poly) 
          
          stopCluster(cl)
          
          cor_boot_delete <- unlist(cor_boot_delete)
          
          a0 <- sum((mean(cor_boot_delete)-cor_boot_delete)^3)/(6*(sum((mean(cor_boot_delete)-cor_boot_delete)^2))^(3/2))
          
          #Obtendo os quantis que serão utilizados para a construção do intervalo BCa
          alpha1 <- pnorm(z0 + (z0 + qnorm(0.05/2))/(1-a0*(z0 + qnorm(0.05/2))))
          alpha2 <- pnorm(z0 + (z0 + qnorm(1-(0.05/2)))/(1-a0*(z0 + qnorm(1-(0.05/2)))))
          
          #Obtendo os limites do intervalo
          as.numeric(c(quantile(estimativas_boot, prob = alpha1), quantile(estimativas_boot, prob = alpha2)))
        }
        
        #Gerando as estimativas bootstrap
        lista_estimativas_boot <- gera_est_boot(X_trans = X_trans, n = n, n.iter_boot = n.iter_boot)
        
        #Calculando as estimativas das correlações policórica e de Pearson
        M1 <- polychoric(x = X_trans)
        M2 <- cor(x = X_trans)
        
        #Calculando o IC da correlação policórica
        ic_boot_poli <- gera_bca(
          X_trans = X_trans, 
          estimativa_original = M1$rho[1, 2],
          estimativas_boot = lista_estimativas_boot$cor_poli_boot,
          n.iter_boot = n.iter_boot,
          poly = TRUE
        )
        
        #Calculando o IC da correlação de Pearson
        ic_boot_pearson <- gera_bca(
          X_trans = X_trans, 
          estimativa_original = M2[1, 2],
          estimativas_boot = lista_estimativas_boot$cor_pearson_boot,
          n.iter_boot = n.iter_boot,
          poly = FALSE
        )
        
        #Calculando os ICs pela teoria normal (não utilizamos)
        ic_normal_poli <- cor.ci(M1$rho, poly = TRUE, n = n)
        ic_normal_pearson <- cor.ci(M2, method = "pearson", n = n)
        
        Cor_poli[k, 1] <- M1$rho[1, 2]
        Cor_poli[k, 2] <- ifelse(ic_normal_poli$ci$lower <= rho & ic_normal_poli$ci$upper >= rho, 1, 0)
        Cor_poli[k, 3] <- ic_normal_poli$ci$upper - ic_normal_poli$ci$lower
        Cor_poli[k, 4] <- ifelse(ic_boot_poli[1] <= rho & ic_boot_poli[2] >= rho, 1, 0)
        Cor_poli[k, 5] <- ic_boot_poli[2] - ic_boot_poli[1]
        
        Cor_pearson[k, 1] <- M2[1, 2]
        Cor_pearson[k, 2] <- ifelse(ic_normal_pearson$ci$lower <= rho & ic_normal_pearson$ci$upper >= rho, 1, 0)
        Cor_pearson[k, 3] <- ic_normal_pearson$ci$upper - ic_normal_pearson$ci$lower
        Cor_pearson[k, 4] <- ifelse(ic_boot_pearson[1] <= rho & ic_boot_pearson[2] >= rho, 1, 0)
        Cor_pearson[k, 5] <- ic_boot_pearson[2] - ic_boot_pearson[1]
        
        if (!file.exists(glue::glue("matrizes_estimativas/matriz_cor_pearson_{rho}_{num_pontos}_{n}.csv"))) {
          write.table(Cor_pearson[k, ], file = glue::glue("matrizes_estimativas/matriz_cor_pearson_{rho}_{num_pontos}_{n}.csv"), row.names = FALSE, sep = ",")
        } else {
          write.table(Cor_pearson[k, ], file = glue::glue("matrizes_estimativas/matriz_cor_pearson_{rho}_{num_pontos}_{n}.csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
        }
        
        if (!file.exists(glue::glue("matrizes_estimativas/matriz_cor_poli_{rho}_{num_pontos}_{n}.csv"))) {
          write.table(Cor_poli[k, ], file = glue::glue("matrizes_estimativas/matriz_cor_poli_{rho}_{num_pontos}_{n}.csv"), row.names = FALSE, sep = ",")
        } else {
          write.table(Cor_poli[k, ], file = glue::glue("matrizes_estimativas/matriz_cor_poli_{rho}_{num_pontos}_{n}.csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
        }
        
        print(k)
      }
      
    }))
    
  }))
})) 

df_completo <- data.frame()

for (rho in c(0, 0.2, 0.5, 0.8)) {
  for (n in c(100, 350, 500)) {
    for (pontos in c(3, 5, 7)) {
      Cor_poli <- read.csv(glue::glue("matrizes_estimativas/matriz_cor_poli_{rho}_{pontos}_{n}.csv"))
      Cor_pearson <- read.csv(glue::glue("matrizes_estimativas/matriz_cor_pearson_{rho}_{pontos}_{n}.csv"))
      
      M <- nrow(Cor_poli)
      
      df_estatisticas_iteracao <- data.frame(
        "rho" = c(rho, rho),
        "n" = c(n, n),
        "num_pontos" = c(pontos, pontos),
        "cor" = c("policorica", "pearson"),
        "vicio" = round(c(mean(Cor_poli$estimativa) - rho, mean(Cor_pearson$estimativa - rho)), 3),
        "eqm" = round(c(sum((Cor_poli$estimativa - rho)^2)/M, sum((Cor_pearson$estimativa - rho)^2)/M), 3),
        "ic_boot_tx_cobertura_amplitude" = paste0(
          round(c(sum(Cor_poli$ic_boot_cobre)/M, sum(Cor_pearson$ic_boot_cobre)/M), 3),
          " (", round(c(mean(Cor_poli$amplitude_ic_boot), mean(Cor_pearson$amplitude_ic_boot)), 3), ")"
        )
      )
      
      df_completo <- rbind(df_completo, df_estatisticas_iteracao)
    }
  }
}


write.csv(df_completo, "resultados_simulacao_correlacao.csv", row.names = FALSE)

