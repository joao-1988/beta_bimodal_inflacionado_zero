#--------- Funções: densidade, acumulada, quantil e aleatoria para BBZ ------------#

dbbz <- function (x, alpha=0.2, Pi=0.3,  mu_1=0.2, phi_1=10, mu_2=0.7, phi_2=10) 
{
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(x < 0) | any(x >= 1)) 
    stop(paste("x must be 0<=x<1, i.e. [0 to 1)", "\n", ""))
  a1 <- mu_1 * phi_1
  b1 <- (1 - mu_1) * phi_1
  a2 <- mu_2 * phi_2
  b2 <- (1 - mu_2) * phi_2
  d <- rep(0, length(x)) 
  d <- ifelse((x > 0), (1-alpha)*(Pi*dbeta(x, shape1 = a1, shape2 = b1) + 
                                    (1-Pi)*dbeta(x, shape1 = a2, shape2 = b2)), 0)
  d <- ifelse((x == 0), alpha, d)
  d
}

pbbz <- function (q, alpha, Pi,  mu_1, phi_1, mu_2, phi_2) 
{
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(q < 0) | any(q >= 1)) 
    stop(paste("y must be 0<=y<1, i.e. [0 to 1)", "\n", ""))
  a1 <- mu_1 * phi_1
  b1 <- (1 - mu_1) * phi_1
  a2 <- mu_2 * phi_2
  b2 <- (1 - mu_2) * phi_2
  p <- ifelse((q > 0), alpha + (1-alpha)*(Pi*pbeta(q, shape1 = a1, shape2 = b1) + (1-Pi)*pbeta(q, shape1 = a2, shape2 = b2)),0)
  p <- ifelse((q == 0), alpha, p)
  p
}

#Quantil resolvendo a inversa da acumulada da mistura de duas beta
# proporciona um quantil aproximado
# a funcao encontra apenas uma raiz, ou seja, isto provoca um problema ao
# utilizar de forma vetorial, o que seria ideal para encontrar varios quantis
# e isto obriga a rbbz a utilizar um "for" ao inves de operacao com vetor
# e assim pode ser mais demorado no R.

qbbz <- function (p, alpha, Pi,  mu_1, phi_1, mu_2, phi_2) 
{
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(p <= 0) | any(p >= 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  a1 <- mu_1 * phi_1
  b1 <- (1 - mu_1) * phi_1
  a2 <- mu_2 * phi_2
  b2 <- (1 - mu_2) * phi_2
  q <- base::sapply(p,function(p){
    f = function(x) { alpha + (1-alpha)*( Pi*pbeta(x, shape1 = a1, shape2 = b1) + (1-Pi)*pbeta(x, shape1 = a2, shape2 = b2) ) - p}
    suppressWarnings(q <- ifelse((p <= (alpha)), 0, ( uniroot(f,c(0.0001,0.9999))$root )))
    q
  })
  q
}

rbbz <- function (n, alpha, Pi,  mu_1, phi_1, mu_2, phi_2) 
{
  if (any(alpha <= 0) | any(alpha >= 1))
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  a1 <- mu_1 * phi_1
  b1 <- (1 - mu_1) * phi_1
  a2 <- mu_2 * phi_2
  b2 <- (1 - mu_2) * phi_2
  n <- ceiling(n)
  y <- vector(length = n)
  for(i in 1:n){
    y[i] <- rbinom(1,1,1-alpha) 
    if( y[i] == 1 ){ y[i] <- rbinom( 1, 1, Pi) 
    ifelse( y[i] == 1 , 
            y[i] <- rbeta(1, shape1 = a1, shape2 = b1) , 
            y[i] <- rbeta(1, shape1 = a2, shape2 = b2) )
    }
  }
  y
}



rbbz_p <- function (n, alpha, Pi,  mu_1, phi_1, mu_2, phi_2) 
{
  if (any(alpha <= 0) | any(alpha >= 1))
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  a1 <- mu_1 * phi_1
  b1 <- (1 - mu_1) * phi_1
  a2 <- mu_2 * phi_2
  b2 <- (1 - mu_2) * phi_2
  n <- ceiling(n)
  y <- matrix(nrow = n, ncol = 2)
  for(i in 1:n){
    y[i,] <- c( rbinom(1,1,1-alpha), 0 ) 
    if( y[i,1] == 1 ){ y[i,1] <- rbinom( 1, 1, Pi) 
    ifelse( y[i,1] == 1 , 
            y[i,] <- c( rbeta(1, shape1 = a1, shape2 = b1) , 1 ), 
            y[i,] <- c( rbeta(1, shape1 = a2, shape2 = b2) , 2 ) )
    }
  }
  y
}

rbbz_aprox <- function (n, alpha, Pi,  mu_1, phi_1, mu_2, phi_2) 
{
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(Pi <= 0) | any(Pi >= 1)) 
    stop(paste("Pi must be between 0 and 1", "\n", ""))
  if (any(mu_1 <= 0) | any(mu_1 >= 1)) 
    stop(paste("mu_1 must be between 0 and 1", "\n", ""))
  if (any(mu_2 <= 0) | any(mu_2 >= 1)) 
    stop(paste("mu_2 must be between 0 and 1", "\n", ""))
  if (any(phi_1 <= 0)) 
    stop(paste("phi_1 must be greated than 0", "\n", ""))
  if (any(phi_2 <= 0)) 
    stop(paste("phi_2 must be greated than 0", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qbbz(p, alpha, Pi, mu_1, phi_1, mu_2, phi_2) 
  r
}


#------------------------------Funcao Esperança--------------------------------------#

esp <- function(alpha,Pi,mu_1,mu_2){
  mu_1*Pi*(1-alpha) + mu_2*(1-Pi)*(1-alpha)
}

#---------------------------Funcao Variancia-----------------------------------------#

varian <- function(alpha,Pi,mu_1,phi_1,mu_2,phi_2){
  (( mu_1*(1-mu_1) )/(phi_1+1) + mu_1^2)*Pi*(1-alpha) +
    (( mu_2*(1-mu_2) )/(phi_2+1) + mu_2^2)*(1-Pi)*(1-alpha) - 
    (mu_1*Pi*(1-alpha) + mu_2*(1-Pi)*(1-alpha))^2
}


#------------------Ajusta um modelo BBZ - Modelo BBZ-----------------------------#
# Não contempla covariáveis

mbbz <- function(y, grupo = NULL){
  
  ## Ajuste do modelo relacionado ao parametro alpha
  
  # y_ast (y* na dissertacao) recebe 1 se y for zero e 0 se y for diferente de zero
  y_ast <- as.numeric(y == 0)
  
  mod_alpha <- stats::glm(y_ast ~ 1,family = binomial(link = "logit"))
  
  ## Ajuste do modelo relacionado a mistura de duas beta
  
  # y_01 recebe apenas as observacoes y tal que y em (0,1)
  y_01 <- y[y != 0]
  
  mod_2beta <- betareg::betamix(formula = y_01 ~ 1, k=2, cluster = grupo, FLXcontrol = list(minprior=0))
  
  list(mod_alpha,mod_2beta)
}

# Exemplo usando kmeans
# n = 1000
# alpha=0.25
# Pi=0.4
# mu_1 = 0.3
# phi_1 = 10
# mu_2 = 0.9
# phi_2 = 12
# y <- rbbz_p( n, alpha, Pi, mu_1, phi_1, mu_2, phi_2)
# mod <- mbbz(y[,1],grupo = kmeans(y[y[,1] != 0,1],2)$cluster)
# mod[[1]]
# mod[[2]]


#------------------Ajusta um modelo regressao BBZ - regressao BBZ ----------------------#

# o argumento dados deve receber um data.frame em que a primeira coluna é a variavel resposta

regbbz <- function(formula_alpha, formula_2beta, 
                   link_alpha="logit", link_2beta="logit",
                   cluster = NULL ,FLXcontrol = list(minprior=0), 
                   dados){
  ## Ajuste do modelo relacionado ao parametro alpha
  dados_alpha = data.frame(y = as.numeric(dados[,1] == 0), dados[,-1]  )
  mod_alpha <- stats::glm(formula = formula_alpha, family = binomial(link = link_alpha), data = dados_alpha)
  
  ## Ajuste do modelo relacionado a mistura de duas beta
  dados_2beta = subset(dados, dados$y != 0)
  mod_2beta <- betareg::betamix(formula = formula_2beta, 
                                k=2, 
                                cluster = cluster, 
                                link = link_2beta,
                                FLXcontrol = FLXcontrol, data = dados_2beta)
  return (list(mod_alpha,mod_2beta))
}  

# Exemplo uso
#fit <- regbbz(formula_alpha = y ~ 1 + x_0, 
#              formula_2beta = y ~ 1 + x_1 + x_2 | 1 | 1 + x_3,
#              cluster = kmeans(dados$y[dados$y != 0], 2)$cluster,
#              dados = dados)
# Duas formas de selecionar os grupos com kmeans
# kmeans(dados$y[dados$y != 0], 2)$cluster 
# kmeans(dados[dados[,1] != 0,1], 2)$cluster neste a variavel resposta so precisa estar na primeira coluna
# Funcoes aplicaveis no modelo
#summary(fit[[1]])
#summary(fit[[2]])
#coef(fit[[1]])
#coef(fit[[2]])
#coef(fit[[2]],which = "concomitant")
#logLik(fit[[1]])
#logLik(fit[[2]])



fitted.regbbz <- function(mod,dados){
  
  k0 <- 0  
  eta0 <- mod[[1]]$coefficients[1]
  for(i in 2:length(mod[[1]]$coefficients)){
    for(j in 1:ncol(dados)){
      if(names(mod[[1]]$coefficients)[i] == names(dados)[j]){
        eta0 <- eta0 + mod[[1]]$coefficients[i]*dados[,j]
        k0 <- k0 + 1
      }
    }
  }
  if( (length(mod[[1]]$coefficients) - 1) != k0){print("Faltou/Sobrou alguma variável para calcular o preditor de alpha")}
  
  k3 <- 0
  eta3 <- coef(mod[[2]],which = "concomitant")[2,1]
  for(i in 2:ncol(coef(mod[[2]],which = "concomitant")) ){
    for(j in 1:ncol(dados)){
      if(names(coef(mod[[2]],which = "concomitant")[2,])[i] == names(dados)[j]){
        eta3 <- eta3 + coef(mod[[2]],which = "concomitant")[2,i]*dados[,j]
        k3 <- k3 + 1
      }
    }
  }
  if( (ncol(coef(mod[[2]],which = "concomitant")) - 1) != k3){print("Faltou/Sobrou alguma variável para calcular o preditor de Pi")}
  
  k1 <- 0
  eta1 <- coef(mod[[2]])[1,1]
  for(i in 2:ncol(coef(mod[[2]])) ){
    for(j in 1:ncol(dados)){
      if(names(coef(mod[[2]])[1,])[i] == names(dados)[j]){
        eta1 <- eta1 + coef(mod[[2]])[1,i]*dados[,j]
        k1 <- k1 + 1 
      }
    }
  }
  if( (ncol(coef(mod[[2]])) - 2) != k1){print("Faltou/Sobrou alguma variável para calcular o preditor de mu_1")}
  
  k2 <- 0
  eta2 <- coef(mod[[2]])[2,1]
  for(i in 2:ncol(coef(mod[[2]])) ){
    for(j in 1:ncol(dados)){
      if(names(coef(mod[[2]])[2,])[i] == names(dados)[j]){
        eta2 <- eta2 + coef(mod[[2]])[2,i]*dados[,j]
        k2 <- k2 +1
      }
    }
  }
  if( (ncol(coef(mod[[2]])) - 2) != k2){print("Faltou/Sobrou alguma variável para calcular o preditor de mu_2")}
  
  cbind(alpha = inv_logit(eta0),
        Pi = 1 - inv_logit(eta3), # a saida do modelo retorna ligacao com o segundo componente densidade, entao fazemos 1 - ele
        mu_1 = inv_logit(eta1),
        mu_2 = inv_logit(eta2))
  
}



# Função Liklihood com estimativa dos parametros


logLikelihood <- function(modelo, dados){
  parametros_ajustados <- fitted.regbbz(modelo,dados)
  dados_est <- data.frame(y = dados$y, 
                          alpha = parametros_ajustados[,1],
                          Pi = parametros_ajustados[,2],
                          mu_1 = parametros_ajustados[,3],
                          phi_1 = rep( exp(coef(modelo[[2]])[1, ncol(coef(modelo[[2]])) ]), nrow(dados) ) ,
                          mu_2 = parametros_ajustados[,4],
                          phi_2 = rep( exp(coef(modelo[[2]])[2, ncol(coef(modelo[[2]])) ]), nrow(dados)  ) )
  
  sum(apply(dados_est, MARGIN = 1, FUN = function(x) {
    ifelse(x[1] == 0, log(x[2]) , log(1-x[2]) + 
             log(x[3]*dbeta(x[1], shape1 = x[4] * x[5] , shape2 = (1 - x[4]) * x[5] ) + 
                   (1-x[3])*dbeta(x[1], shape1 = x[6] * x[7], shape2 = (1 - x[6]) * x[7] ))
    )  
  }) 
  )
  
}

