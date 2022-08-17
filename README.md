# beta_bimodal_inflacionado_zero
# Introdução

Frequentemente observamos dados de proporção com presença de zeros e/ou uns. Para representar as altas concentrações de observações nos extremos 0 e 1, Martinez (2008) propôs alternativas para modelagem de dados em que a variável resposta pode variar nos intervalos [0,1], [0,1) e (0,1], os modelos assumem que os dados são originários de uma distribuição de mistura entre uma beta e uma Bernoulli, degenerada em zero e degenerada em 1, respectivamente. 

Neste sentido, propomos um modelo de regressão em que as variáveis respostas possuem inflações em zero e bimodalidade no intervalo $(0,1)$, a partir de uma distribuição de mistura de duas distribuições beta e uma degenerada em zero. Desenvolvemos a distribuição beta bimodal inflacionada em zero, definimos estimadores via máxima verossimilhança e construímos o modelo de regressão para este modelo probabilístico, realizamos estudo de simulação para avaliar o desempenho dos estimadores de máxima verossimilhança (obtidos por meio do algoritmo EM inicializado com método K-médias) dos parâmetros do modelo de regressão beta bimodal inflacionado em zero.


# Modelo de regressão beta bimodal inflacionada em zero

Considere, inicialmente, uma variável aleatória $W$ proveniente de uma mistura de duas distribuições beta, a distribuição desta variável aleatória possui função densidade de probabilidade da forma

$$bb(w;\pi,\mu_1,\phi_1,\mu_2,\phi_2)= \pi f(w;\mu_1,\phi_1) + (1-\pi)f(w;\mu_2,\phi_2)$$

em que $0 < w$,$\pi < 1$, $f(w;\mu_1,\phi_1)$ (Parametrização alternativa da distribuição beta em
termos da média $\mu_1$, e parâmetro de precisão $\phi_1$) e $f(w;\mu_2,\phi_2)$ (densidades componentes) são funções densidade de probabilidade da distribuição beta referentes às duas subpopulações misturadas aleatoriamente com proporções $\pi$ e $(1-\pi)$, respectivamente, estes são chamados pesos das componentes. Assim dizemos que $W$ segue distribuição beta bimodal com parâmetros $\pi,\mu_1,\phi_1,\mu_2$ e $\phi_2$.

Neste momento, considere $Y$ uma variável aleatória, que assume valores no intervalo $[0,1)$, originária da mistura de uma distribuição degenerada em zero e uma distribuição beta bimodal \eqref{dens_beta_bimodal}, desta forma sua função densidade de probabilidade é dada por

$$bbz(y;\alpha,\pi,\mu_1,\phi_1,\mu_2,\phi_2) = 
\alpha I(y=0) + (1-\alpha)bb(y;\pi,\mu_1,\phi_1,\mu_2,\phi_2) I(y \in (0,1))$$

com $\alpha,\pi,\mu_1, \mu_2 \in (0,1)$ e $\phi_1,\phi_2 > 0$. A função $bb(y;\pi,\mu_1,\phi_1,\mu_2)$ refere-se a densidade beta bimodal \eqref{dens_beta_bimodal}. Observe que $\alpha=P(Y=0)$, configura a probabilidade de se observar o valor zero, e com probabilidade $(1-\alpha)$ a variável aleatória tem origem beta bimodal. Então $Y \sim BBZ(\alpha,\pi,\mu_1,\phi_1,\mu_2$ e $\phi_2)$, isto é, $Y$ é uma variável aleatória com distribuição beta bimodal inflacionada em zero (BBZ) com parâmetros $\alpha,\pi,\mu_1,\phi_1,\mu_2$ e $\phi_2$. 


A esperança e a variância de $Y \sim BBZ(\alpha,\pi,\mu_1,\phi_1,\mu_2$ e $\phi_2)$ são dadas, respectivamente, por 
$$E[Y] = \mu_{1} \pi(1-\alpha) + \mu_{2} (1-\pi) (1-\alpha),$$
$$Var[Y] = \left(\frac{\mu_1(1-\mu_1)}{\phi_1+1}+\mu_1^2\right)\pi(1-\alpha)+\left(\frac{\mu_2(1-\mu_2)}{\phi_2+1} +\mu_2^2\right)(1-\pi) (1-\alpha)-\mu^2,$$

em que $\mu=E[Y]$.

Sejam $Y_1,...,Y_n$ variáveis aleatórias independentes, em que cada $Y_i$, $i=1,2,...,n$, possui função densidade de probabilidade beta bimodal inflacionada em zero da forma \eqref{dens_beta_bimodal_zero}, com parâmetros $\alpha_i,\pi_i,\mu_{1i},\phi_{1i},\mu_{2i}$ $\phi_{2i}$, respectivamente. Os modelos de regressão beta bimodal inflacionados em zero (RBBZ) são definidos pelos seguintes componentes sistemáticos:
$$g_{0}({\alpha}) = {\eta}_{0} = X_{0} {\beta}_{0},$$
$$g_{1}({\mu}_{1}) = {\eta}_{1} = X_{1} {\beta}_{1},$$

$$g_{2}({\mu}_{2}) = {\eta}_{2} = X_{2} {\beta}_{2},$$

$$g_{3}({\pi}) = {\eta}_{3} = X_{3} {\beta}_{3},$$

em que ${\alpha},{\mu}_1,{\mu}_2,{\pi}$ e ${\eta}_{k}$, $k=0,1,2$ e $3$, são vetores de tamanho $n$, ${\beta}_{k}^{T}=(\beta_{k1},\beta_{k2},...,\beta_{kd_{k}})$ é um vetor de tamanho $d_k$, $X_{k}$ é uma matriz de valores conhecidos da ordem $n \times d_{k}$. As funções $g_{k}(\cdot)$, são denominadas funções de ligação, relacionam os vetores de parâmetros ${\alpha},{\mu}_1,{\mu}_2$ e ${\pi}$ às variáveis explanatórias em $X_{0},X_{1},X_{2}$ e $X_{3}$, respectivamente. As funções de ligação são conhecidas e devem ser estritamente monótonas e $g_{k}(\cdot):(0,1) \longrightarrow \mathbb{R}$. 

O ajuste do modelo RBBZ é realizado pela estimação de máxima verossimilhança. 
