# TODO

-   [ OK! ] Arrumar estimação

-   Refactor do código

    -   [ OK ] Pegar centróides para sortloc
    -   [ OK ] Trabalhar com `sf` no `locs` para suporte de dados de área
    -   [ OK ] Separar geração do modelo `createHGP(sf, n.partition)`
        -   [ OK ] Separar a blocagem
        -   [ OK ] Separar matriz de precisão
        -   [ OK ] fazer uma função generatemodel()

-   Implementações

    -   [ ] Implementação de dados de área
        -   [ OK ] rodar com dados do CARbayesdata
        -   [ OK ] Alterar a variável resposta para os dados de área -> A variavel resposta deve ser respiratorydata$observed e deve ser introduzido um offset no modelo. no inla, vc bota o offset com E = <offset>
            -   https://becarioprecario.bitbucket.io/inla-gitbook/ch-INLA.html#sec:GLM

            -   The goals of the analyses are to provide smooth estimates of the SMR, adjusted for the percentage of income deprivation. In this application, we have:

                -   Y(s_i) is the observed number of admissions in the i-th IZ (s_i) `$observed`
                -   x_i as the percentage of people who are defined to be income deprived `$incomedep`
                -   E_i is the expected number of admissions `$expected`

            -   The basic model for area-level count data is given by $y_i∼Poisson(μ_i)$

                -   $μi​=Ei​+ηi​$
                -   $η_i$ is the linear predictor: $$β0+β1xi+f(i) = β0+β1 * text{incomedep}_i+f(idx, model = blockNNGP.model)$$



            -   Configurar o offset no INLA

-   Ajustes:

    -   [ OK ] exigir sf no hdist_sf
    -   [ OK ] usar dist euclidiana na blocagem

# Commits:

### Implementado função generateHGP para gerar o modelo

-   Criadas as funções `get_precMatrixData()`, `get_blocksdata()` e uma função que junta as duas, `get_HGPdata()`

### 10/2/25 Implementação do sf no ajuste do modelo

-   Implementado centróides na criação dos blocos
-   A matriz ordenada `orderedLoc` agora é um objeto `sf`, que foi passada para fazer a matriz de distância usando a função sf `hdist_sf`

### 3/2/25 SETBACK: hdist funcionando (sem dados de área)

-   Baixei a versão `"Ajustes do uso dos parâmetros do ajuste (Funcionando direito agora)"`
-   implementado o `hdist`
-   Ajustadas as funções dos parâmetros no INLA (`[[blockNNGPfunctionIRREGULAR.R:188]]`)
-   Arrumei o runblock para isso