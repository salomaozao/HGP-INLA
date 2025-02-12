# TODO
- [ OK! ] Arrumar estimação

- Refactor do código
    - [ OK ] Pegar centróides para sortloc
    - [ OK ] Trabalhar com `sf` no `locs` para suporte de dados de área
    - [ OK ] Separar a blocagem
    - [ ] Separar geração do modelo 
    - [ ] Juntar todos pontos acima em uma função `createHGP(sf, n.partition)`
    
    - [ ] Separar a parte da geração dos dados e ajuste do modelo, implementar pensando na geração dos dados de área
    


- Implementações
    - [ ] Implementar regressão Poisson
    - [ ] Implementação de dados de área
    - [ ] Ajustar os dados do CAR para regressão Poisson ( dados discretos ) (?)
        -> Os dados da variável dependente SMR são uma taxa, mas a variável dependente de um processo poisson seria contagem?
        -> A variavel resposta deve ser respiratorydata$observed e deve ser introduzido um offset no modelo. no inla, vc bota o offset com E = <offset>
        - https://becarioprecario.bitbucket.io/inla-gitbook/ch-INLA.html#sec:GLM


- Ajustes:
    - [ ] exigir sf no hdist_sf
    - [ ] usar dist euclidiana na blocagem

- ideias:
    - [ ] mudar hdist para makeDistMatrix() e ter opção haussdorff ou euclidiana

-shapefile
# Problemas no Código:
    passar coords.D no NNGP erra a estimação (porque?)


# Commits:

### 10/2/25 Implementação do sf no ajuste do modelo
- Implementado centróides na criação dos blocos
- A matriz ordenada `orderedLoc` agora é um objeto `sf`, que foi passada para fazer a matriz de distância usando a função sf `hdist_sf`


### 3/2/25 SETBACK: hdist funcionando (sem dados de área)
- Baixei a versão `"Ajustes do uso dos parâmetros do ajuste (Funcionando direito agora)"`
- implementado o `hdist`
- Ajustadas as funções dos parâmetros no INLA (`[[blockNNGPfunctionIRREGULAR.R:188]]`)
- Arrumei o runblock para isso
 

