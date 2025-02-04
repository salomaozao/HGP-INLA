# O que precisa ser feito:

- [ OK ] MELHORAR ESTIMAÇÃO

- [ ] Ajustar os dados do CAR para regressão Poisson ( dados discretos ) (?)
    -> Os dados da variável dependente SMR são uma taxa, mas a variável dependente de um processo poisson seria contagem?
    -> A variavel resposta deve ser respiratorydata$observed e deve ser introduzido um offset no modelo. no inla, vc bota o offset com E = <offset>
    - https://becarioprecario.bitbucket.io/inla-gitbook/ch-INLA.html#sec:GLM

- [ ] Separar o block apenas para o rgeneric (?)
    -> Separar a blocagem em uma outra função?






NOTAS:
- Apenas 1 hdist é possível? Se não, melhorar nomeclatura

=========================
O que quero fazer:
1. Arrumar estimação.

dados reais:                    0.1       1         1/3
caso pointdata antigo:          0.1154132 1.4536849 0.5149386  \\ 0.09881827 1.20887559 0.33474234
caso pointdata antigo ajustado: 0.1015100 1.7776850 0.4391706
caso atual:                     0.2732637 0.4066729 0.3088849





2. Refactor
    1. runblock -> DRY, calcular matriz D