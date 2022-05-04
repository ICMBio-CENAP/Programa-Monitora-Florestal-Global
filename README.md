# Programa-Monitora-Florestal-Global

Analise de dados de mamíferos e aves ("Mastoaves") do protocolo florestal global do [**Programa Monitora-ICMBio**](https://www.icmbio.gov.br/portal/monitoramento-2016/programas-de-monitoramento-da-biodiversidade-em-ucs)

Prepara os dados do Monitora-ICMBio no formato exigido pelo pacote *rlpi* desenvolvido pela [**Zoological Society of London**](https://github.com/Zoological-Society-of-London/rlpi);

Cria subsets dos dados com UC e grupo taxonômico de interesse (e.g., mamíferos, aves);

Calcula o **LPI (Living Planet Index)** a partir de dados do protocolo florestal básico (mamíferos e aves)

## ESTRUTURA DO DIRETÓRIO

``` bash
├───01_dadosDeEntrada
├───02_script
│   ├───01_limpezaDeDados
│   ├───02_dataPrepLPI
│   ├───03_calculoDeLPI
│   │   └───lpi_temp
│   └───funcoes
└───03_dadosDeSaida
    ├───dados
    └───plots
```

## WORKFLOW

[01_dadosDeEntrada](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/01_dadosDeEntrada "01_dadosDeEntrada") contém os dados brutos do Programa Florestal Global.

[02_script](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/02_script "02_script") contém scripts e funções a serem usados na seguintes ordem:

1.  **LIMPEZA** - [01_limpezaDeDados](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/02_script/01_limpezaDeDados "01_limpezaDeDados")

    -   o arquivo [parte01-limpezaDeDados.Rmd](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/blob/master/02_script/01_limpezaDeDados/parte01-limpezaDeDados.Rmd "parte01-limpezaDeDados.Rmd") corrige dados de taxonomia, ajusta formatos de data e ajusta nomes de colunas

2.  **PREPARAÇÃO** - [02_dataPrepLPI](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/02_script/02_dataPrepLPI "02_dataPrepLPI")

    -   após a limpeza, o arquivo [part02-dataPrepLPI.Rmd](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/blob/master/02_script/02_dataPrepLPI/part02-dataPrepLPI.Rmd "part02-dataPrepLPI.Rmd") usa os arquivos de saída gerados na etapa anterior para calcular e reorganizar os dados no formato necessário para a análise de LPI.

3.  **CÁLCULO** - [03_calculoDeLPI](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/02_script/03_calculoDeLPI "03_calculoDeLPI")

    -   o arquivo [part03-calculateLPI.Rmd](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/blob/master/02_script/03_calculoDeLPI/part03-calculateLPI.Rmd "part03-calculateLPI.Rmd") usa os arquivos de saída preparados na etapa anterior para o cálculo básico de LPI.

        -   as funções localizadas em [funcoes](https://github.com/ICMBio-CENAP/Programa-Monitora-Florestal-Global/tree/master/02_script/funcoes "funcoes") descartam a necessidade de usar o pacote `rLPI`

`03_dadosDeSaida` contém os arquivos de saída gerados nas etapas anteriores.

## Instruções para instalar o repositório localmente no R

Se você pretende usar o R mas nunca fez isso antes, veja abaixo as instruções para instalar o R, RStudio e Git.

1.  Instalar o R: Baixe a versão mais atualizada do R [aqui](https://cran.rstudio.com).
2.  Instalar o Rstudio: [link para download do Rstudio](https://www.rstudio.com/products/rstudio/download/)
3.  Instalar o Git: [Veja essas notas para instalação](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN) - *somente se você quiser compartilhar seu código com esse repositório e/ou quiser manter sua cópia local atualizada na medida em que atualizamos e melhoramos nosso código*.
4.  Clonar ou baixar esse repositório - *botão verde à direita no alto da página*.

# Instruções para instalar o *RLPI*

1 - Instalar o pacote devtools do R

``` r
install.packages("devtools")
```

2 - Instalar o pacote **rlpi** da Zoological-Society-of-London

``` r
library(devtools)

install_github("Zoological-Society-of-London/rlpi", dependencies=TRUE)
```

# Contato

[elildojr\@gmail.com](mailto:elildojr@gmail.com)

[pardalismitis\@gmail.com](mailto:pardalismitis@gmail.com)
