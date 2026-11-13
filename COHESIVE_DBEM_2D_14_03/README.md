# DBEM 2D Coesivo Multidomínio (COHESIVE_DBEM_2D)

Código computacional em Fortran moderno (OpenMP + LAPACK/BLAS) para análise de fratura e destacamento interfacial (debonding) em juntas coladas e estruturas compósitas bidimensionais, utilizando o **Método dos Elementos de Contorno Multidomínio com Zona Coesiva (Multizone Cohesive DBEM 2D)**.

Este repositório e formulação estão diretamente associados ao artigo científico:

> **Comparative analysis of adhesive debonding in three-point bending tests using boundary element methods: Coupled criterion versus interface damage model**  
> *Sérgio Gustavo Ferreira Cordeiro, Thiago V. Birro, René Quispe Rodriguez, Eder Lima de Albuquerque, Lucas Silveira Campos, Jairo Francisco Useche Vivero, Éric Paroissien, Maëlenn Aufray, Frédéric Lachaud.*  
> **Engineering Analysis with Boundary Elements (EABE)**, 2025.  
> DOI / Referência: [ArtigoEABE2025.pdf](../ArtigoEABE2025.pdf)

---

## 1. Visão Geral da Formulação Teórica

### 1.1 Método dos Elementos de Contorno Multidomínio (Multizone BEM)
Diferente do Método dos Elementos Finitos (MEF), que exige a discretização volumétrica de todo o domínio e sofre com elevadas concentrações de tensão nas bordas e singularidades de interface, o BEM discretiza exclusivamente as fronteiras e interfaces. 

Para corpos heterogêneos (ex.: substrato metálico colado a uma camada adesiva polimérica), o domínio é decomposto em sub-regiões $\Omega_1, \Omega_2, \dots, \Omega_N$. Em cada sub-região homogênea e isotrópica, a Equação Integral de Contorno de Deslocamento (BIE) é dada por:

$$c_{ij}(\mathbf{x}') u_j(\mathbf{x}') + \int_{\Gamma} t_{ij}^*(\mathbf{x}', \mathbf{x}) u_j(\mathbf{x}) \, \mathrm{d}\Gamma = \int_{\Gamma} u_{ij}^*(\mathbf{x}', \mathbf{x}) t_j(\mathbf{x}) \, \mathrm{d}\Gamma$$

onde:
- $\mathbf{x}'$ e $\mathbf{x}$ são os pontos fonte e campo;
- $u_{ij}^*$ e $t_{ij}^*$ são as soluções fundamentais de Kelvin para elasticidade 2D (estado plano de tensão ou de deformação);
- $u_j$ e $t_j$ são os campos de deslocamento e tração no contorno $\Gamma$.

### 1.2 Elementos Descontínuos de Contorno (DBEM)
Para evitar indefinições do vetor normal unitário nos cantos geométricos e nas pontas de interface, o código emprega **elementos lineares descontínuos**:
- Os nós geométricos definem os extremos do elemento.
- Nós duplos são alocados nas quinas e transições de interface.
- Os **pontos de colocação** onde as equações integrais são efetivamente avaliadas são deslocados para o interior do elemento paramétrico ($\xi = -1/3$ e $\xi = +1/3$ para elementos lineares), garantindo continuidade de trações e eliminando singularidades geométricas.

### 1.3 Modelo de Dano Progressivo de Interface (PID / Leis Coesivas)
Nas interfaces entre sub-regiões, as superfícies são modeladas como pares coesivos:
- Superfície **$I^-$** (Interface minus, associada à sub-região 1);
- Superfície **$I^+$** (Interface plus, associada à sub-região 2).

A relação entre as trações normais ($t_n$) e cisalhantes ($t_t$) e os saltos de deslocamento ($\Delta u_n, \Delta u_t$) segue uma lei de tração-separação bilinear com dano progressivo:
- **Modo I Puro (Abertura):** rigidez elástica inicial $K_{n0}$, tensão crítica $\sigma_{max}$, e abertura crítica de falha $\delta_{nf}$ baseada na taxa crítica de liberação de energia $\mathcal{G}_{Ic}$.
- **Modo II Puro (Cisalhamento):** rigidez elástica inicial $K_{t0}$, tensão crítica de cisalhamento $\tau_{max}$, e salto de deslocamento de falha $\delta_{tf}$ com base em $\mathcal{G}_{IIc}$.
- **Modo Misto:** critério de interação de tensões e degradação acoplada via parâmetro escalar de dano $\beta \in [0, 1]$.

O equilíbrio não linear e a evolução do dano ao longo dos incrementos de carga são resolvidos através do algoritmo iterativo de **Newton-Raphson**, atualizando a matriz de rigidez tangente interfacial $K_{coh}$ e os resíduos de forças.

---

## 2. Arquitetura dos Módulos Fortran

O código está estruturado de forma modular em Fortran 90/95/2003, compilado com paralelismo OpenMP:

| Módulo | Arquivo | Responsabilidade |
|---|---|---|
| **Principal** | `MAIN_COHESIVE_DBEM_2D.f90` | Ponto de entrada, inicialização das *threads* OpenMP, controle de tempo de CPU e fluxo de execução. |
| **Dados Globais** | `MODULE_DATA.f90` | Definição de tipos (`real64`), matrizes de influência $G$ e $H$, vetores nodais, graus de liberdade e propriedades. |
| **Entrada de Dados** | `MODULE_INPUT_DATA.f90` | Leitura e interpretação dos 5 arquivos de entrada, pontos de quadratura de Gauss e alocações. |
| **Tratamento de Dados** | `MODULE_DATA_TREATMENT.f90` | Mapeamento dos pontos de colocação descontínuos, cálculo de normais, pareamento coesivo $I^- \leftrightarrow I^+$ e valores prescritos $u, t$. |
| **Mapeamento Geométrico** | `MODULE_GEOMETRIC_MAPPINGS.f90` | Funções de forma de Lagrange para a geometria, cálculo do Jacobiano e vetores normais unitários. |
| **Funções de Aproximação** | `MODULE_APPROXIMATION_FUNCTIONS.f90` | Funções de interpolação para deslocamentos e trações nos pontos de colocação. |
| **Integração Numérica** | `MODULE_NUMERICAL_INTEGRATION.f90` | Quadratura de Gauss-Legendre para integrais regulares e tratamento analítico/numérico de singularidades. |
| **Solver Não Linear** | `MODULE_SOLVER.f90` | Montagem das matrizes globais do sistema, aplicação das condições coesivas, iterações de Newton-Raphson e resolução direta via LAPACK (`dgetrf`, `dgetrs`). |
| **Análise Incremental** | `MODULE_ANALYSIS.f90` | Gestão dos passos de propagação coesiva ($N_s$), critérios de convergência e atualização do dano. |
| **Pós-processamento** | `MODULE_POSTPROCESSING.f90` | Geração dos arquivos de saída com as respostas estruturais, coesivas e tensões nos pontos de contorno e internos. |

---

## 3. Arquivos de Entrada e Saída

### 3.1 Arquivos de Entrada (Input)
1. **`1_BOUNDARY_MESH.txt`**:
   - Número de sub-regiões, elementos de contorno e nós.
   - Lista de nós duplicados para tratamento de quinas e descontinuidade.
   - Coordenadas nodais $(X, Y, Z)$.
   - Ordem polinomial de cada elemento (linear = 1) e identificador da sub-região.
   - Conectividade nodal dos elementos.
2. **`2_MATERIALS.txt`**:
   - Propriedades elásticas isotrópicas por sub-região: Módulo de Young ($E$) e Coeficiente de Poisson ($\nu$).
   - Propriedades coesivas em Modo I: $K_{n0}$, $\Phi_n$ ($\mathcal{G}_{Ic}$), $\sigma_{max}$.
   - Propriedades coesivas em Modo II: $K_{t0}$, $\Phi_t$ ($\mathcal{G}_{IIc}$), $\tau_{max}$.
3. **`3_BOUNDARY_CONDITIONS.txt`**:
   - Prescrições de $U_x$, $U_y$, $T_x$ e $T_y$ expressas através de polinômios ao longo do comprimento $S$ do elemento:
     $$v(S) = A \cdot S^2 + B \cdot S + C$$
   - Definição do grau de liberdade monitorado/controlado (tipo $U$, nó, direção e valor).
4. **`4_ANALYSIS_PARAMETERS.txt`**:
   - Número de pontos de Gauss (ex.: 400).
   - Estado plano: Estado Plano de Deformação (`PT` - Plane Strain) ou Tensão (`PS` - Plane Stress).
   - Tipo de análise: Elástica Linear (`LE`) ou Propagação Coesiva (`CP`).
   - Número de passos de carga ($N_s$), tolerância de Newton-Raphson e número máximo de iterações.
   - Classificação das superfícies dos elementos: contorno convencional (`B`), interface menos (`Im`) ou interface mais (`Ip`).
5. **`5_POSTPROCESSING_MESH.txt`**:
   - Coordenadas de pontos internos para avaliação de deslocamentos e tensões no domínio.

### 3.2 Arquivos de Saída (Output)
- **`6_OUTPUT_COLLOCATION_POINTS_COORDINATES.txt`**: Coordenadas físicas $(X, Y)$ de todos os pontos de colocação.
- **`7_OUTPUT_ANALYSIS_TIME.txt`**: Tempos computacionais e convergência por passo de carga.
- **`8_OUTPUT_NODAL_STRUCTURAL_RESPONSE.txt`**: Histórico de deslocamentos $(U_x, U_y)$ e trações $(T_x, T_y)$ no nó monitorado.
- **`9_OUTPUT_COHESIVE_RESPONSE.txt`**: Salto de deslocamento relativo coesivo e tração de interface ao longo dos incrementos.
- **`10_OUTPUT_COHESIVE_STIFFNESS.txt`**: Evolução da matriz de rigidez interfacial degradada.
- **`11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt`**: Solução nodal completa $(X, Y, U_x, U_y, T_x, T_y)$ em todos os nós do contorno no estágio final.

---

## 4. Instruções de Compilação e Execução

### 4.1 Pré-requisitos
- Compilador **GNU Fortran (`gfortran`)** com suporte a OpenMP.
- **`mingw32-make`** (ambiente MSYS2 UCRT64).
- Bibliotecas **LAPACK** e **BLAS** (`liblapack.a`, `libblas.a`).
- Python 3 com `numpy` e `matplotlib` para visualização gráfica.

### 4.2 Compilação via Makefile
No terminal PowerShell, defina o caminho das ferramentas UCRT64 e execute o makefile:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;" + $env:PATH
mingw32-make
```

O comando gerará o executável otimizado `MAIN_COHESIVE_DBEM_2D.exe`.

Para recompilar do zero:
```powershell
mingw32-make veryclean
mingw32-make
```

### 4.3 Execução da Simulação
Execute o binário diretamente no diretório do projeto:

```powershell
.\MAIN_COHESIVE_DBEM_2D.exe
```

A saída no terminal indicará o número de *threads* OpenMP utilizadas e o tempo total de CPU gasto na solução.

---

## 5. Visualização da Malha e Condições de Contorno (`plot_mesh_bc.py`)

Para inspecionar visualmente a geometria, a malha de elementos de contorno, as sub-regiões, a interface coesiva e as condições de contorno (Dirichlet e Neumann), foi desenvolvido o script Python `plot_mesh_bc.py`, baseado no estilo gráfico do `preproc.py` da biblioteca `BEM2D_elastostatic`.

### 5.1 Execução do Script
```powershell
python plot_mesh_bc.py
```

### 5.2 Elementos Visuais Gerados
- **Sub-regiões ($\Omega_1, \Omega_2$):** Preenchimentos translúcidos e contornos com cores distintas para cada material/corpo.
- **Interface Coesiva ($I^- / I^+$):** Linha tracejada em destaque vermelho representando as superfícies emparelhadas do adesivo.
- **Elementos e Nós:** Numeração dos elementos de contorno (`e1, e2, ...`) e nós geométricos.
- **Pontos de Colocação:** Pontos descontínuos internos aos elementos onde as equações integrais são montadas.
- **Condições de Dirichlet:** Marcadores vermelhos (`rx` para $u_x$ prescrito, `rP` para $u_y$ prescrito) com vetores de direção.
- **Condições de Neumann:** Marcadores azuis (`bx` para $t_x$ prescrito, `bP` para $t_y$ prescrito) com setas indicando o sentido das trações aplicadas.
- **Nó de Controle:** Estrela verde destacada com anotação do grau de liberdade controlado pelo método de Newton-Raphson.

A imagem em alta resolução é salva automaticamente em:
`mesh_and_boundary_conditions.png`

---

## 6. Casos de Teste e Validação

Além do modelo simplificado na raiz, a pasta `EXEMPLOS/Artigo em conjunto/Teste - Tese Thiago` contém as malhas completas dos ensaios de flexão em três pontos (3PBT) analisados no artigo EABE 2025:
- `M 4 elements-mm`: Malha preliminar com 4 elementos/mm.
- `M 8 elements-mm`: Malha intermediária com 8 elementos/mm.
- `M 16 elements-mm`: Malha de referência recomendada para análise de fratura.
- `M 32 elements-mm`: Malha altamente refinada para convergência e validação fina da zona coesiva.

---

## 7. Geração de Modelos com Gmsh e Conversor Automatizado (`gmsh_to_dbem.py`)

Para automatizar e facilitar a criação de novos modelos para o DBEM 2D Coesivo, foi implementado um fluxo integrado baseado no Gmsh e no módulo Python `meshio`, seguindo a metodologia de pré-processamento do notebook de referência `Contato.ipynb`.

### 7.1 Arquivo de Geometria Gmsh (`exemplo_small.geo`)
O arquivo `exemplo_small.geo` define as superfícies das sub-regiões e os grupos físicos de contorno:
- **Curvas Físicas (1D)**:
  - `Physical Curve("bot")`: borda inferior (condição de apoio $U_y = 0$).
  - `Physical Curve("interface_m")`: interface coesiva $I^-$ associada à sub-região 1.
  - `Physical Curve("left")`: borda esquerda (engastamento $U_x = 0$).
  - `Physical Curve("right")`: borda direita (tração prescrita $T_x = -1$).
  - `Physical Curve("top")`: borda superior (tração prescrita $T_y = 1$ e grau de controle de deslocamento).
  - `Physical Curve("interface_p")`: interface coesiva $I^+$ associada à sub-região 2.
- **Superfícies Físicas (2D)**:
  - `Physical Surface("subregion_1")`: domínio da sub-região 1 (triângulo inferior).
  - `Physical Surface("subregion_2")`: domínio da sub-região 2 (triângulo superior).

### 7.2 Conversor `gmsh_to_dbem.py`
O script realiza automaticamente o pipeline completo:
1. Gera a malha executando o Gmsh (`.geo` $\to$ `.msh`).
2. Lê o arquivo `.msh` utilizando o `meshio`.
3. Extrai as coordenadas nodais, elementos de linha 1D, identificadores de sub-região e mapeamento das curvas físicas.
4. Exporta e formata com precisão todos os 5 arquivos de entrada necessários para o solver Fortran:
   - `1_BOUNDARY_MESH.txt`
   - `2_MATERIALS.txt`
   - `3_BOUNDARY_CONDITIONS.txt`
   - `4_ANALYSIS_PARAMETERS.txt`
   - `5_POSTPROCESSING_MESH.txt`

### 7.3 Execução do Pré-processador
Para gerar os arquivos de entrada a partir da geometria Gmsh:
```powershell
python gmsh_to_dbem.py
```

Em seguida, basta rodar o solver Fortran:
```powershell
.\MAIN_COHESIVE_DBEM_2D.exe
python plot_mesh_bc.py
```

