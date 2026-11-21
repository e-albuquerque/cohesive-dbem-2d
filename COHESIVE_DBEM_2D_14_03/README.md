# DBEM 2D Coesivo Multidomínio (`COHESIVE_DBEM_2D`)

Código computacional em Fortran moderno (**OpenMP Multithread + LAPACK/BLAS**) com pré e pós-processamento em Python (**Gmsh + `meshio` + `matplotlib`**) para análise elástica linear e propagação de fratura/descolamento interfacial (*debonding*) em juntas coladas e estruturas compósitas bidimensionais, utilizando o **Método dos Elementos de Contorno Dual Multidomínio com Zona Coesiva (Multizone Cohesive DBEM 2D)**.

Este repositório e sua formulação estão diretamente associados ao artigo científico:

> **Comparative analysis of adhesive debonding in three-point bending tests using boundary element methods: Coupled criterion versus interface damage model**  
> *Sérgio Gustavo Ferreira Cordeiro, Thiago V. Birro, René Quispe Rodriguez, Eder Lima de Albuquerque, Lucas Silveira Campos, Jairo Francisco Useche Vivero, Éric Paroissien, Maëlenn Aufray, Frédéric Lachaud.*  
> **Engineering Analysis with Boundary Elements (EABE)**, vol. 180, 106474, 2025.  
> Arquivo de referência: [ArtigoEABE2025.pdf](../ArtigoEABE2025.pdf)

---

## 1. Visão Geral da Formulação Teórica

### 1.1 Método dos Elementos de Contorno Multidomínio (Multizone DBEM)
Diferente do Método dos Elementos Finitos (MEF), que exige a discretização volumétrica de todo o domínio e sofre com severa distorção de malha em camadas adesivas finas, o BEM discretiza exclusivamente as fronteiras e interfaces dos corpos.

Para sistemas heterogêneos (ex.: substrato metálico colado a um bloco de adesivo polimérico), o domínio é decomposto em sub-regiões $\Omega_1, \Omega_2, \dots, \Omega_N$. Em cada sub-região homogênea e isotrópica, a Equação Integral de Contorno de Deslocamento (DBIE) e a Equação Integral de Contorno de Tração (TBIE) são dadas por:

$$c_{ij}(\mathbf{x}') u_j(\mathbf{x}') + \text{CPV} \int_{\Gamma} T_{ij}^*(\mathbf{x}', \mathbf{x}) u_j(\mathbf{x}) \, \mathrm{d}\Gamma = \int_{\Gamma} U_{ij}^*(\mathbf{x}', \mathbf{x}) t_j(\mathbf{x}) \, \mathrm{d}\Gamma$$

$$c_{ij}(\mathbf{x}') t_j(\mathbf{x}') + n_k(\mathbf{x}') \, \text{HFP} \int_{\Gamma} S_{kij}^*(\mathbf{x}', \mathbf{x}) u_j(\mathbf{x}) \, \mathrm{d}\Gamma = n_k(\mathbf{x}') \, \text{CPV} \int_{\Gamma} D_{kij}^*(\mathbf{x}', \mathbf{x}) t_j(\mathbf{x}) \, \mathrm{d}\Gamma$$

onde:
- $\mathbf{x}'$ e $\mathbf{x}$ são os pontos fonte (colocação) e campo;
- $U_{ij}^*, T_{ij}^*, D_{kij}^*, S_{kij}^*$ são os tensores fundamentais de Kelvin para elasticidade 2D (estado plano de tensão `PT` ou estado plano de deformação `PD`);
- $\text{CPV}$ e $\text{HFP}$ indicam integrais avaliadas no sentido do Valor Principal de Cauchy e Parte Finita de Hadamard.

### 1.2 Elementos Descontínuos e Semidescontínuos
Para evitar indefinições do vetor normal unitário em quinas geométricas e transições de condições de contorno/interface:
- Nós duplos são alocados nos vértices e extremidades de interface;
- Nos elementos descontínuos (`D`, `LD`, `RD`), os **pontos de colocação** são posicionados no interior do elemento paramétrico ($\xi = \pm 1/3$ para elementos lineares), garantindo avaliação suave e precisa das trações interfaciais.

### 1.3 Acoplamento Interfacial e Modelo de Dano Coesivo (PID)
As interfaces entre sub-regiões são divididas em pares de superfícies:
- Superfície **$I^-$ (`Im`)**: associada à sub-região 1 (ex.: substrato);
- Superfície **$I^+$ (`Ip`)**: associada à sub-região 2 (ex.: adesivo).

O solver suporta três modos de análise configuráveis em `4_ANALYSIS_PARAMETERS.txt`:
1. **Análise Elástica Linear (`LE`)**: avalia as tensões normais de descolamento (*Peel Stress*, $\sigma$) e de cisalhamento (*Shear Stress*, $\tau$) ao longo da interface íntegra com rigidez penalizada ($K_{n0}, K_{t0}$), usada na aplicação do **Critério Acoplado (Coupled Criterion - CC)**.
2. **Propagação Coesiva Não Linear Simétrica (`CP`)**: utiliza uma lei de tração-separação bilinear em Modo Misto (Modo I: $K_{n0}, \Phi_n, \sigma_{\max}$; Modo II: $K_{t0}, \Phi_t, \tau_{\max}$) com variável escalar de dano $\beta \in [0, 1]$ resolvida incrementalmente via **Newton-Raphson** com controle de deslocamento ou salto relativo.
3. **Propagação Coesiva Não Linear Assimétrica (`CA`)**: reproduz o descolamento assimétrico observado experimentalmente no ensaio 3PBT (Seção 6.3.2 do artigo EABE 2025), bloqueando a iniciação espúria na extremidade direita ($\sigma_c = 10^{50}\text{ MPa}$) para que uma única trinca nuclee e se propague a partir da extremidade esquerda.

---

## 2. Arquitetura Modular e Otimizações de Performance (OpenMP)

O solver Fortran foi revisado e otimizado para alta performance mantendo a clareza matemática da formulação:

| Módulo | Arquivo | Responsabilidade e Otimizações |
|---|---|---|
| **Principal** | `MAIN_COHESIVE_DBEM_2D.f90` | Inicialização do ambiente OpenMP e medição de tempo real (*wall-clock*) via `omp_get_wtime()`. |
| **Dados Globais** | `MODULE_DATA.f90` | Estruturas globais em precisão dupla (`real64`), matrizes $G$ e $H$, conectividades e propriedades. |
| **Entrada de Dados** | `MODULE_INPUT_DATA.f90` | Leitura dos 5 arquivos `.txt`, quadratura de Gauss-Legendre e mapeamento topológico sem vazamentos de memória. |
| **Tratamento de Dados** | `MODULE_DATA_TREATMENT.f90` | Cálculo dos pontos de colocação descontínuos e pareamento paralelo (`!$OMP PARALLEL DO`) dos nós coesivos $I^- \leftrightarrow I^+$. |
| **Mapeamento Geométrico** | `MODULE_GEOMETRIC_MAPPINGS.f90` | Interpolação geométrica de Lagrange, Jacobiano (`DSQRT`) e vetor normal exterior. |
| **Funções de Forma** | `MODULE_APPROXIMATION_FUNCTIONS.f90` | Funções de base e derivadas nos pontos nodais e de colocação. |
| **Integração Numérica** | `MODULE_NUMERICAL_INTEGRATION.f90` | Integrações regulares e singulares (CPV/HFP) para DBIE e TBIE utilizando *arrays* automáticos na *stack* (eliminando $\sim 400.000$ alocações dinâmicas no *heap*) e constantes elásticas pré-calculadas fora do laço de Gauss. |
| **Análise e Montagem** | `MODULE_ANALYSIS.f90` | Montagem das matrizes $H$ e $G$ paralelizada com **`!$OMP PARALLEL DO SCHEDULE(DYNAMIC)`** e tabela de busca rápida $O(1)$ de elementos por nó (`NODE_ELEM`). |
| **Solver Linear/Não Linear** | `MODULE_SOLVER.f90` | Montagem paralela da matriz do sistema acoplado, produto esparso de condições de contorno prescritas e solução direta via LAPACK (`DGESV`). |
| **Pós-processamento** | `MODULE_POSTPROCESSING.f90` | Avaliação paralelizada com OpenMP de deslocamentos e tensões em pontos internos e de contorno, exportação de `11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt` e arquivo `OUTPUT.vtk` portável. |

### Ganho de Desempenho Alcançado
No benchmark do artigo EABE 2025 ($N_{\text{nodes}} = 445$ nós de contorno, $N_{\text{elem}} = 445$ elementos, $N_{\text{Gauss}} = 200$ pontos de quadratura por elemento, processador de 8 threads):
- **Versão original (`-O0 -fbounds-check`, sequencial)**: $\approx 38,05\text{ s}$ por simulação.
- **Com otimização de compilador (`-O3`, sequencial)**: $\approx 11,75\text{ s}$ por simulação ($3,2\times$ mais rápido).
- **Com paralelização OpenMP + eliminação de alocações no laço interno + `-O3`**: **$\approx 2,97\text{ s}$ por simulação** (**$12,8\times$ mais rápido**, mantendo precisão numérica idêntica até a 9ª casa decimal).

---

## 3. Notebooks Interativos (Execução Local e Google Colab)

O projeto conta com três Jupyter Notebooks completos que cobrem tanto a análise elástica linear quanto a análise não linear incremental de dano coesivo (PID):

1. **`Run_cohesive_dbem_2d.ipynb`** (*Análise Linear Elástica — Local / Multiplataforma*):
   - Executa o fluxo completo desde a criação do arquivo `.geo` no Gmsh, geração da malha `.msh`, pré-processamento via `meshio`, compilação do Fortran e resolução elástica para $h_s = 1.08\text{ mm}$ e $h_s = 1.62\text{ mm}$, reproduzindo a **Fig. 10** e as **Figs. 12, 13, 14 e 15** do artigo.
2. **`Run_cohesive_dbem_2d_on_Colab.ipynb`** (*Análise Linear Elástica — Google Colab Linux*):
   - Projetado especificamente para execução no **Google Colab (Ubuntu Linux)**, seguindo o padrão de comandos diretos (`!apt-get`, `!make`, `!gmsh`, `!./MAIN_COHESIVE_DBEM_2D`).
3. **`Run_cohesive_dbem_2d_nlinear.ipynb`** (*Análise Não Linear de Dano Progressivo PID — Passo a Passo*):
   - Detalha passo a passo a teoria e a implementação da análise não linear incremental por controle de deslocamento (Newton-Raphson com lei coesiva bilinear em Modo Misto) para descolamento assimétrico (`ANALYSIS_TYPE = 'CA'`), reproduzindo a **Fig. 26** (*Force–displacement results for asymmetric debonding: $h_s = 1.08\text{ mm}$*) e a evolução do comprimento da trinca (análogo à **Fig. 30**).

---

## 4. Benchmark do Artigo EABE 2025 (Ensaio de Flexão em 3 Pontos - Fig. 10 e Fig. 26)

### 4.1 Parâmetros Geométricos e Físicos
- **Substrato ($\Omega_1$)**: comprimento $L_s = 50,0\text{ mm}$, espessura $h_s = 1,08\text{ mm}$ ou $h_s = 1,62\text{ mm}$, $E_1 = 136\,000\text{ MPa}$ ($2 \times 68\,000\text{ MPa}$ pela razão de larguras $b_1/b_2$), $\nu_1 = 0,33$.
- **Adesivo ($\Omega_2$)**: comprimento de colagem $L_a = 25,0\text{ mm}$ ($x \in [-12,5; +12,5]\text{ mm}$), espessura $h_a = 4,0\text{ mm}$, $E_2 = 3\,000\text{ MPa}$, $\nu_2 = 0,35$.
- **Apoios e Carregamento**: vão entre apoios $d_{\text{supp}} = 35,90\text{ mm}$ (rolete à esquerda $u_y = 0$, apoio fixo à direita $u_x = u_y = 0$, largura $0,12\text{ mm}$) e deslocamento vertical prescrito no topo central ($w_{\text{load}} = 2,0\text{ mm}$):
  - Para $h_s = 1,08\text{ mm}$: $U_c = 0,576\text{ mm}$ ($F_c = 178\text{ N}$, $\sigma_c = 44,654\text{ MPa}$, $\mathcal{G}_c = 0,520\text{ N/mm}$, $k_{n0} = k_{t0} = 2000\text{ MPa/mm}$);
  - Para $h_s = 1,62\text{ mm}$: $U_c = 0,6093\text{ mm}$ ($F_c = 326\text{ N}$, $\sigma_c = 45,906\text{ MPa}$, $\mathcal{G}_c = 0,5371\text{ N/mm}$, $k_{n0} = k_{t0} = 948\text{ MPa/mm}$).

### 4.2 Figuras Geradas Automaticamente
- `mesh_and_boundary_conditions.png`: Reprodução da **Fig. 10** (malha e condições de contorno para análise linear $h_s = 1,08\text{ mm}$).
- `mesh_and_boundary_conditions_nlinear.png`: Malha, condições de contorno e nó de controle para a análise não linear por controle de deslocamento ($h_s = 1,08\text{ mm}$).
- `fig12_peel_hs108.png`: **Fig. 12** — *Peel stress over the undamaged interface: $h_s = 1.08\text{ mm}$* (pico de $\approx 43,45\text{ MPa}$ nas bordas $x = \pm 12,5\text{ mm}$).
- `fig13_shear_hs108.png`: **Fig. 13** — *Shear stress over the undamaged interface: $h_s = 1.08\text{ mm}$* (extremos antissimétricos de $\approx \mp 13,96\text{ MPa}$).
- `fig14_peel_hs162.png`: **Fig. 14** — *Comparison of peel Stress over the undamaged interface: $h_s = 1.62\text{ mm}$* (pico de $\approx 46,15\text{ MPa}$).
- `fig15_shear_hs162.png`: **Fig. 15** — *Comparison of shear stress over the undamaged interface: $h_s = 1.62\text{ mm}$* (extremos de $\approx \mp 15,09\text{ MPa}$).
- `fig26_force_displacement_hs108.png`: **Fig. 26** — *Force–displacement results for asymmetric debonding: $h_s = 1.08\text{ mm}$* (curva não linear **2D BEM - PID**).
- `fig26_cohesive_damage_evolution_hs108.png`: Evolução do comprimento da trinca assimétrica $a(u)$ (análogo à **Fig. 30**) e perfil de degradação da rigidez coesiva $K_{nn}/k_{n0}$.
- `debonding_animation_hs108.mp4` (gerado por `animate_debonding.py`): Animação em vídeo `.mp4` mostrando a malha DBEM 2D se deformando em escala $1:1$ (**Fig. 28**) em sincronia com o gráfico *Force × Central displacement* (**Fig. 26**).

---

## 5. Compilação e Execução via Linha de Comando

### No Windows (MSYS2 MinGW-w64 UCRT64)
```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
mingw32-make clean
mingw32-make
.\MAIN_COHESIVE_DBEM_2D.exe
```

### No Linux / Google Colab
```bash
sudo apt-get update && sudo apt-get install -y gfortran liblapack-dev libblas-dev gmsh
pip install numpy matplotlib meshio gmsh
make clean && make
./MAIN_COHESIVE_DBEM_2D
```
