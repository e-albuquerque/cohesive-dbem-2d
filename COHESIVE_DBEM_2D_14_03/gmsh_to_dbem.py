"""
gmsh_to_dbem.py
---------------
Script de pré-processamento para geração de malhas com Gmsh e conversão
para os arquivos de entrada do código Fortran DBEM 2D Coesivo:
  - 1_BOUNDARY_MESH.txt
  - 2_MATERIALS.txt
  - 3_BOUNDARY_CONDITIONS.txt
  - 4_ANALYSIS_PARAMETERS.txt
  - 5_POSTPROCESSING_MESH.txt

Utiliza como referência o fluxo do Jupyter Notebook 'Contato.ipynb'.
"""

import os
import sys
import subprocess
import numpy as np
import meshio


# =============================================================================
# 1. GERAÇÃO DA GEOMETRIA E EXECUÇÃO DO GMSH (.geo -> .msh)
# =============================================================================

def generate_geometry_and_mesh(file_name, geo_code=None):
    """
    Gera o arquivo .geo e executa o Gmsh para criar a malha .msh.
    Suporta execução via CLI gmsh/gmsh.bat com fallback para o módulo python gmsh.
    """
    geo_file = f"{file_name}.geo"
    msh_file = f"{file_name}.msh"

    # Se um código .geo for fornecido, grava o arquivo
    if geo_code is not None:
        with open(geo_file, "w", encoding="utf-8") as f:
            f.write(geo_code.strip() + "\n")
        print(f"[Gmsh] Arquivo de geometria gravado: '{geo_file}'")

    if not os.path.exists(geo_file):
        raise FileNotFoundError(f"Arquivo de geometria '{geo_file}' não encontrado.")

    print(f"[Gmsh] Gerando malha 2D a partir de '{geo_file}'...")

    # Tentativa 1: Execução via executável de linha de comando ou script gmsh.bat
    gmsh_cmds = ['gmsh', 'gmsh.bat']
    # Adiciona diretório comum de scripts do Python se necessário
    python_scripts = os.path.join(os.path.dirname(sys.executable), "Scripts")
    if os.path.exists(python_scripts):
        gmsh_cmds.append(os.path.join(python_scripts, "gmsh.bat"))

    mesh_ok = False
    for cmd in gmsh_cmds:
        try:
            res = subprocess.run(
                [cmd, '-2', geo_file, '-o', msh_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if res.returncode == 0 and os.path.exists(msh_file):
                mesh_ok = True
                print(f"[Gmsh] Malha gerada com sucesso via '{cmd}': '{msh_file}'")
                break
        except FileNotFoundError:
            continue

    # Tentativa 2: Fallback utilizando a biblioteca oficial Python `gmsh`
    if not mesh_ok:
        try:
            import gmsh
        except ImportError:
            print("[Gmsh] Módulo python 'gmsh' não encontrado. Instalando automaticamente via pip...")
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", "gmsh"])
                import importlib
                importlib.invalidate_caches()
                import gmsh
            except Exception as err_install:
                raise RuntimeError(
                    f"O executável Gmsh não foi encontrado e a instalação do pacote 'gmsh' falhou: {err_install}"
                )

        try:
            gmsh.initialize()
            gmsh.open(geo_file)
            gmsh.model.mesh.generate(2)
            gmsh.write(msh_file)
            gmsh.finalize()
            if os.path.exists(msh_file):
                mesh_ok = True
                print(f"[Gmsh] Malha gerada com sucesso via API Python gmsh: '{msh_file}'")
        except Exception as e:
            raise RuntimeError(f"Falha ao executar API do Gmsh. Erro: {e}")

    if not mesh_ok:
        raise RuntimeError("Não foi possível gerar a malha .msh.")


# =============================================================================
# 2. DEFINIÇÃO DOS DADOS DE ENTRADA DO PROBLEMA
# =============================================================================

def input_data():
    """
    Define parâmetros geométricos, materiais, condições de contorno e
    parâmetros de análise para o modelo de referência (Exemplo Small - Placa Bipartida).
    """
    file_name = "exemplo_small"

    # Código .geo correspondente à geometria de referência
    geo_code = """// ==============================================================================
// Exemplo DBEM 2D Coesivo: Placa Bipartida Diagonal (Exemplo Small)
// ==============================================================================
ref = 0.1;

// Nós geométricos da Sub-região 1 (Triângulo Inferior)
Point(1) = {0.000, 0.000, 0, ref};
Point(2) = {0.005, 0.000, 0, ref};
Point(3) = {0.095, 0.000, 0, ref};
Point(4) = {0.100, 0.000, 0, ref};
Point(5) = {0.000, 0.100, 0, ref};
Point(6) = {0.000, 0.095, 0, ref};

// Linhas da Sub-região 1 (Sentido anti-horário)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5}; // Interface Im
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Nós geométricos da Sub-região 2 (Triângulo Superior)
Point(7) = {0.100, 0.000, 0, ref};
Point(8) = {0.100, 0.100, 0, ref};
Point(9) = {0.000, 0.100, 0, ref};

// Linhas da Sub-região 2 (Sentido anti-horário)
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 7}; // Interface Ip

Curve Loop(2) = {7, 8, 9};
Plane Surface(2) = {2};

// Manter exatamente 1 elemento por linha para o exemplo small
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9} = 2 Using Progression 1;

// Grupos Físicos de Contorno (1D)
Physical Curve("bot") = {1, 2, 3};
Physical Curve("interface_m") = {4};
Physical Curve("left") = {5, 6};
Physical Curve("right") = {7};
Physical Curve("top") = {8};
Physical Curve("interface_p") = {9};

// Grupos Físicos de Superfície (2D)
Physical Surface("subregion_1") = {1};
Physical Surface("subregion_2") = {2};
"""

    # Propriedades dos Materiais (Sub-regiões e Interface Coesiva)
    materials = {
        'subregions': {
            1: {'E': 32.0e3, 'nu': 0.2},
            2: {'E': 32.0e3, 'nu': 0.2},
        },
        'cohesive_mode_I': {
            'Kn0': 1.0e9,
            'phin': 0.0001,
            'smax': 4.0,
        },
        'cohesive_mode_II': {
            'Kto': 1.0e9,
            'phit': 0.0002,
            'tmax': 3.0,
        }
    }

    # Condições de Contorno por Physical Curve (conforme Contato.ipynb)
    # Formato: (A, B, C) onde prescrição é A*S^2 + B*S + C
    bound_cond = {
        'bot': {
            'surface_type': 'B',
            'Ux': None,
            'Uy': (0.0, 0.0, 0.0),
        },
        'interface_m': {
            'surface_type': 'Im',
            'Ux': None, 'Uy': None, 'Tx': None, 'Ty': None,
        },
        'left': {
            'surface_type': 'B',
            'Ux': (0.0, 0.0, 0.0),
            'Uy': None,
        },
        'right': {
            'surface_type': 'B',
            'Tx': (0.0, 0.0, -1.0),
            'Ty': None,
        },
        'top': {
            'surface_type': 'B',
            'Tx': None,
            'Ty': (0.0, 0.0, 1.0),
        },
        'interface_p': {
            'surface_type': 'Ip',
            'Ux': None, 'Uy': None, 'Tx': None, 'Ty': None,
        },
    }

    # Condição específica de contorno adicional (Elemento 1 tem Ux=0 no canto inferior esquerdo)
    element_specific_bc = {
        1: {'Ux': (0.0, 0.0, 0.0)}
    }

    # Parâmetros Numéricos e de Análise (4_ANALYSIS_PARAMETERS.txt)
    analysis_params = {
        'N_GAUSS': 400,
        'CP_DISTANCING': 0.0,
        'PLANE_STATE': 'PT',      # PT = Plane Strain, PS = Plane Stress
        'ANALYSIS_TYPE': 'CP',    # CP = Cohesive Propagation, LE = Linear Elastic
        'Ns': 50,                 # Número de passos de carga
        'tol_Newton': 1.0e-2,     # Tolerância Newton-Raphson
        'itemax': 100000,         # Número máximo de iterações
        'structural_node': 9,     # Nó monitorado da resposta estrutural
        'cohesive_node': 18,      # Nó monitorado da resposta coesiva
        'control': {
            'type': 'U',
            'node': 17,           # Nó com deslocamento prescrito/controlado
            'direction': 2,       # 1 = X, 2 = Y
            'value': 0.00012,     # Valor do incremento
        }
    }

    return {
        'file_name': file_name,
        'geo_code': geo_code,
        'materials': materials,
        'bound_cond': bound_cond,
        'element_specific_bc': element_specific_bc,
        'analysis_params': analysis_params,
    }


# =============================================================================
# 3. EXTRAÇÃO DE DADOS DA MALHA (.msh via meshio)
# =============================================================================

def extract_surface_data(mesh):
    """
    Coleta as sub-regiões (superfícies) e seus respectivos elementos e nós a partir da malha.
    Adaptado diretamente de extract_surface_data de 'Contato.ipynb'.
    """
    surface_data = {}
    for key, value in mesh.field_data.items():
        if key.startswith("surf") or "subregion" in key:
            surface_data[key] = value[0]

    triangles = mesh.cell_data_dict.get('gmsh:physical', {}).get('triangle', [])
    triangles_surf = {surface: [] for surface in surface_data}

    for i, triangle in enumerate(triangles):
        for surface, element in surface_data.items():
            if triangle == element:
                triangles_surf[surface].append(i)
                break

    tri_elem = mesh.cells_dict.get('triangle', np.array([]))
    lines_surf = {surface: [] for surface in surface_data}
    nodes_surf = {surface: [] for surface in surface_data}

    ind_elem = 0
    line_elements = mesh.cells_dict.get('line', [])
    for elem in line_elements:
        for surface, tri_list in triangles_surf.items():
            found = False
            for triangle in tri_list:
                t = tri_elem[triangle]
                if (elem[0] == t[0] and elem[1] == t[1]) or \
                   (elem[0] == t[1] and elem[1] == t[2]) or \
                   (elem[0] == t[2] and elem[1] == t[0]) or \
                   (elem[1] == t[0] and elem[0] == t[1]) or \
                   (elem[1] == t[1] and elem[0] == t[2]) or \
                   (elem[1] == t[2] and elem[0] == t[0]):
                    lines_surf[surface].append(ind_elem)
                    nodes_surf[surface].extend([elem[0], elem[1]])
                    found = True
                    break
            if found:
                break
        ind_elem += 1

    return triangles_surf, nodes_surf, lines_surf


def compute_inodes(file_name, bound_cond):
    """
    Lê o arquivo .msh via meshio e extrai as coordenadas, conectividade 1D,
    associação com sub-regiões e mapeamento de condições de contorno.
    Baseado na função compute_inodes de 'Contato.ipynb'.
    """
    msh_file = f"{file_name}.msh"
    if not os.path.exists(msh_file):
        raise FileNotFoundError(f"Arquivo de malha '{msh_file}' não encontrado.")

    mesh = meshio.read(msh_file)
    coordinates = mesh.points
    line_elements = mesh.cells_dict['line']
    segments = mesh.cell_data_dict['gmsh:physical']['line']

    # Mapeamento de tags para nomes de curvas físicas
    tag_to_curve = {v[0]: k for k, v in mesh.field_data.items() if v[1] == 1}
    tag_to_surface = {v[0]: k for k, v in mesh.field_data.items() if v[1] == 2}

    # Identificar superfícies (sub-regiões)
    triangles_surf, nodes_surf, lines_surf = extract_surface_data(mesh)

    # Identificar sub-região de cada elemento 1D
    elem_subregion = np.ones(len(line_elements), dtype=int)
    surf_names_ordered = sorted(list(lines_surf.keys()))
    for s_idx, sname in enumerate(surf_names_ordered, start=1):
        for line_idx in lines_surf[sname]:
            elem_subregion[line_idx] = s_idx

    # Caso fallback: se nenhuma superfície 2D foi capturada pelo teste de triângulos,
    # agrupa por conectividade sequencial
    n_subreg = len(surf_names_ordered) if len(surf_names_ordered) > 0 else 1

    # Filtra e renumera apenas os nós pertencentes aos elementos de contorno 1D (BEM)
    unique_nodes, inv = np.unique(line_elements, return_inverse=True)
    line_elements_renumbered = inv.reshape(line_elements.shape)
    coordinates_boundary = coordinates[unique_nodes]

    return {
        'coordinates': coordinates_boundary,
        'line_elements': line_elements_renumbered,
        'segments': segments,
        'tag_to_curve': tag_to_curve,
        'tag_to_surface': tag_to_surface,
        'elem_subregion': elem_subregion,
        'n_subregions': n_subreg,
        'boundary_nodes': np.arange(len(coordinates_boundary)),
        'mesh': mesh,
    }


# =============================================================================
# 4. ESCRITA DOS ARQUIVOS DE ENTRADA DO FORTRAN
# =============================================================================

def format_fortran_real(val):
    """Formata número real em formato aceito pelo compilador Fortran (ex: 0.D0, 1.D9)."""
    if abs(val) < 1e-15:
        return "0.D0"
    exp_str = f"{val:E}".replace('E', 'D')
    # Simplificações visuais comuns
    if exp_str.endswith("D+00"):
        return f"{val:.1f}D0" if float(val).is_integer() else f"{val}D0"
    return exp_str


def create_fortran_input_files(output_dir, comp_data, inp_data):
    """
    Gera todos os arquivos de entrada estruturados para o código Fortran DBEM 2D:
      - 1_BOUNDARY_MESH.txt
      - 2_MATERIALS.txt
      - 3_BOUNDARY_CONDITIONS.txt
      - 4_ANALYSIS_PARAMETERS.txt
      - 5_POSTPROCESSING_MESH.txt
    """
    os.makedirs(output_dir, exist_ok=True)

    coords = comp_data['coordinates']
    line_elems = comp_data['line_elements']
    segments = comp_data['segments']
    tag_to_curve = comp_data['tag_to_curve']
    elem_subreg = comp_data['elem_subregion']
    n_subregions = comp_data['n_subregions']

    n_elem = len(line_elems)
    n_nodes = len(coords)
    n_dup = n_nodes

    bound_cond = inp_data['bound_cond']
    elem_spec_bc = inp_data.get('element_specific_bc', {})
    materials = inp_data['materials']
    params = inp_data['analysis_params']

    # --- 1. Escrever 1_BOUNDARY_MESH.txt ---
    mesh_path = os.path.join(output_dir, '1_BOUNDARY_MESH.txt')
    with open(mesh_path, 'w', encoding='utf-8') as f:
        f.write("**************** BOUNDARY MESH DEFINITIONS *****************\n")
        f.write("WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA\n")
        f.write("AND THE NUMBER FORMAT\n\n")
        f.write("NUMBER OF SUBREGIONS\n")
        f.write(f"{n_subregions}\n\n")
        f.write("NUMBER OF BOUNDARY ELEMENTS\n")
        f.write(f"{n_elem}\n\n")
        f.write("NUMBER OF NODES \n")
        f.write(f"{n_nodes}\n\n")
        f.write("NUMBER OF DUPLICATE NODES \n")
        f.write(f"{n_dup}\n\n")
        f.write("DUPLICATE NODES\n")
        f.write("NODE\n")
        for i in range(1, n_dup + 1):
            f.write(f"{i}\n")
        f.write("\n")
        f.write("NODAL COORDINATES\n")
        f.write("NODE\tX\t     Y\t          Z\n")
        for i in range(n_nodes):
            x, y = coords[i, 0], coords[i, 1]
            z = coords[i, 2] if coords.shape[1] > 2 else 0.0
            f.write(f"{i+1}\t{format_fortran_real(x)}\t{format_fortran_real(y)}\t{z:.2f}\n")
        f.write("\n")
        f.write("INFORM THE POLYNOMIAL ORDER FOR ALL ELEMENTS\n")
        f.write("ELEMENT\t        ORDER\t        SUBREGION         \t        \n")
        for i in range(n_elem):
            f.write(f" {i+1}\t          1\t          {elem_subreg[i]}\n")
        f.write("               \t\t\n")
        f.write("INFORM THE CONECTIVITY FOR ALL ELEMENTS\n")
        f.write("NODES\n")
        for i in range(n_elem):
            n1 = line_elems[i, 0] + 1
            n2 = line_elems[i, 1] + 1
            f.write(f"{n1}\t{n2}\n")

    print(f"[OK] Gerado: {mesh_path}")

    # --- 2. Escrever 2_MATERIALS.txt ---
    mat_path = os.path.join(output_dir, '2_MATERIALS.txt')
    with open(mat_path, 'w', encoding='utf-8') as f:
        f.write("**************** MATERIALS DEFINITIONS *****************\n")
        f.write("WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA\n")
        f.write("AND THE NUMBER FORMAT\n\n")
        f.write("INFORM THE MATERIAL PROPERTIES\n")
        f.write("SUBREGION   Young_Modulus       Poisson \n")
        for s in range(1, n_subregions + 1):
            s_mat = materials['subregions'].get(s, materials['subregions'][1])
            e_val = format_fortran_real(s_mat['E'])
            nu_val = format_fortran_real(s_mat['nu'])
            f.write(f"   {s}            {e_val}            {nu_val}\n")
        f.write("\n")
        f.write("INFORM THE COHESIVE PROPERTIES FOR MODE I\n")
        f.write("Kn0          phin          smax             \n")
        m1 = materials['cohesive_mode_I']
        f.write(f"{format_fortran_real(m1['Kn0'])}       {format_fortran_real(m1['phin'])}        {format_fortran_real(m1['smax'])}                         \n\n")
        f.write("INFORM THE COHESIVE PROPERTIES FOR MODE II\n")
        f.write("Kto          phit          tmax           \n")
        m2 = materials['cohesive_mode_II']
        f.write(f"{format_fortran_real(m2['Kto'])}       {format_fortran_real(m2['phit'])}        {format_fortran_real(m2['tmax'])}                       \n")

    print(f"[OK] Gerado: {mat_path}")

    # --- 3. Escrever 3_BOUNDARY_CONDITIONS.txt ---
    # Coletar elementos com condições prescritas
    ux_list, uy_list, tx_list, ty_list = [], [], [], []

    for i in range(n_elem):
        eid = i + 1
        tag = segments[i]
        cname = tag_to_curve.get(tag, '')
        bc_rule = bound_cond.get(cname, {})

        # Prescrição em Ux
        ux_val = bc_rule.get('Ux', None)
        if eid in elem_spec_bc and 'Ux' in elem_spec_bc[eid]:
            ux_val = elem_spec_bc[eid]['Ux']
        if ux_val is not None:
            a, b, c = ux_val if isinstance(ux_val, (tuple, list)) else (0.0, 0.0, ux_val)
            ux_list.append((eid, a, b, c))

        # Prescrição em Uy
        uy_val = bc_rule.get('Uy', None)
        if eid in elem_spec_bc and 'Uy' in elem_spec_bc[eid]:
            uy_val = elem_spec_bc[eid]['Uy']
        if uy_val is not None:
            a, b, c = uy_val if isinstance(uy_val, (tuple, list)) else (0.0, 0.0, uy_val)
            uy_list.append((eid, a, b, c))

        # Prescrição em Tx
        tx_val = bc_rule.get('Tx', None)
        if eid in elem_spec_bc and 'Tx' in elem_spec_bc[eid]:
            tx_val = elem_spec_bc[eid]['Tx']
        if tx_val is not None:
            a, b, c = tx_val if isinstance(tx_val, (tuple, list)) else (0.0, 0.0, tx_val)
            tx_list.append((eid, a, b, c))

        # Prescrição em Ty
        ty_val = bc_rule.get('Ty', None)
        if eid in elem_spec_bc and 'Ty' in elem_spec_bc[eid]:
            ty_val = elem_spec_bc[eid]['Ty']
        if ty_val is not None:
            a, b, c = ty_val if isinstance(ty_val, (tuple, list)) else (0.0, 0.0, ty_val)
            ty_list.append((eid, a, b, c))

    bc_path = os.path.join(output_dir, '3_BOUNDARY_CONDITIONS.txt')
    with open(bc_path, 'w', encoding='utf-8') as f:
        f.write("**************** BOUNDARY CONDITIONS DEFINITIONS *****************\n")
        f.write("WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA\n")
        f.write("AND THE NUMBER FORMAT\n\n")
        f.write("NUMBER OF ELEMENTS WITH Ux PRECRIBED\n")
        f.write(f"{len(ux_list)}\n\n")
        f.write("NUMBER OF ELEMENTS WITH Uy PRECRIBED\n")
        f.write(f"{len(uy_list)}\n\n")
        f.write("NUMBER OF ELEMENTS WITH Tx PRESCRIBED\n")
        f.write(f"{len(tx_list)}\n\n")
        f.write("NUMBER OF ELEMENTS WITH Ty PRESCRIBED\n")
        f.write(f"{len(ty_list)}\n\n")

        f.write("BOUNDARY CONDITIONS\n")
        f.write("II\tAux\t Bux\t Cux  \n")
        for (eid, a, b, c) in ux_list:
            f.write(f"{eid}\t{format_fortran_real(a)}\t{format_fortran_real(b)}\t{format_fortran_real(c)}\n")
        f.write("\n")

        f.write("II\tAuy\t Buy\t Cuy\n")
        for (eid, a, b, c) in uy_list:
            f.write(f"{eid}\t{format_fortran_real(a)}\t{format_fortran_real(b)}\t{format_fortran_real(c)}\n")
        f.write("\n")

        f.write("II\tAtx\t Btx \t Ctx\n")
        for (eid, a, b, c) in tx_list:
            f.write(f"{eid}\t{format_fortran_real(a)}\t{format_fortran_real(b)}\t{format_fortran_real(c)}\n")
        f.write("\n")

        f.write("II\tAty\t Bty \t Cty \n")
        for (eid, a, b, c) in ty_list:
            f.write(f"{eid}\t{format_fortran_real(a)}\t{format_fortran_real(b)}\t{format_fortran_real(c)}\n")
        f.write("\n")

        f.write("CONTROLED PARAMETER (U=DISPLACEMENT, T=TRACTION, DU=CRACK_RELATIVE_DISPLACEMENT)  \n")
        f.write("TYPE    NODE   DIRECTION      VALUE \n")
        ctrl = params['control']
        c_val = format_fortran_real(ctrl['value'])
        f.write(f" {ctrl['type']}        {ctrl['node']}        {ctrl['direction']}         {c_val} \n")

    print(f"[OK] Gerado: {bc_path}")

    # --- 4. Escrever 4_ANALYSIS_PARAMETERS.txt ---
    param_path = os.path.join(output_dir, '4_ANALYSIS_PARAMETERS.txt')
    with open(param_path, 'w', encoding='utf-8') as f:
        f.write("**************** ANALYSIS PARAMETERS DEFINITIONS *****************\n")
        f.write("WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA\n")
        f.write("AND THE NUMBER FORMAT\n\n")
        f.write("NUMBER OF GAUSS POINTS\n")
        f.write(f"{params['N_GAUSS']}\n\n")
        f.write("COLLOCATION POINT DISTANCING\n")
        f.write(f"{format_fortran_real(params['CP_DISTANCING'])}\n\n")
        f.write("PLANE STRAIN OR PLANE STRESS CONDITION\n")
        f.write(f"{params['PLANE_STATE']}\n\n")
        f.write("ANALYSIS TYPE: LINEAR ELASTIC (LE) OR COHESIVE PROPAGATION (CP)\n")
        f.write(f"{params['ANALYSIS_TYPE']}\n\n")
        f.write("NUMBER OF STEPS IN COHESIVE PROPAGATION ANALYSIS\n")
        f.write(f"{params['Ns']}\n\n")
        f.write("TOLERANCE FOR CONVERGENCE OF THE NEWTON-RAPHSON\n")
        f.write(f"{format_fortran_real(params['tol_Newton'])}\n\n")
        f.write("MAXIMUM NUMBER OF ITERACTIONS\n")
        f.write(f"{params['itemax']}\n\n")
        f.write("NODE FOR THE EVALUATION OF THE STRUCTURAL RESPONSE\n")
        f.write(f"{params['structural_node']}\n\n")
        f.write("NODE FOR THE EVALUATION OF THE COHESIVE RESPONSE\n")
        f.write(f"{params['cohesive_node']}\n\n")
        f.write("ELEMENT AND SURFACE TYPE \n")
        f.write("ELEMENT      SURFACE TYPE  (B = BOUNDARY Ip = INTERFACE plus Im = INTERFACE minus)    IE TYPE  (D = Displacement BIE    T = Traction BIE)   \t\n")
        for i in range(n_elem):
            tag = segments[i]
            cname = tag_to_curve.get(tag, '')
            stype = bound_cond.get(cname, {}).get('surface_type', 'B')
            f.write(f"  {i+1}\t\t   {stype}                                           \t\t                 D\n")

    print(f"[OK] Gerado: {param_path}")

    # --- 5. Escrever 5_POSTPROCESSING_MESH.txt ---
    post_path = os.path.join(output_dir, '5_POSTPROCESSING_MESH.txt')
    with open(post_path, 'w', encoding='utf-8') as f:
        f.write("**************** POSTPROCESSING MESH DEFINITIONS DEFINITIONS *****************\n")
        f.write("WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA\n")
        f.write("AND THE NUMBER FORMAT\n\n")
        f.write("NUMBER OF ELEMENTS\n0\n\n")
        f.write("NUMBER OF POINTS \n0\n\n")
        f.write("NUMBER OF BOUNDARY POINTS\n0\n\n")
        f.write("POINTS COORDINATES\n")
        f.write("II    \t X         Y         Z         SUBREGION         \n\n")
        f.write("INFORM THE CONECTIVITY FOR ALL POSTPROCESSING ELEMENTS\n")
        f.write("POINTS\n")

    print(f"[OK] Gerado: {post_path}")


# =============================================================================
# 5. EXECUÇÃO PRINCIPAL
# =============================================================================

if __name__ == '__main__':
    print("=================================================================")
    print("  Pré-processador Gmsh -> DBEM 2D Coesivo (meshio pipeline)")
    print("=================================================================")

    # 1. Obter dados de entrada e código .geo
    inp = input_data()

    # 2. Gerar malha com o Gmsh
    generate_geometry_and_mesh(inp['file_name'], inp['geo_code'])

    # 3. Ler arquivo .msh via meshio
    print(f"\n[meshio] Lendo malha '{inp['file_name']}.msh'...")
    comp_data = compute_inodes(inp['file_name'], inp['bound_cond'])

    print(f"[meshio] Nós extraídos: {len(comp_data['coordinates'])}")
    print(f"[meshio] Elementos 1D de contorno: {len(comp_data['line_elements'])}")
    print(f"[meshio] Sub-regiões identificadas: {comp_data['n_subregions']}")

    # 4. Escrever arquivos de entrada para o código Fortran
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    print(f"\n[Export] Escrevendo arquivos de entrada em '{os.path.abspath(output_dir)}'...")
    create_fortran_input_files(output_dir, comp_data, inp)

    print("\n[Sucesso] Todos os 5 arquivos de entrada foram gerados com sucesso!")
    print("          Você já pode rodar o solver Fortran executando '.\\MAIN_COHESIVE_DBEM_2D.exe'")
