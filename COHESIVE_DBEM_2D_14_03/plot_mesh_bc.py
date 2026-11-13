"""
plot_mesh_bc.py
---------------
Script de pós/pré-processamento para visualização da malha de elementos de
contorno (DBEM 2D), sub-regiões, interface coesiva e condições de contorno
(Dirichlet e Neumann), baseado na referência 'preproc.py' da biblioteca
BEM2D_elastostatic.

Suporta leitura automática de:
  - 1_BOUNDARY_MESH.txt
  - 2_MATERIALS.txt
  - 3_BOUNDARY_CONDITIONS.txt
  - 4_ANALYSIS_PARAMETERS.txt
  - 6_OUTPUT_COLLOCATION_POINTS_COORDINATES.txt (se disponível)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import LineCollection


def parse_fortran_float(val_str):
    """Converte representações numéricas Fortran (ex: 1.D-2, 0.D0) para float."""
    clean = val_str.strip().upper().replace('D', 'E')
    return float(clean)


def read_mesh_and_bcs(dir_path='.'):
    """Lê todos os arquivos de definição do modelo DBEM Coesivo 2D."""
    mesh_file = os.path.join(dir_path, '1_BOUNDARY_MESH.txt')
    bc_file = os.path.join(dir_path, '3_BOUNDARY_CONDITIONS.txt')
    param_file = os.path.join(dir_path, '4_ANALYSIS_PARAMETERS.txt')
    colloc_file = os.path.join(dir_path, '6_OUTPUT_COLLOCATION_POINTS_COORDINATES.txt')

    if not os.path.exists(mesh_file):
        raise FileNotFoundError(f"Arquivo de malha não encontrado: {mesh_file}")

    # --- 1. Leitura de 1_BOUNDARY_MESH.txt ---
    with open(mesh_file, 'r', encoding='utf-8', errors='ignore') as f:
        mesh_lines = [l.strip() for l in f if l.strip()]

    idx = 0
    while idx < len(mesh_lines) and 'NUMBER OF SUBREGIONS' not in mesh_lines[idx]:
        idx += 1
    idx += 1
    n_subregions = int(mesh_lines[idx])

    while idx < len(mesh_lines) and 'NUMBER OF BOUNDARY ELEMENTS' not in mesh_lines[idx]:
        idx += 1
    idx += 1
    n_elem = int(mesh_lines[idx])

    while idx < len(mesh_lines) and 'NUMBER OF NODES' not in mesh_lines[idx]:
        idx += 1
    idx += 1
    n_nodes_raw = int(mesh_lines[idx])

    while idx < len(mesh_lines) and 'NUMBER OF DUPLICATE NODES' not in mesh_lines[idx]:
        idx += 1
    idx += 1
    n_dup = int(mesh_lines[idx])

    while idx < len(mesh_lines) and 'DUPLICATE NODES' not in mesh_lines[idx]:
        idx += 1
    idx += 2
    dup_nodes = []
    for _ in range(n_dup):
        dup_nodes.append(int(mesh_lines[idx]))
        idx += 1

    while idx < len(mesh_lines) and 'NODAL COORDINATES' not in mesh_lines[idx]:
        idx += 1
    idx += 2
    raw_coords = {}
    for _ in range(n_nodes_raw):
        parts = mesh_lines[idx].split()
        nid = int(parts[0])
        x = parse_fortran_float(parts[1])
        y = parse_fortran_float(parts[2])
        raw_coords[nid] = (x, y)
        idx += 1

    # Construção de coordenadas com nós duplicados (replicando o algoritmo do módulo Fortran)
    node_coords = dict(raw_coords)
    cont = 0
    for nid in range(1, n_nodes_raw + 1):
        if nid in dup_nodes:
            cont += 1
            new_id = n_nodes_raw + cont
            node_coords[new_id] = raw_coords[nid]

    while idx < len(mesh_lines) and 'INFORM THE POLYNOMIAL ORDER' not in mesh_lines[idx]:
        idx += 1
    idx += 2
    elem_order = {}
    elem_subregion = {}
    for _ in range(n_elem):
        parts = mesh_lines[idx].split()
        eid = int(parts[0])
        order = int(parts[1])
        subreg = int(parts[2])
        elem_order[eid] = order
        elem_subregion[eid] = subreg
        idx += 1

    while idx < len(mesh_lines) and 'INFORM THE CONECTIVITY' not in mesh_lines[idx]:
        idx += 1
    idx += 2
    elem_conec = {}
    for eid in range(1, n_elem + 1):
        parts = [int(p) for p in mesh_lines[idx].split()]
        elem_conec[eid] = parts
        idx += 1

    # --- 2. Leitura de 4_ANALYSIS_PARAMETERS.txt ---
    surf_types = {}
    ie_types = {}
    analysis_type = "CP"
    control_info = {}
    if os.path.exists(param_file):
        with open(param_file, 'r', encoding='utf-8', errors='ignore') as f:
            param_lines = [l.strip() for l in f if l.strip()]

        pidx = 0
        while pidx < len(param_lines) and 'ANALYSIS TYPE:' not in param_lines[pidx]:
            pidx += 1
        if pidx < len(param_lines):
            pidx += 1
            analysis_type = param_lines[pidx]

        while pidx < len(param_lines) and 'ELEMENT AND SURFACE TYPE' not in param_lines[pidx]:
            pidx += 1
        if pidx < len(param_lines):
            pidx += 2
            for _ in range(n_elem):
                if pidx < len(param_lines):
                    parts = param_lines[pidx].split()
                    eid = int(parts[0])
                    stype = parts[1]
                    ietype = parts[2] if len(parts) > 2 else 'D'
                    surf_types[eid] = stype
                    ie_types[eid] = ietype
                    pidx += 1

    # --- 3. Leitura de 3_BOUNDARY_CONDITIONS.txt ---
    bc_ux, bc_uy, bc_tx, bc_ty = {}, {}, {}, {}
    if os.path.exists(bc_file):
        with open(bc_file, 'r', encoding='utf-8', errors='ignore') as f:
            bc_lines = [l.strip() for l in f if l.strip()]

        bidx = 0
        while bidx < len(bc_lines) and 'NUMBER OF ELEMENTS WITH Ux PRECRIBED' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        nu1 = int(bc_lines[bidx])

        while bidx < len(bc_lines) and 'NUMBER OF ELEMENTS WITH Uy PRECRIBED' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        nu2 = int(bc_lines[bidx])

        while bidx < len(bc_lines) and 'NUMBER OF ELEMENTS WITH Tx PRESCRIBED' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        nt1 = int(bc_lines[bidx])

        while bidx < len(bc_lines) and 'NUMBER OF ELEMENTS WITH Ty PRESCRIBED' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        nt2 = int(bc_lines[bidx])

        while bidx < len(bc_lines) and 'II' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        for _ in range(nu1):
            parts = bc_lines[bidx].split()
            eid = int(parts[0])
            bc_ux[eid] = (parse_fortran_float(parts[1]), parse_fortran_float(parts[2]), parse_fortran_float(parts[3]))
            bidx += 1

        while bidx < len(bc_lines) and 'Auy' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        for _ in range(nu2):
            parts = bc_lines[bidx].split()
            eid = int(parts[0])
            bc_uy[eid] = (parse_fortran_float(parts[1]), parse_fortran_float(parts[2]), parse_fortran_float(parts[3]))
            bidx += 1

        while bidx < len(bc_lines) and 'Atx' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        for _ in range(nt1):
            parts = bc_lines[bidx].split()
            eid = int(parts[0])
            bc_tx[eid] = (parse_fortran_float(parts[1]), parse_fortran_float(parts[2]), parse_fortran_float(parts[3]))
            bidx += 1

        while bidx < len(bc_lines) and 'Aty' not in bc_lines[bidx]:
            bidx += 1
        bidx += 1
        for _ in range(nt2):
            parts = bc_lines[bidx].split()
            eid = int(parts[0])
            bc_ty[eid] = (parse_fortran_float(parts[1]), parse_fortran_float(parts[2]), parse_fortran_float(parts[3]))
            bidx += 1

        while bidx < len(bc_lines) and 'CONTROLED PARAMETER' not in bc_lines[bidx]:
            bidx += 1
        if bidx < len(bc_lines):
            bidx += 2
            parts = bc_lines[bidx].split()
            control_info = {
                'type': parts[0],
                'node': int(parts[1]),
                'direction': int(parts[2]),
                'value': parse_fortran_float(parts[3])
            }

    # --- 4. Leitura opcional dos pontos de colocação gerados ---
    colloc_pts = []
    if os.path.exists(colloc_file):
        with open(colloc_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    colloc_pts.append((float(parts[0]), float(parts[1])))

    return {
        'n_subregions': n_subregions,
        'n_elem': n_elem,
        'n_nodes_raw': n_nodes_raw,
        'node_coords': node_coords,
        'raw_coords': raw_coords,
        'elem_order': elem_order,
        'elem_subregion': elem_subregion,
        'elem_conec': elem_conec,
        'surf_types': surf_types,
        'bc_ux': bc_ux,
        'bc_uy': bc_uy,
        'bc_tx': bc_tx,
        'bc_ty': bc_ty,
        'control_info': control_info,
        'colloc_pts': colloc_pts,
        'analysis_type': analysis_type,
    }


def find_subregion_polygon(subreg_id, data):
    """Gera a lista ordenada de vértices para uma sub-região."""
    elems = [eid for eid, s in data['elem_subregion'].items() if s == subreg_id]
    if not elems:
        return None
    poly_nodes = []
    for eid in elems:
        c = data['elem_conec'][eid]
        poly_nodes.append(c[0])
    pts = [data['node_coords'][n] for n in poly_nodes]
    return pts


def show_dbem_geometry_and_bcs(data, save_path='mesh_and_boundary_conditions.png'):
    """
    Desenha a malha e condições de contorno de acordo com o padrão preproc.py:
      - Linhas e nós de contorno
      - Destaque da interface coesiva (I- / I+)
      - Condições de Dirichlet (ux: 'rx', uy: 'rP')
      - Condições de Neumann (tx: 'bx', ty: 'bP')
      - Vetores de tração/deslocamento (quivers)
      - Ponto de controle de propagação
      - Pontos de colocação descontínuos
    """
    coords = data['node_coords']
    raw_coords = data['raw_coords']
    elem_conec = data['elem_conec']
    elem_subreg = data['elem_subregion']
    surf_types = data['surf_types']
    bc_ux = data['bc_ux']
    bc_uy = data['bc_uy']
    bc_tx = data['bc_tx']
    bc_ty = data['bc_ty']
    control_info = data.get('control_info', {})
    colloc_pts = data.get('colloc_pts', [])

    fig, ax = plt.subplots(figsize=(10, 8.5), dpi=300)
    ax.grid(True, linestyle='--', alpha=0.4)

    # Coordenadas globais para delimitação
    all_x = [raw_coords[n][0] for n in raw_coords]
    all_y = [raw_coords[n][1] for n in raw_coords]
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    dx = x_max - x_min if (x_max - x_min) > 0 else 1.0
    dy = y_max - y_min if (y_max - y_min) > 0 else 1.0
    max_delta = max(dx, dy)

    # Preenchimento suave de cada sub-região
    subreg_colors = {1: '#1f77b4', 2: '#2ca02c', 3: '#ff7f0e'}
    for s_id in range(1, data['n_subregions'] + 1):
        poly_pts = find_subregion_polygon(s_id, data)
        if poly_pts and len(poly_pts) >= 3:
            poly = Polygon(poly_pts, closed=True, facecolor=subreg_colors.get(s_id, '#888888'),
                           alpha=0.08, edgecolor='none', zorder=1)
            ax.add_patch(poly)
            cx = np.mean([p[0] for p in poly_pts])
            cy = np.mean([p[1] for p in poly_pts])
            ax.text(cx, cy, fr'Subregion $\Omega_{s_id}$',
                    fontsize=10.5, color=subreg_colors.get(s_id, '#333333'),
                    fontweight='bold', ha='center', va='center', zorder=2)

    # 1. Traçado dos elementos de contorno
    for eid, con in elem_conec.items():
        n1, n2 = con[0], con[1]
        x1, y1 = coords[n1]
        x2, y2 = coords[n2]
        stype = surf_types.get(eid, 'B')
        subreg = elem_subreg.get(eid, 1)

        if stype in ['Im', 'Ip']:
            # Interface Coesiva (vermelho tracejado)
            lbl = 'Cohesive Interface ($I^- / I^+$)' if eid == 4 else ""
            ax.plot([x1, x2], [y1, y2], color='#d62728', linestyle='--', linewidth=2.6,
                    label=lbl, zorder=5)
        else:
            col = subreg_colors.get(subreg, 'black')
            lbl = fr'Subregion $\Omega_{subreg}$ Boundary' if eid in [1, 7] else ""
            ax.plot([x1, x2], [y1, y2], color=col, linewidth=1.8, label=lbl, zorder=3)

        # Rótulo do elemento e1, e2, ...
        xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        de_x, de_y = (x2 - x1), (y2 - y1)
        L = np.hypot(de_x, de_y)
        nx = -de_y / L if L > 1e-12 else 0.0
        ny = de_x / L if L > 1e-12 else 1.0
        ax.text(xm + 0.03 * max_delta * nx, ym + 0.03 * max_delta * ny, f'e{eid}',
                color='#111111', fontsize=8.5, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='#ffffff', edgecolor='#bbbbbb', alpha=0.85),
                zorder=6)

    # 2. Nós geométricos
    for nid, (x, y) in raw_coords.items():
        ax.plot(x, y, 'ko', markersize=4, zorder=7)
        ax.text(x - 0.02 * max_delta, y - 0.02 * max_delta, f'{nid}',
                color='#222222', fontsize=8, ha='right', va='top', zorder=8)

    # 3. Pontos de colocação (se disponíveis)
    if colloc_pts:
        c_arr = np.array(colloc_pts)
        ax.plot(c_arr[:, 0], c_arr[:, 1], 'm.', markersize=4.5,
                label='Collocation Points (DBEM)', zorder=6)

    # 4. Condições de contorno (Dirichlet & Neumann)
    arrow_len = 0.08 * max_delta

    # Marcadores coletados
    known_ux_pts, known_uy_pts = [], []
    known_tx_pts, known_ty_pts = [], []

    # Ux prescrito
    for eid, (a, b, c) in bc_ux.items():
        for nid in elem_conec[eid]:
            pt = coords[nid]
            known_ux_pts.append(pt)
            val = c
            vx = -arrow_len if abs(val) < 1e-12 else (arrow_len * np.sign(val))
            ax.quiver(pt[0], pt[1], vx, 0.0, color='red', angles='xy', scale_units='xy',
                      scale=1, width=0.004, headwidth=4, headlength=5, zorder=9, alpha=0.8)

    # Uy prescrito
    for eid, (a, b, c) in bc_uy.items():
        for nid in elem_conec[eid]:
            pt = coords[nid]
            known_uy_pts.append(pt)
            val = c
            vy = -arrow_len if abs(val) < 1e-12 else (arrow_len * np.sign(val))
            ax.quiver(pt[0], pt[1], 0.0, vy, color='red', angles='xy', scale_units='xy',
                      scale=1, width=0.004, headwidth=4, headlength=5, zorder=9, alpha=0.8)

    # Tx prescrito (Neumann)
    for eid, (a, b, c) in bc_tx.items():
        val = c
        for nid in elem_conec[eid]:
            pt = coords[nid]
            known_tx_pts.append(pt)
            vx = arrow_len * (np.sign(val) if abs(val) > 1e-12 else 1.0)
            ax.quiver(pt[0], pt[1], vx, 0.0, color='blue', angles='xy', scale_units='xy',
                      scale=1, width=0.005, headwidth=4, headlength=5, zorder=9)

    # Ty prescrito (Neumann)
    for eid, (a, b, c) in bc_ty.items():
        val = c
        for nid in elem_conec[eid]:
            pt = coords[nid]
            known_ty_pts.append(pt)
            vy = arrow_len * (np.sign(val) if abs(val) > 1e-12 else 1.0)
            ax.quiver(pt[0], pt[1], 0.0, vy, color='blue', angles='xy', scale_units='xy',
                      scale=1, width=0.005, headwidth=4, headlength=5, zorder=9)

    # Plot dos marcadores característicos do preproc.py:
    # 'rx' para Ux conhecido, 'rP' para Uy conhecido
    # 'bx' para Tx conhecido, 'bP' para Ty conhecido
    if known_tx_pts:
        arr = np.array(known_tx_pts)
        ax.plot(arr[:, 0], arr[:, 1], 'bx', markersize=7, markeredgewidth=1.8,
                label='Known $t_x$ (Neumann)', zorder=10)
    if known_ty_pts:
        arr = np.array(known_ty_pts)
        ax.plot(arr[:, 0], arr[:, 1], 'bP', markersize=6,
                label='Known $t_y$ (Neumann)', zorder=10)
    if known_ux_pts:
        arr = np.array(known_ux_pts)
        ax.plot(arr[:, 0], arr[:, 1], 'rx', markersize=7, markeredgewidth=1.8,
                label='Known $u_x$ (Dirichlet)', zorder=11)
    if known_uy_pts:
        arr = np.array(known_uy_pts)
        ax.plot(arr[:, 0], arr[:, 1], 'rP', markersize=6,
                label='Known $u_y$ (Dirichlet)', zorder=11)

    # 5. Marcador para Ponto Controlado (Newton-Raphson control DOF)
    if control_info:
        cnid = control_info['node']
        cval = control_info['value']
        cdir = control_info['direction']
        cpt = coords.get(cnid, None)
        if cpt is not None:
            dir_str = "x" if cdir == 1 else "y"
            ax.plot(cpt[0], cpt[1], '*', color='#006400', markersize=14, zorder=12,
                    label=f'Control Node {cnid} ($u_{dir_str} = {cval}$)')
            ax.annotate(f'Control Node {cnid}\n$u_{dir_str} = {cval}$',
                        xy=cpt, xytext=(cpt[0] + 0.1 * max_delta, cpt[1] + 0.1 * max_delta),
                        arrowprops=dict(facecolor='#006400', shrink=0.08, width=1.4, headwidth=6),
                        fontsize=9, fontweight='bold', color='#006400',
                        bbox=dict(boxstyle='round,pad=0.25', facecolor='#eafaf1', edgecolor='#006400'),
                        zorder=13)

    # Ajustes finais da figura
    margin = 0.20 * max_delta
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin, y_max + margin)
    ax.set_aspect('equal')
    ax.set_xlabel('Coordenada $X$ (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Coordenada $Y$ (m)', fontsize=11, fontweight='bold')
    ax.set_title('Multidomain Cohesive DBEM 2D: Malha e Condições de Contorno',
                 fontsize=13, fontweight='bold', pad=14)

    # Legenda estruturada
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper center',
              bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5,
              frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[OK] Visualização gerada e salva com sucesso em '{save_path}'.")
    plt.close(fig)


if __name__ == '__main__':
    data_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    output_img = sys.argv[2] if len(sys.argv) > 2 else 'mesh_and_boundary_conditions.png'

    print(f"Lendo dados de entrada da pasta: {os.path.abspath(data_dir)}")
    parsed_data = read_mesh_and_bcs(data_dir)
    print(f"Sub-regiões: {parsed_data['n_subregions']} | Elementos: {parsed_data['n_elem']} | Nós: {parsed_data['n_nodes_raw']}")
    show_dbem_geometry_and_bcs(parsed_data, output_img)
