"""
animate_debonding.py
--------------------
Gera uma animação em vídeo (.mp4) sincronizando:
  1. A configuração deformada da malha de elementos de contorno (DBEM 2D) durante
     o descolamento progressivo assimétrico (Substrato Omega_1 e Adesivo Omega_2
     em escala real 1:1, nos moldes da Fig. 28 do artigo EABE 2025);
  2. A evolução em tempo real da curva Força x Deslocamento Central
     (Fig. 26 - Force-displacement results for asymmetric debonding: hs = 1.08 mm).

Uso via linha de comando:
  python animate_debonding.py [pasta_simulacao] [arquivo_saida.mp4]
"""

import os
import sys
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Polygon
from matplotlib.collections import LineCollection

# Garante que o diretório do módulo plot_mesh_bc esteja no sys.path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from plot_mesh_bc import read_mesh_and_bcs


def configure_ffmpeg():
    """Configura o caminho do executável ffmpeg (do sistema ou via imageio-ffmpeg)."""
    sys_ffmpeg = shutil.which("ffmpeg")
    if sys_ffmpeg:
        mpl.rcParams["animation.ffmpeg_path"] = sys_ffmpeg
        return sys_ffmpeg

    try:
        import imageio_ffmpeg
        exe = imageio_ffmpeg.get_ffmpeg_exe()
        mpl.rcParams["animation.ffmpeg_path"] = exe
        return exe
    except ImportError:
        raise RuntimeError(
            "FFmpeg não encontrado. Instale o pacote 'imageio-ffmpeg' (pip install imageio-ffmpeg) "
            "ou adicione o ffmpeg ao PATH do sistema."
        )


def load_simulation_steps(sim_dir):
    """
    Carrega a malha, a curva Força-Deslocamento, a evolução de rigidez coesiva
    e os deslocamentos de contorno em cada passo de carga.
    """
    mesh_data = read_mesh_and_bcs(sim_dir)
    n_elem = mesh_data["n_elem"]
    n_nodes_raw = mesh_data["n_nodes_raw"]
    n_nodes_total = 2 * n_nodes_raw

    # 1. Resposta estrutural nodal (Força x Deslocamento Central)
    resp_file = os.path.join(sim_dir, "8_OUTPUT_NODAL_STRUCTURAL_RESPONSE.txt")
    if not os.path.exists(resp_file):
        raise FileNotFoundError(f"Arquivo não encontrado: {resp_file}")

    s_data = np.loadtxt(resp_file, skiprows=1)
    n_steps = len(s_data)
    u_mm = np.concatenate([[0.0], -s_data[:, 1]])
    # Área carregada: w_load (2.0 mm) * b (5.0 mm) = 10.0 mm^2
    F_N = np.concatenate([[0.0], -10.0 * s_data[:, 3]])

    # 2. Rigidez coesiva em cada passo (10_OUTPUT_COHESIVE_STIFFNESS.txt)
    stiff_file = os.path.join(sim_dir, "10_OUTPUT_COHESIVE_STIFFNESS.txt")
    if os.path.exists(stiff_file):
        stiff_raw = np.loadtxt(stiff_file, skiprows=1)
        steps_knn = stiff_raw[:, 1].reshape(n_steps, n_nodes_total)
        # Adiciona o passo 0 (interface 100% íntegra)
        knn_0 = steps_knn[0].copy()
        knn_0[knn_0 > 0] = knn_0.max()
        steps_knn = np.vstack([knn_0[None, :], steps_knn])
    else:
        steps_knn = np.ones((n_steps + 1, n_nodes_total)) * 2.0e6

    # 3. Coordenadas iniciais dos vértices dos elementos (em mm)
    # No Fortran (INPUT_BOUNDARY_MESH), para cada elemento J:
    #   CONECTIVITY(J, 1) = n_nodes_raw + raw_n1 (ponto de colocação xi = -1/3)
    #   CONECTIVITY(J, 2) = raw_n2               (ponto de colocação xi = +1/3)
    raw_coords = mesh_data["raw_coords"]
    elem_conec = mesh_data["elem_conec"]
    elem_subreg = mesh_data["elem_subregion"]
    surf_types = mesh_data["surf_types"]

    x0_elem = np.zeros((n_elem, 2))
    y0_elem = np.zeros((n_elem, 2))
    fn1_idx = np.zeros(n_elem, dtype=int)
    fn2_idx = np.zeros(n_elem, dtype=int)

    for eid in range(1, n_elem + 1):
        n1_raw, n2_raw = elem_conec[eid][0], elem_conec[eid][1]
        x0_elem[eid - 1, 0] = raw_coords[n1_raw][0] * 1e3
        y0_elem[eid - 1, 0] = raw_coords[n1_raw][1] * 1e3
        x0_elem[eid - 1, 1] = raw_coords[n2_raw][0] * 1e3
        y0_elem[eid - 1, 1] = raw_coords[n2_raw][1] * 1e3
        fn1_idx[eid - 1] = n_nodes_raw + n1_raw - 1
        fn2_idx[eid - 1] = n2_raw - 1

    # 4. Deslocamentos de todos os elementos em cada passo (em mm)
    steps_disp_file = os.path.join(sim_dir, "12_OUTPUT_NODAL_DISPLACEMENTS_STEPS.txt")
    ux_steps = np.zeros((n_steps + 1, n_elem, 2))
    uy_steps = np.zeros((n_steps + 1, n_elem, 2))

    if os.path.exists(steps_disp_file):
        disp_raw = np.loadtxt(steps_disp_file, skiprows=1)
        u_colloc_x = disp_raw[:, 2].reshape(n_steps, n_nodes_total) * 1e3
        u_colloc_y = disp_raw[:, 3].reshape(n_steps, n_nodes_total) * 1e3

        # Extrapolação linear idêntica à sub-rotina BOUNDARY_OUTPUT (MODULE_POSTPROCESSING.f90):
        #   xi_colloc = (-1/3, +1/3) -> xi_vertices = (-1, +1)
        #   u(xi=-1) =  2.0 * u(xi=-1/3) - 1.0 * u(xi=+1/3)
        #   u(xi=+1) = -1.0 * u(xi=-1/3) + 2.0 * u(xi=+1/3)
        for s in range(1, n_steps + 1):
            ux1 = u_colloc_x[s - 1, fn1_idx]
            ux2 = u_colloc_x[s - 1, fn2_idx]
            uy1 = u_colloc_y[s - 1, fn1_idx]
            uy2 = u_colloc_y[s - 1, fn2_idx]

            ux_steps[s, :, 0] = 2.0 * ux1 - 1.0 * ux2
            ux_steps[s, :, 1] = -1.0 * ux1 + 2.0 * ux2
            uy_steps[s, :, 0] = 2.0 * uy1 - 1.0 * uy2
            uy_steps[s, :, 1] = -1.0 * uy1 + 2.0 * uy2
    else:
        # Fallback caso 12_OUTPUT_NODAL_DISPLACEMENTS_STEPS.txt ainda não exista:
        # utiliza 11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt escalonado pelo deslocamento
        b_file = os.path.join(sim_dir, "11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt")
        b_resp = np.loadtxt(b_file, skiprows=1, usecols=(0, 1, 2, 3))
        ux_final = b_resp[:, 2] * 1e3
        uy_final = b_resp[:, 3] * 1e3
        for s in range(1, n_steps + 1):
            ratio = u_mm[s] / max(u_mm[-1], 1e-12)
            ux_steps[s, :, 0] = ratio * ux_final[fn1_idx]
            ux_steps[s, :, 1] = ratio * ux_final[fn2_idx]
            uy_steps[s, :, 0] = ratio * uy_final[fn1_idx]
            uy_steps[s, :, 1] = ratio * uy_final[fn2_idx]

    # 5. Comprimento da trinca assimétrica a(u) em cada passo (em mm)
    colloc_file = os.path.join(sim_dir, "6_OUTPUT_COLLOCATION_POINTS_COORDINATES.txt")
    crack_len_mm = np.zeros(n_steps + 1)
    if os.path.exists(colloc_file):
        fcoords = np.loadtxt(colloc_file)
        x_colloc_mm = fcoords[:, 0] * 1e3
        for s in range(1, n_steps + 1):
            st_k = steps_knn[s]
            broken = np.where((st_k > 0) & (st_k < 1.0e4))[0]
            if len(broken) > 0:
                crack_len_mm[s] = max(0.0, x_colloc_mm[broken].max() - (-12.5))

    return {
        "n_steps": n_steps,
        "n_elem": n_elem,
        "u_mm": u_mm,
        "F_N": F_N,
        "crack_len_mm": crack_len_mm,
        "x0_elem": x0_elem,
        "y0_elem": y0_elem,
        "ux_steps": ux_steps,
        "uy_steps": uy_steps,
        "steps_knn": steps_knn,
        "fn1_idx": fn1_idx,
        "fn2_idx": fn2_idx,
        "elem_subreg": elem_subreg,
        "surf_types": surf_types,
    }


def build_interpolated_timeline(sim_data, substeps_elastic=2, substeps_post=3, pause_peak=8, pause_end=15):
    """
    Cria uma linha do tempo suavizada por interpolação linear entre passos de carga
    consecutivos, com breves pausas no início, no pico de força e no final.
    """
    n_steps = sim_data["n_steps"]
    peak_step = int(np.argmax(sim_data["F_N"]))

    timeline = []
    # Pausa curta no estado inicial (u = 0)
    for _ in range(5):
        timeline.append((0, 0, 0.0))

    for s in range(n_steps):
        # Se o passo seguinte ultrapassar 1.60 mm (limite do eixo da Fig. 26), encerra exatamente em u = 1.600 mm
        u0, u1 = sim_data["u_mm"][s], sim_data["u_mm"][s + 1]
        if u1 > 1.6001 and u0 < 1.600:
            alpha_max = (1.600 - u0) / (u1 - u0)
            for k in range(1, substeps_post + 1):
                alpha = (k / float(substeps_post)) * alpha_max
                timeline.append((s, s + 1, alpha))
            break

        n_sub = substeps_elastic if (s + 1) <= peak_step else substeps_post
        # No salto brusco logo após o pico, adiciona mais quadros intermediários para visualizar a abertura
        if s == peak_step:
            n_sub = 8
        for k in range(1, n_sub + 1):
            alpha = k / float(n_sub)
            timeline.append((s, s + 1, alpha))
        if (s + 1) == peak_step:
            for _ in range(pause_peak):
                timeline.append((peak_step, peak_step, 0.0))

    # Pausa no último quadro
    last_frame = timeline[-1]
    for _ in range(pause_end):
        timeline.append(last_frame)

    return timeline


def create_debonding_animation(sim_dir="sim_fig26_nlinear", output_mp4="debonding_animation_hs108.mp4",
                               scale_factor=1.0, fps=20):
    """
    Gera e salva o arquivo .mp4 com a malha se deformando em sincronia com a curva
    Força x Deslocamento Central (Fig. 26 & Fig. 28 do artigo EABE 2025).
    """
    ffmpeg_exe = configure_ffmpeg()
    print(f"[*] Usando FFmpeg: {ffmpeg_exe}")
    print(f"[*] Carregando resultados não lineares de '{sim_dir}'...")
    sim = load_simulation_steps(sim_dir)

    n_elem = sim["n_elem"]
    u_mm = sim["u_mm"]
    F_N = sim["F_N"]
    crack_len_mm = sim["crack_len_mm"]
    x0 = sim["x0_elem"]
    y0 = sim["y0_elem"]
    ux_s = sim["ux_steps"]
    uy_s = sim["uy_steps"]
    steps_knn = sim["steps_knn"]
    fn1_idx = sim["fn1_idx"]
    fn2_idx = sim["fn2_idx"]
    elem_subreg = sim["elem_subreg"]
    surf_types = sim["surf_types"]

    sub1_indices = [i for i in range(n_elem) if elem_subreg[i + 1] == 1]
    sub2_indices = [i for i in range(n_elem) if elem_subreg[i + 1] == 2]
    int_indices = [i for i in range(n_elem) if surf_types.get(i + 1, "B") in ("Im", "Ip")]

    timeline = build_interpolated_timeline(sim)
    print(f"[*] Total de quadros da animação: {len(timeline)} ({len(timeline)/fps:.1f}s a {fps} fps)")

    # Configuração da figura com 2 painéis (Superior: Malha Deformando | Inferior: Fig. 26 Força x Deslocamento)
    fig = plt.figure(figsize=(12.0, 9.2), dpi=160, facecolor="white")
    gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.25], hspace=0.30)
    ax_mesh = fig.add_subplot(gs[0])
    ax_plot = fig.add_subplot(gs[1])

    # =========================================================================
    # PAINEL SUPERIOR: Malha DBEM 2D se deformando (Configuração Cn vs C0)
    # =========================================================================
    ax_mesh.set_xlim(-27.5, 27.5)
    ax_mesh.set_ylim(-10.8, 5.4)
    ax_mesh.set_aspect("equal")
    ax_mesh.grid(True, linestyle="--", alpha=0.35)
    ax_mesh.set_xlabel("Coordenada $x$ (mm)", fontsize=11, fontfamily="serif")
    ax_mesh.set_ylabel("Coordenada $y$ (mm)", fontsize=11, fontfamily="serif")
    ax_mesh.set_title(
        r"Descolamento Assimétrico no Ensaio 3PBT — Configuração Deformada $C_n$ ($h_s = 1.08$ mm, Escala 1:1)",
        fontsize=12, fontweight="bold", pad=10
    )

    # Contorno indeformado de referência C0 (tracejado cinza claro)
    c0_segs = [np.column_stack([x0[i], y0[i]]) for i in range(n_elem)]
    lc_c0 = LineCollection(c0_segs, colors="#bbbbbb", linewidths=0.9, linestyles="--", zorder=1, label=r"Configuração Inicial $C_0$")
    ax_mesh.add_collection(lc_c0)

    # Apoios fixo (direita) e móvel/rolete (esquerda) em y = 0, x = +-17.95 mm
    ax_mesh.plot([-17.95], [-0.35], marker="^", color="#d62728", markersize=8, zorder=6, label=r"Apoios ($x = \pm 17.95$ mm)")
    ax_mesh.plot([17.95], [-0.35], marker="^", color="#d62728", markersize=8, zorder=6)

    # Polígonos preenchidos para Substrato (Omega_1) e Adesivo (Omega_2)
    poly_sub1 = Polygon(np.column_stack([x0[sub1_indices, 0], y0[sub1_indices, 0]]),
                        closed=True, facecolor="#1f77b4", alpha=0.18, edgecolor="none", zorder=2,
                        label=r"Substrato $\Omega_1$ (Alumínio 2024-T3)")
    poly_sub2 = Polygon(np.column_stack([x0[sub2_indices, 0], y0[sub2_indices, 0]]),
                        closed=True, facecolor="#2ca02c", alpha=0.22, edgecolor="none", zorder=2,
                        label=r"Adesivo $\Omega_2$ (DGEBA/DETA)")
    ax_mesh.add_patch(poly_sub1)
    ax_mesh.add_patch(poly_sub2)

    # Coleção de segmentos dos elementos de contorno deformados
    lc_elems = LineCollection(c0_segs, colors="black", linewidths=1.35, zorder=4)
    ax_mesh.add_collection(lc_elems)

    # Coleção destacando os elementos da interface coesiva (vermelho onde descolou, verde onde íntegro)
    int_segs_init = [np.column_stack([x0[i], y0[i]]) for i in int_indices]
    lc_interface = LineCollection(int_segs_init, colors="#2ca02c", linewidths=2.2, zorder=5,
                                  label="Interface Coesiva (Verde: Íntegra | Vermelho: Trinca)")
    ax_mesh.add_collection(lc_interface)

    # Nós dos elementos de contorno (pontos pretos como na Fig. 28 do artigo)
    nodes_scatter, = ax_mesh.plot(x0[:, 0], y0[:, 0], "ko", markersize=1.9, zorder=6)

    # Seta dinâmica indicando o rolete central aplicando o deslocamento u_y
    load_quiver = ax_mesh.quiver([0.0], [1.08 + 2.2], [0.0], [-1.9], color="blue",
                                 angles="xy", scale_units="xy", scale=1, width=0.0045,
                                 headwidth=4.2, headlength=5.0, zorder=8)

    # Caixa de informações em tempo real no painel da malha (canto inferior esquerdo, 3 linhas compactas)
    info_text = ax_mesh.text(
        0.012, 0.045, "", transform=ax_mesh.transAxes, fontsize=8.4, fontweight="bold",
        color="#111111", bbox=dict(boxstyle="round,pad=0.30", facecolor="#f8f9fa", edgecolor="#888888", alpha=0.92),
        zorder=10
    )
    ax_mesh.legend(loc="lower right", fontsize=7.8, ncol=2, frameon=True, facecolor="white", framealpha=0.92)

    # =========================================================================
    # PAINEL INFERIOR: Gráfico Força x Deslocamento Central (Fig. 26)
    # =========================================================================
    ax_plot.set_xlim(0.0, 1.6)
    ax_plot.set_ylim(0.0, 205.0)
    ax_plot.set_xticks([0.0, 0.4, 0.8, 1.2, 1.6])
    ax_plot.set_yticks([0, 50, 100, 150, 200])
    ax_plot.grid(True, linestyle="-", color="#cccccc", linewidth=0.7)
    ax_plot.set_xlabel("Central displacement (mm)", fontsize=11.5, fontfamily="serif")
    ax_plot.set_ylabel("Force (N)", fontsize=11.5, fontfamily="serif")
    ax_plot.set_title(
        r"Fig. 26. Force–displacement results for asymmetric debonding: $h_s = 1.08$ mm",
        fontsize=12, fontweight="bold", pad=10
    )

    # Curva completa de fundo (tracejada suave como guia visual)
    ax_plot.plot(u_mm, F_N, "--", color="#b0c4de", linewidth=1.3, zorder=2, label="Trajetória Completa (2D BEM - PID)")

    # Curva progressiva azul com marcadores circulares abertos (idêntica à Fig. 26)
    line_pid_curve, = ax_plot.plot([], [], "-", color="blue", linewidth=2.0, zorder=4)
    line_pid_markers, = ax_plot.plot([], [], "o", color="blue", markersize=5.2,
                                     markerfacecolor="none", markeredgecolor="blue",
                                     markeredgewidth=1.25, zorder=5, label="2D BEM - PID")

    # Marcador de destaque do estado atual (ponto vermelho/dourado sincronizado)
    current_pt, = ax_plot.plot([], [], "o", color="#d62728", markersize=8.5,
                               markerfacecolor="#ffdd57", markeredgecolor="#d62728",
                               markeredgewidth=2.0, zorder=7, label="Estado Atual")

    # Linha vertical guia sincronizando o deslocamento atual no gráfico
    vline_u = ax_plot.axvline(0.0, color="#d62728", linestyle=":", linewidth=1.2, alpha=0.7, zorder=3)

    # Anotação móvel (u, F) junto ao cursor atual
    cursor_annot = ax_plot.annotate(
        "", xy=(0, 0), xytext=(12, -18), textcoords="offset points",
        fontsize=9.5, fontweight="bold", color="#b2182b",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff8dc", edgecolor="#d62728", alpha=0.92),
        zorder=8
    )

    ax_plot.legend(loc="lower right", frameon=True, edgecolor="black", fancybox=False, fontsize=9.5)

    # =========================================================================
    # FUNÇÃO DE ATUALIZAÇÃO QUADRO A QUADRO
    # =========================================================================
    def update_frame(frame_idx):
        s0, s1, alpha = timeline[frame_idx]

        # Interpolação linear dos deslocamentos, força e comprimento de trinca
        u_curr = (1.0 - alpha) * u_mm[s0] + alpha * u_mm[s1]
        F_curr = (1.0 - alpha) * F_N[s0] + alpha * F_N[s1]
        a_curr = (1.0 - alpha) * crack_len_mm[s0] + alpha * crack_len_mm[s1]

        ux_curr = (1.0 - alpha) * ux_s[s0] + alpha * ux_s[s1]
        uy_curr = (1.0 - alpha) * uy_s[s0] + alpha * uy_s[s1]
        knn_curr = (1.0 - alpha) * steps_knn[s0] + alpha * steps_knn[s1]

        x_def = x0 + scale_factor * ux_curr
        y_def = y0 + scale_factor * uy_curr

        # 1. Atualiza segmentos de todos os elementos de contorno
        segs = [np.column_stack([x_def[i], y_def[i]]) for i in range(n_elem)]
        lc_elems.set_segments(segs)

        # 2. Atualiza segmentos e cores da interface coesiva (verde = íntegro, laranja = dano, vermelho = trinca)
        int_segs = []
        int_colors = []
        for i in int_indices:
            int_segs.append(np.column_stack([x_def[i], y_def[i]]))
            k_avg = 0.5 * (knn_curr[fn1_idx[i]] + knn_curr[fn2_idx[i]])
            if k_avg < 1.0e4:
                int_colors.append("#d62728")  # Trinca aberta (vermelho)
            elif k_avg < 1.95e6:
                int_colors.append("#ff7f0e")  # Zona de processo / amolecimento (laranja)
            else:
                int_colors.append("#2ca02c")  # Interface íntegra (verde)
        lc_interface.set_segments(int_segs)
        lc_interface.set_color(int_colors)

        # 3. Atualiza polígonos sombreados das sub-regiões Omega_1 e Omega_2
        poly_sub1.set_xy(np.column_stack([x_def[sub1_indices, 0], y_def[sub1_indices, 0]]))
        poly_sub2.set_xy(np.column_stack([x_def[sub2_indices, 0], y_def[sub2_indices, 0]]))

        # 4. Atualiza posição dos nós de contorno e da seta de carga central
        nodes_scatter.set_data(x_def[:, 0], y_def[:, 0])
        top_center_y = 1.08 - scale_factor * u_curr
        load_quiver.set_offsets(np.array([[0.0, top_center_y + 2.1]]))

        # 5. Atualiza texto de diagnóstico na malha
        step_display = s1 if alpha > 0.5 else s0
        regime_str = "Elástico Íntegro" if step_display <= 13 else ("Nucleação / Salto Brusco" if step_display == 14 else "Propagação Coesiva (PID)")
        info_text.set_text(
            f"Passo: {step_display:02d}/{sim['n_steps']} ({regime_str})\n"
            f"u = {u_curr:.3f} mm  |  F = {F_curr:.1f} N\n"
            f"Trinca assimétrica: a = {a_curr:.2f} mm"
        )

        # 6. Atualiza curva Força x Deslocamento no painel inferior
        u_hist = np.append(u_mm[: s0 + 1], u_curr)
        F_hist = np.append(F_N[: s0 + 1], F_curr)
        line_pid_curve.set_data(u_hist, F_hist)

        done_idx = s1 if alpha > 0.99 else s0
        line_pid_markers.set_data(u_mm[: done_idx + 1], F_N[: done_idx + 1])

        current_pt.set_data([u_curr], [F_curr])
        vline_u.set_xdata([u_curr, u_curr])

        cursor_annot.xy = (u_curr, F_curr)
        offset_x = -95 if u_curr > 1.22 else 12
        offset_y = -26 if F_curr > 165 else 10
        cursor_annot.xyann = (offset_x, offset_y)
        cursor_annot.set_text(f"u = {u_curr:.3f} mm\nF = {F_curr:.1f} N")

        return (lc_elems, lc_interface, poly_sub1, poly_sub2, nodes_scatter,
                load_quiver, info_text, line_pid_curve, line_pid_markers,
                current_pt, vline_u, cursor_annot)

    anim = animation.FuncAnimation(
        fig, update_frame, frames=len(timeline), interval=1000.0 / fps, blit=False
    )

    writer = animation.FFMpegWriter(
        fps=fps, metadata=dict(title="DBEM 2D Asymmetric Cohesive Debonding", artist="Antigravity DBEM 2D"),
        bitrate=3200, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"]
    )

    print(f"[*] Renderizando e codificando vídeo MP4 em '{output_mp4}'...")
    anim.save(output_mp4, writer=writer)
    plt.close(fig)
    print(f"[OK] Animação MP4 gerada com sucesso: {os.path.abspath(output_mp4)}")
    return output_mp4


if __name__ == "__main__":
    default_sim_dir = os.path.join(SCRIPT_DIR, "sim_fig26_nlinear")
    default_out_mp4 = os.path.join(SCRIPT_DIR, "debonding_animation_hs108.mp4")

    sim_folder = sys.argv[1] if len(sys.argv) > 1 else default_sim_dir
    out_video = sys.argv[2] if len(sys.argv) > 2 else default_out_mp4

    create_debonding_animation(sim_folder, out_video)
