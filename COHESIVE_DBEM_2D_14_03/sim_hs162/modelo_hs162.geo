SetFactory("OpenCASCADE");

// Parâmetros geométricos (em metros)
hs = 0.00162;
ha = 4.0e-3;
L_sub = 50.0e-3;
L_adh = 25.0e-3;
d_supp = 35.90e-3;
w_supp = 0.12e-3;
w_load = 2.0e-3;

x_sub = L_sub / 2.0;
x_adh = L_adh / 2.0;
x_s_c = d_supp / 2.0;
x_s1 = x_s_c - w_supp / 2.0;
x_s2 = x_s_c + w_supp / 2.0;
x_l = w_load / 2.0;

// Tamanhos de elemento característicos
lc_int = 0.00025;
lc_ext = 1.0e-3;
lc_supp = 0.06e-3;
lc_load = 0.5e-3;

// --- Sub-região 1: Substrato (y = 0 até hs) ---
p1  = newp; Point(p1)  = {-x_sub, 0, 0, lc_ext};
p2  = newp; Point(p2)  = {-x_s2,  0, 0, lc_supp};
p3  = newp; Point(p3)  = {-x_s1,  0, 0, lc_supp};
p4  = newp; Point(p4)  = {-x_adh, 0, 0, lc_int};
p5  = newp; Point(p5)  = { x_adh, 0, 0, lc_int};
p6  = newp; Point(p6)  = { x_s1,  0, 0, lc_supp};
p7  = newp; Point(p7)  = { x_s2,  0, 0, lc_supp};
p8  = newp; Point(p8)  = { x_sub, 0, 0, lc_ext};

p9  = newp; Point(p9)  = { x_sub, hs, 0, lc_ext};
p10 = newp; Point(p10) = { x_l,   hs, 0, lc_load};
p11 = newp; Point(p11) = {-x_l,   hs, 0, lc_load};
p12 = newp; Point(p12) = {-x_sub, hs, 0, lc_ext};

l_sub_bot1 = newl; Line(l_sub_bot1) = {p1, p2};
l_supp_esq = newl; Line(l_supp_esq) = {p2, p3};
l_sub_bot2 = newl; Line(l_sub_bot2) = {p3, p4};
l_int_m    = newl; Line(l_int_m)    = {p4, p5}; // Interface Im (Substrato)
l_sub_bot3 = newl; Line(l_sub_bot3) = {p5, p6};
l_supp_dir = newl; Line(l_supp_dir) = {p6, p7};
l_sub_bot4 = newl; Line(l_sub_bot4) = {p7, p8};

l_sub_dir  = newl; Line(l_sub_dir)  = {p8, p9};
l_sub_top1 = newl; Line(l_sub_top1) = {p9, p10};
l_load     = newl; Line(l_load)     = {p10, p11}; // Carga aplicada
l_sub_top2 = newl; Line(l_sub_top2) = {p11, p12};
l_sub_esq  = newl; Line(l_sub_esq)  = {p12, p1};

Curve Loop(1) = {l_sub_bot1, l_supp_esq, l_sub_bot2, l_int_m, l_sub_bot3, l_supp_dir, l_sub_bot4, l_sub_dir, l_sub_top1, l_load, l_sub_top2, l_sub_esq};
Plane Surface(1) = {1};

// --- Sub-região 2: Adesivo (y = -ha até 0) ---
p13 = newp; Point(p13) = { x_adh, 0, 0, lc_int};
p14 = newp; Point(p14) = {-x_adh, 0, 0, lc_int};
p15 = newp; Point(p15) = {-x_adh, -ha, 0, lc_ext};
p16 = newp; Point(p16) = { x_adh, -ha, 0, lc_ext};

l_int_p    = newl; Line(l_int_p)    = {p13, p14}; // Interface Ip (Adesivo)
l_adh_esq  = newl; Line(l_adh_esq)  = {p14, p15};
l_adh_bot  = newl; Line(l_adh_bot)  = {p15, p16};
l_adh_dir  = newl; Line(l_adh_dir)  = {p16, p13};

Curve Loop(2) = {l_int_p, l_adh_esq, l_adh_bot, l_adh_dir};
Plane Surface(2) = {2};

// Discretização transfinita ao longo da interface
n_int = Ceil(L_adh / lc_int) + 1;
Transfinite Curve {l_int_m, l_int_p} = n_int Using Progression 1;

// Grupos Físicos de Contorno (1D)
Physical Curve("contorno_livre_sub") = {l_sub_bot1, l_sub_bot2, l_sub_bot3, l_sub_bot4, l_sub_dir, l_sub_top1, l_sub_top2, l_sub_esq};
Physical Curve("apoio_esq")           = {l_supp_esq};
Physical Curve("apoio_dir")           = {l_supp_dir};
Physical Curve("carga")               = {l_load};
Physical Curve("interface_m")         = {l_int_m};
Physical Curve("interface_p")         = {l_int_p};
Physical Curve("contorno_livre_adh") = {l_adh_esq, l_adh_bot, l_adh_dir};

// Grupos Físicos de Superfície (2D)
Physical Surface("subregion_1") = {1};
Physical Surface("subregion_2") = {2};
