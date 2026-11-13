// ==============================================================================
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
