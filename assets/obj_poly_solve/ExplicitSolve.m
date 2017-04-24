(* ::Package:: *)
\[Alpha] = ( {
    {\[Alpha]1},
    {\[Alpha]2},
    {\[Alpha]3},
    {\[Alpha]4}
   } );
c0 = ( {
    {c01},
    {c02},
    {c03},
    {c04}
   } );
c1 = ( {
    {c11},
    {c12},
    {c13},
    {c14}
   } );
\[Beta]1 = ( {
    {\[Beta]11},
    {\[Beta]12},
    {\[Beta]13},
    {\[Beta]14}
   } );
\[Beta]2 = ( {
    {\[Beta]21},
    {\[Beta]22},
    {\[Beta]23},
    {\[Beta]24}
   } );
\[Beta]3 = ( {
    {\[Beta]31},
    {\[Beta]32},
    {\[Beta]33},
    {\[Beta]34}
   } );
\[Beta]4 = ( {
    {\[Beta]41},
    {\[Beta]42},
    {\[Beta]43},
    {\[Beta]44}
   } );
   A0 = ( {
   {0, 0, 0, 0},
   {A021, 0, 0, 0},
   {A031, A032, 0, 0},
   {A041, A042, A043, 0}
  } ); (* A0=(A011	A012	A013	A014
A021	A022	A023	A024
A031	A032	A033	A034
A041	A042	A043	A044

); *)
A1 = ( {
    {0, 0, 0, 0},
    {A121, 0, 0, 0},
    {A131, A132, 0, 0},
    {A141, A142, A143, 0}
   } ); (* Subscript[A, 1]=(A111	A112	A113	A114
A121	A122	A123	A124
A131	A132	A133	A134
A141	A142	A143	A144

); *)
B0 = ( {
    {0, 0, 0, 0},
    {B021, 0, 0, 0},
    {B031, B032, 0, 0},
    {B041, B042, B043, 0}
   } ); (* Subscript[B, 0]=(B011	B012	B013	B014
B021	B022	B023	B024
B031	B032	B033	B034
B041	B042	B043	B044

); *)
B1 = ( {
    {0, 0, 0, 0},
    {B121, 0, 0, 0},
    {B131, B132, 0, 0},
    {B141, B142, B143, 0}
   } ); (* Subscript[B, 1]=(B111	B112	B113	B114
B121	B122	B123	B124
B131	B132	B133	B134
B141	B142	B143	B144

); *)
H0 = ( {
    {H01},
    {H02},
    {H03},
    {H04}
   } );
H1 = ( {
    {H11},
    {H12},
    {H13},
    {H14}
   } );
VU = ( {
    {U},
    {U},
    {U},
    {U}
   } );
e = ( {
    {1},
    {1},
    {1},
    {1}
   } );

   K1 = IdentityMatrix[4] + (\[Sigma] I10)/
    h B0 . Inverse[IdentityMatrix[4] - \[Sigma] Sqrt[h] B1];
K2 = Inverse[
   IdentityMatrix[
     4] - \[Mu] h A0 - \[Mu] \[Sigma] I10 A1. B0.
      Inverse[IdentityMatrix[4] - \[Sigma] Sqrt[h] B1]];
H0 = K2.K1 .VU;
K1 = IdentityMatrix[4] + (\[Sigma] I10)/
    h A1 . Inverse[IdentityMatrix[4] - \[Sigma] Sqrt[h] A0];
K2 = Inverse[
   IdentityMatrix[
     4] - \[Mu] h A1 - \[Mu] \[Sigma] I10 A1. B0.
      Inverse[IdentityMatrix[4] - \[Sigma] Sqrt[h] B0]];
H1 = K2.K1.VU;
U2 = U + \[Mu] h \[Alpha]\[Transpose].H0 + \[Sigma] I1 \[Beta]1\
\[Transpose].H1 + (\[Sigma] I11)/Sqrt[
    h] \[Beta]2\[Transpose].H1 + \[Sigma] I10/
    h \[Beta]3\[Transpose].H1 + (\[Sigma] I111)/
    h \[Beta]4\[Transpose].H1;
U2 = U2[[1]][[1]];
k = Expand[U2^2/U^2];

I10 = 1/2 h (W + Z/Sqrt[3]);
I1 = W;
I11 = 1/2 (W^2 - h);
I111 = 1/6 (W^3 - 3 h W);

kW = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(5\)]\(Coefficient[k, W,
     2  i]
\*SuperscriptBox[\(h\), \(i\)]\ Factorial2[2  i - 1]\)\);
kW = Coefficient[k, W, 0] + kW;
kDet = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(3\)]\(Coefficient[kW,
      Z, 2  i]
\*SuperscriptBox[\(h\), \(i\)]\ Factorial2[2  i - 1]\)\) -
   1; (* subtract off the constant *)


   stableRegionsEq = kDet /.  {\[Mu] -> z/h, \[Sigma] -> w/Sqrt[h]};
stableRegionsEq = Simplify[stableRegionsEq];
DumpSave["stableRegionsEq.mx", stableRegionsEq];

Get["stableRegionsEq.mx"]
stableRegionsEq >> "stableeqs.txt"
Exit[]
