/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label  16.192.5.cg.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, has 4 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*z + y*t - z*w + w*t, 2*x^2 - y*z + w*t, 3*y^2 - 2*y*z + 2*y*w + 2*y*t + 3*z^2 - 2*z*w - 2*z*t + w^2 - 2*w*t + t^2]);

S := AutomorphismGroup(C); 

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotient by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);prj;
if Genus(CG) eq 1 then
try
l := Append(l,CG);
catch e
m := Append(m,CG);
end try;
end if;
CG; Genus(CG);
end if;
print ".........";
end for;

l;
/*
[
    Curve over Rational Field defined by
    13*x[1]^2 - 21*x[4]^2 + 50*x[4]*x[6] - 49*x[6]^2 + 2*x[6]*x[7] + 
        34*x[4]*x[8] - 2*x[5]*x[8] - 22*x[6]*x[8] + 2*x[7]*x[8] + x[8]^2,
    x[1]*x[2] - x[8]^2,
    x[1]*x[3] - x[4]*x[6],
    x[1]*x[4] + x[4]^2 + x[4]*x[5] - x[3]*x[6] - x[4]*x[6] - x[5]*x[6],
    x[1]*x[5] - x[4]*x[8],
    -8*x[4]^2 + 13*x[4]*x[5] + 13*x[1]*x[6] - 13*x[3]*x[6] + 37*x[4]*x[6] - 
        13*x[5]*x[6] - 36*x[6]^2 + 2*x[6]*x[7] + 21*x[4]*x[8] - 2*x[5]*x[8] - 
        9*x[6]*x[8] + 2*x[7]*x[8] + x[8]^2,
    x[1]*x[7] - x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8] + x[4]*x[8] + x[5]*x[8] - x[6]*x[8] - x[7]*x[8],
    x[2]^2 + 252*x[3]*x[5] - 558*x[4]*x[5] - 621*x[5]^2 + 1088*x[5]*x[6] + 
        806*x[5]*x[7] - 516*x[6]*x[7] - 313*x[7]^2 - 574*x[4]*x[8] - 
        82*x[5]*x[8] + 406*x[6]*x[8] + 34*x[7]*x[8] - 27*x[8]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[5] - 21*x[3]*x[5] + 50*x[4]*x[5] + 50*x[5]^2 - 99*x[5]*x[6] - 
        63*x[5]*x[7] + 49*x[6]*x[7] + 24*x[7]^2 + 50*x[4]*x[8] + x[5]*x[8] - 
        36*x[6]*x[8] + x[7]*x[8] + x[8]^2,
    x[2]*x[6] - x[7]*x[8],
    -21*x[3]*x[5] + 50*x[4]*x[5] + 50*x[5]^2 - 99*x[5]*x[6] + x[2]*x[7] - 
        63*x[5]*x[7] + 49*x[6]*x[7] + 25*x[7]^2 + 50*x[4]*x[8] - 36*x[6]*x[8] + 
        2*x[7]*x[8],
    -21*x[4]*x[5] + 50*x[5]*x[6] - 36*x[6]*x[7] + x[2]*x[8] - 13*x[4]*x[8] + 
        23*x[5]*x[8] + 13*x[6]*x[8] - 11*x[7]*x[8] + 2*x[8]^2,
    5733*x[3]^2 + 6168*x[4]^2 - 13650*x[3]*x[5] + 7176*x[4]*x[5] + 17199*x[5]^2 
        + 4550*x[3]*x[6] - 23704*x[4]*x[6] - 18928*x[5]*x[6] + 19111*x[6]^2 - 
        24024*x[5]*x[7] + 1864*x[6]*x[7] + 6552*x[7]^2 + 1008*x[4]*x[8] + 
        1542*x[5]*x[8] - 1186*x[6]*x[8] - 2192*x[7]*x[8] - 771*x[8]^2,
    273*x[3]*x[4] + 288*x[4]^2 - 468*x[4]*x[5] - 182*x[3]*x[6] - 695*x[4]*x[6] +
        637*x[5]*x[6] + 659*x[6]^2 - 397*x[6]*x[7] - 756*x[4]*x[8] + 
        72*x[5]*x[8] + 298*x[6]*x[8] - 85*x[7]*x[8] - 36*x[8]^2,
    -x[4]*x[5] - x[5]^2 + x[5]*x[6] + x[3]*x[7] + x[5]*x[7] - x[4]*x[8],
    -x[5]*x[6] + x[3]*x[8],
    -x[5]*x[6] + x[4]*x[7]
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
