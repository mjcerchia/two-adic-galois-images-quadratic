64.96.5.d.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 64.96.5.d.2
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Over F7, the first and third
have 4 and 8 points, respectively, while the single rank one elliptic curve factor of the Jacobian has 12 points over F7. Over F5,
the second has 12 points, while the single rank one elliptic curve factor of the Jacobian has 8 points.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[x*z + y*w, 4*x*w + y*z - 2*t^2, 4*x^2 + y^2 - z*w]);
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are three genus one quotient by an involution
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
    4*x[1]^2 + x[4]^2 - 64*x[6]*x[7],
    x[1]*x[2] - x[6]^2,
    x[1]*x[3] - x[7]^2,
    x[1]*x[5] - x[4]*x[6],
    x[1]*x[6] - x[4]*x[7],
    x[4]*x[6] + 4*x[1]*x[7] - 64*x[7]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] + 4*x[6]*x[7] - 64*x[8]^2,
    x[2]*x[6] - x[5]*x[8],
    x[2]*x[7] - x[6]*x[8],
    -x[5]^2 - 4*x[6]^2 + 64*x[2]*x[8],
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    -x[6]^2 - 4*x[7]^2 + 64*x[3]*x[8],
    x[4]*x[5] + 4*x[4]*x[7] - 64*x[6]*x[8],
    -x[6]^2 + x[4]*x[8],
    x[5]*x[6] + 4*x[6]*x[7] - 64*x[8]^2,
    -x[6]^2 + x[5]*x[7],
    Curve over Rational Field defined by
    3*x[1]^2 + 36*x[1]*x[5] + 120*x[5]^2 + 336*x[5]*x[6] - 48*x[6]^2 + 
        6*x[4]*x[7] + 100*x[7]^2 - 27*x[4]*x[8] + 12*x[7]*x[8] - 12*x[8]^2,
    3*x[1]*x[2] + 12*x[1]*x[5] - 72*x[5]*x[6] - 144*x[6]^2 + 6*x[4]*x[7] + 
        76*x[7]^2 - 6*x[4]*x[8] + 44*x[7]*x[8],
    x[1]*x[3] + 4*x[6]^2 - x[7]*x[8] - 2*x[8]^2,
    9*x[1]*x[4] - 132*x[3]*x[4] + 60*x[4]*x[5] + 132*x[4]*x[6] + 40*x[5]*x[7] + 
        352*x[6]*x[7] - 144*x[5]*x[8] - 304*x[6]*x[8],
    12*x[1]*x[6] + 24*x[5]*x[6] + 48*x[6]^2 - 3*x[4]*x[7] - 20*x[7]^2 - 
        12*x[7]*x[8],
    3*x[3]*x[4] - 3*x[4]*x[5] - 12*x[4]*x[6] + 3*x[1]*x[7] + 16*x[5]*x[7] + 
        40*x[6]*x[7] - 16*x[6]*x[8],
    -6*x[3]*x[4] - 3*x[4]*x[6] + 16*x[6]*x[7] + 3*x[1]*x[8] + 8*x[6]*x[8],
    3*x[2]^2 - 48*x[5]^2 - 288*x[5]*x[6] - 432*x[6]^2 + 24*x[7]^2 - 3*x[4]*x[8] 
        + 28*x[7]*x[8] + 24*x[8]^2,
    12*x[2]*x[3] + 24*x[5]*x[6] + 144*x[6]^2 + 3*x[4]*x[7] + 8*x[7]^2 + 
        6*x[4]*x[8] + 4*x[7]*x[8] - 36*x[8]^2,
    9*x[2]*x[4] - 36*x[3]*x[4] + 36*x[4]*x[5] + 132*x[4]*x[6] + 112*x[6]*x[7] - 
        72*x[5]*x[8] - 96*x[6]*x[8],
    x[2]*x[5] + 4*x[5]^2 + 12*x[5]*x[6] - x[7]^2 + 2*x[7]*x[8] - x[8]^2,
    x[2]*x[6] + 4*x[5]*x[6] + 12*x[6]^2 - x[7]*x[8],
    -3*x[4]*x[6] + 3*x[2]*x[7] + 12*x[5]*x[7] + 40*x[6]*x[7],
    -3*x[3]*x[4] + 3*x[2]*x[8] + 12*x[5]*x[8] + 40*x[6]*x[8],
    24*x[3]^2 - 24*x[6]^2 - 3*x[4]*x[8] - 8*x[7]*x[8] + 12*x[8]^2,
    24*x[3]*x[5] + 24*x[5]*x[6] + 3*x[4]*x[7] + 8*x[7]^2 - 3*x[4]*x[8] - 
        20*x[7]*x[8] + 12*x[8]^2,
    48*x[3]*x[6] - 24*x[5]*x[6] - 48*x[6]^2 - 3*x[4]*x[7] - 8*x[7]^2 + 
        12*x[7]*x[8],
    x[3]*x[7] - x[6]*x[8],
    x[6]*x[7] + x[3]*x[8] - x[5]*x[8] - 2*x[6]*x[8],
    9*x[4]^2 + 12*x[4]*x[7] + 40*x[7]^2 - 36*x[4]*x[8] + 48*x[7]*x[8] - 
        72*x[8]^2,
    Curve over Rational Field defined by
    3*x[1]^2 + 36*x[1]*x[5] + 120*x[5]^2 - 336*x[5]*x[6] - 48*x[6]^2 - 
        6*x[4]*x[7] - 100*x[7]^2 - 27*x[4]*x[8] + 12*x[7]*x[8] + 12*x[8]^2,
    3*x[1]*x[2] - 12*x[1]*x[5] - 72*x[5]*x[6] + 144*x[6]^2 + 6*x[4]*x[7] + 
        76*x[7]^2 + 6*x[4]*x[8] - 44*x[7]*x[8],
    x[1]*x[3] + 4*x[6]^2 - x[7]*x[8] + 2*x[8]^2,
    9*x[1]*x[4] - 132*x[3]*x[4] + 60*x[4]*x[5] - 132*x[4]*x[6] + 40*x[5]*x[7] - 
        352*x[6]*x[7] + 144*x[5]*x[8] - 304*x[6]*x[8],
    12*x[1]*x[6] + 24*x[5]*x[6] - 48*x[6]^2 - 3*x[4]*x[7] - 20*x[7]^2 + 
        12*x[7]*x[8],
    3*x[3]*x[4] - 3*x[4]*x[5] + 12*x[4]*x[6] + 3*x[1]*x[7] + 16*x[5]*x[7] - 
        40*x[6]*x[7] - 16*x[6]*x[8],
    6*x[3]*x[4] - 3*x[4]*x[6] + 16*x[6]*x[7] + 3*x[1]*x[8] - 8*x[6]*x[8],
    3*x[2]^2 - 48*x[5]^2 + 288*x[5]*x[6] - 432*x[6]^2 - 24*x[7]^2 - 3*x[4]*x[8] 
        + 28*x[7]*x[8] - 24*x[8]^2,
    12*x[2]*x[3] + 24*x[5]*x[6] - 144*x[6]^2 + 3*x[4]*x[7] + 8*x[7]^2 - 
        6*x[4]*x[8] - 4*x[7]*x[8] - 36*x[8]^2,
    9*x[2]*x[4] + 36*x[3]*x[4] - 36*x[4]*x[5] + 132*x[4]*x[6] + 112*x[6]*x[7] - 
        72*x[5]*x[8] + 96*x[6]*x[8],
    x[2]*x[5] - 4*x[5]^2 + 12*x[5]*x[6] - x[7]^2 - 2*x[7]*x[8] - x[8]^2,
    x[2]*x[6] - 4*x[5]*x[6] + 12*x[6]^2 - x[7]*x[8],
    -3*x[4]*x[6] + 3*x[2]*x[7] - 12*x[5]*x[7] + 40*x[6]*x[7],
    -3*x[3]*x[4] + 3*x[2]*x[8] - 12*x[5]*x[8] + 40*x[6]*x[8],
    24*x[3]^2 - 24*x[6]^2 - 3*x[4]*x[8] - 8*x[7]*x[8] - 12*x[8]^2,
    24*x[3]*x[5] - 24*x[5]*x[6] - 3*x[4]*x[7] - 8*x[7]^2 - 3*x[4]*x[8] - 
        20*x[7]*x[8] - 12*x[8]^2,
    48*x[3]*x[6] - 24*x[5]*x[6] + 48*x[6]^2 - 3*x[4]*x[7] - 8*x[7]^2 - 
        12*x[7]*x[8],
    x[3]*x[7] - x[6]*x[8],
    x[6]*x[7] + x[3]*x[8] - x[5]*x[8] + 2*x[6]*x[8],
    9*x[4]^2 + 12*x[4]*x[7] + 40*x[7]^2 + 36*x[4]*x[8] - 48*x[7]*x[8] - 
        72*x[8]^2
]
*/

#EllipticCurve(Curve(Reduction((l[1]),7))); //8
#EllipticCurve(Curve(Reduction((l[2]),5))); //8
#EllipticCurve(Curve(Reduction((l[3]),7))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,7))); //12
#EllipticCurve(Curve(Reduction(E,5))); //4
