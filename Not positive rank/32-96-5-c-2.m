32.96.5.c.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.c.2.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Over F5, the second and third each
have 4 points, while the single rank one elliptic curve factor of the Jacobian has 8 points over F5. Over F7, the first has 8 points, while
the single rank one elliptic curve factor of the Jacobian has 12 points.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[-y*t + z*w, 2*x^2 - y*w + z*t, 2*y^2 + 2*z^2 + w*t]);
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are three genus one quotients by an involution
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
    x[1]^2 + x[4]^2 + 2*x[6]*x[7],
    x[1]*x[2] - x[6]^2,
    x[1]*x[3] - x[7]^2,
    x[1]*x[5] - x[4]*x[6],
    x[1]*x[6] - x[4]*x[7],
    x[4]*x[6] + x[1]*x[7] + 2*x[7]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] + x[6]*x[7] + 2*x[8]^2,
    x[2]*x[6] - x[5]*x[8],
    x[2]*x[7] - x[6]*x[8],
    x[5]^2 + x[6]^2 + 2*x[2]*x[8],
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    x[6]^2 + x[7]^2 + 2*x[3]*x[8],
    x[4]*x[5] + x[4]*x[7] + 2*x[6]*x[8],
    -x[6]^2 + x[4]*x[8],
    x[5]*x[6] + x[6]*x[7] + 2*x[8]^2,
    -x[6]^2 + x[5]*x[7],
    Curve over Rational Field defined by
    48*x[1]^2 - 288*x[1]*x[4] + 480*x[4]^2 - 672*x[4]*x[5] - 48*x[5]^2 - 
        12*x[6]*x[7] + 100*x[7]^2 - 27*x[6]*x[8] - 6*x[7]*x[8] - 3*x[8]^2,
    3*x[1]*x[2] - 3*x[6]*x[7] - 8*x[7]^2 + 3*x[6]*x[8] + 2*x[7]*x[8],
    4*x[1]*x[3] - 8*x[5]^2 - x[7]*x[8] + x[8]^2,
    24*x[1]*x[5] - 24*x[4]*x[5] + 24*x[5]^2 + 3*x[6]*x[7] - 10*x[7]^2 + 
        3*x[7]*x[8],
    18*x[1]*x[6] + 33*x[3]*x[6] - 60*x[4]*x[6] + 66*x[5]*x[6] + 20*x[4]*x[7] - 
        88*x[5]*x[7] + 36*x[4]*x[8] - 38*x[5]*x[8],
    3*x[3]*x[6] - 12*x[4]*x[6] + 24*x[5]*x[6] + 12*x[1]*x[7] - 32*x[4]*x[7] + 
        40*x[5]*x[7] + 8*x[5]*x[8],
    3*x[3]*x[6] - 3*x[5]*x[6] - 8*x[5]*x[7] + 3*x[1]*x[8] + 2*x[5]*x[8],
    3*x[2]^2 - 12*x[6]*x[8] - 8*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    4*x[2]*x[4] - 4*x[7]^2 - 4*x[7]*x[8] - x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    9*x[2]*x[6] - 36*x[3]*x[6] - 48*x[5]*x[6] + 112*x[5]*x[7] - 72*x[4]*x[8] + 
        48*x[5]*x[8],
    -12*x[5]*x[6] + 3*x[2]*x[7] - 8*x[5]*x[7],
    -12*x[3]*x[6] + 3*x[2]*x[8] - 8*x[5]*x[8],
    6*x[3]^2 - 24*x[5]^2 - 3*x[6]*x[8] + 4*x[7]*x[8] + 3*x[8]^2,
    24*x[3]*x[4] - 48*x[4]*x[5] - 6*x[6]*x[7] + 8*x[7]^2 - 3*x[6]*x[8] + 
        10*x[7]*x[8] + 3*x[8]^2,
    12*x[3]*x[5] - 24*x[4]*x[5] + 24*x[5]^2 - 3*x[6]*x[7] + 4*x[7]^2 + 
        3*x[7]*x[8],
    x[3]*x[7] - x[5]*x[8],
    4*x[5]*x[7] + x[3]*x[8] - 4*x[4]*x[8] + 4*x[5]*x[8],
    18*x[6]^2 - 12*x[6]*x[7] + 20*x[7]^2 - 18*x[6]*x[8] - 12*x[7]*x[8] - 
        9*x[8]^2,
    Curve over Rational Field defined by
    48*x[1]^2 - 288*x[1]*x[4] + 480*x[4]^2 + 672*x[4]*x[5] - 48*x[5]^2 + 
        12*x[6]*x[7] - 100*x[7]^2 - 27*x[6]*x[8] - 6*x[7]*x[8] + 3*x[8]^2,
    3*x[1]*x[2] - 3*x[6]*x[7] - 8*x[7]^2 - 3*x[6]*x[8] - 2*x[7]*x[8],
    4*x[1]*x[3] - 8*x[5]^2 - x[7]*x[8] - x[8]^2,
    24*x[1]*x[5] - 24*x[4]*x[5] - 24*x[5]^2 + 3*x[6]*x[7] - 10*x[7]^2 - 
        3*x[7]*x[8],
    18*x[1]*x[6] + 33*x[3]*x[6] - 60*x[4]*x[6] - 66*x[5]*x[6] + 20*x[4]*x[7] + 
        88*x[5]*x[7] - 36*x[4]*x[8] - 38*x[5]*x[8],
    3*x[3]*x[6] - 12*x[4]*x[6] - 24*x[5]*x[6] + 12*x[1]*x[7] - 32*x[4]*x[7] - 
        40*x[5]*x[7] + 8*x[5]*x[8],
    -3*x[3]*x[6] - 3*x[5]*x[6] - 8*x[5]*x[7] + 3*x[1]*x[8] - 2*x[5]*x[8],
    3*x[2]^2 - 12*x[6]*x[8] - 8*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    4*x[2]*x[4] - 4*x[7]^2 + 4*x[7]*x[8] - x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    9*x[2]*x[6] + 36*x[3]*x[6] - 48*x[5]*x[6] + 112*x[5]*x[7] - 72*x[4]*x[8] - 
        48*x[5]*x[8],
    -12*x[5]*x[6] + 3*x[2]*x[7] - 8*x[5]*x[7],
    -12*x[3]*x[6] + 3*x[2]*x[8] - 8*x[5]*x[8],
    6*x[3]^2 - 24*x[5]^2 - 3*x[6]*x[8] + 4*x[7]*x[8] - 3*x[8]^2,
    24*x[3]*x[4] + 48*x[4]*x[5] + 6*x[6]*x[7] - 8*x[7]^2 - 3*x[6]*x[8] + 
        10*x[7]*x[8] - 3*x[8]^2,
    12*x[3]*x[5] - 24*x[4]*x[5] - 24*x[5]^2 - 3*x[6]*x[7] + 4*x[7]^2 - 
        3*x[7]*x[8],
    x[3]*x[7] - x[5]*x[8],
    4*x[5]*x[7] + x[3]*x[8] - 4*x[4]*x[8] - 4*x[5]*x[8],
    18*x[6]^2 - 12*x[6]*x[7] + 20*x[7]^2 + 18*x[6]*x[8] + 12*x[7]*x[8] - 
        9*x[8]^2
]
*/
#EllipticCurve(Curve(Reduction((l[1]),7))); //8
#EllipticCurve(Curve(Reduction((l[2]),5))); //4
#EllipticCurve(Curve(Reduction((l[3]),5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,5))); //8
#EllipticCurve(Curve(Reduction(E,7))); //12
