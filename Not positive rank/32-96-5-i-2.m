32.96.5.i.2
32.96.5.c.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.c.2.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Over F7, the first and second have 8 and 4 
points, while the single rank one elliptic curve factor of the Jacobian has 12 points over F7. Over F11, the third has 14 points, while
the single rank one elliptic curve factor of the Jacobian has 10 points.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*w + z*t, 2*x^2 - y*t - z*w, 8*y*z - w^2 - t^2]);
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
    x[1]*x[2] + 8*x[5]*x[6] + 8*x[6]*x[7],
    x[1]*x[3] - x[6]^2,
    x[1]*x[4] + 8*x[5]^2 + 8*x[6]^2,
    -x[4]*x[5] + x[1]*x[6],
    x[1]*x[7] + 8*x[5]*x[8] + 8*x[7]*x[8],
    -x[5]*x[6] + x[1]*x[8],
    x[2]*x[3] - x[7]^2,
    x[2]*x[4] + 8*x[6]^2 + 8*x[7]^2,
    x[2]*x[5] + 8*x[5]*x[8] + 8*x[7]*x[8],
    x[2]*x[6] - x[4]*x[7],
    -x[6]*x[7] + x[2]*x[8],
    8*x[3]^2 + x[6]*x[7] + 8*x[8]^2,
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    8*x[3]*x[7] + x[4]*x[7] + 8*x[6]*x[8],
    x[4]^2 + 8*x[5]*x[6] + 8*x[6]*x[7],
    x[4]*x[6] + 8*x[5]*x[8] + 8*x[7]*x[8],
    -x[6]^2 + x[4]*x[8],
    -x[6]^2 + x[5]*x[7],
    Curve over Rational Field defined by
    45*x[1]*x[2] - 45*x[6]^2 - 4*x[7]^2 - 48*x[6]*x[8] + 32*x[7]*x[8] + 
        48*x[8]^2,
    240*x[1]*x[3] - 192*x[5]^2 - 15*x[6]*x[8] - 14*x[7]*x[8] - 24*x[8]^2,
    5*x[1]*x[4] - x[4]^2 + 64*x[3]*x[5] + 16*x[5]^2 - x[7]^2 - 2*x[7]*x[8],
    240*x[1]*x[5] + 768*x[3]*x[5] + 3*x[7]^2 + 30*x[6]*x[8] - 20*x[7]*x[8] + 
        60*x[8]^2,
    45*x[2]*x[6] - 1344*x[5]*x[6] + 360*x[1]*x[7] - 72*x[4]*x[7] + 
        1184*x[5]*x[7] - 720*x[4]*x[8] + 4032*x[5]*x[8],
    -45*x[2]*x[6] - 96*x[5]*x[6] - 512*x[5]*x[7] + 720*x[1]*x[8] + 576*x[4]*x[8]
        - 2304*x[5]*x[8],
    3*x[2]^2 - 48*x[6]*x[8] + 32*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[7]^2 - 4*x[7]*x[8] - 4*x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    -48*x[5]*x[6] + 3*x[2]*x[7] + 32*x[5]*x[7],
    24*x[5]*x[6] - 16*x[5]*x[7] + 3*x[2]*x[8] + 24*x[4]*x[8] - 96*x[5]*x[8],
    64*x[3]^2 + 16*x[5]^2 + x[7]*x[8] + 2*x[8]^2,
    16*x[3]*x[4] - 64*x[3]*x[5] + x[7]*x[8] + 2*x[8]^2,
    6*x[3]*x[6] + 3*x[5]*x[6] - 2*x[5]*x[7] + 3*x[4]*x[8] - 16*x[5]*x[8],
    x[3]*x[7] - x[5]*x[8],
    x[5]*x[7] + 4*x[3]*x[8] - x[4]*x[8] + 4*x[5]*x[8],
    16*x[4]*x[5] - 64*x[5]^2 + x[7]^2 + 2*x[7]*x[8],
    3*x[4]*x[6] + x[4]*x[7] - 12*x[5]*x[7] + 6*x[4]*x[8] - 24*x[5]*x[8],
    3*x[6]*x[7] + x[7]^2 + 6*x[6]*x[8] - 4*x[7]*x[8] + 12*x[8]^2,
    Curve over Rational Field defined by
    45*x[1]*x[2] - 45*x[6]^2 - 4*x[7]^2 + 48*x[6]*x[8] - 32*x[7]*x[8] + 
        48*x[8]^2,
    240*x[1]*x[3] - 192*x[5]^2 - 15*x[6]*x[8] - 14*x[7]*x[8] + 24*x[8]^2,
    5*x[1]*x[4] - x[4]^2 - 64*x[3]*x[5] + 16*x[5]^2 + x[7]^2 - 2*x[7]*x[8],
    240*x[1]*x[5] + 768*x[3]*x[5] + 3*x[7]^2 - 30*x[6]*x[8] + 20*x[7]*x[8] + 
        60*x[8]^2,
    -45*x[2]*x[6] + 1344*x[5]*x[6] + 360*x[1]*x[7] - 72*x[4]*x[7] - 
        1184*x[5]*x[7] + 720*x[4]*x[8] + 4032*x[5]*x[8],
    -45*x[2]*x[6] - 96*x[5]*x[6] - 512*x[5]*x[7] + 720*x[1]*x[8] + 576*x[4]*x[8]
        + 2304*x[5]*x[8],
    3*x[2]^2 - 48*x[6]*x[8] + 32*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[7]^2 + 4*x[7]*x[8] - 4*x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    -48*x[5]*x[6] + 3*x[2]*x[7] + 32*x[5]*x[7],
    -24*x[5]*x[6] + 16*x[5]*x[7] + 3*x[2]*x[8] - 24*x[4]*x[8] - 96*x[5]*x[8],
    64*x[3]^2 + 16*x[5]^2 + x[7]*x[8] - 2*x[8]^2,
    16*x[3]*x[4] + 64*x[3]*x[5] + x[7]*x[8] - 2*x[8]^2,
    6*x[3]*x[6] - 3*x[5]*x[6] + 2*x[5]*x[7] - 3*x[4]*x[8] - 16*x[5]*x[8],
    x[3]*x[7] - x[5]*x[8],
    x[5]*x[7] + 4*x[3]*x[8] - x[4]*x[8] - 4*x[5]*x[8],
    16*x[4]*x[5] + 64*x[5]^2 + x[7]^2 - 2*x[7]*x[8],
    3*x[4]*x[6] + x[4]*x[7] + 12*x[5]*x[7] - 6*x[4]*x[8] - 24*x[5]*x[8],
    3*x[6]*x[7] + x[7]^2 - 6*x[6]*x[8] + 4*x[7]*x[8] + 12*x[8]^2
]
*/
#EllipticCurve(Curve(Reduction((l[1]),7))); //8
#EllipticCurve(Curve(Reduction((l[2]),7))); //4
#EllipticCurve(Curve(Reduction((l[3]),11))); //14

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,7))); //12
#EllipticCurve(Curve(Reduction(E,11))); //10
