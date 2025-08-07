16.96.5.l.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.96.5.l.1.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. The first and third have 8 and 4 points mod 7,
respectively, while the single rank one elliptic curve factor of the Jacobian has 12 points mod 7. The second curve has 14 points mod 11 while the
while the single rank one elliptic curve factor of the Jacobian has 10 points mod 11.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*z + w*t, 2*y*w + z^2 + t^2, 4*x^2 + y*t + z*w]);

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
    x[1]*x[2] - 2*x[5]*x[6] - 2*x[6]*x[7],
    x[1]*x[3] - x[6]^2,
    x[1]*x[4] - 2*x[5]^2 - 2*x[6]^2,
    -x[4]*x[5] + x[1]*x[6],
    x[1]*x[7] - 2*x[5]*x[8] - 2*x[7]*x[8],
    -x[5]*x[6] + x[1]*x[8],
    x[2]*x[3] - x[7]^2,
    x[2]*x[4] - 2*x[6]^2 - 2*x[7]^2,
    x[2]*x[5] - 2*x[5]*x[8] - 2*x[7]*x[8],
    x[2]*x[6] - x[4]*x[7],
    -x[6]*x[7] + x[2]*x[8],
    2*x[3]^2 - x[6]*x[7] + 2*x[8]^2,
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    2*x[3]*x[7] - x[4]*x[7] + 2*x[6]*x[8],
    x[4]^2 - 2*x[5]*x[6] - 2*x[6]*x[7],
    x[4]*x[6] - 2*x[5]*x[8] - 2*x[7]*x[8],
    -x[6]^2 + x[4]*x[8],
    -x[6]^2 + x[5]*x[7],
    Curve over Rational Field defined by
    45*x[1]*x[2] - 45*x[6]^2 - 4*x[7]^2 - 24*x[6]*x[8] - 16*x[7]*x[8] + 
        12*x[8]^2,
    120*x[1]*x[3] - 96*x[5]^2 + 15*x[6]*x[8] - 14*x[7]*x[8] + 12*x[8]^2,
    5*x[1]*x[4] - x[4]^2 - 8*x[3]*x[5] + 4*x[5]^2 + x[7]^2 - x[7]*x[8],
    120*x[1]*x[5] + 96*x[3]*x[5] + 3*x[7]^2 + 15*x[6]*x[8] + 10*x[7]*x[8] + 
        15*x[8]^2,
    45*x[2]*x[6] - 672*x[5]*x[6] + 360*x[1]*x[7] - 72*x[4]*x[7] - 592*x[5]*x[7] 
        + 360*x[4]*x[8] + 1008*x[5]*x[8],
    45*x[2]*x[6] + 48*x[5]*x[6] - 256*x[5]*x[7] + 360*x[1]*x[8] + 288*x[4]*x[8] 
        + 576*x[5]*x[8],
    3*x[2]^2 + 24*x[6]*x[8] + 16*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[7]^2 + 2*x[7]*x[8] - x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    24*x[5]*x[6] + 3*x[2]*x[7] + 16*x[5]*x[7],
    24*x[5]*x[6] + 16*x[5]*x[7] + 3*x[2]*x[8] - 24*x[4]*x[8] - 48*x[5]*x[8],
    8*x[3]^2 + 8*x[5]^2 + x[7]*x[8] - x[8]^2,
    8*x[3]*x[4] + 16*x[3]*x[5] + x[7]*x[8] - x[8]^2,
    3*x[3]*x[6] - 3*x[5]*x[6] - 2*x[5]*x[7] + 3*x[4]*x[8] + 8*x[5]*x[8],
    x[3]*x[7] - x[5]*x[8],
    x[5]*x[7] + x[3]*x[8] - x[4]*x[8] - 2*x[5]*x[8],
    8*x[4]*x[5] + 16*x[5]^2 + x[7]^2 - x[7]*x[8],
    3*x[4]*x[6] - x[4]*x[7] - 6*x[5]*x[7] + 3*x[4]*x[8] + 6*x[5]*x[8],
    3*x[6]*x[7] - x[7]^2 - 3*x[6]*x[8] - 2*x[7]*x[8] - 3*x[8]^2,
    Curve over Rational Field defined by
    45*x[1]*x[2] - 45*x[6]^2 - 4*x[7]^2 + 24*x[6]*x[8] + 16*x[7]*x[8] + 
        12*x[8]^2,
    120*x[1]*x[3] - 96*x[5]^2 + 15*x[6]*x[8] - 14*x[7]*x[8] - 12*x[8]^2,
    5*x[1]*x[4] - x[4]^2 + 8*x[3]*x[5] + 4*x[5]^2 - x[7]^2 - x[7]*x[8],
    120*x[1]*x[5] + 96*x[3]*x[5] + 3*x[7]^2 - 15*x[6]*x[8] - 10*x[7]*x[8] + 
        15*x[8]^2,
    -45*x[2]*x[6] + 672*x[5]*x[6] + 360*x[1]*x[7] - 72*x[4]*x[7] + 592*x[5]*x[7]
        - 360*x[4]*x[8] + 1008*x[5]*x[8],
    45*x[2]*x[6] + 48*x[5]*x[6] - 256*x[5]*x[7] + 360*x[1]*x[8] + 288*x[4]*x[8] 
        - 576*x[5]*x[8],
    3*x[2]^2 + 24*x[6]*x[8] + 16*x[7]*x[8],
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[7]^2 - 2*x[7]*x[8] - x[8]^2,
    x[2]*x[5] - x[7]*x[8],
    24*x[5]*x[6] + 3*x[2]*x[7] + 16*x[5]*x[7],
    -24*x[5]*x[6] - 16*x[5]*x[7] + 3*x[2]*x[8] + 24*x[4]*x[8] - 48*x[5]*x[8],
    8*x[3]^2 + 8*x[5]^2 + x[7]*x[8] + x[8]^2,
    8*x[3]*x[4] - 16*x[3]*x[5] + x[7]*x[8] + x[8]^2,
    3*x[3]*x[6] + 3*x[5]*x[6] + 2*x[5]*x[7] - 3*x[4]*x[8] + 8*x[5]*x[8],
    x[3]*x[7] - x[5]*x[8],
    x[5]*x[7] + x[3]*x[8] - x[4]*x[8] + 2*x[5]*x[8],
    8*x[4]*x[5] - 16*x[5]^2 + x[7]^2 + x[7]*x[8],
    3*x[4]*x[6] - x[4]*x[7] + 6*x[5]*x[7] - 3*x[4]*x[8] + 6*x[5]*x[8],
    3*x[6]*x[7] - x[7]^2 + 3*x[6]*x[8] + 2*x[7]*x[8] - 3*x[8]^2
]
*/
#EllipticCurve(Curve(Reduction(l[1],7))); //8
#EllipticCurve(Curve(Reduction(l[3],7))); //4

#EllipticCurve(Curve(Reduction(l[2],11))); //14

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,7))); //12
#EllipticCurve(Curve(Reduction(E,11))); //10
