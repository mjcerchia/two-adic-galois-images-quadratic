16.192.5.cf.1

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.5.cf.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, has 4 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*z - y*t + z*w + w*t, 4*x^2 - y^2 - z^2 - w^2 - t^2, 2*x^2 - y^2 + 3*y*z + y*w + y*t - z^2 - z*w - z*t - w*t]);

S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotients by an involution
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
    3*x[1]^2 - 27*x[4]^2 + 20*x[4]*x[5] + 20*x[3]*x[6] - 46*x[4]*x[6] + 
        20*x[5]*x[6] - 43*x[6]^2 + 6*x[6]*x[7] - 10*x[4]*x[8] + 6*x[5]*x[8] - 
        10*x[6]*x[8] + 6*x[7]*x[8] - x[8]^2,
    x[1]*x[2] - x[8]^2,
    x[1]*x[3] - x[4]*x[6],
    x[1]*x[4] + x[4]^2 - x[4]*x[5] - x[3]*x[6] + x[4]*x[6] - x[5]*x[6],
    x[1]*x[5] - x[4]*x[8],
    24*x[4]^2 - 17*x[4]*x[5] + 3*x[1]*x[6] - 17*x[3]*x[6] + 43*x[4]*x[6] - 
        17*x[5]*x[6] + 40*x[6]^2 - 6*x[6]*x[7] + 7*x[4]*x[8] - 6*x[5]*x[8] + 
        7*x[6]*x[8] - 6*x[7]*x[8] + x[8]^2,
    x[1]*x[7] - x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8] + x[4]*x[8] - x[5]*x[8] + x[6]*x[8] - x[7]*x[8],
    x[2]^2 + 42*x[4]*x[5] + 7*x[5]^2 + 156*x[5]*x[6] + 26*x[5]*x[7] + 
        240*x[6]*x[7] + 39*x[7]^2 - 102*x[4]*x[8] - 10*x[5]*x[8] + 18*x[6]*x[8] 
        + 10*x[7]*x[8] - 35*x[8]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[5] - 7*x[3]*x[5] - 26*x[4]*x[5] + 26*x[5]^2 - 49*x[5]*x[6] + 
        65*x[5]*x[7] - 43*x[6]*x[7] + 40*x[7]^2 - 26*x[4]*x[8] - x[5]*x[8] - 
        40*x[6]*x[8] + 5*x[7]*x[8] - x[8]^2,
    x[2]*x[6] - x[7]*x[8],
    7*x[3]*x[5] + 26*x[4]*x[5] - 26*x[5]^2 + 49*x[5]*x[6] + x[2]*x[7] - 
        65*x[5]*x[7] + 43*x[6]*x[7] - 39*x[7]^2 + 26*x[4]*x[8] + 40*x[6]*x[8] - 
        6*x[7]*x[8],
    7*x[4]*x[5] + 26*x[5]*x[6] + 40*x[6]*x[7] + x[2]*x[8] - 17*x[4]*x[8] + 
        x[5]*x[8] + 3*x[6]*x[8] + x[7]*x[8] - 6*x[8]^2,
    147*x[3]^2 + 17736*x[4]^2 - 546*x[3]*x[5] - 13928*x[4]*x[5] + 1365*x[5]^2 - 
        13562*x[3]*x[6] + 31432*x[4]*x[6] - 12608*x[5]*x[6] + 27949*x[6]^2 + 
        2184*x[5]*x[7] - 3240*x[6]*x[7] + 840*x[7]^2 + 3808*x[4]*x[8] - 
        4434*x[5]*x[8] + 4822*x[6]*x[8] - 4512*x[7]*x[8] + 739*x[8]^2,
    21*x[3]*x[4] - 960*x[4]^2 + 680*x[4]*x[5] + 758*x[3]*x[6] - 1651*x[4]*x[6] +
        563*x[5]*x[6] - 1471*x[6]^2 + 123*x[6]*x[7] - 280*x[4]*x[8] + 
        240*x[5]*x[8] - 298*x[6]*x[8] + 243*x[7]*x[8] - 40*x[8]^2,
    -x[4]*x[5] + x[5]^2 - x[5]*x[6] + x[3]*x[7] + x[5]*x[7] - x[4]*x[8],
    -x[5]*x[6] + x[3]*x[8],
    -x[5]*x[6] + x[4]*x[7]
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
