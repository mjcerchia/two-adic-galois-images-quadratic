16.96.5.c.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.96.5.c.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F3, has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*t + z*w, y^2 + 2*y*w - z^2 - 2*z*t - w^2 + t^2, 4*x^2 + y*z - w*t]);

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
\*
[
    Curve over Rational Field defined by
    x[1]^2 - x[4]^2 - 4*x[4]*x[6] + 13*x[6]^2 + 42*x[3]*x[7] - 29*x[7]^2 + 
        6*x[6]*x[8],
    x[1]*x[2] - x[6]^2,
    x[1]*x[3] - x[7]^2,
    x[1]*x[5] - x[4]*x[6],
    x[1]*x[6] - x[4]*x[7],
    -x[4]*x[6] + 2*x[6]^2 + x[1]*x[7] + 7*x[3]*x[7] - 6*x[7]^2 + x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    x[2]^2 + 2*x[2]*x[5] - x[5]^2 + x[6]^2 - 6*x[6]*x[8] + 7*x[8]^2,
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[5]*x[6],
    x[2]*x[6] - x[5]*x[8],
    x[2]*x[7] - x[6]*x[8],
    -x[5]*x[6] + x[6]*x[7] + x[2]*x[8] + 7*x[3]*x[8] + 2*x[5]*x[8] - 
        6*x[7]*x[8],
    7*x[3]^2 - x[6]^2 - 6*x[3]*x[7] + x[7]^2 + 2*x[6]*x[8] + x[8]^2,
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    x[4]*x[5] - 2*x[5]*x[6] - x[4]*x[7] + 6*x[6]*x[7] - x[5]*x[8] - 7*x[7]*x[8],
    -x[6]^2 + x[4]*x[8],
    -x[6]^2 + x[5]*x[7]
]
*/

#EllipticCurve(Curve(Reduction(l[1],3))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
