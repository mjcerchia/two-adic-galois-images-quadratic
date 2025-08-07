64.192.9.q.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 64.192.9.q.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, it has 4 points, 
while the single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<[x]> := ProjectiveSpace(Rationals(),4);


C := Curve(P,[y*v - u*r, y*v + z*w, x*z - t*v, x*z - y*r, x*w - t*u, x*v + w\
*r, x*v + y*z, x*w - y^2, z*u - w*v, x*y + w*t, y*u + w^2, x*r + z*t, z^2 + v*\
r, x^2 + y*t, x*u + y*w, x*w + y^2 + 2*t*u - u^2 - v^2, 2*x*y + z*v - 2*w*t + \
w*u, 2*x*r - 2*z*t - z*u + s^2, 2*x^2 - 2*y*t + y*u + v*r, x*w - 4*t^2 + r^2, \
4*x*t - x*u + z*r]);

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
     x[1]*x[2] - x[7]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - x[4]*x[6],
    -x[3]^2 + 64*x[1]*x[6] + 64*x[6]*x[7],
    -x[6]^2 + x[1]*x[7],
    -x[6]*x[7] + x[1]*x[8],
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    x[2]*x[7] - x[8]^2,
    -x[5]^2 + 64*x[2]*x[8] + 64*x[7]*x[8],
    x[3]*x[4] - 64*x[6]^2 - 64*x[7]^2,
    x[3]*x[5] - 64*x[6]*x[7] - 64*x[7]*x[8],
    -x[4]*x[6] + x[3]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    x[4]^2 - 64*x[6]*x[7] - 64*x[7]*x[8],
    x[4]*x[5] - 64*x[7]^2 - 64*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    -x[5]*x[7] + x[4]*x[8],
    -x[7]^2 + x[6]*x[8]
]
*/

#EllipticCurve(Curve(Reduction((l[1]),5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
