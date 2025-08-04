16.192.9.bf.1

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.9.bf.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, has 4 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,v,r,s> := ProjectiveSpace(Rationals(),8);
C := Curve(P,[v^2 + r*s, t*r + u*v, z*s + w*s - u*v, y*v + z*s - t*s, t*v - u*s, y*v - w*r, y*s + w*v, y*t + w*u, y*u + z*t - t^2, y*u - z*t + t^2 - r*s, y*r - y*s + u*r + u*s, z*t + w*t - u^2, y^2 - z*w + w*t, y*u + z*w + w^2, y*r - z*v + u*s, y*t - w*u + v*s, z*t - w*t - u^2 + s^2, 2*x^2 + y*r + u*s, y*z + y*w + z*u - t*u, y*z + y*w - z*u + t*u - v*r, y^2 + 2*z^2 - w^2 - t^2 - u^2 - r^2]);

S := AutomorphismGroup(C); 

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
\*
[
    Curve over Rational Field defined by
    2*x[1]^2 - x[3]^2 - 2*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - x[4]*x[6],
    -x[3]*x[4] + 2*x[1]*x[6] - 2*x[7]*x[8],
    -x[6]^2 + x[1]*x[7],
    -x[6]*x[7] + x[1]*x[8],
    2*x[2]^2 + x[5]^2 - 2*x[7]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    x[2]*x[7] - x[8]^2,
    x[4]*x[5] - 2*x[6]*x[7] + 2*x[2]*x[8],
    x[3]*x[5] - 2*x[6]^2 + 2*x[8]^2,
    -x[4]*x[6] + x[3]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    x[4]^2 - 2*x[6]^2 + 2*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    -x[5]*x[7] + x[4]*x[8],
    -x[7]^2 + x[6]*x[8]
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
