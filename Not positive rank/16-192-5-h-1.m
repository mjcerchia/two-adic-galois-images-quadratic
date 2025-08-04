16.192.5.h.1

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.5.h.1
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, it has 4 points.
But the single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x*y + x*z - y*z + z^2, 2*x*y - 2*x*z - 2*y^2 + 2*y*z + t^2, x*t + z*t + 4*w^2]);

S := AutomorphismGroup(C); 

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are four genus one quotients by an involution
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
    x[1]^2 - 65536*x[4]*x[6] - 512*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[3] - 256*x[6]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - 256*x[6]*x[8],
    -256*x[3]*x[4] + x[1]*x[6] - 512*x[7]*x[8],
    x[1]*x[7] - 512*x[2]*x[7] - 65536*x[4]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    512*x[2]^2 + 256*x[5]^2 - x[7]^2,
    x[2]*x[3] - 256*x[8]^2,
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    256*x[4]*x[5] - x[6]*x[7] + 512*x[2]*x[8],
    x[3]^2 - 256*x[4]*x[6],
    x[3]*x[5] - 256*x[4]*x[8],
    x[3]*x[7] - 256*x[6]*x[8],
    -x[5]*x[6] + x[3]*x[8],
    256*x[4]^2 - x[6]^2 + 512*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    x[5]*x[7] - 256*x[8]^2
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
