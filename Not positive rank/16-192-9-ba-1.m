16.192.9.ba.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.9.ba.1.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F3, it has 4 points.
But the single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,r,u,s,v> := ProjectiveSpace(Rationals(),8);
C := Curve(P,[t*r - u*s, x*r + w*u, x*s + w*t, x*t - z*s, t^2 - v*s + s^2, t*u - v*s - r*s, t*u - v*r + r*s, u^2 - v*r - r^2, x*t + w*v - w*s, x*u + w*v + w*r, x*u - z*r, x*v + x*r - z*u, x*v - x*s - z*t, x^2 + z*w, t*v + t*r - u*v + u*s, x^2 - z*w - 2*w^2 + t*u - s^2, x^2 - z*w + 2*w^2 + t*u + r^2, x^2 - 2*z^2 - z*w + t^2 - v^2 - r*s, x*v - 2*y^2, 2*x*z - 2*x*w + t*v + t*r + u*r, 2*x*z + 2*x*w + t*v + t*r - t*s]);


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
    x[3]*x[5] + 2*x[1]*x[11] + 2*x[1]*x[15],
    x[3]*x[6] + 2*x[1]*x[12] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[3]*x[7] + 2*x[1]*x[13],
    x[3]*x[8] + 2*x[2]*x[10] + 2*x[1]*x[14],
    x[3]*x[9] + 2*x[2]*x[11] + 2*x[1]*x[15],
    -x[1]*x[4] + x[3]*x[10],
    -x[1]*x[5] + x[3]*x[11],
    -x[1]*x[6] + x[3]*x[12],
    -x[1]*x[7] + x[3]*x[13],
    -x[1]*x[8] + x[3]*x[14],
    -x[1]*x[9] + x[3]*x[15],
    -x[9]*x[10] + x[3]*x[16],
    x[4]^2 + 2*x[1]*x[11] + 2*x[1]*x[15],
    x[4]*x[5] + 2*x[1]*x[12] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[4]*x[6] + 2*x[1]*x[13],
    x[4]*x[7] + 2*x[2]*x[10] + 2*x[1]*x[14],
    x[4]*x[8] + 2*x[2]*x[11] + 2*x[1]*x[15],
    x[4]*x[9] + 2*x[2]*x[12] + 2*x[1]*x[16],
    -x[1]*x[5] + x[4]*x[10],
    -x[1]*x[6] + x[4]*x[11],
    -x[1]*x[7] + x[4]*x[12],
    -x[1]*x[8] + x[4]*x[13],
    -x[1]*x[9] + x[4]*x[14],
    -x[9]*x[10] + x[4]*x[15],
    -x[2]*x[3] + x[4]*x[16],
    2*x[1]*x[2] + x[5]^2 + 2*x[1]*x[13],
    x[5]*x[6] + 2*x[2]*x[10] + 2*x[1]*x[14],
    x[5]*x[7] + 2*x[2]*x[11] + 2*x[1]*x[15],
    x[5]*x[8] + 2*x[2]*x[12] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[5]*x[9] + 2*x[2]*x[13],
    -x[1]*x[6] + x[5]*x[10],
    -x[1]*x[7] + x[5]*x[11],
    -x[1]*x[8] + x[5]*x[12],
    -x[1]*x[9] + x[5]*x[13],
    -x[9]*x[10] + x[5]*x[14],
    -x[2]*x[3] + x[5]*x[15],
    -x[2]*x[4] + x[5]*x[16],
    x[6]^2 + 2*x[2]*x[11] + 2*x[1]*x[15],
    x[6]*x[7] + 2*x[2]*x[12] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[6]*x[8] + 2*x[2]*x[13],
    x[6]*x[9] + 2*x[2]*x[10] + 2*x[2]*x[14],
    -x[1]*x[7] + x[6]*x[10],
    -x[1]*x[8] + x[6]*x[11],
    -x[1]*x[9] + x[6]*x[12],
    -x[9]*x[10] + x[6]*x[13],
    -x[2]*x[3] + x[6]*x[14],
    -x[2]*x[4] + x[6]*x[15],
    -x[2]*x[5] + x[6]*x[16],
    2*x[1]*x[2] + x[7]^2 + 2*x[2]*x[13],
    x[7]*x[8] + 2*x[2]*x[10] + 2*x[2]*x[14],
    x[7]*x[9] + 2*x[2]*x[11] + 2*x[2]*x[15],
    -x[1]*x[8] + x[7]*x[10],
    -x[1]*x[9] + x[7]*x[11],
    -x[9]*x[10] + x[7]*x[12],
    -x[2]*x[3] + x[7]*x[13],
    -x[2]*x[4] + x[7]*x[14],
    -x[2]*x[5] + x[7]*x[15],
    -x[2]*x[6] + x[7]*x[16],
    x[8]^2 + 2*x[2]*x[11] + 2*x[2]*x[15],
    x[8]*x[9] + 2*x[2]*x[12] + 2*x[2]*x[16],
    -x[1]*x[9] + x[8]*x[10],
    -x[9]*x[10] + x[8]*x[11],
    -x[2]*x[3] + x[8]*x[12],
    -x[2]*x[4] + x[8]*x[13],
    -x[2]*x[5] + x[8]*x[14],
    -x[2]*x[6] + x[8]*x[15],
    -x[2]*x[7] + x[8]*x[16],
    2*x[2]^2 + x[9]^2 + 2*x[2]*x[13],
    x[3]*x[4] + 2*x[1]*x[10] + 2*x[1]*x[14],
    -x[2]*x[3] + x[9]*x[11],
    -x[2]*x[4] + x[9]*x[12],
    -x[2]*x[5] + x[9]*x[13],
    -x[2]*x[6] + x[9]*x[14],
    -x[2]*x[7] + x[9]*x[15],
    -x[2]*x[8] + x[9]*x[16],
    x[10]^2 - x[1]*x[11],
    x[10]*x[11] - x[1]*x[12],
    x[10]*x[12] - x[1]*x[13],
    x[10]*x[13] - x[1]*x[14],
    x[10]*x[14] - x[1]*x[15],
    x[10]*x[15] - x[1]*x[16],
    -x[1]*x[2] + x[10]*x[16],
    x[11]^2 - x[1]*x[13],
    x[11]*x[12] - x[1]*x[14],
    x[11]*x[13] - x[1]*x[15],
    x[11]*x[14] - x[1]*x[16],
    -x[1]*x[2] + x[11]*x[15],
    -x[2]*x[10] + x[11]*x[16],
    x[12]^2 - x[1]*x[15],
    x[12]*x[13] - x[1]*x[16],
    -x[1]*x[2] + x[12]*x[14],
    -x[2]*x[10] + x[12]*x[15],
    -x[2]*x[11] + x[12]*x[16],
    -x[1]*x[2] + x[13]^2,
    -x[2]*x[10] + x[13]*x[14],
    -x[2]*x[11] + x[13]*x[15],
    -x[2]*x[12] + x[13]*x[16],
    -x[2]*x[11] + x[14]^2,
    -x[2]*x[12] + x[14]*x[15],
    -x[2]*x[13] + x[14]*x[16],
    -x[2]*x[13] + x[15]^2,
    -x[2]*x[14] + x[15]*x[16],
    -x[2]*x[15] + x[16]^2,
    2*x[1]^2 + x[3]^2 + 2*x[1]*x[13]
]
*/

#EllipticCurve(Curve(Reduction(l[1],3))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
