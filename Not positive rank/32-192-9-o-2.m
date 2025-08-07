32.192.9.o.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.o.2.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[u*r - v*s, x^2 + r*s - s^2, x^2 + v^2 - v*r, y^2 + u^2 - u*s, \
x*z - u*s + v*s, x*y - u*v + v*s, x*y - t*s, x^2 + t*r, x*z + t*v, y*z - t*u, \
x*y + y^2 - t^2, x*u - z*s, x*v - z*r, z^2 - u^2 + u*v, x*t + x*v - z*v, x*u -\
 y*t - z*u, x*u + y*u - z*t, x*t - x*s - y*s, x*s + y*r, x*u + y*v, x*u + x*r \
- y*u - y*v + y*r + y*s + z*t - z*v - z*r + z*s + w^2]);
 
C3 := Curve(Reduction(C,3)); 
S := AutomorphismGroup(C3); 

auts := [];
Stemp := Automorphisms(C3);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotient by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C3,[g]);
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
    Curve over GF(3) defined by
    x[3]*x[5] + 2*x[1]*x[12],
    x[3]*x[6] + 2*x[1]*x[13] + x[1]*x[15],
    x[3]*x[7] + 2*x[2]*x[10] + 2*x[1]*x[14] + x[1]*x[16],
    x[3]*x[8] + 2*x[1]*x[15],
    x[3]*x[9] + x[2]*x[10] + 2*x[1]*x[16],
    2*x[1]*x[4] + x[1]*x[6] + x[3]*x[10],
    2*x[1]*x[5] + x[1]*x[7] + x[3]*x[11],
    2*x[1]*x[6] + x[1]*x[8] + 2*x[9]*x[10] + x[3]*x[12],
    2*x[1]*x[7] + x[3]*x[13],
    2*x[1]*x[8] + x[9]*x[10] + x[3]*x[14],
    x[2]*x[3] + 2*x[1]*x[9] + x[3]*x[15],
    2*x[9]*x[10] + x[3]*x[16],
    x[4]^2 + 2*x[1]*x[12] + 2*x[1]*x[14],
    x[4]*x[5] + 2*x[1]*x[13],
    x[4]*x[6] + 2*x[1]*x[14],
    x[4]*x[7] + 2*x[1]*x[15],
    x[4]*x[8] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[4]*x[9],
    2*x[1]*x[5] + x[4]*x[10],
    2*x[1]*x[6] + x[4]*x[11],
    2*x[1]*x[7] + x[4]*x[12],
    2*x[1]*x[8] + x[4]*x[13],
    2*x[1]*x[9] + x[4]*x[14],
    2*x[9]*x[10] + x[4]*x[15],
    2*x[2]*x[3] + 2*x[2]*x[5] + x[4]*x[16],
    x[5]^2 + 2*x[1]*x[14],
    x[5]*x[6] + 2*x[1]*x[15],
    x[5]*x[7] + x[2]*x[10] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[5]*x[8],
    x[5]*x[9] + 2*x[2]*x[10],
    2*x[1]*x[6] + x[5]*x[10],
    2*x[1]*x[7] + x[5]*x[11],
    2*x[1]*x[8] + x[9]*x[10] + x[5]*x[12],
    2*x[1]*x[9] + x[5]*x[13],
    2*x[9]*x[10] + x[5]*x[14],
    2*x[2]*x[3] + x[5]*x[15],
    2*x[2]*x[4] + x[5]*x[16],
    x[6]^2 + x[2]*x[10] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[6]*x[7] + x[2]*x[11],
    x[6]*x[8] + 2*x[2]*x[10],
    x[6]*x[9] + 2*x[2]*x[11],
    2*x[1]*x[7] + x[6]*x[10],
    2*x[1]*x[8] + x[9]*x[10] + x[6]*x[11],
    x[2]*x[3] + 2*x[1]*x[9] + x[6]*x[12],
    2*x[9]*x[10] + x[6]*x[13],
    2*x[2]*x[3] + x[6]*x[14],
    2*x[2]*x[4] + x[2]*x[6] + x[6]*x[15],
    2*x[2]*x[5] + x[6]*x[16],
    x[7]^2 + 2*x[2]*x[10] + x[2]*x[12],
    x[7]*x[8] + 2*x[2]*x[11],
    x[7]*x[9] + 2*x[2]*x[12],
    2*x[1]*x[8] + x[7]*x[10] + x[9]*x[10],
    x[2]*x[3] + 2*x[1]*x[9] + x[7]*x[11],
    x[2]*x[4] + 2*x[2]*x[6] + 2*x[9]*x[10] + x[7]*x[12],
    2*x[2]*x[3] + x[7]*x[13],
    2*x[2]*x[4] + x[2]*x[6] + x[7]*x[14],
    2*x[2]*x[5] + x[2]*x[7] + x[7]*x[15],
    2*x[2]*x[6] + x[7]*x[16],
    x[8]^2 + 2*x[2]*x[12] + 2*x[2]*x[14],
    x[8]*x[9] + 2*x[2]*x[13],
    2*x[1]*x[9] + x[8]*x[10],
    2*x[9]*x[10] + x[8]*x[11],
    2*x[2]*x[3] + x[8]*x[12],
    2*x[2]*x[4] + x[8]*x[13],
    2*x[2]*x[5] + x[8]*x[14],
    2*x[2]*x[6] + x[8]*x[15],
    2*x[2]*x[7] + 2*x[2]*x[9] + x[8]*x[16],
    x[9]^2 + 2*x[2]*x[14],
    x[3]*x[4] + 2*x[1]*x[11],
    2*x[2]*x[3] + x[9]*x[11],
    2*x[2]*x[4] + x[2]*x[6] + x[9]*x[12],
    2*x[2]*x[5] + x[9]*x[13],
    2*x[2]*x[6] + x[9]*x[14],
    2*x[2]*x[7] + x[9]*x[15],
    2*x[2]*x[8] + x[9]*x[16],
    x[10]^2 + 2*x[1]*x[11],
    x[10]*x[11] + 2*x[1]*x[12],
    x[10]*x[12] + 2*x[1]*x[13] + x[1]*x[15],
    x[10]*x[13] + 2*x[1]*x[14],
    x[10]*x[14] + 2*x[1]*x[15],
    x[2]*x[10] + x[10]*x[15] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[10]*x[16],
    x[11]^2 + 2*x[1]*x[13] + x[1]*x[15],
    2*x[2]*x[10] + x[11]*x[12] + 2*x[1]*x[14] + x[1]*x[16],
    x[11]*x[13] + 2*x[1]*x[15],
    x[2]*x[10] + x[11]*x[14] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[2]*x[11] + x[11]*x[15],
    2*x[2]*x[10] + x[11]*x[16],
    x[1]*x[2] + 2*x[2]*x[11] + x[12]^2 + 2*x[1]*x[15],
    x[2]*x[10] + x[12]*x[13] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[2]*x[11] + x[12]*x[14],
    2*x[2]*x[10] + x[2]*x[12] + x[12]*x[15],
    2*x[2]*x[11] + x[12]*x[16],
    2*x[1]*x[2] + x[13]^2,
    2*x[2]*x[10] + x[13]*x[14],
    2*x[2]*x[11] + x[13]*x[15],
    2*x[2]*x[12] + 2*x[2]*x[14] + x[13]*x[16],
    2*x[2]*x[11] + x[14]^2,
    2*x[2]*x[12] + x[14]*x[15],
    2*x[2]*x[13] + x[14]*x[16],
    2*x[2]*x[13] + x[2]*x[15] + x[15]^2,
    2*x[2]*x[14] + x[15]*x[16],
    2*x[2]^2 + 2*x[2]*x[15] + x[16]^2,
    x[3]^2 + 2*x[1]*x[10] + x[1]*x[12]
]
*/
#EllipticCurve(Curve(l[1])); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
