16.96.5.cy.1

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.96.5.cy.1
We compute the automorphism group over F5 and find there are two genus one quotient by an involution. They have 8 and 6 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[2*x^2 + y*t - z*t - t^2, 2*x*y - 2*x*z - y^2 + y*t - z^2 - z*t, 2*x^2 + x*y - x*z - 2*x*t - y*t + z*t - 2*w^2 + t^2]);
C5 := Curve(Reduction(C,5)); 
S := AutomorphismGroup(C5); 

auts := [];
Stemp := Automorphisms(C5);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are two genus one quotients by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C5,[g]);
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
    Curve over GF(5) defined by
    x[3]*x[5] + 4*x[1]*x[7],
    2*x[1]*x[2] + 3*x[2]^2 + 2*x[1]*x[3] + 2*x[1]*x[4] + 2*x[2]*x[4] + 
        4*x[1]*x[5] + x[2]*x[5] + 4*x[3]*x[6] + x[1]*x[7] + 2*x[2]*x[7] + 
        2*x[1]*x[8] + 4*x[2]*x[8] + 2*x[6]*x[8],
    4*x[1]*x[2] + 3*x[2]^2 + x[2]*x[3] + 3*x[2]*x[5] + 3*x[1]*x[7] + 2*x[2]*x[7]
        + x[3]*x[7] + x[1]*x[8] + 4*x[2]*x[8] + 2*x[6]*x[8],
    3*x[2]^2 + 4*x[2]*x[4] + 2*x[2]*x[5] + 3*x[1]*x[7] + 2*x[2]*x[7] + 
        2*x[1]*x[8] + 4*x[2]*x[8] + 4*x[3]*x[8] + 2*x[6]*x[8],
    3*x[1]*x[3] + x[1]*x[4] + 4*x[4]^2 + 4*x[1]*x[5] + 3*x[1]*x[6] + x[1]*x[7] +
        3*x[1]*x[8],
    2*x[4]*x[5] + 3*x[1]*x[8],
    3*x[1]*x[2] + 2*x[1]*x[3] + 3*x[2]*x[3] + 3*x[1]*x[4] + x[2]*x[4] + 
        2*x[2]*x[5] + x[1]*x[6] + 3*x[4]*x[6] + x[1]*x[8],
    3*x[2]^2 + 4*x[2]*x[4] + 2*x[2]*x[5] + 3*x[1]*x[7] + 2*x[2]*x[7] + 
        4*x[4]*x[7] + 2*x[1]*x[8] + 4*x[2]*x[8] + 2*x[6]*x[8],
    3*x[1]*x[2] + 3*x[2]^2 + 2*x[2]*x[3] + 2*x[2]*x[5] + 4*x[1]*x[7] + 
        2*x[2]*x[7] + 4*x[1]*x[8] + 4*x[2]*x[8] + 3*x[4]*x[8] + 2*x[6]*x[8],
    2*x[1]*x[2] + 3*x[5]^2,
    3*x[2]^2 + 4*x[2]*x[4] + 2*x[2]*x[5] + 4*x[5]*x[6] + 3*x[1]*x[7] + 
        2*x[2]*x[7] + 2*x[1]*x[8] + 4*x[2]*x[8] + 2*x[6]*x[8],
    3*x[2]*x[3] + 2*x[5]*x[7],
    2*x[2]*x[4] + 3*x[5]*x[8],
    4*x[2]^2 + x[1]*x[3] + 4*x[2]*x[3] + 2*x[1]*x[4] + 2*x[2]*x[4] + 3*x[1]*x[5]
        + 3*x[2]*x[5] + 3*x[1]*x[6] + 4*x[2]*x[6] + 4*x[6]^2 + x[1]*x[7] + 
        4*x[2]*x[7] + x[1]*x[8] + 2*x[6]*x[8],
    x[1]*x[2] + x[2]^2 + 2*x[2]*x[3] + 3*x[2]*x[5] + 4*x[2]*x[6] + 3*x[1]*x[7] +
        x[6]*x[7] + 3*x[1]*x[8] + 2*x[2]*x[8],
    3*x[3]*x[4] + 2*x[1]*x[6],
    2*x[2]^2 + 3*x[2]*x[4] + 3*x[2]*x[5] + 2*x[2]*x[6] + 2*x[2]*x[7] + 2*x[7]^2 
        + 2*x[2]*x[8],
    x[2]*x[6] + 4*x[7]*x[8],
    2*x[2]*x[3] + 4*x[2]*x[4] + x[2]*x[5] + 2*x[2]*x[6] + 4*x[2]*x[7] + 
        2*x[2]*x[8] + x[8]^2,
    3*x[1]*x[2] + 3*x[3]^2 + 2*x[1]*x[4] + 2*x[1]*x[5] + 3*x[1]*x[6] + 
        3*x[1]*x[7] + 3*x[1]*x[8],
    Curve over GF(5) defined by
    4*x[2]^2 + 2*x[2]*x[3] + x[1]*x[5] + x[3]*x[5] + 4*x[1]*x[6],
    3*x[1]*x[2] + x[2]^2 + 4*x[1]*x[3] + 2*x[2]*x[3] + 3*x[1]*x[5] + 
        2*x[3]*x[6],
    x[2]*x[4] + 2*x[1]*x[7] + 3*x[2]*x[7] + 2*x[3]*x[7] + 3*x[1]*x[8] + 
        x[2]*x[8],
    x[1]*x[7] + 3*x[2]*x[7] + x[2]*x[8] + x[3]*x[8],
    3*x[1]*x[2] + 2*x[4]^2,
    2*x[2]*x[4] + 4*x[4]*x[5] + 3*x[2]*x[7] + x[1]*x[8] + x[2]*x[8],
    x[2]*x[4] + 2*x[4]*x[6] + 2*x[1]*x[7] + 3*x[2]*x[7] + 3*x[1]*x[8] + 
        x[2]*x[8],
    4*x[1]*x[2] + 2*x[1]*x[3] + 2*x[2]*x[3] + 3*x[1]*x[5] + x[4]*x[7],
    2*x[2]*x[3] + 3*x[4]*x[8],
    x[1]*x[2] + 2*x[2]^2 + 3*x[1]*x[3] + 4*x[2]*x[3] + 3*x[1]*x[5] + 3*x[2]*x[5]
        + x[5]^2 + x[6]^2,
    4*x[1]*x[2] + 4*x[2]^2 + 2*x[1]*x[3] + 4*x[2]*x[3] + x[1]*x[5] + 4*x[2]*x[5]
        + x[1]*x[6] + 3*x[5]*x[6] + x[6]^2,
    4*x[1]*x[7] + 4*x[5]*x[7] + x[2]*x[8] + 4*x[6]*x[8],
    4*x[2]*x[4] + x[1]*x[7] + x[2]*x[7] + 2*x[1]*x[8] + x[2]*x[8] + x[5]*x[8] + 
        2*x[6]*x[8],
    4*x[3]*x[4] + x[1]*x[8],
    4*x[2]*x[4] + x[1]*x[7] + x[2]*x[7] + x[6]*x[7] + 2*x[1]*x[8] + x[2]*x[8] + 
        2*x[6]*x[8],
    3*x[1]*x[2] + 4*x[1]*x[3] + 3*x[3]^2 + 4*x[1]*x[5] + 2*x[1]*x[6],
    3*x[2]*x[5] + 2*x[7]^2,
    x[1]*x[2] + 2*x[2]^2 + 3*x[1]*x[3] + 4*x[2]*x[3] + 2*x[1]*x[5] + 3*x[2]*x[5]
        + x[7]*x[8],
    3*x[1]*x[2] + 4*x[2]^2 + 4*x[1]*x[3] + x[2]*x[3] + x[1]*x[5] + 3*x[2]*x[5] +
        3*x[8]^2,
    x[1]*x[2] + 2*x[2]^2 + 3*x[1]*x[3] + 4*x[2]*x[3] + 2*x[1]*x[5] + 3*x[2]*x[5]
        + x[2]*x[6]
*/
#EllipticCurve(Curve(l[1])); //8
#EllipticCurve(Curve(l[2])); //6

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
