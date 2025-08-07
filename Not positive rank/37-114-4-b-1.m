/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 37.114.4.b.1
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 3 points.
The single rank one elliptic curve factor of the Jacobian has 7 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[x^2 - x*z - x*w - y*z + z^2 + z*w, x^3 - x^2*y + x*y^2 - x*z*w\
 - y^2*z - y^2*w + 2*y*w^2 + z^2*w - w^3]);
C3 := Curve(Reduction(C,3)); 
S := AutomorphismGroup(C3); 

auts := [];
Stemp := Automorphisms(C3);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotients by an involution
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
/*
[
    Curve over GF(3) defined by
    2*x[1]*x[2] + x[2]^2 + 2*x[2]*x[4] + x[2]*x[5] + x[3]*x[5] + x[1]*x[6] + 
        x[2]*x[6] + x[5]*x[6],
    x[1]*x[2] + 2*x[2]^2 + 2*x[2]*x[3] + 2*x[1]*x[4] + 2*x[2]*x[4] + x[1]*x[5] +
        x[1]*x[6] + x[2]*x[6] + x[3]*x[6],
    x[1]*x[2] + 2*x[2]^2 + x[2]*x[3] + x[2]*x[4] + x[4]^2 + x[1]*x[5] + 
        2*x[2]*x[5] + x[1]*x[6] + x[5]*x[6],
    2*x[1]*x[2] + x[2]^2 + 2*x[2]*x[3] + 2*x[2]*x[4] + 2*x[1]*x[5] + x[4]*x[5] +
        x[1]*x[6] + x[2]*x[6],
    x[1]*x[2] + 2*x[2]^2 + x[1]*x[4] + x[1]*x[5] + x[2]*x[5] + 2*x[1]*x[6] + 
        x[4]*x[6],
    x[1]*x[2] + 2*x[2]^2 + 2*x[2]*x[3] + 2*x[1]*x[4] + 2*x[2]*x[4] + x[1]*x[5] +
        x[5]^2 + x[1]*x[6] + x[2]*x[6],
    x[1]*x[2] + 2*x[2]^2 + x[2]*x[3] + x[2]*x[4] + x[3]*x[4] + x[2]*x[5] + 
        2*x[1]*x[6] + 2*x[2]*x[6] + x[5]*x[6],
    x[2]*x[3] + 2*x[1]*x[4] + x[2]*x[4] + 2*x[2]*x[5] + 2*x[1]*x[6] + x[6]^2,
    x[1]^2 + x[1]*x[2] + x[2]^2 + x[1]*x[3] + x[2]*x[3] + x[3]^2 + x[2]*x[5] + 
        2*x[2]*x[6]
]
*/
#EllipticCurve(Curve(l[1])); //6

E := EllipticCurve([0,0,1,-1,0]);
#EllipticCurve(Curve(Reduction(E,3))); //7
