/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 37.114.4.b.2
We compute the automorphism group over F5 and find there is one genus one quotient by an involution. It has 6 points.
The single rank one elliptic curve factor of the Jacobian has 8 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w>:=ProjectiveSpace(Rationals(),3);
C := Curve(P,[12*x^2 - 13*x*y - x*z + x*w + 12*y^2 - y*z + y*w - z^2 + z*w - w^2, 160*x^2*y - 21*x^2*z + 44*x^2*w + 172*x*y^2 + 31*x*y*z - 22*x*y*w - 23*x*z^2 + 16*x*z*w - 55*x*w^2 + 49*y^3 + 89*y^2*z - 66*y^2*w + 87*y*z^2 - 94*y*z*w + 55*y*w^2 + 10*z^3 - 39*z^2*w + 6*z*w^2 + 11*w^3]);
C5 := Curve(Reduction(C,5)); 
S := AutomorphismGroup(C5); 

auts := [];
Stemp := Automorphisms(C5);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;


//There is one genus one quotient by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C5,[g]);
CG,prj := CurveQuotient(AG);
if Genus(CG) eq 1 then
try
l := Append(l,CG);
catch e
m := Append(m,CG);
end try;
end if;

end if;

end for;

#l; //1

#EllipticCurve(Curve(l[1])); //6

E := EllipticCurve([0,0,1,-1,0]);
#EllipticCurve(Curve(Reduction(E,5))); //8
