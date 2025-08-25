/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 37.114.4.b.1
We compute the automorphism group over F5 and find there is one genus one quotient by an involution. It has 6 points.
The single rank one elliptic curve factor of the Jacobian has 8 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w>:=ProjectiveSpace(Rationals(),3);
C := Curve(P,[x^2 - x*z - x*w - y*z + z^2 + z*w, x^3 - x^2*y + x*y^2 - x*z*w - y^2*z - y^2*w + 2*y*w^2 + z^2*w - w^3]);
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

#l;//1

#EllipticCurve(Curve(l[1])); //6

E := EllipticCurve([0,0,1,-1,0]);
#EllipticCurve(Curve(Reduction(E,5))); //8
