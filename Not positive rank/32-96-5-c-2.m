/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.c.2.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Over F5, the second and third each
have 4 points, while the single rank one elliptic curve factor of the Jacobian has 8 points over F5. Over F7, the first has 8 points, while
the single rank one elliptic curve factor of the Jacobian has 12 points.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[-y*t + z*w, 2*x^2 - y*w + z*t, 2*y^2 + 2*z^2 + w*t]);
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;


//There are three genus one quotients by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
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

#l; //3

#EllipticCurve(Curve(Reduction((l[1]),7))); //8
#EllipticCurve(Curve(Reduction((l[2]),5))); //4
#EllipticCurve(Curve(Reduction((l[3]),5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,5))); //8
#EllipticCurve(Curve(Reduction(E,7))); //12
