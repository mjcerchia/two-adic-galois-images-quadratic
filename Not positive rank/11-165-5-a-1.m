11.165.5.a.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 11.165.5.a.1.
We compute the automorphism group over F5 and find there are no genus one quotient by an involution. 
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[3*x^2 + 3*x*y - x*z - 3*x*w + 3*x*t + 3*y^2 + 2*y*z - y*w - 2*\
z^2 - 2*z*w + 2*z*t + w^2, 2*x*y - 2*x*z + 6*x*w - x*t - 3*y^2 + 3*y*z - y*t +\
 2*z^2 - 7*z*w - z*t - t^2, 8*x^2 - 10*x*y + 5*x*z - 3*x*w - x*t - 3*y^2 + 3*y\
*z + 4*y*w + 3*y*t - 3*z^2 - 7*z*w - z*t + 3*w^2 + 3*w*t - t^2]);

C5 := Curve(Reduction(C,5));
S := AutomorphismGroup(C5); 

auts := [];
Stemp := Automorphisms(C5);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are no genus one quotients by an involution
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

