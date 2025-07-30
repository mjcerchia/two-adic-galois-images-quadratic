
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.l.1. 
We find there are three genus one quotients by an involution. 
Over F3, each of these have a different number of points, so they can't be isogenous.
There are only three factors of the Jacobian up to isogeny.
One of the quotients must be isogenous to the rank 1 factor of the Jacobian. 
******************************************************************************/
P<x,y,z,u,w,t> := ProjectiveSpace(Rationals(),5);
C := Curve(P,[2*x*t + z^2 + w*u, -2*x*w - x*u + y*z, -2*x*t + 2*y^2 + w*u, 4*w^2 + 2*t^2 + u^2, 4*x^2 - x*t + y^2, 4*x*z + 2*y*w + y*u - z*t, 8*x*y + 2*y*t + 2*z*w + z*u]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are three genus one quotients by an involution
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

//These are different from each other:
#EllipticCurve(Curve(Reduction(l[1],3))); 
#EllipticCurve(Curve(Reduction(l[2],3))); 
#EllipticCurve(Curve(Reduction(l[3],3))); 
