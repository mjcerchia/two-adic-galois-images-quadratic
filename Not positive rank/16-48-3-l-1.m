
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.l.1. 
Using a singular model, we find there are three genus one quotients by an involution. Apriori there could be more quotients from canonical model but these three quotients correspond to three elliptic curves in Jacobian decomposition.
All of these are pointless so this is not positive rank bielliptic.

******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[8*x^6 - 8*x^4*y^2 + 12*x^4*z^2 + 2*x^2*y^4 + 8*x^2*y^2*z^2 + 6*x^2*z^4 + y^4*z^2 - 2*y^2*z^4 + z^6]);

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
l; // It is easy to observe that all genus 1 quotients are pointless

Jacobian(l[1]); //corresponds to 128.2.a.d
Jacobian(l[2]); //corresponds to 32.2.a.a
Jacobian(l[3]); //corresponds to 128.2.a.a
