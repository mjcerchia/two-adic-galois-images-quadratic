
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.l.1. 
Using a singular model, we find there are three genus one quotients by an involution. 

******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[8*x^6 - 8*x^4*y^2 + 12*x^4*z^2 + 2*x^2*y^4 + 8*x^2*y^2*z^2 + 6*x^2*z^4 + y^4\
*z^2 - 2*y^2*z^4 + z^6]);

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

/*
[
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    1/32768*x[1]^2 + 1/4096*x[2]^2 + 1/1024*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    1/256*x[1]^2 + 1/1024*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    1/32768*x[1]^2 + 1/2048*x[2]^2 + 1/256*x[3]^2 + x[4]^2
]
*/


