/* Summary : We compute the number of genus 1 quotients mod 5 and there are three. 
We know that there are at least 3 quotients over Q, since we only have three quotients over F_5 we know these are all. Further, LMFDB gives all those quotients
they are either pointless or rank 0*/

P<x,y,z>:=ProjectiveSpace(Rationals(),2);
C:=Curve(P,[x^4 + 4*x^2*z^2 + 2*y^4 + 2*z^4]);
p:=5;
Cp:=Curve(Reduction(C,p));
G:=AutomorphismGroup(Cp);
S:=Automorphisms(Cp);
for s in S do;
s1:=G!s;
if Order(s1) eq 2 then
AG:=AutomorphismGroup(Cp,[s1]);
CG:=CurveQuotient(AG);
if Genus(CG) eq 1 then
CG;
end if;
end if;
print "........";
end for;
