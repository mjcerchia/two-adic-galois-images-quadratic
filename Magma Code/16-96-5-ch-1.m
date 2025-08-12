/* Summary : We compute the number of genus 1 quotients mod 5 and there are three. 
We know that there are at least 3 quotients over Q, since we only have three quotients over F_5 we know these are all. Further, LMFDB gives all those quotients
they are either pointless or rank 0*/

P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C:=Curve(P,[z^2 + w^2 - t^2, 2*x^2 - y*z, 8*y^2 + w^2 + t^2]);
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

