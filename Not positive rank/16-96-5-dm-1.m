/* Summary : We compute the number of genus 1 quotients mod 5 and there are two. 
No quotient corresponds to 256.2.a.b*/

P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C:=Curve(P,[y^2 + z^2 - z*w, 2*x*z + 2*x*w + t^2, 8*x^2 - z^2 + 2*w^2]);
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
E:=EllipticCurve(CG);
#Points(E); // 6 and 4
end if;
end if;
print "........";
end for;

// For 256.2.a.b

E := EllipticCurve([0, 0, 0, -2, 0]);
Cp:=Curve(Reduction(E,5));
Ep:=EllipticCurve(Cp);
#Points(Ep); // 10
