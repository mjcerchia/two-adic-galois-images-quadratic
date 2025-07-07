//////////////////////////////////////////////////////////////////////////////////////////
//  General code to find all quotients by involutions of a modular curve associated to 
//  subgroup G of GL2(Z_hat). 
//////////////////////////////////////////////////////////////////////////////////////////

GL2:=GL(2,Integers(16));
//G:=sub<GL2|[[],[],[]]>
G:=sub<GL2|[[3,2,4,1],[13,5,10,3],[15,14,12,13]]>;
Grec:=CreateModularCurveRec(16,Generators(G));
Model:=FindModelOfXG(Grec,67)`psi;

//x:=AssociativeArray();
P2:=ProjectiveSpace(Rationals(),8);
//for i in [1..6] do;
//    x[i]:=P2.i;
//end for;
C:=Curve(P2,Model);
// C5 := Curve(Reduction(C,5)); 
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C5);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

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
