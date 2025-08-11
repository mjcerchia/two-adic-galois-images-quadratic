/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.x.1. 
We find there are five genus one quotients by an involution. 
We are able to find a point on the second curve. From this, we are able to construct an elliptic curve, 
which we find to be rank 1. 
******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[2*x^3*y + x^3*z + 3*x^2*y*z + 2*x*y^3 + 3*x*y^2*z - 2*x*z^3 + y^3*z - 2*y*z^3 - z^4]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are five genus one quotients by an involution
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

l;
/*
[
    Curve over Rational Field defined by
    1156*x[1]^2 - 2312*x[1]*x[2] + 1156*x[2]^2 - 6800*x[2]*x[3] - 3600*x[3]^2 + 
        4300*x[1]*x[4] - 1700*x[2]*x[4] + 15400*x[3]*x[4] - 225*x[4]^2,
    8*x[1]*x[3] + 16*x[3]^2 - 6*x[1]*x[4] - 24*x[3]*x[4] + x[4]^2,
    Curve over Rational Field defined by
    872*x[1]^2 - 1744*x[1]*x[2] + 872*x[2]^2 - 126420*x[2]*x[3] - 457905*x[3]^2 
        - 10816*x[1]*x[4] + 32768*x[2]*x[4] + 246960*x[3]*x[4] - 13184*x[4]^2,
    872*x[1]*x[3] + 2544*x[2]*x[3] + 24031*x[3]^2 - 128*x[1]*x[4] - 
        768*x[2]*x[4] - 17056*x[3]*x[4] + 2816*x[4]^2,
    Curve over Rational Field defined by
    35836*x[1]^2 - 71672*x[1]*x[2] + 35836*x[2]^2 - 14688*x[2]*x[3] - 
        44928*x[3]^2 - 21672*x[1]*x[4] + 22644*x[2]*x[4] - 3888*x[3]*x[4] + 
        675*x[4]^2,
    248*x[1]*x[3] - 272*x[2]*x[3] - 336*x[3]^2 - 50*x[1]*x[4] + 68*x[2]*x[4] - 
        72*x[3]*x[4] - 3*x[4]^2,
    Curve over Rational Field defined by
    4*x[1]^2 + 16*x[1]*x[3] - 8*x[2]*x[3] - 17*x[3]^2 - 16*x[1]*x[4] + 
        4*x[2]*x[4] + 49*x[3]*x[4] - 16*x[4]^2,
    4*x[1]*x[2] + 25*x[1]*x[3] - 4*x[2]*x[3] - 40*x[1]*x[4] + 16*x[4]^2,
    Curve over Rational Field defined by
    4624*x[1]^2 - 4624*x[2]^2 + 11560*x[1]*x[3] - 3400*x[2]*x[3] + 2265*x[3]^2 +
        71961*x[1]*x[4] - 2890*x[2]*x[4] + 14110*x[3]*x[4] + 11849*x[4]^2,
    1156*x[1]*x[2] - 1156*x[2]^2 + 1292*x[1]*x[3] - 272*x[2]*x[3] + 273*x[3]^2 +
        8381*x[1]*x[4] + 578*x[2]*x[4] + 1802*x[3]*x[4] + 1445*x[4]^2
]
*/

//One of these quotients has the following model:
P<[x]> := ProjectiveSpace(Rationals(),3);
C1 := Curve(P,[872*x[1]^2 - 1744*x[1]*x[2] + 872*x[2]^2 - 126420*x[2]*x[3] - 457905*x[3]^2 - 
    10816*x[1]*x[4] + 32768*x[2]*x[4] + 246960*x[3]*x[4] - 13184*x[4]^2,
872*x[1]*x[3] + 2544*x[2]*x[3] + 24031*x[3]^2 - 128*x[1]*x[4] - 768*x[2]*x[4] - 
    17056*x[3]*x[4] + 2816*x[4]^2]);

//We search for a point
  rationalPoints := function(D : Bound := 1)
    return {@D![t :t in tup]
              : tup in CartesianPower([-Bound..Bound],Dimension(AmbientSpace(D)\
)+1)
              | not {i : i in tup} eq {0}
                and {Evaluate(eqns,[t : t in tup]) : eqns in
DefiningEquations(D)} eq {0}
            @};
  end function;
  rationalPoints(C1:Bound := 1);
pt :=   rationalPoints(C1:Bound := 1)[1];

E:=EllipticCurve(C1,pt);
Rank(E); //1
