
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 64.96.5.a.1, 
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve by intersecting with all hyperplanes with coefficients 
that are either 1,-1, or 0. From this, we are able to construct an elliptic curve, 
which we find to have rank 1. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x^2 + z*w, 4*y*z - y*w + 2*t^2, 4*x*z + x*w + y^2]);

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

//One of these quotients has the following model:
P<[x]> := ProjectiveSpace(Rationals(),7);
C1 := Curve(P,[50*x[1]^2 + x[3]^2 - 25*x[7]^2,
x[1]*x[2] - x[7]^2,
x[1]*x[4] - x[3]*x[6],
x[1]*x[5] - x[4]*x[6],
x[3]*x[4] + 50*x[1]*x[6] - 25*x[7]*x[8],
-x[6]^2 + x[1]*x[7],
-x[6]*x[7] + x[1]*x[8],
25*x[2]^2 - x[5]^2 - 50*x[7]^2,
x[2]*x[3] - x[5]*x[7],
x[2]*x[4] - x[5]*x[8],
x[2]*x[6] - x[7]*x[8],
x[2]*x[7] - x[8]^2,
-x[4]*x[5] - 50*x[6]*x[7] + 25*x[2]*x[8],
x[3]*x[5] + 50*x[6]^2 - 25*x[8]^2,
-x[4]*x[6] + x[3]*x[7],
-x[5]*x[6] + x[3]*x[8],
x[4]^2 + 50*x[6]^2 - 25*x[8]^2,
-x[5]*x[6] + x[4]*x[7],
-x[5]*x[7] + x[4]*x[8],
-x[7]^2 + x[6]*x[8] ]);

//We can't immediately find a point, so we intersect with hyperplanes to find a point
P<[X]> := AmbientSpace(C1);
deg := 1;
tups := {};
d := Dimension(P);
ideals := {@@};
tupsOld := tups;
coeffs := [ 0, 1 ];
    tups   :=
        [ P![ tup[i] : i in [1..#tup] ] :
          tup in CartesianPower(coeffs, d+1) |
          not tup eq < 0 : i in [1..d+1] >  
          and not P![tup[i] : i in [1..#tup]] in tupsOld
        ];
    for c in tups  do  
            L := Scheme(P,&+[c[i]*X[i] : i in [1..Dimension(AmbientSpace(C1)) + 1]]);
              irr := IrreducibleComponents(L meet C1);
              degs := [Degree(ReducedSubscheme(i)) : i in irr];
              if  deg in degs then
        for cpt in irr do  
         if  Degree(ReducedSubscheme(cpt)) eq deg
             and not Basis(Ideal(ReducedSubscheme(cpt))) in ideals then      
               ideals := ideals join {@Basis(Ideal(ReducedSubscheme(cpt)))@}; 
        end if;
    end for;
    end if;
    end for;
    ideals;

//This produces a point
pt := C1!Points(Scheme(AmbientSpace(C1),ideals[1]))[1];
E := EllipticCurve(C1,pt);
Rank(E); //1
