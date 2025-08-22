
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.l.1. 
Using a singular model, we find there are three genus one quotients by an involution. 
We construct a simpler model of the third quotient and find it is a rank one elliptic curve.

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

C1 := l[1];

//Searching for degree two divisors
P<[X]> := AmbientSpace(C1);
deg := 2;
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

/*
{@
    [
        X[1]^2 + 32768*X[4]^2,
        X[2],
        X[3]
    ],
    [
        X[3]^2 + 256*X[4]^2,
        X[1],
        X[2]
    ],
    [
        X[3]^2 + 32768/145*X[4]^2,
        X[1] - X[3],
        X[2] + X[3]
    ]
@}
*/

//building a simpler model
gens := ideals[1];
P7 := AmbientSpace(C1);
R  := CoordinateRing(P7);
I := ideal< R | gens >;
Pl := Place(C1, I);
D  := Divisor(Pl);
a,b := RiemannRochSpace(D);
x1 := b(a.1); 
x2 := b(a.2); 
a,b := RiemannRochSpace(2*D);
y := b(a.1);
F := BaseField(C1);
P := PolynomialRing( F, [2,1,1], "grevlexw", [2,1,1] );
vars := ["yv","x1v", "x2v"]; //(calling them yv and xv so magma doesn't confuse the functions with x and y)
P := ProjectiveSpace(P);
AssignNames(~P, vars);
phi := map<C1 -> P | [y,x1,x2]>;
C2:=Image(phi); 

C2;
/*
Curve over Rational Field defined by
yv^2 + 2*x1v^4 + 32*x1v^2*x2v^2 + 256*x2v^4
*/

P<x> := PolynomialRing(Rationals());
f := -(2*x^4 + 32*x^2*x + 256*x);
H := HyperellipticCurve(f);
pt := H!RationalPoints(H : Bound := 10000)[1];

E := EllipticCurve(H,pt);
Rank(E); //1

