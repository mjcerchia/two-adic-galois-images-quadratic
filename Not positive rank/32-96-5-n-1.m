32.96.5.n.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.n.1.
We compute the automorphism group over Q and find there are three genus one quotients by an involution. Mod 5, the first has 4 points and the third
has 2 points, while the single rank one newform -- call it E -- has 10 points. 
This leaves us with the second curve. We construct a new model of this quotient curve using Riemann-Roch and find that it 
fails local solubility. 
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[x^2 + z*w, x*t + y^2, z^2 + 16*w^2 + t^2]);
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

l;
\*
[
    Curve over Rational Field defined by
    x[1]^2 - 32768*x[4]*x[6] + 1024*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[3] + 512*x[6]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] + 512*x[6]*x[8],
    64*x[3]*x[4] + x[1]*x[6] + 1024*x[7]*x[8],
    x[1]*x[7] + 1024*x[2]*x[7] - 32768*x[4]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    1024*x[2]^2 + 64*x[5]^2 + x[7]^2,
    x[2]*x[3] + 512*x[8]^2,
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    64*x[4]*x[5] + x[6]*x[7] + 1024*x[2]*x[8],
    x[3]^2 + 512*x[4]*x[6],
    x[3]*x[5] + 512*x[4]*x[8],
    x[3]*x[7] + 512*x[6]*x[8],
    -x[5]*x[6] + x[3]*x[8],
    64*x[4]^2 + x[6]^2 + 1024*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    x[5]*x[7] + 512*x[8]^2,
    Curve over Rational Field defined by
    128*x[1]^2 + x[3]^2 + 16*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - x[4]*x[6],
    x[3]*x[4] + 128*x[1]*x[6] + 16*x[7]*x[8],
    -x[6]^2 + x[1]*x[7],
    -x[6]*x[7] + x[1]*x[8],
    16*x[2]^2 + x[5]^2 + 128*x[7]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    x[2]*x[7] - x[8]^2,
    x[4]*x[5] + 128*x[6]*x[7] + 16*x[2]*x[8],
    x[3]*x[5] + 128*x[6]^2 + 16*x[8]^2,
    -x[4]*x[6] + x[3]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    x[4]^2 + 128*x[6]^2 + 16*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    -x[5]*x[7] + x[4]*x[8],
    -x[7]^2 + x[6]*x[8],
    Curve over Rational Field defined by
    128*x[1]^2 - x[3]^2 - 16*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - x[4]*x[6],
    -x[3]*x[4] + 128*x[1]*x[6] - 16*x[7]*x[8],
    -x[6]^2 + x[1]*x[7],
    -x[6]*x[7] + x[1]*x[8],
    16*x[2]^2 + x[5]^2 - 128*x[7]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    x[2]*x[7] - x[8]^2,
    x[4]*x[5] - 128*x[6]*x[7] + 16*x[2]*x[8],
    x[3]*x[5] - 128*x[6]^2 + 16*x[8]^2,
    -x[4]*x[6] + x[3]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    x[4]^2 - 128*x[6]^2 + 16*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    -x[5]*x[7] + x[4]*x[8],
    -x[7]^2 + x[6]*x[8]
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4
#EllipticCurve(Curve(Reduction(l[3],5))); //2

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10


C1 := l[2];

//Search for degree 2 divisors on C1
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
    
ideals[1];
/*
[
    X[1]^2 + 1/128*X[3]^2,
    X[2],
    X[4],
    X[5],
    X[6],
    X[7],
    X[8]
]
*/


//build a new model
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
yv^2 + 128*x1v^4 + 16*x2v^4
*/

P<x> := PolynomialRing(Rationals());
f := -(128*x^4 + 16);
H := HyperellipticCurve(f);
HasPointsEverywhereLocally(f,2); // false


    
