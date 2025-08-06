32.96.5.bf.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.bf.1.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Mod 97, the first has 80 points,
while the single rank one newform -- call it E -- has 100 points. Mod 83, the third has 90 points, while E has 78 points.
This leaves us with the second curve. We construct a new model of this quotient curve using Riemann-Roch and find that it 
fails local solubility. 
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*w - y*t - z*t, 2*x^2 - y*w - z*t, 6*y^2 + 8*y*z + 8*z^2 + 2*w^2 - 2*w*t + t^2]);
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
    3*x[1]^2 - 4*x[1]*x[4] + 4*x[4]^2 + 64*x[6]^2 - 192*x[6]*x[7] + 160*x[7]^2,
    x[1]*x[2] - x[6]^2,
    x[1]*x[3] - x[7]^2,
    x[1]*x[5] - x[4]*x[6],
    4*x[4]*x[5] + 3*x[1]*x[6] - 4*x[4]*x[6] + 64*x[5]*x[8] - 128*x[6]*x[8] + 
        160*x[7]*x[8],
    4*x[4]*x[5] - 4*x[4]*x[6] + 3*x[1]*x[7] + 3*x[4]*x[7] + 64*x[5]*x[8] - 
        128*x[6]*x[8] + 160*x[7]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    64*x[2]^2 + 4*x[5]^2 + 8*x[5]*x[6] - 21*x[6]^2 + 33*x[6]*x[7] + 
        480*x[3]*x[8] - 416*x[8]^2,
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[5]*x[6],
    x[2]*x[6] - x[5]*x[8] - x[6]*x[8],
    x[2]*x[7] - x[6]*x[8],
    4*x[5]*x[6] - 8*x[6]^2 + 11*x[6]*x[7] + 64*x[2]*x[8] + 160*x[3]*x[8] - 
        192*x[8]^2,
    160*x[3]^2 + 4*x[6]^2 - 12*x[6]*x[7] + 11*x[7]^2 - 192*x[3]*x[8] + 
        64*x[8]^2,
    x[3]*x[4] - x[6]*x[7] + x[7]^2,
    x[3]*x[5] - x[6]*x[8] + x[7]*x[8],
    x[3]*x[6] - x[7]*x[8],
    -4*x[4]*x[5] + 8*x[4]*x[6] + 160*x[3]*x[7] - 11*x[4]*x[7] - 64*x[5]*x[8] + 
        192*x[6]*x[8] - 352*x[7]*x[8],
    -x[6]^2 + x[6]*x[7] + x[4]*x[8],
    -x[6]^2 + x[5]*x[7] + x[6]*x[7],
    Curve over Rational Field defined by
    162*x[1]^2 + 20736*x[3]^2 - 45568*x[6]^2 - 351*x[7]^2 + 648*x[5]*x[8] - 
        5076*x[7]*x[8] - 3072*x[8]^2,
    54*x[1]*x[2] - 432*x[2]*x[3] + 10368*x[3]^2 + 1056*x[2]*x[6] - 2048*x[6]^2 +
        81*x[5]*x[7] - 225*x[7]^2 + 216*x[5]*x[8] - 2424*x[7]*x[8] - 
        1504*x[8]^2,
    1296*x[1]*x[3] - 1296*x[2]*x[3] + 10368*x[3]^2 + 576*x[2]*x[6] + 7168*x[6]^2
        + 27*x[5]*x[7] - 135*x[7]^2 + 504*x[5]*x[8] - 2754*x[7]*x[8] + 
        1152*x[8]^2,
    24*x[1]*x[4] + 512*x[6]^2 - 3*x[5]*x[7] + 3*x[7]^2 - 4*x[5]*x[8] + 
        4*x[7]*x[8] + 48*x[8]^2,
    3*x[1]*x[5] - 48*x[3]*x[5] + 64*x[4]*x[5] - 16*x[5]*x[6] + 72*x[3]*x[7] - 
        84*x[4]*x[7] + 128*x[6]*x[7],
    8*x[1]*x[6] + 128*x[6]^2 - x[5]*x[8] + x[7]*x[8],
    -48*x[3]*x[5] + 64*x[4]*x[5] - 64*x[5]*x[6] + 3*x[1]*x[7] + 48*x[3]*x[7] - 
        64*x[4]*x[7] + 112*x[6]*x[7],
    -3*x[4]*x[5] + 4*x[5]*x[6] + 3*x[4]*x[7] - 4*x[6]*x[7] + 3*x[1]*x[8] + 
        48*x[6]*x[8],
    162*x[2]^2 - 2592*x[2]*x[3] + 31104*x[3]^2 + 2304*x[2]*x[6] + 14336*x[6]^2 +
        432*x[5]*x[7] - 729*x[7]^2 + 1152*x[5]*x[8] - 7524*x[7]*x[8] - 
        3840*x[8]^2,
    144*x[2]*x[4] - 192*x[2]*x[6] - 18*x[5]*x[7] + 9*x[7]^2 - 48*x[5]*x[8] + 
        48*x[7]*x[8] + 64*x[8]^2,
    9*x[2]*x[5] - 216*x[3]*x[5] + 264*x[4]*x[5] - 224*x[5]*x[6] + 216*x[3]*x[7] 
        - 216*x[4]*x[7] + 432*x[6]*x[7] + 384*x[6]*x[8],
    -144*x[3]*x[5] + 168*x[4]*x[5] - 160*x[5]*x[6] + 9*x[2]*x[7] + 72*x[3]*x[7] 
        - 72*x[4]*x[7] + 96*x[6]*x[7],
    -18*x[4]*x[5] - 24*x[5]*x[6] + 9*x[4]*x[7] + 36*x[6]*x[7] + 18*x[2]*x[8] + 
        64*x[6]*x[8],
    3456*x[3]*x[4] - 2048*x[6]^2 - 27*x[7]^2 - 612*x[7]*x[8] - 960*x[8]^2,
    1152*x[3]*x[6] - 512*x[6]^2 - 9*x[7]*x[8] - 192*x[8]^2,
    -9*x[4]*x[7] - 180*x[6]*x[7] + 144*x[3]*x[8] - 64*x[6]*x[8],
    72*x[4]^2 - 128*x[6]^2 - 9*x[7]*x[8] - 24*x[8]^2,
    24*x[4]*x[6] - 32*x[6]^2 - 3*x[8]^2,
    -3*x[6]*x[7] + 3*x[4]*x[8] - 4*x[6]*x[8],
    2*x[5]^2 - 4*x[5]*x[7] + 3*x[7]^2 + 8*x[7]*x[8] + 32*x[8]^2,
    Curve over Rational Field defined by
    162*x[1]^2 + 20736*x[3]^2 - 45568*x[6]^2 + 351*x[7]^2 + 648*x[5]*x[8] - 
        5076*x[7]*x[8] + 3072*x[8]^2,
    54*x[1]*x[2] + 432*x[2]*x[3] + 10368*x[3]^2 + 1056*x[2]*x[6] - 2048*x[6]^2 -
        81*x[5]*x[7] + 225*x[7]^2 + 216*x[5]*x[8] - 2424*x[7]*x[8] + 
        1504*x[8]^2,
    1296*x[1]*x[3] - 1296*x[2]*x[3] - 10368*x[3]^2 - 576*x[2]*x[6] - 7168*x[6]^2
        + 27*x[5]*x[7] - 135*x[7]^2 - 504*x[5]*x[8] + 2754*x[7]*x[8] + 
        1152*x[8]^2,
    24*x[1]*x[4] - 512*x[6]^2 - 3*x[5]*x[7] + 3*x[7]^2 + 4*x[5]*x[8] - 
        4*x[7]*x[8] + 48*x[8]^2,
    3*x[1]*x[5] + 48*x[3]*x[5] - 64*x[4]*x[5] - 16*x[5]*x[6] - 72*x[3]*x[7] + 
        84*x[4]*x[7] + 128*x[6]*x[7],
    8*x[1]*x[6] + 128*x[6]^2 - x[5]*x[8] + x[7]*x[8],
    48*x[3]*x[5] - 64*x[4]*x[5] - 64*x[5]*x[6] + 3*x[1]*x[7] - 48*x[3]*x[7] + 
        64*x[4]*x[7] + 112*x[6]*x[7],
    -3*x[4]*x[5] - 4*x[5]*x[6] + 3*x[4]*x[7] + 4*x[6]*x[7] + 3*x[1]*x[8] + 
        48*x[6]*x[8],
    162*x[2]^2 + 2592*x[2]*x[3] + 31104*x[3]^2 + 2304*x[2]*x[6] + 14336*x[6]^2 -
        432*x[5]*x[7] + 729*x[7]^2 + 1152*x[5]*x[8] - 7524*x[7]*x[8] + 
        3840*x[8]^2,
    144*x[2]*x[4] + 192*x[2]*x[6] - 18*x[5]*x[7] + 9*x[7]^2 + 48*x[5]*x[8] - 
        48*x[7]*x[8] + 64*x[8]^2,
    9*x[2]*x[5] + 216*x[3]*x[5] - 264*x[4]*x[5] - 224*x[5]*x[6] - 216*x[3]*x[7] 
        + 216*x[4]*x[7] + 432*x[6]*x[7] - 384*x[6]*x[8],
    144*x[3]*x[5] - 168*x[4]*x[5] - 160*x[5]*x[6] + 9*x[2]*x[7] - 72*x[3]*x[7] +
        72*x[4]*x[7] + 96*x[6]*x[7],
    -18*x[4]*x[5] + 24*x[5]*x[6] + 9*x[4]*x[7] - 36*x[6]*x[7] + 18*x[2]*x[8] + 
        64*x[6]*x[8],
    3456*x[3]*x[4] - 2048*x[6]^2 + 27*x[7]^2 - 612*x[7]*x[8] + 960*x[8]^2,
    1152*x[3]*x[6] + 512*x[6]^2 + 9*x[7]*x[8] - 192*x[8]^2,
    9*x[4]*x[7] - 180*x[6]*x[7] + 144*x[3]*x[8] + 64*x[6]*x[8],
    72*x[4]^2 - 128*x[6]^2 - 9*x[7]*x[8] + 24*x[8]^2,
    24*x[4]*x[6] + 32*x[6]^2 - 3*x[8]^2,
    -3*x[6]*x[7] + 3*x[4]*x[8] + 4*x[6]*x[8],
    2*x[5]^2 - 4*x[5]*x[7] + 3*x[7]^2 - 8*x[7]*x[8] + 32*x[8]^2
]

*/

#EllipticCurve(Curve(Reduction(l[1],97))); //80
#EllipticCurve(Curve(Reduction(l[3],83))); //90

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,83))); //78
#EllipticCurve(Curve(Reduction(E,97))); //100

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
yv^2 + 4096*x1v^4 - 4*yv*x2v^2 + 128*x1v^2*x2v^2 + 6*x2v^4
*/

//We complete the square to get the model:
P<x> := PolynomialRing(Rationals());
f := -(4096*x^4+128*x^2+10);
H := HyperellipticCurve(f);
HasPointsEverywhereLocally(f,2); // false


    
