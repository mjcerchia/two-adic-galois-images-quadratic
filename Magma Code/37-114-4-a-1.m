/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 37.114.4.a.1. 
We find there are four genus one quotients by an involution. 
One of these is isomorphic to 1369.d2, which is a rank 1 elliptic curve. 

******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[778*x^6 - 2294*x^5*y + 547*x^5*z + 814*x^4*y^2 + 3367*x^4*y*z - 3513*x^4*z^2 - 962*x^3*y^3 - 14023*x^3*y^2*z + 18907*x^3*y*z^2 - 8863*x^3*z^3 - 259*x^2*y^3*z - 24235*x^2*y^2*z^2 + 31413*x^2*y*z^3 - 8564*x^2*z^4 - 11063*x*y^3*z^2 - 24013*x*y^2*z^3 + 19832*x*y*z^4 - 4944*x*z^5 - 3145*y^3*z^3 - 4292*y^2*z^4 + 5328*y*z^5 - 1664*z^6]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are four genus one quotients by an involution
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
C1:=l[3];

//This is quite complicated, so we search for divisors to construct a simpler model.
P<[X]> := AmbientSpace(C1);
deg := 3;
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

//Using a degree three divisor found above, we build a new model using Riemann-Roch.
gens := ideals[1];
P7 := AmbientSpace(C1);
R  := CoordinateRing(P7);
I := ideal< R | gens >;
Pl := Place(C1, I);
D  := Divisor(Pl);
a,b := RiemannRochSpace(D);
x1 := b(a.1); 
x2 := b(a.2);
x3:=b(a.3); 
F := BaseField(C1);
P := PolynomialRing( F, [1,1,1], "grevlexw", [1,1,1] );
vars := ["x", "y","z"]; 
P := ProjectiveSpace(P);
AssignNames(~P, vars);
phi := map<C1 -> P | [x1,x2,x3]>;
C2:=Image(phi); 
C2;
/*
Curve over Rational Field defined by
x^3 - 14224201279478138869386111025040519028063317452888544517*x^2*y - 
    166056007867155068360946856005971698912073216968525233271647680597257199683\
    849883336460184366314037670945069184*x*y^2 - 
    539446075020546776921821746979826630279898922993215885245938082631991584859\
    549472794373687108857495543773928974613811274896474606671941884467487504298\
    229549441339056*y^3 + 11314449924730929673793457923983490787070679872150760\
    724*x^2*z + 294139694211657937359813252965782349776038668475917254742462394\
    499268748354413029645022194196826201759437286432*x*y*z + 
    180709352001489078451630110508684207892604064115490953805061776099863324536\
    381272307864275691527202581168773109634115675577817421424193110449390260213\
    9669751578314048*y^2*z - 12876462370162805303884530692771972383064451390445\
    9687507634367989499331001754131838138298654933136466930037632*x*z^2 - 
    188836816028363905211866688336005565802967907579554156303716381136008761226\
    097116135907426348196248347670706886245594433863164779296012340943665296379\
    4132513311143424*y*z^2 + 63234308121841064732829426151863027949988747119466\
    785472671237377533668044197887134413025619140501078557860641887720577551361\
    4111483400190167849026813279145775648768*z^3
*/

E := EllipticCurve(C2);

//Magma realizes C2 as an elliptic curve, but it is difficult to find a point or compute the rank. 
//Instead, we show that it is isomorphic to a rank 1 elliptic curve.

E1 := EllipticCurve([1, -1, 0, 3166, -59359]);
IsIsomorphic(E,E1); //true
    
