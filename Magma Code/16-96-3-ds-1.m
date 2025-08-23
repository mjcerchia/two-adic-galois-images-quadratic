
//We observe an automorphism for 16.96.3.ds.1 that gives a rank 1 elliptic curve as a quotient. This is positive rank bielliptic.

P<x,y,z,w,t,u> := ProjectiveSpace(Rationals(),5);
C := Curve(P,[y^2 + 2*y*u + w*t + t^2, 
2*y^2 + w^2 + t^2,
-2*x*u + z^2, 
4*x*y + y*u + w*t,
4*x*w + 4*x*t - y*w - w*u + t*u, 
4*x*w - 4*x*t - y*t + w*u + t*u, 
32*x^2 - y^2 + 2*u^2]);




s:=iso<C->C|[u/4,y,z,-t,-w,4*x],[1/4*u,y,z,-t,-w,4*x]>;

AG := AutomorphismGroup(C,[s]);
CG,prj := CurveQuotient(AG);

Pts:=PointSearch(CG,10000);

E:=EllipticCurve(CG,Pts[1]);
Rank(E); // 1 true
