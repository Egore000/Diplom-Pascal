unit TwoBody;

interface
const
      k = pi/180;
      g = 180/pi;
      
      N_round = 100;
      N_ONE_STEP = 30;
      N_STEP = N_round*N_ONE_STEP;
      
      // ///orbit elements
      // M0 = 282.3536014404697;
      ecc = 1e-3;
      // i = 0.2291831180523293;
      // a = 42166.473;
      // Omega = 0;
      // w = 0;
      
      ///constants
      mu = 3.986004418e+5;
      // mu = 1.32712442099e+11;
      pi = 3.1415926535897932;
      eps = 1e-12;
      
      ///time
      t0 = 0; 
      dt = N_round*43077.61;
     
      
type matrix = array[1..3,1..3] of extended;
     mas = array[1..3] of extended; 
     
var t, h, n: extended;
    a0, e, Omega0, w0, M, i0: extended;
    X, V, Orbit: mas;
    data: text;

function Sign(n:extended):Shortint;
function ArcTg(x,y: extended):extended;   
function Arctg2(x, y: extended):extended;
function f(x, M: extended): extended;
function df(x: extended): extended;
procedure TwoPoints(t, n, M0, w, Omega, i, a:extended; var X, V: mas);
procedure CoordsToElements(Coords, Velocities: mas; mu: extended; 
                           var a, e, i, Omega, w, M: extended);

implementation

function Sign(n:extended):Shortint;
  begin
    if (n>0) then 
      Sign := 1
    else if (n<0) then 
      Sign := -1
      else Sign := 0;
  end; {Sign}

function ArcTg(x,y: extended):extended;
{Arctang alpha [0,2*Pi];
x = c*sin(alpha)
y = c*cos(alpha),
}
var a: extended;
begin
  if (abs(y)<1e-18) then 
    ArcTg := sign(x)*0.5*Pi
  else
  begin
    a := ArcTan(x/y);
    if (y<0) then
      a := a+Pi;
    if (a<0) then 
      a := a+2*Pi;
    ArcTg := a;
  end;
end; {ArcTg}    


function Arctg2(x, y: extended):extended;
var a: extended;
begin
   if (abs(y)<1e-18) then
      if (x>0) then  
        Arctg2 := sign(x)*0.5*pi
      else
        Arctg2 := -sign(x)*0.5*pi
  else
  begin
    a := ArcTan(x/y);
    if (y>0) then
      
    else if (x>=0) then 
            a := a + pi
         else
            a := a - pi;
   
    Arctg2 := a;
  end;
end; {Arctg2}



function f(x, M: extended): extended;
begin
  f := x - ecc*sin(x) - M;
end; 

function df(x: extended): extended;
begin
  df := 1 - ecc*cos(x);
end;

procedure TwoPoints(t, n, M0, w, Omega, i, a:extended; var X, V: mas);
var M, E, E0, v0, u, dif, sum: extended;
    parametr, alpha0, beta0, gamma0, alpha, beta, gamma: extended;
    i1,i2,i3,i4: integer;
    Z1, Z2, Xi, prod: matrix;
begin
  
  for i1 := 1 to 3 do
  begin  
    X[i1] := 0;
    for i2 := 1 to 3 do
    begin
      Z1[i1,i2] := 0;
      Z2[i1,i2] := 0;
      Xi[i1,i2] := 0;
      prod[i1,i2] := 0;
    end;
  end;
  
  M := n*(t - t0) + M0*k;
  
  E0 := M;
  dif := 1;
  
  while (abs(dif)>eps) do
  begin
    E := E0 - f(E0,M)/df(E0);
    dif := E - E0;
    E0 := E;
  end;
  
  v0 := ArcTg((sqrt(1-sqr(ecc))*sin(E)),(cos(E)-ecc));
  
  Orbit[1] := a*(cos(E) - ecc);
  Orbit[2] := a*sqrt(1-sqr(ecc))*sin(E);
  Orbit[3] := 0;
  
  Z1[1,1] := cos(Omega*k);  Z2[1,1] := cos(w*k);
  Z1[1,2] := -sin(Omega*k); Z2[1,2] := -sin(w*k);
  Z1[2,1] := sin(Omega*k);  Z2[2,1] := sin(w*k);
  Z1[2,2] := cos(Omega*k);  Z2[2,2] := cos(w*k);
  Z1[3,3] := 1;             Z2[3,3] := 1; 
      
  Xi[1,1] := 1;
  Xi[2,2] := cos(i*k);
  Xi[2,3] := -sin(i*k);
  Xi[3,2] := sin(i*k);
  Xi[3,3] := cos(i*k);
  
  for i1 := 1 to 3 do
    for i2 := 1 to 3 do
      for i3 := 1 to 3 do
        for i4 := 1 to 3 do
          prod[i1, i4] := prod[i1, i4] + Z1[i1, i2]*Xi[i2, i3]*Z2[i3, i4];      
        
  for i1 := 1 to 3 do
  begin
    sum := 0;
    for i2 := 1 to 3 do  
      sum := sum + prod[i1,i2] * orbit[i2];
    X[i1] := sum;  
  end;
 
  parametr := a*(1 - sqr(ecc));
  
  u := v0 + w*k;
  
  alpha0 := cos(u)*cos(Omega*k) - sin(u)*sin(Omega*k)*cos(i*k);
  alpha := -sin(u)*cos(Omega*k) - cos(u)*sin(Omega*k)*cos(i*k);
  
  beta0 := cos(u)*sin(Omega*k) + sin(u)*cos(Omega*k)*cos(i*k);
  beta := -sin(u)*sin(Omega*k) + cos(u)*cos(Omega*k)*cos(i*k);
  
  gamma0 := sin(u)*sin(i*k);
  gamma := cos(u)*sin(i*k);
  
  V[1] := sqrt(mu/parametr)*(ecc*sin(v0)*alpha0 + (1 + ecc*cos(v0))*alpha);
  V[2] := sqrt(mu/parametr)*(ecc*sin(v0)*beta0 + (1 + ecc*cos(v0))*beta);
  V[3] := sqrt(mu/parametr)*(ecc*sin(v0)*gamma0 + (1 + ecc*cos(v0))*gamma);
end;



procedure CoordsToElements(Coords, Velocities: mas; mu: extended;
                           var a, e, i, Omega, w, M: extended);
var x, y, z, Vx, Vy, Vz: extended;
    r, V2, h, c1, c2, c3, l1, l2, l3, E0: extended;
    c, l: extended;
begin
    x := Coords[1];
    y := Coords[2];
    z := Coords[3];
    Vx := Velocities[1];
    Vy := Velocities[2];
    Vz := Velocities[3];

    r := sqrt(x*x + y*y + z*z);
    V2 := Vx*Vx + Vy*Vy + Vz*Vz;
    
    h := V2/2 - mu/r;
    
    c1 := y*Vz - z*Vy;
    c2 := z*Vx - x*Vz;
    c3 := x*Vy - y*Vx;
    
    l1 := -mu*x/r + Vy*c3 - Vz*c2;
    l2 := -mu*y/r + Vz*c1 - Vx*c3;
    l3 := -mu*z/r + Vx*c2 - Vy*c1;
    
    c := sqrt(c1*c1 + c2*c2 + c3*c3);
    l := sqrt(l1*l1 + l2*l2 + l3*l3);
    
    a := -mu/(2*h);
    e := l/mu;
    i := ArcTg(sqrt(1 - (c3*c3/(c*c))), c3/c);

    if (i = 0) then i := i + 1e-12;

    Omega := Arctg2(c1/(c*sin(i)), -c2/(c*sin(i)));
    w := Arctg2(l3/(l*sin(i)), l1*cos(Omega)/l + l2*sin(Omega)/l);
    
    E0 := Arctg2((x*Vx + y*Vy + z*Vz)/(e*sqrt(mu*a)), (1-r/a)/e);

    M := E0 - e*sin(E0);
      
end;  

begin {Main}
end.