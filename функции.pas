 type mas = array[1..3] of real;
 var jd: real;

procedure REDUCE(A:real;VAR B:real);
      CONST PI2=180*2;
      BEGIN
      B:=A-trunc(A/PI2)*PI2;
      IF A<0 THEN B:=B+PI2;
      END;

function Sign(n:real):Shortint;
begin
if n>0 then Sign:=1
else if n<0 then Sign:=-1
else Sign:=0;
end; {Sign}

function ArcTg(x,y: real):real;
var a:real;
begin

if abs(y)<1e-18 then ArcTg:=sign(x)*0.5*Pi
else
begin
a:=ArcTan(x/y);
if y<0 then a:=a+Pi;
if a<0 then a:=a+2*Pi;
ArcTg:=a;
end;
end; {ArcTg}

function ArcSin(x: real):real;
{ arcsin in [Pi/2,-Pi/2] }
begin
 if (x<1.0)and(x>-1.0) then
    Arcsin:=ArcTan(x/sqrt(1-sqr(x)))
 else if (abs(x)-1.0)<1e-18 then ArcSin:=Sign(x)*Pi/2
 else begin Writeln('Error: Argument of ArcSin >1.0 or <-1.0',x); Halt; end;
end; {Arcsin}

function date_jd(year,month:integer;day:real):real;
  label 1;
  const han=100;
  var m1,d,date,jd:real;
      i,me,ja,jb:integer;
  begin
  date:=year+month/100+day/1e4;
  i:=trunc(date);
  m1:=(date-i)*han;
  me:=trunc(m1);
  d:=(m1-me)*han;
  if me>2 then goto 1;
  i:=i-1;
  me:=me+12;
1:jd:=trunc(365.25*i)+trunc(30.6001*(me+1))+d+1720994.5;
  if date<1582.1015 then
    begin
    date_jd:=jd;
    exit;
    end;
  ja:=trunc(i/100.0);
  jb:=2-ja+trunc(ja/4);
  jd:=jd+jb;
  date_jd:=jd;
  end;

function sid2000(jd:real):real; {v radianah}
const jd2000=2451545;
      jdyear=36525;
var m,mm,d,t,s,sr :real;
  begin
  m:=frac(jd)-0.5;
  d:=jd-m-jd2000;

  t:=(d+m)/jdyear;
  mm:=m*86400;
  s:=(24110.54841+mm+236.555367908*(d+m)+(0.093104*t-6.21E-6*sqr(t))*t)/86400*2*pi;
  sid2000:=s;
  end;

procedure perehod(HH:real; xx:mas; var lambda,phi:real);
 var A:array[1..3,1..3]of real;
     i,j:integer;
     s,r,argum:real;
     yy:mas;
 begin

  A[1,1]:=cos(HH);     A[1,2]:=sin(HH);    A[1,3]:=0;
  A[2,1]:=-sin(HH);    A[2,2]:=cos(HH);    A[2,3]:=0;
  A[3,1]:=0;           A[3,2]:=0;          A[3,3]:=1;

  for i:=1 to 3 do
   begin
    s:=0;
    for j:=1 to 3 do
      s:=s+A[i,j]*xx[j];
    yy[i]:=s;
  end;
  r:=sqrt(yy[1]*yy[1]+yy[2]*yy[2]+yy[3]*yy[3]);

  lambda:=arctg(yy[2],yy[1]);

  argum:=yy[3]/r;
  phi:=arcsin(argum);
 end;

begin
  jd := date_jd(2023, 10, 20);
  writeln(sid2000(jd));
  writeln(jd);
end.