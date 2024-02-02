unit readfond;
{$J+}
interface
type masc=array[1..66] of extended;
const mno=180/pi;
      e0=84381.448/3600/mno;

procedure read200(ecl,vel:integer;t:extended;var x:masc);
procedure read406(ecl,vel:integer;t:extended;var x:masc);
procedure read405(ecl,vel:integer;t:extended;var x:masc);

implementation

type mas826=array[1..826] of double;
var  de200:file of mas826;
     buf1:mas826;
     sie0,coe0:extended;

procedure read200(ecl,vel:integer;t:extended;var x:masc);
type vec=array[1..11] of integer;
     mas15=array[0..15] of extended;
const ratl:extended=1/(81.300587+1);
      nper:vec=(3,147,183,273,303,330,354,378,396,414,702);
      pow:vec=(11, 11, 14,  9,  8,  7,  7,  5,  5, 11, 14);
 {             12  12  15  10   9   8   8   6   6  12  15    10     0
                4   1   2   1   1   1   1   1   1   8   1     4     0}
      ae=149597870.66;
      tmin=2305424.5;
      tmax=2513392.5;
      nrc:longint=-1;
      ntp200:boolean=true;

label 1;
var tc,tcp:mas15;
    a,b,c,d,t0,s,s1,s2:extended;
    i,j,k,ii,jj,ik,ist,n1,nr:integer;
    xecl:masc;

procedure cheb(vel:integer; a,b,t:extended;st:integer;var tc,tcp:mas15);
var i:integer;
    tau,tau2,dlin,rat:extended;
begin dlin:=b-a;
tau:=2*(t-a)/dlin-1; tau2:=tau*2;
tc[0]:=1;
tc[1]:=tau;
for i:=2 to st do tc[i]:=tau2*tc[i-1]-tc[i-2];
if vel<>0 then
  begin
  rat:=4/dlin;
  tcp[0]:=0;
  tcp[1]:=2/dlin;
  for i:=2 to st do tcp[i]:=rat*tc[i-1]+tau2*tcp[i-1]-tcp[i-2];
  end;
end;
procedure coor(vel:integer; i,n1:integer;tc,tcp:mas15;var x:masc);
var k,j,ii,jj:integer;
    s,s1:extended;
begin ii:=6*i-6;ist:=pow[i]+1;
for k:=1 to 3 do
  begin s:=0; ik:=ii+k;
  jj:=nper[i]+(k-1)*ist+n1;
  for j:=0 to pow[i] do s:=s+tc[j]*buf1[jj+j];
  x[ik]:=s;
  end;
if vel<>0 then
for k:=1 to 3 do
  begin s:=0; ik:=ii+k;
  jj:=nper[i]+(k-1)*ist+n1;
  s:=0;
  for j:=1 to pow[i] do s:=s+tcp[j]*buf1[jj+j];
  x[ik+3]:=s;
  end
end;
begin
if (t<tmin) or (t>tmax) then
  begin writeln('��� �� �।����� 䮭��!');
  writeln('tmin=',tmin:10:1);
  writeln('tmax=',tmax:10:1);readln;
  halt;
  end;

if ntp200 then
  begin assign(de200,'d:\users\de200\16002169.fpd');
  reset(de200);
  ntp200:=false;
  end;

nr:=trunc((t-tmin)/32);

if nr<>nrc then
  begin
1:seek(de200,nr);
  read(de200,buf1);
  if t<buf1[1] then begin nr:=nr-1; goto 1 end;
  if t>buf1[2] then begin nr:=nr+1; goto 1 end;
  nrc:=nr;
  end;
a:=buf1[1]; b:=buf1[2]; n1:=0;
cheb(vel,a,b,t,14,tc,tcp);
for i:=2 to 11 do if(i<>3)and(i<>10)then
coor(vel,i,n1,tc,tcp,x);

c:=a-8;d:=a;n1:=-36;
repeat c:=c+8;d:=d+8;n1:=n1+36; until (t>=c)and(t<=d);
cheb(vel,c,d,t,pow[1],tc,tcp); coor(vel,1,n1,tc,tcp,x);

c:=a-16;d:=a;n1:=-45;
repeat c:=c+16;d:=d+16;n1:=n1+45 until (t>=c)and(t<=d);
cheb(vel,c,d,t,pow[3],tc,tcp); coor(vel,3,n1,tc,tcp,x);

c:=a-4;d:=a;n1:=-36;
repeat c:=c+4;d:=d+4;n1:=n1+36 until (t>=c)and(t<=d);
cheb(vel,c,d,t,pow[10],tc,tcp); coor(vel,10,n1,tc,tcp,x);

for k:=1 to 3 do x[12+k]:=x[12+k]-ratl*x[54+k];
for k:=1 to 3 do x[54+k]:=x[12+k]+x[54+k];
if vel<>0 then
  begin
  for k:=4 to 6 do x[12+k]:=x[12+k]-ratl*x[54+k];
  for k:=4 to 6 do x[54+k]:=x[12+k]+x[54+k];
  end;
for i:=1 to 10 do
 begin ii:=6*i-6;
 for k:=1 to 3 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
 if vel<>0 then for k:=4 to 6 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
 end;
if ecl=1 then
  begin
  for i:=1 to 10 do
     begin ii:=6*i-6;
     xecl[ii+1]:=x[ii+1];
     xecl[ii+2]:=x[ii+2]*coe0+x[ii+3]*sie0;
     xecl[ii+3]:=-x[ii+2]*sie0+x[ii+3]*coe0;
     if vel=1 then
       begin
       xecl[ii+4]:=x[ii+4];
       xecl[ii+5]:=x[ii+5]*coe0+x[ii+6]*sie0;
       xecl[ii+6]:=-x[ii+5]*sie0+x[ii+6]*coe0;
       end;
     end;
    x:=xecl;
  end;
end;

type mas728=array[1..728] of double;
var de406:file of mas728;
    buf:mas728;

procedure read406(ecl,vel:integer;t:extended;var x:masc);
type vec=array[1..11] of integer;
     mas15=array[0..15] of extended;
const ratl:extended=1/(81.30056+1);
      nper:vec=( 3, 171, 207, 261, 291, 309, 327, 345, 363, 381, 693);
      pow:vec=( 13,  11,   8,   9,   5,   5,   5,   5,   5,  12,  11);
{               14   12    9   10    6    6    6    6    6   13   12     0     0
                 4    1    2    1    1    1    1    1    1    8    1     0     0}
      ae=149597870.691;
      tmin=625360.5;
      tmax=2816912.5;
      nrc:longint=-1;
      ntp406:boolean=true;

label 1;
var tc,tcp:mas15;
    a,b,c,d,t0,s,s1,s2,TT:extended;
    i,j,k,ii,jj,ik,ist,n1:integer;
    nr:longint;
    xecl:masc;
  procedure cheb(vel:integer; a,b,t:extended;st:integer;var tc,tcp:mas15);
  var i:integer;
      tau,tau2,dlin,rat:extended;
   begin
     dlin:=b-a;
     tau:=2*(t-a)/dlin-1; tau2:=tau*2;
     tc[0]:=1;           tc[1]:=tau;
     for i:=2 to st do tc[i]:=tau2*tc[i-1]-tc[i-2];
     if vel<>0 then
       begin
       rat:=4/dlin;
       tcp[0]:=0;         tcp[1]:=2/dlin;
       for i:=2 to st do tcp[i]:=rat*tc[i-1]+tau2*tcp[i-1]-tcp[i-2];
       end;
   end;
  procedure coor(vel:integer; i,n1:integer;tc,tcp:mas15;var x:masc);
  var k,j,ii,jj:integer;
      s,s1:extended;
   begin
     ii:=6*i-6;
     ist:=pow[i]+1;
     for k:=1 to 3 do
       begin s:=0; ik:=ii+k;
       jj:=nper[i]+(k-1)*ist+n1;
       for j:=0 to pow[i] do s:=s+tc[j]*buf[jj+j];
       x[ik]:=s;
       end;
    if vel<>0 then
    for k:=1 to 3 do
       begin s:=0; ik:=ii+k;
       jj:=nper[i]+(k-1)*ist+n1;
       s:=0;
       for j:=1 to pow[i] do s:=s+tcp[j]*buf[jj+j];
       x[ik+3]:=s;
       end
   end;
 begin
if (t<tmin) or (t>tmax) then
  begin writeln('��� �� �।����� 䮭��!');
  writeln('tmin=',tmin:10:1);
  writeln('tmax=',tmax:10:1);readln;
  halt;
  end;

if ntp406 then
  begin assign(de406,'d:\users\de406\windows.406');
  reset(de406);
  ntp406:=false;
  end;

   nr:=trunc((t-tmin)/64);

if nr<>nrc then
  begin
1:seek(de406,nr);
  read(de406,buf);
  if t<buf[1] then begin nr:=nr-1; goto 1 end;
  if t>buf[2] then begin nr:=nr+1; goto 1 end;
  nrc:=nr;
  end;
a:=buf[1]; b:=buf[2]; n1:=0;
   cheb(vel,a,b,t,11,tc,tcp);
   for i:=2 to 11 do  if(i<>3)and(i<>10) then coor(vel,i,n1,tc,tcp,x);

   c:=a-16; d:=a; n1:=-42;
   repeat c:=c+16; d:=d+16; n1:=n1+42; until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[1],tc,tcp); coor(vel,1,n1,tc,tcp,x);

   c:=a-32; d:=a; n1:=-27;
   repeat c:=c+32; d:=d+32; n1:=n1+27  until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[3],tc,tcp); coor(vel,3,n1,tc,tcp,x);

   c:=a-8; d:=a; n1:=-39;
   repeat c:=c+8; d:=d+8; n1:=n1+39 until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[10],tc,tcp); coor(vel,10,n1,tc,tcp,x);

   for k:=1 to 3 do x[12+k]:=x[12+k]-ratl*x[54+k];
   for k:=1 to 3 do x[54+k]:=x[12+k]+x[54+k];
   if vel<>0 then
     begin
     for k:=4 to 6 do x[12+k]:=x[12+k]-ratl*x[54+k];
     for k:=4 to 6 do x[54+k]:=x[12+k]+x[54+k];
     end;
   for i:=1 to 10 do
    begin
      ii:=6*i-6;
      for k:=1 to 3 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
      if vel<>0 then for k:=4 to 6 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
     end;

   if ecl=1 then
    begin
      for i:=1 to 10 do
       begin
         ii:=6*i-6;
         xecl[ii+1]:=x[ii+1];
         xecl[ii+2]:=x[ii+2]*coe0+x[ii+3]*sie0;
         xecl[ii+3]:=-x[ii+2]*sie0+x[ii+3]*coe0;
         if vel=1 then
           begin
           xecl[ii+4]:=x[ii+4];
           xecl[ii+5]:=x[ii+5]*coe0+x[ii+6]*sie0;
           xecl[ii+6]:=-x[ii+5]*sie0+x[ii+6]*coe0;
           end;
       end;
      x:=xecl;
    end;
 for k:=1 to 6 do x[60+k]:=0;
 end;

type mas1018=array[1..1018] of double;
var de405:file of mas1018;
    buf2:mas1018;

procedure read405(ecl,vel:integer;t:extended;var x:masc);
type vec=array[1..11] of integer;
     mas15=array[0..15] of extended;
const ratl:extended=1/(81.30056+1);
      nper:vec=( 3, 171, 231, 309, 342, 366, 387, 405, 423, 441, 753);{819   899}
      pow:vec=( 13,   9,  12,  10,   7,   6,   5,   5,   5,  12,  10);
{               14   10   13   11    8    7    6    6    6   13   11    10    10
                 4    2    2    1    1    1    1    1    1    8    2     4     4}
      ae=149597870.691;{!!!!!!!!!!!!!}
      tmin=2305424.5;
      tmax=2524624.5;
      nrc:longint=-1;
      ntp405:boolean=true;

label 1;
var tc,tcp:mas15;
    a,b,c,d,t0,s,s1,s2,TT:extended;
    i,j,k,ii,jj,ik,ist,n1:integer;
    nr:longint;
    xecl:masc;
  procedure cheb(vel:integer; a,b,t:extended;st:integer;var tc,tcp:mas15);
  var i:integer;
      tau,tau2,dlin,rat:extended;
   begin
     dlin:=b-a;
     tau:=2*(t-a)/dlin-1; tau2:=tau*2;
     tc[0]:=1;           tc[1]:=tau;
     for i:=2 to st do tc[i]:=tau2*tc[i-1]-tc[i-2];
     if vel<>0 then
       begin
       rat:=4/dlin;
       tcp[0]:=0;         tcp[1]:=2/dlin;
       for i:=2 to st do tcp[i]:=rat*tc[i-1]+tau2*tcp[i-1]-tcp[i-2];
       end;
   end;
  procedure coor(vel:integer; i,n1:integer;tc,tcp:mas15;var x:masc);
  var k,j,ii,jj:integer;
      s,s1:extended;
   begin
     ii:=6*i-6;
     ist:=pow[i]+1;
     for k:=1 to 3 do
       begin
       s:=0;
       ik:=ii+k;
       jj:=nper[i]+(k-1)*ist+n1;
       for j:=0 to pow[i] do s:=s+tc[j]*buf2[jj+j];
       x[ik]:=s;
       end;
    if vel<>0 then
    for k:=1 to 3 do
       begin
       s:=0;
       ik:=ii+k;
       jj:=nper[i]+(k-1)*ist+n1;
       s:=0;
       for j:=1 to pow[i] do s:=s+tcp[j]*buf2[jj+j];
       x[ik+3]:=s;
       end
   end;
 begin
if (t<tmin) or (t>tmax) then
  begin writeln('��� �� �।����� 䮭��!');
  writeln('tmin=',tmin:10:1);
  writeln('tmax=',tmax:10:1);readln;
  halt;
  end;

if ntp405 then
  begin assign(de405,'16002200.405');
  reset(de405);
  ntp405:=false;
  end;

   nr:=trunc((t-tmin)/32);

if nr<>nrc then
  begin
1:seek(de405,nr);
  read(de405,buf2);
  if t<buf2[1] then begin nr:=nr-1; goto 1 end;
  if t>buf2[2] then begin nr:=nr+1; goto 1 end;
  nrc:=nr;
  end;

   a:=buf2[1]; b:=buf2[2]; n1:=0;
   cheb(vel,a,b,t,10,tc,tcp);
   for i:=4 to 9 do coor(vel,i,n1,tc,tcp,x);

   c:=a-8; d:=a; n1:=-42;
   repeat c:=c+8; d:=d+8; n1:=n1+42; until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[1],tc,tcp); coor(vel,1,n1,tc,tcp,x);

   c:=a-16; d:=a; n1:=-30;
   repeat c:=c+16; d:=d+16; n1:=n1+30 until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[2],tc,tcp); coor(vel,2,n1,tc,tcp,x);

   c:=a-16; d:=a; n1:=-39;
   repeat  c:=c+16; d:=d+16; n1:=n1+39  until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[3],tc,tcp); coor(vel,3,n1,tc,tcp,x);

   c:=a-4; d:=a; n1:=-39;
   repeat  c:=c+4; d:=d+4; n1:=n1+39  until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[10],tc,tcp); coor(vel,10,n1,tc,tcp,x);

   c:=a-16; d:=a; n1:=-33;
   repeat  c:=c+16; d:=d+16; n1:=n1+33  until (t>=c)and(t<=d);
   cheb(vel,c,d,t,pow[11],tc,tcp); coor(vel,11,n1,tc,tcp,x);

   for k:=1 to 3 do x[12+k]:=x[12+k]-ratl*x[54+k];
   for k:=1 to 3 do x[54+k]:=x[12+k]+x[54+k];
   if vel<>0 then
     begin
     for k:=4 to 6 do x[12+k]:=x[12+k]-ratl*x[54+k];
     for k:=4 to 6 do x[54+k]:=x[12+k]+x[54+k];
     end;
   for i:=1 to 10 do
    begin
      ii:=6*i-6;
      for k:=1 to 3 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
      if vel<>0 then for k:=4 to 6 do x[ii+k]:=(x[ii+k]-x[60+k])/ae;
     end;

   if ecl=1 then
    begin
      for i:=1 to 10 do
       begin
         ii:=6*i-6;
         xecl[ii+1]:=x[ii+1];
         xecl[ii+2]:=x[ii+2]*coe0+x[ii+3]*sie0;
         xecl[ii+3]:=-x[ii+2]*sie0+x[ii+3]*coe0;
         if vel=1 then
           begin
           xecl[ii+4]:=x[ii+4];
           xecl[ii+5]:=x[ii+5]*coe0+x[ii+6]*sie0;
           xecl[ii+6]:=-x[ii+5]*sie0+x[ii+6]*coe0;
           end;
       end;
      x:=xecl;
    end;
 for k:=1 to 6 do x[60+k]:=0;
 end;

begin
sie0:=sin(e0); coe0:=cos(e0);
end.
