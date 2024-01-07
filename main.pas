uses SysUtils, 
      readfond, 
      TwoBody, 
      service;

const a_e = 149597870.691; //[km]
      pi = 3.1415926;
      mu = 3.986004418e+5; // Гравитационная постоянная для спутника
      Gmu = 1.32712442099e+11; // Гравитационная постоянная для Солнца

      znak = -1; // Знак у lambda_s при вторичных возмущениях
      toRad = pi/180;
      toDeg = 180/pi;

      rows = 12;
      cols = 30;
      row_step = 360 / rows;
      col_step = 30 / cols;

type 
    arr = array[1..5] of extended;
    mas = array[1..3] of extended;
    NETWORK = array[1..5, 1..rows + 1, 1..cols + 1] of integer;
    angle_data = array[1..5, 1..2000] of extended;
    time_data = array[1..2000] of extended;

var coords, velocities: mas;
    angles, freq: arr; 
    net: NETWORK;
    
    jd, a, e, i, Omega, w, M, megno, mean_megno: extended; 
    tm, time, day: extended;
    year, month, num, x, row, col, class_: integer;
    idx, time_idx, angle_idx, freq_idx: integer;

    data, orbit_res, second_res: text;
    // ss: string[11];
    ss: string[7];

    phi, dot_phi: angle_data;
    t: time_data;

procedure OutNET(net: NETWORK);
// Вывод матриц
var 
  num, row, col: integer;
begin
  for num := 1 to 5 do
  begin
    writeln('num = ', num);
    for row := 1 to rows do
    begin
      for col := 1 to cols do
      begin
        write(net[num, row, col], #9);
      end;
      writeln;
    end;
    writeln;
  end;
end;


procedure perehod(HH: extended; 
                  xx: mas; 
                  var lambda, phi: extended);
// Переход во вращающуюся СК 
var A: array[1..3,1..3] of extended;
    i, j: integer;
    s, r, argum: extended;
    yy: mas;
begin

    A[1,1] := cos(HH);     A[1,2] := sin(HH);    A[1,3] := 0;
    A[2,1] := -sin(HH);    A[2,2] := cos(HH);    A[2,3] := 0;
    A[3,1] := 0;           A[3,2] := 0;          A[3,3] := 1;

    for i:=1 to 3 do
    begin
      s := 0;
      for j:=1 to 3 do
        s := s + A[i,j]*xx[j];
      yy[i] := s;
    end;
    r := sqrt(yy[1]*yy[1] + yy[2]*yy[2] + yy[3]*yy[3]);

    lambda := arctg(yy[2], yy[1]);

    argum := yy[3]/r;
    phi := arcsin(argum);
  end; //perehod

procedure fond405(jd:extended; var xm_,xs_,vm_,vs_:mas);  var x_planet:masc;
begin
    read405(0, 1, jd, x_planet);    
    // процедура выдает координаты и скорости Луны и Солнца в геоцентрической 
    // экваториальной — км и км/с

    //переход в геоцентрическую экваториальную — в км

    // координаты Луны и Солнца
    xm_[1] := (x_planet[55] - x_planet[13]) * a_e;  //
    xm_[2] := (x_planet[56] - x_planet[14]) * a_e;  //
    xm_[3] := (x_planet[57] - x_planet[15]) * a_e;  //   {xle=xl-xe}
        
    xs_[1] := -x_planet[13] * a_e;  //
    xs_[2] := -x_planet[14] * a_e; //
    xs_[3]:= - x_planet[15] * a_e; //    {xs=-xe }

    //скорости Луны и Солнца
    vm_[1] := (x_planet[58] - x_planet[16]) * a_e/86400;
    vm_[2] := (x_planet[59] - x_planet[17]) * a_e/86400;
    vm_[3] := (x_planet[60] - x_planet[18])*a_e/86400;     {vle=vl-ve}

    vs_[1] := -x_planet[16] * a_e/86400;
    vs_[2] := -x_planet[17] * a_e/86400;
    vs_[3] := -x_planet[18] * a_e/86400;     {vs=-ve}

end;

procedure Resonance(res, znak, year, month: integer;
                    day: extended;
                    M, Omega, w, ecc, i, a: extended;
                    var angles, freq: arr);
// Процедура для вычисления орбитальных резонансов
// res - порядок резонанса (1 - без учёта вторичных возмущений)
//                         (2 - с учётом вторичных возмущений)
// znak - знак lambda_s в формулах для вторичных возмущений
// angles - массив значений критического аргумента
// freq - массив частот резонанса
const
      mL = 1/81.3005690699;
      mS = 332946.048166;
      J2 = 1.0826359e-3;
      r0 = 6363.6726;
      i_L = 23.45 * toRad;
      i_S = 23.45 * toRad;
      a_L = 384748;
      a_S = 149597868;
      n_L = 2 * pi/(27.32166 * 86400);
      n_S = 2 * pi/(365.25 * 86400);
      d_theta = 7.292115e-5;

var theta, n, d_OmegaJ2, d_wJ2, d_Omega_L, d_Omega_S, d_w_L, d_w_S, d_Omega, d_w: extended;
    xm_, xs_, vm_, vs_: mas;
    b, ec, i_b, OmegaS, ws, M_s, lmd_s: extended;

begin
    jd := date_jd(year, month, day);
    theta := sid2000(jd);

    // Учёт влияния Солнца при res = 2
    lmd_s := 0;
    if (res = 2) then
    begin
      fond405(jd, xm_, xs_, vm_, vs_);
      CoordsToElements(xs_, vs_, Gmu, b, ec, i_b, OmegaS, ws, M_s);
      lmd_s := OmegaS + ws + M_s;
    end;  

    Reduce(M + Omega + w - theta + znak * lmd_s, angles[1]);
    Reduce(M + Omega + w - theta + znak * lmd_s, angles[2]);
    Reduce(M + Omega + w - theta + znak * lmd_s, angles[3]);
    Reduce(M + w - theta + znak * lmd_s, angles[4]);
    Reduce(M + 2*Omega - w - theta + znak * lmd_s, angles[5]);

    n := sqrt(mu/(a*sqr(a)));
    d_OmegaJ2 := -1.5*J2 * n * sqr(r0/a) * cos(i) / sqr(1 - sqr(ecc));
    d_wJ2 := 0.75*J2 * n * sqr(r0/a) * (5*sqr(cos(i)) - 1)/sqr(1 - sqr(ecc));

    d_Omega_L := -3/16 * n_L * mL * (a/a_L)*sqr(a/a_L) * (2 + 3*sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_L))) * cos(i);
    d_Omega_S := -3/16 * n_S * mS * (a/a_S)*sqr(a/a_S) * (2 + 3*sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_S))) * cos(i);

    d_w_L := 3/16 * n_L * mL * (a/a_L)*sqr(a/a_L) * (4 - 5*sqr(sin(i)) + sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_L)));
    d_w_S := 3/16 * n_S * mS * (a/a_S)*sqr(a/a_S) * (4 - 5*sqr(sin(i)) + sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_S)));

    d_Omega := d_OmegaJ2 + d_Omega_L + d_Omega_S;
    d_w := d_wJ2 + d_w_L + d_w_S;

    freq[1] := n + d_Omega + d_w - d_theta;
    freq[2] := n + d_Omega + d_w - d_theta;
    freq[3] := n + d_Omega + d_w - d_theta;
    freq[4] := n + d_w - d_theta;
    freq[5] := n + 2*d_Omega - d_w - d_theta;

    // Вычисление частот резонанса при res = 2
    if (res = 2) then
    begin
      freq[1] := freq[1] + znak * ( d_Omega_S + d_w_S + n_S );
      freq[2] := freq[2] + znak * ( d_Omega_S + d_w_S + n_S );
      freq[3] := freq[3] + znak * ( d_Omega_S + d_w_S + n_S );
      freq[4] := freq[4] + znak * ( d_Omega_S + d_w_S + n_S );
      freq[5] := freq[5] + znak * ( d_Omega_S + d_w_S + n_S );
    end;
end; {Resonance}


procedure Classification(net: NETWORK;
                        t: time_data;
                        phi, dot_phi: angle_data;
                        num: integer;
                        var class_: integer);
// Классификация резонанса
// Параметры:
// net - сетка графика
// num - номер исследуемого объекта
// class_ - выходное значение класса резонанса

// 1 - циркуляция
// 2 - либрация
// 3 - смешанный тип
var
  res, i, j: integer;
  zero_cols_counter, zero_rows_counter, zero_counter, count: integer;
  inc_count, dec_count: integer;
begin
  writeln('[FILE]', #9, num);

  // Цикл по компонентам резонанса
  for res := 1 to 5 do
  begin  
    zero_counter := 0;
    count := 0;

    for i := 1 to rows do
    begin
      zero_rows_counter := 0;
    
      for j := 1 to cols do
      begin
        count := count + net[res, i, j];

        if (net[res, i, j] = 0) then 
        begin
          inc(zero_counter); 
          inc(zero_rows_counter);
        end;
      end;

      if (zero_rows_counter = cols) then class_ := 2;
      if (zero_rows_counter > 1) then class_ := 3;
    end;

    for j:= 1 to cols do
    begin
      zero_cols_counter := 0;

      for i := 1 to rows do
      begin
        if (net[res, i, j] = 0) then inc(zero_cols_counter);

      end;  

      if (zero_cols_counter > 1) then class_ := 3;
    end;

    if (zero_counter = 0) then class_ := 1;

    inc_count := 0;
    dec_count := 0;
    for i := 1 to count-1 do
    begin
      if (phi[res, i] > phi[res, i+1]) or 
      ((phi[res, i] < 10 * toRad) and (phi[res, i+1] > 350 * toRad)) then
        inc(dec_count); 
      
      if (phi[res, i] < phi[res, i+1]) or 
      ((phi[res, i] > 350 * toRad) and (phi[res, i+1] < 10 * toRad)) then
        inc(inc_count);         
    end;

    if (inc_count = count) then 
    begin
      writeln(inc_count, '=',count);
      class_ := 1;
    end;
    if (dec_count = count) then
    begin
      writeln(dec_count, '=',count);
      class_ := 1;
    end;

    writeln('[RESONANCE]', #9, res, #9, '[CLASS]', #9, class_);
  end;
end;

procedure WriteToFile(var f: text;
                      time: extended;
                      angles, freq: arr);
begin
  writeln(f, time/(86400 * 365), ' ',
              angles[1] * toDeg, ' ',
              angles[2] * toDeg, ' ',
              angles[3] * toDeg, ' ',
              angles[4] * toDeg, ' ',
              angles[5] * toDeg, ' ',
              freq[1], ' ', 
              freq[2], ' ',
              freq[3], ' ',
              freq[4], ' ',
              freq[5]);
end;


procedure Check_all(folder: integer);
var 
  data, outdata: text;
  num: integer;
  file_num: string;

begin
  // assign(outdata, '..\Данные\Классификация.csv');
  // rewrite(outdata);
  
  for num := 1 to 1003 do
  begin
    if (num < 10) then file_num := '000' + inttostr(num);
    if (num >= 10) and (num < 100) then file_num := '00' + inttostr(num);
    if (num >= 100) and (num < 1000) then file_num := '0' + inttostr(num);
    if (num >= 1000) then file_num := inttostr(num);

    assign(data, '..\Данные\' + inttostr(folder) + '\EPH_' + file_num + '.DAT');
    reset(data);

    close(data);
  end;

  // close(outdata);
end;

begin {Main}
    // assign(data, 'C:\Users\egorp\Desktop\диплом\файлы\ЧМ ИСЗ (для ПК) 28.04.20 Lobbie III\EPH_0001.DAT');
    assign(data, '..\Данные\1\EPH_0001.DAT');
    reset(data);
    // assign(orbit_res, 'Орбитальные резонансы.dat');
    // rewrite(orbit_res);
    // assign(second_res, 'Вторичные резонансы.dat');
    // rewrite(second_res);

    // writeln(orbit_res, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');
    // writeln(second_res, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');

    for num := 1 to 5 do 
    begin
      for row := 1 to rows do
        for col := 1 to cols do
          net[num, row, col] := 0;
      
      for row := 1 to 2000 do
      begin
        t[row] := 0;
        phi[num, row] := 0;
        dot_phi[num, row] := 0;
      end;
    end;
    
    idx := 0;

    while not eof(data) do
    begin
        readln(data, tm, time, ss, year, month, day);
        readln(data, x, coords[1], coords[2], coords[3], megno);
        readln(data, velocities[1], velocities[2], velocities[3], mean_megno);

        CoordsToElements(coords, velocities, mu, a, e, i, Omega, w, M);

        Resonance(1, znak, year, month, day, M, Omega, w, ecc, i, a, angles, freq);
        
        t[idx] := time;
        time_idx := trunc(time / (86400 * 365 * col_step)) + 1;
        for num := 1 to 5 do
        begin
          angle_idx := trunc(angles[num] * toDeg / row_step) + 1;
          inc(net[num, angle_idx, time_idx]);

          phi[num, idx] := angles[num];
          dot_phi[num, idx] := freq[num];
        end;
        // WriteToFile(orbit_res, time, angles, freq);

        // Resonance(2, znak, year, month, day, M, Omega, w, ecc, i, a, angles, freq);
        // WriteToFile(second_res, time, angles, freq);
    end;

    // OutNET(net);

    Classification(net, t, phi, dot_phi, 1, class_);

    Check_all(1);
    // close(second_res);
    // close(orbit_res);
    close(data);
end.