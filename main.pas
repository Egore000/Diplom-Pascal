uses SysUtils, 
      readfond, 
      TwoBody, 
      service;

const a_e = 149597870.691; //[km]
      pi = 3.1415926;
      mu = 3.986004418e+5; // Гравитационная постоянная для спутника
      Gmu = 1.32712442099e+11; // Гравитационная постоянная для Солнца

      znak = -1; // Знак у lambda_s при вторичных возмущениях
      toRad = pi/180; // Перевод в радианы
      toDeg = 180/pi; // Перевод в градусы

      rows = 8; // Количество строк в сетке
      cols = 80; // Количество столбцов в сетке
      row_step = 360 / rows; // Шаг по строкам
      col_step = 100 / cols; // Шаг по столбцам

      // Запись в файл
      ORBIT = false;
      SECOND_PLUS = false;
      SECOND_MINUS = false;

type 
    CLS = array[1..5] of integer;
    arr = array[1..5] of extended;
    mas = array[1..3] of extended;
    NETWORK = array[1..5, 1..rows + 1, 1..cols + 1] of integer;
    angle_data = array[1..5, 1..2000] of extended;
    time_data = array[1..2000] of extended;

var coords, velocities: mas; // Массив координат и скоростей
    angles, angles2, angles3: arr; // Массивы резонансных углов Ф
    freq, freq2, freq3: arr; // Массивы резонансных частот Ф'

    net, net2, net3: NETWORK; // Сетки для разных наборов данных
    classes, classes2, classes3: CLS; // Массивы с классификацией резонансов

    jd, a, e, i, Omega, w, M, megno, mean_megno: extended; 
    tm, time, day: extended;
    year, month, num, number, x, row, col: integer;
    idx, time_idx, angle_idx, angle2_idx, angle3_idx: integer; // Индексы

    data, outdata, orbit_res, second_plus, second_minus: text; // Файлы
                                                // data - файл с исходными данными
                                                // outdata - файл для записи элементов
                                                // orbit_res - выходной файл с орбитальными резонансами
                                                // second_plus - выходной файл с вторичными резонансами (+)
                                                // second_minus - выходной файл с вторичными резонансами (-)
    // ss: string[11];
    ss: string[7]; // Служебная строка

    phi, phi2, phi3: angle_data; // Массивы с резонансными углами
    dot_phi, dot_phi2, dot_phi3: angle_data; // Массивы с резонансными частотами
    t: time_data; // Массив с моментами времени

    folder: integer; // Папка с исходными файлами
    file_num: string;

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
                    M, Omega, w, ecc, i, a: extended; // Элементы орбиты
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
                        var classes: CLS);
// Классификация резонанса
// Параметры:
// net - сетка графика
// t - массив времени
// phi, dot_phi - массивы углов и частот
// num - номер исследуемого объекта
// classes - выходной массив классов резонанса

// 0 - циркуляция
// 1 - либрация
// 2 - смешанный тип
var
  res, i, j: integer;
  zero_cols_counter, zero_rows_counter, zero_counter, count: integer;
  inc_count, dec_count, class_: integer;

begin
  // writeln('[FILE]', #9, num);

  // Цикл по компонентам резонанса
  for res := 1 to 5 do
  begin  
    class_ := 0;
    zero_counter := 0;
    count := 0;

    // Цикл по строчкам сетки
    for i := 1 to rows do
    begin
      zero_rows_counter := 0; {Счётчик нулевых ячеек в строке}
    
      // Цикл по столбцам
      for j := 1 to cols do
      begin
        count := count + net[res, i, j];

        if (net[res, i, j] = 0) then 
        begin
          inc(zero_counter); 
          inc(zero_rows_counter);
        end;
      end;

      if (zero_rows_counter <> 0) then class_ := 2;
      if (zero_rows_counter = cols) then class_ := 1;
    end;

    for j:= 1 to cols do
    begin
      zero_cols_counter := 0;

      for i := 1 to rows do
        if (net[res, i, j] = 0) then inc(zero_cols_counter);

      if (zero_cols_counter > 1) then class_ := 2;
    end;

    if (zero_counter = 0) then class_ := 0;

    inc_count := 0;
    dec_count := 0;
    for i := 1 to count-1 do
    begin
      // Подсчёт убывющих точек
      if (phi[res, i] > phi[res, i+1]) or 
      ((phi[res, i] < 10 * toRad) and (phi[res, i+1] > 350 * toRad)) then inc(dec_count);

      // Подсчёт возрастающих точек
      if (phi[res, i] < phi[res, i+1]) or 
      ((phi[res, i] > 350 * toRad) and (phi[res, i+1] < 10 * toRad)) then inc(inc_count);         
    end;

    if (inc_count = count) then 
    begin
      writeln(inc_count, '=',count);
      class_ := 0;
    end;
    if (dec_count = count) then
    begin
      writeln(dec_count, '=',count);
      class_ := 0;
    end;

    // writeln('[RESONANCE]', #9, res, #9, '[CLASS]', #9, class_);
    classes[res] := class_;
  end;
end;

procedure WriteToFile(var f: text;
                      time: extended;
                      angles, freq: arr);
// Запись данных о резонансах в файл f
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

begin {Main}
  assign(outdata, '..\Данные\Выход\Классификация.csv');
  rewrite(outdata);
  folder := 1;

  writeln(outdata, 'folder,file,F1,F2,F3,F4,F5,dF1(+),dF2(+),dF3(+),dF4(+),dF5(+),dF1(-),dF2(-),dF3(-),dF4(-),dF5(-)');

  { Цикл по файлам в папке folder }
  for number := 1 to 9000 do
  begin
    if (number < 10) then file_num := '000' + inttostr(number);
    if (number >= 10) and (number < 100) then file_num := '00' + inttostr(number);
    if (number >= 100) and (number < 1000) then file_num := '0' + inttostr(number);
    if (number >= 1000) then file_num := inttostr(number);

    assign(data, '..\Данные\' + inttostr(folder) + '\EPH_' + file_num + '.DAT');
    reset(data);

    if ORBIT then
    begin
      assign(orbit_res, '..\Данные\Выход\Орбитальные\' + file_num + '.dat');
      rewrite(orbit_res);
      writeln(orbit_res, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');
    end;

    if SECOND_PLUS then
    begin
      assign(second_plus, '..\Данные\Выход\Вторичные\плюс\' + file_num + '.dat');
      rewrite(second_plus);
      writeln(second_plus, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');
    end;

    if SECOND_MINUS then
    begin
      assign(second_minus, '..\Данные\Выход\Вторичные\минус\' + file_num + '.dat');
      rewrite(second_minus);
      writeln(second_minus, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');
    end;

    for num := 1 to 5 do 
    begin
      for row := 1 to rows do
        for col := 1 to cols do
        begin
          net[num, row, col] := 0;
          net2[num, row, col] := 0;
          net3[num, row, col] := 0;
        end;
      
      for row := 1 to 2000 do
      begin
        t[row] := 0;

        phi[num, row] := 0;
        phi2[num, row] := 0;
        phi3[num, row] := 0;

        dot_phi[num, row] := 0;
        dot_phi2[num, row] := 0;
        dot_phi3[num, row] := 0;
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
      Resonance(2, -1, year, month, day, M, Omega, w, ecc, i, a, angles2, freq2);
      Resonance(2, 1, year, month, day, M, Omega, w, ecc, i, a, angles3, freq3);
      
      t[idx] := time;
      time_idx := trunc(time / (86400 * 365 * col_step)) + 1;
      for num := 1 to 5 do
      begin
        {Заполнение массивов для орбитального резонанса}
        angle_idx := trunc(angles[num] * toDeg / row_step) + 1;
        inc(net[num, angle_idx, time_idx]);
        phi[num, idx] := angles[num];
        dot_phi[num, idx] := freq[num];
        
        {Заполнение массивов для вторичных резонансов (знак -)}
        angle2_idx := trunc(angles2[num] * toDeg / row_step) + 1;
        inc(net2[num, angle2_idx, time_idx]);
        phi2[num, idx] := angles2[num];
        dot_phi2[num, idx] := freq2[num];

        {Заполнение массивов для вторичных резонансов (знак +)}
        angle3_idx := trunc(angles3[num] * toDeg / row_step) + 1;
        inc(net3[num, angle3_idx, time_idx]);
        phi3[num, idx] := angles3[num];
        dot_phi3[num, idx] := freq3[num];
      end;
      
      if ORBIT then WriteToFile(orbit_res, time, angles, freq);
      if SECOND_PLUS then WriteToFile(second_minus, time, angles2, freq2);
      if SECOND_PLUS then WriteToFile(second_plus, time, angles3, freq3);

      inc(idx);
    end;
    
    // OutNET(net);
    Classification(net, t, phi, dot_phi, number, classes);
    Classification(net2, t, phi2, dot_phi2, number, classes2);
    Classification(net3, t, phi3, dot_phi3, number, classes3);

    writeln(outdata, folder,',',
                      number,',',
                      classes[1],',',
                      classes[2],',',
                      classes[3],',',
                      classes[4],',',
                      classes[5],',',
                      classes2[1],',',
                      classes2[2],',',
                      classes2[3],',',
                      classes2[4],',',
                      classes2[5],',',
                      classes3[1],',',
                      classes3[2],',',
                      classes3[3],',',
                      classes3[4],',',
                      classes3[5]);

    if SECOND_PLUS then close(second_plus);
    if SECOND_MINUS then close(second_minus);
    if ORBIT then close(orbit_res);
    close(data);
  end;
  close(outdata);
end.