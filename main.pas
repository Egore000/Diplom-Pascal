uses SysUtils, 
      readfond, 
      TwoBody, 
      service;

const 
      // Порядок резонанса (u:v)
      u = 1;
      v = 2;

      // Запись в файлы
      WRITE_ORBIT = false;
      WRITE_SECOND_PLUS = false;
      WRITE_SECOND_MINUS = false;

      // Пути к файлам и директориям
      // TARGER_FOLDER = 'Без светового давления';
      TARGER_FOLDER = 'Со световым давлением';

      PATH_CLASSIFICATION = '..\Данные\' + TARGER_FOLDER +'\Классификация.csv'; // Путь к файлу с классификацией
      PATH_DATA = '..\Данные\' + TARGER_FOLDER + '\'; // Путь к папке с исходными данными
      PATH_ORBITAL = '..\Данные\Выход\Орбитальные\'; // Путь к папке с данными об орбитальных резонансах
      PATH_SECOND_PLUS = '..\Данные\Выход\Вторичные\плюс\'; // Путь к папке с данными о вторичных резонансах (+)
      PATH_SECOND_MINUS = '..\Данные\Выход\Вторичные\минус\'; // Путь к папке с данными о вторичных резонансах (-)

var coords, velocities: mas; // Массивы координат и скоростей
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

    // Учёт влияния Солнца (при res = 2)
    lmd_s := 0;
    if (res = 2) then
    begin
      fond405(jd, xm_, xs_, vm_, vs_);
      CoordsToElements(xs_, vs_, Gmu, b, ec, i_b, OmegaS, ws, M_s);
      lmd_s := OmegaS + ws + M_s;
    end;  

    Reduce(u * (M + Omega + w) - v * theta + znak * lmd_s, angles[1]);
    Reduce(u * (M + w) + v * (Omega - theta) + znak * lmd_s, angles[2]);
    Reduce(u * M + v * (Omega + w - theta) + znak * lmd_s, angles[3]);
    Reduce(angles[1] - v * Omega + znak * lmd_s, angles[4]);
    Reduce(angles[3] + v * Omega - 2 * v * w + znak * lmd_s, angles[5]);

    n := sqrt(mu/(a*sqr(a)));
    d_OmegaJ2 := -1.5*J2 * n * sqr(r0/a) * cos(i) / sqr(1 - sqr(ecc));
    d_wJ2 := 0.75*J2 * n * sqr(r0/a) * (5*sqr(cos(i)) - 1)/sqr(1 - sqr(ecc));

    d_Omega_L := -3/16 * n_L * mL * (a/a_L)*sqr(a/a_L) * (2 + 3*sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_L))) * cos(i);
    d_Omega_S := -3/16 * n_S * mS * (a/a_S)*sqr(a/a_S) * (2 + 3*sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_S))) * cos(i);

    d_w_L := 3/16 * n_L * mL * (a/a_L)*sqr(a/a_L) * (4 - 5*sqr(sin(i)) + sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_L)));
    d_w_S := 3/16 * n_S * mS * (a/a_S)*sqr(a/a_S) * (4 - 5*sqr(sin(i)) + sqr(ecc))/sqrt(1 - sqr(ecc)) * (2 - 3*sqr(sin(i_S)));

    d_Omega := d_OmegaJ2 + d_Omega_L + d_Omega_S;
    d_w := d_wJ2 + d_w_L + d_w_S;

    freq[1] := u * (n + d_Omega + d_w) - v * d_theta;
    freq[2] := u * (n + d_w) + v * (d_Omega - d_theta);
    freq[3] := u * n + v * (d_Omega + d_w - d_theta);
    freq[4] := freq[1] - v * d_Omega;
    freq[5] := freq[3] + v * d_Omega - 2 * v * d_w;

    // Вычисление частот вторичного резонанса (при res = 2)
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
  inc_count, dec_count, perehod_count, class_: integer;

begin
  // writeln('[FILE]', #9, num);

  // Цикл по компонентам резонанса
  for res := res_start to res_end do
  begin  
    class_ := 0;
    zero_counter := 0;
    count := 0;
    perehod_count := 0;

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
      (phi[res, i]/2 > phi[res, i+1]) then inc(dec_count);
      // ((phi[res, i] < 10 * toRad) and (phi[res, i+1] > 350 * toRad)) then inc(dec_count);


      // Подсчёт возрастающих точек
      if (phi[res, i] < phi[res, i+1]) or 
      (phi[res, i]/2 < phi[res, i+1]) then inc(inc_count);
      // ((phi[res, i] > 350 * toRad) and (phi[res, i+1] < 10 * toRad)) then inc(inc_count);         
    
      // Подсчёт переходов частоты через 0
      if (dot_phi[res, i] * dot_phi[res, i+1]) < 0 then inc(perehod_count);
    end;

    // writeln('[COUNT]    ', count);
    // writeln('[INCREASE]   ', inc_count);
    // writeln('[DECREASE]   ', dec_count);

    // writeln('[ZERO FREQUENCE TRANSITION]   ', perehod_count);
    if (perehod_count > count * coef) then class_ := 2;

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



begin {Main}
  assign(outdata, PATH_CLASSIFICATION);
  rewrite(outdata);

  WriteHeader(outdata, res_start, res_end);

  for folder := start_folder to finish_folder do
    { Цикл по файлам в папке folder }
    for number := start to finish do
    begin
      if (number < 10) then file_num := '000' + inttostr(number);
      if (number >= 10) and (number < 100) then file_num := '00' + inttostr(number);
      if (number >= 100) and (number < 1000) then file_num := '0' + inttostr(number);
      if (number >= 1000) then file_num := inttostr(number);

      assign(data, PATH_DATA + inttostr(folder) + '\EPH_' + file_num + '.DAT');
      reset(data);

      if WRITE_ORBIT then
        Create_File(orbit_res, PATH_ORBITAL + file_num + '.dat');

      if WRITE_SECOND_PLUS then
        Create_File(second_plus, PATH_SECOND_PLUS + file_num + '.dat');

      if WRITE_SECOND_MINUS then
        Create_File(second_minus, PATH_SECOND_MINUS + file_num + '.dat');

      for num := res_start to res_end do 
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

        Resonance(1, 0, year, month, day, M, Omega, w, ecc, i, a, angles, freq);
        Resonance(2, -1, year, month, day, M, Omega, w, ecc, i, a, angles2, freq2);
        Resonance(2, 1, year, month, day, M, Omega, w, ecc, i, a, angles3, freq3);
        
        t[idx] := time;
        time_idx := trunc(time / (86400 * 365 * col_step)) + 1;
        for num := res_start to res_end do
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
        
        if WRITE_ORBIT then WriteToFile(orbit_res, time, angles, freq);
        if WRITE_SECOND_PLUS then WriteToFile(second_minus, time, angles2, freq2);
        if WRITE_SECOND_PLUS then WriteToFile(second_plus, time, angles3, freq3);

        inc(idx);
      end;
      
      // OutNET(net);
      Classification(net, t, phi, dot_phi, number, classes);
      Classification(net2, t, phi2, dot_phi2, number, classes2);
      Classification(net3, t, phi3, dot_phi3, number, classes3);

      WriteClassification(outdata, folder, number, classes, classes2, classes3);

      if WRITE_SECOND_PLUS then close(second_plus);
      if WRITE_SECOND_MINUS then close(second_minus);
      if WRITE_ORBIT then close(orbit_res);
      close(data);
    end;
  close(outdata);
end.