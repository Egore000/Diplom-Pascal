uses SysUtils, 
      Classifier,
      ResonanceUnit,
      readfond, 
      TwoBody,
      service;

const 
      // Запись в файлы
      WRITE_ORBIT = false;
      WRITE_SECOND_PLUS = false;
      WRITE_SECOND_MINUS = false;

      // Пути к файлам и директориям
      // TARGER_FOLDER = 'Без светового давления';
      TARGER_FOLDER = 'Со световым давлением';

      PATH_DATA = '..\Исходные данные\' + TARGER_FOLDER + '\'; // Путь к папке с исходными данными
      PATH_CLASSIFICATION = '..\Выходные данные\' + TARGER_FOLDER + '\Классификация1.csv'; // Путь к файлу с классификацией
      PATH_ORBITAL = '..\Выходные данные\' + TARGER_FOLDER + '\Орбитальные\'; // Путь к папке с данными об орбитальных резонансах
      PATH_SECOND_PLUS = '..\Выходные данные\' + TARGER_FOLDER + '\Вторичные\плюс\'; // Путь к папке с данными о вторичных резонансах (+)
      PATH_SECOND_MINUS = '..\Выходные данные\' + TARGER_FOLDER + '\Вторичные\минус\'; // Путь к папке с данными о вторичных резонансах (-)

var coords, velocities: mas; // Массивы координат и скоростей
    angles, angles2, angles3: arr; // Массивы резонансных углов Ф
    freq, freq2, freq3: arr; // Массивы резонансных частот Ф'

    net, net2, net3: NETWORK; // Сетки для разных наборов данных
    classes, classes2, classes3: CLS; // Массивы с классификацией резонансов

    a, e, i, Omega, w, M, megno, mean_megno: extended; 
    tm, time, day: extended;
    year, month, num, number, x: integer;
    idx, time_idx, angle_idx, angle2_idx, angle3_idx: integer; // Индексы

    data, outdata, orbit_res, second_plus, second_minus: text; // Файлы
                                                // data - файл с исходными данными
                                                // outdata - файл для записи элементов
                                                // orbit_res - выходной файл с орбитальными резонансами
                                                // second_plus - выходной файл с вторичными резонансами (+)
                                                // second_minus - выходной файл с вторичными резонансами (-)
    ss: string[7]; // Служебная строка

    phi, phi2, phi3: angle_data; // Массивы с резонансными углами
    dot_phi, dot_phi2, dot_phi3: angle_data; // Массивы с резонансными частотами
    t: time_data; // Массив с моментами времени

    folder: integer; // Папка с исходными файлами
    file_num: string;


begin {Main}
  assign(outdata, PATH_CLASSIFICATION);

  if FileExists(PATH_CLASSIFICATION) then
    append(outdata)
  else
    rewrite(outdata);

  {Заполнение заголовка в файле классификации}
  WriteHeader(outdata, res_start, res_end);

  {Цикл по папкам}
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

      {Связь с файлами, в случае, если осуществляется запись}
      if WRITE_ORBIT then
        Create_File(orbit_res, PATH_ORBITAL + inttostr(folder) + '\' + file_num + '.dat');

      if WRITE_SECOND_PLUS then
        Create_File(second_plus, PATH_SECOND_PLUS + inttostr(folder) + '\' + file_num + '.dat');

      if WRITE_SECOND_MINUS then
        Create_File(second_minus, PATH_SECOND_MINUS + inttostr(folder) + '\' + file_num + '.dat');

      {Заполнение массивов нулями}
      FillZero(net, net2, net3, 
              t,
              phi, phi2, phi3,
              dot_phi, dot_phi2, dot_phi3);
        
      idx := 0;

      while not eof(data) do
      begin
        {Считывание данных из файла}
        readln(data, tm, time, ss, year, month, day);
        readln(data, x, coords[1], coords[2], coords[3], megno);
        readln(data, velocities[1], velocities[2], velocities[3], mean_megno);

        {Расчёт элментов орбиты}
        CoordsToElements(coords, velocities, mu, a, e, i, Omega, w, M);

        {Вычисление аргументов резонансов}
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
        
        {Запись в файлы}
        if WRITE_ORBIT then WriteToFile(orbit_res, time, angles, freq);
        if WRITE_SECOND_PLUS then WriteToFile(second_minus, time, angles2, freq2);
        if WRITE_SECOND_PLUS then WriteToFile(second_plus, time, angles3, freq3);

        inc(idx);
      end;
      
      // OutNET(net);
      writeln('[FILE]', #9, number);

      {Классификация резонансов}
      Classification(net, t, phi, dot_phi, number, classes);
      Classification(net2, t, phi2, dot_phi2, number, classes2);
      Classification(net3, t, phi3, dot_phi3, number, classes3);

      {Запись классификации в файл}
      WriteClassification(outdata, folder, number, classes, classes2, classes3);

      {Закрытие файлов, если они были открыты на запись}
      if WRITE_SECOND_PLUS then close(second_plus);
      if WRITE_SECOND_MINUS then close(second_minus);
      if WRITE_ORBIT then close(orbit_res);
      close(data);
    end;
  close(outdata);
end.