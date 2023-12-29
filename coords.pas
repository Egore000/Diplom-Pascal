uses TwoBody, SysUtils;

const Omega = 240;
    w = 0;
    M0 = 0;

    mu = 3.986004418e+5;
    pi = 3.1415926;
    toRad = pi/180;
    toDeg = 180/pi;
    ecc = 1e-3;

type 
    mas = array[1..3] of extended;

var 
    X, V: mas;
    count, file_number: integer;
    a, i, da, di, n: extended;
    a1, e1, i1, Omega1, w1, M1: extended;
    data, elements: text;

begin
    assign(data, '../coords/coords - 1.dat');
    assign(elements, '../coords/elements.dat');
    rewrite(elements);
    rewrite(data);
    a := 26500;
    da := 0.5;
    i := 0;
    di := 2;
    count := 0;
    n := 0;
    file_number := 1;

    while (a <= 26600) do
    begin
        // writeln(count);
        n := sqrt(mu/(a * sqr(a)));
        while (i <= 180) do
        begin
            TwoPoints(0, n, M0, w, Omega, i*toRad, a, X, V);
            inc(count);
            writeln(elements, file_number, ' ', count, a, i);
            writeln(data, '1.0000');
            writeln(data, '5.0000');
            writeln(data, X[1], X[2], X[3]);
            writeln(data, V[1], V[2], V[3]);
            // writeln(data, a, ecc, i);
            if (count = 9000) then
            begin
                count := 0;
                close(data);
                inc(file_number);
                assign(data, '../coords/coords - ' + inttostr(file_number) + '.dat');
                rewrite(data);
            end;
            i := i + di;
        end;
        i := 0;
        a := a + da;
    end;
    CoordsToElements(X, V, mu, a1, e1, i1, Omega1, w1, M1);
    writeln(data, a1, e1, i1*toDeg);

    writeln('finished');
    close(elements);
    close(data);
end.