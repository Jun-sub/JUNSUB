function Z = Impedance_input()

    % Input parameters
    f_min = input('Enter the minimum frequency (Hz): ');
    f_max = input('Enter the maximum frequency (Hz): ');
    points = input('Enter the number of points: ');
    R = input('Enter the resistance (ohms): ');
    C = input('Enter the capacitance (farads): ');
    A = input('Enter the amplitude factor: ');

    % Frequency vector
    f = logspace(f_min, f_max, points);
    w = 2 * pi * f;

    % Warburg impedance
    Z_W = A .* (1 - 1i) ./ sqrt(w);

    % Impedance components
    Z_RW = R + Z_W;
    Z_C = 1 ./ (1i * w * C);

    % Overall impedance
    Z = (Z_RW .* Z_C) ./ (Z_RW + Z_C);

    % Plot
    figure;
    plot(real(Z), -imag(Z));
    xlim([0, 3]);
    xlabel('Re(Z) / Ohm');
    ylabel('-Im(Z) / Ohm');
    title('Impedance');
    grid on;
end