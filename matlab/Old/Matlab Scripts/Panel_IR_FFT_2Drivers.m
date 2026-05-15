clear all
close all

% Panel Properties

% Plate Dimensions
Lx = 0.38;
Ly = 0.28;
H = 0.006;

% divide panel into 100x100 grid points
points = 100;

% 0 to Lx
x = linspace(0,Lx,points);  % Define a coordinate grid on the panel
y = linspace(0,Ly,points);

% Plate Material Properties
% Acrylic
% Youngs Modulus
E = 3.2e9;
% Density
p = 1180;
% Poisson's Ratio
v = 0.35;
% Quality Factor
Q = 10;

% Drivers: two points, equal amplitude, in-phase
D_x = [0.25, 0.75];
D_y = [0.50, 0.50];
amp = [1, 1];
phi = [0, 0];        % radians; try phi(2)=pi for anti-phase

% possible modes
m = 15;
n = 15;
counter = 1;
for i = 1:m
    for j = 1:n
    Mode(counter,1) = i;
    Mode(counter,2) = j;
    counter = counter+1;
    end
end
     
% Calculate Resonant Frequencies
Res = zeros(1,m*n);
for i = 1:(counter-1)
    Res(i) = single_clamped_mode_freqs(Lx,Ly,H,E,v,p,Mode(i,1),Mode(i,2)); 
end

% Order Modes by resonant frequencies
[wmn,I] = sort(Res);

for i = 1:length(I)
    Modes(i,1) = Mode(I(i),1);
    Modes(i,2) = Mode(I(i),2);
end


% so ifft is a 1 second IR at 48kHz
f = 0:1:24000;
w = 2*pi*f;
displacement_sensor = zeros(1,length(f));
velocity_sensor = zeros(1,length(f));
average_velocity = zeros(1,length(f));

x_pos = 0;  y_pos = 0;  z_pos = 3.0;
pressure_far_pt = zeros(1, length(f));

pressure_far_pt_mode12 = zeros(1, length(f));
pressure_far_pt_mode21 = zeros(1, length(f));

for i = 1:length(f)
    displacement = zeros(points);

    displacement_mode12 = zeros(points);
    displacement_mode21 = zeros(points);
    % velocity = zeros(100);
    for j = 1:length(Modes)
        % alpha1 = 4*(sin(Modes(j,1)*pi*D_x)*sin(Modes(j,2)*pi*D_y))/Lx/Ly/p/H;
        coupling=0;
        for k = 1:numel(D_x)
            coupling = coupling + amp(k)*exp(1i*phi(k)) * sin(Modes(j,1)*pi*D_x(k)) * sin(Modes(j,2)*pi*D_y(k));
        end
        alpha1 = 4 * coupling / (Lx * Ly * p * H);
        displacement = displacement + alpha1*(sin(Modes(j,1)*pi*x/Lx)'*sin(Modes(j,2)*pi*y/Ly))'/(wmn(j)^2 - w(i)^2 + 1i*w(i)*wmn(j)/Q);

        % add a probe to check if mode(1,2) and mode(2,1) are all 0 in the
        % displacement matrix
        if(Modes(j,1) == 1 && Modes(j,2) == 2)
            displacement_mode12 = displacement_mode12 + alpha1*(sin(Modes(j,1)*pi*x/Lx)'*sin(Modes(j,2)*pi*y/Ly))'/(wmn(j)^2 - w(i)^2 + 1i*w(i)*wmn(j)/Q);
        % % elseif (Modes(j,1) == 2 && Modes(j,2) == 1)
        % %     displacement_mode21 = displacement_mode21 + alpha1*(sin(Modes(j,1)*pi*x/Lx)'*sin(Modes(j,2)*pi*y/Ly))'/(wmn(j)^2 - w(i)^2 + 1i*w(i)*wmn(j)/Q);
        end
    end
    pressure_far_pt_temp = discrete_Rayleigh(displacement, Lx, Ly, f(i),x_pos, y_pos, z_pos);
    pressure_far_pt(i) = pressure_far_pt_temp;

    pressure_far_pt_mode12(i) = discrete_Rayleigh(displacement_mode12, Lx, Ly, f(i),x_pos, y_pos, z_pos);
    % pressure_far_pt_mode21(i) = discrete_Rayleigh(displacement_mode21, Lx, Ly, f(i), x_pos, y_pos, z_pos);
end

figure
loglog(f,abs(pressure_far_pt));grid on;
title(sprintf('Pressure Magnitude  (x=%.1f m, y=%.1f m, z=%.1f m)',x_pos, y_pos, z_pos))
xlabel('Frequency (Hz)')
ylabel('Magnitude (Pa)')
%% 
IRdata_pressure = ifft([pressure_far_pt conj(fliplr(pressure_far_pt(2:end)))]);
IRdata_pressure_truncated = IRdata_pressure(1:10000);

figure
plot(IRdata_pressure);
title(sprintf('Impulse Response — point simulation (z = %.1f m)', z_pos))
xlabel('Index n')
ylabel('Pressure (Pa)')

figure
plot(IRdata_pressure_truncated);
title(sprintf('Impulse Response truncated — point simulation (z = %.1f m)', z_pos))
xlabel('Index n')
ylabel('Pressure (Pa) ')

%% check the magnitude of mode(1,2)
IR_mode12 = ifft([pressure_far_pt_mode12 conj(fliplr(pressure_far_pt_mode12(2:end)))]);
IR_mode12_truncated = IR_mode12(1:10000);
figure
plot(IR_mode12_truncated);
xlabel('Index n')
ylabel('Pressure (Pa)')

IR_mode21 = ifft([pressure_far_pt_mode21 conj(fliplr(pressure_far_pt_mode21(2:end)))]);
IR_mode21_truncated = IR_mode21(1:10000);
figure
plot(IR_mode21_truncated);
xlabel('Index n')
ylabel('Pressure (Pa)')