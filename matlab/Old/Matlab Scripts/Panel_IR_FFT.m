% MH 03/13/23
% This program computes the velocity at a point on a panel due to an
% impulse being played at a different point on the panel. Note that there
% is a small artifact at the end of the IR that can be chopped off. 

% This program takes a long time to run. If you load the included .mat
% file, you can just run the last section
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

% Normalized Driver Location
D_x = 0.14;
D_y = 0.5;

% Normalized Sensor Location
S_x = 0.5;
S_y = 0.5;

%% compute resonant modes and sort from lowest to highest frequency
% possible modes
m = 30;
n = 30;
counter = 1;
for i = 1:m
    for j = 1:n
    Mode(counter,1) = i;
    Mode(counter,2) = j;
    counter = counter+1;
    end
end
     
% Calculate Resonant Frequencies
Res = zeros(1,900);
for i = 1:(counter-1)
    Res(i) = single_clamped_mode_freqs(Lx,Ly,H,E,v,p,Mode(i,1),Mode(i,2)); 
end

% Order Modes by resonant frequencies
[wmn,I] = sort(Res);

for i = 1:length(I)
    Modes(i,1) = Mode(I(i),1);
    Modes(i,2) = Mode(I(i),2);
end

%% Simulate Time Response
% so ifft is a 1 second IR at 48kHz
f = 0:1:24000;
%f = logspace(log10(20),log10(20000),500);
w = 2*pi*f;
displacement_sensor = zeros(1,length(f));
velocity_sensor = zeros(1,length(f));
average_velocity = zeros(1,length(f));

x_pos = 0;  y_pos = 0;  z_pos = 1.0;
pressure_far_pt = zeros(1, length(f));

for i = 1:length(f)
    displacement = zeros(100);
    velocity = zeros(100);
    for j = 1:length(Modes)
        alpha1 = 4*(sin(Modes(j,1)*pi*D_x)*sin(Modes(j,2)*pi*D_y))/Lx/Ly/p/H;
        displacement = displacement + alpha1*(sin(Modes(j,1)*pi*x/Lx)'*sin(Modes(j,2)*pi*y/Ly))'/(wmn(j)^2 - w(i)^2 + 1i*w(i)*wmn(j)/Q);
        displacement_sensor(i) = displacement_sensor(i) + alpha1*sin(Modes(j,1)*pi*S_x)*sin(Modes(j,2)*pi*S_y)/(wmn(j)^2 - w(i)^2 + 1i*w(i)*wmn(j)/Q);
        %velocity = 1i*w(i)*displacement;   
    end
    velocity_sensor(i) = 1i*w(i)*displacement_sensor(i);
    average_velocity(i) = mean(mean(abs(1i*w(i)*displacement)));

    % compute the sound pressure using rayleigh intergral function
    % out : pressure_far_pt_temp: total sound pressure at point(x,y,z)
    % v_map = 1i * w(i) * displacement;
    pressure_far_pt_temp = discrete_Rayleigh(displacement, Lx, Ly, f(i),x_pos, y_pos, z_pos);
    pressure_far_pt(i) = pressure_far_pt_temp;

%     figure(1)
%     surf(x,y,real(displacement))
%     axis equal
%     colorbar
%     caxis([-1 1])
%     title(['f = ' num2str(f(i)) ' Hz'])
%     view(0,90)
%     pause(0.001)
end

figure
loglog(f,abs(average_velocity));grid on;
title('Panel Average Normal Velocity')
xlabel('Frequency (Hz)')
ylabel('Magnitude (m·s^{-1})')

figure
loglog(f,abs(pressure_far_pt));grid on;
title(sprintf('Pressure Magnitude  (x=%.1f m, y=%.1f m, z=%.1f m)',x_pos, y_pos, z_pos))
xlabel('Frequency (Hz)')
ylabel('Magnitude (Pa)')
%% compute IR

IRdata = ifft([velocity_sensor conj(fliplr(velocity_sensor(2:end)))]);

% truncate IR 
IRdata_truncated = IRdata(1:10000);

figure
plot(IRdata_truncated);
title('Impulse Response — Panel Velocity at Sensor (S_x = 0.5, S_y = 0.5)')
xlabel('Index n')
ylabel('Velocity')

% test IR on speech data
% [sig,FS] = audioread('Voice_Test.wav');
% out = conv(IRdata_truncated,sig);
% soundsc(out,FS)

%% compute ir use pressure

IRdata_pressure = ifft([pressure_far_pt conj(fliplr(pressure_far_pt(2:end)))]);

IRdata_pressure_truncated = IRdata_pressure(1:10000);

figure
plot(IRdata_pressure_truncated);
title(sprintf('Impulse Response — point simulation (z = %.1f m)', z_pos))
xlabel('Index n')
ylabel('Pressure (Pa)')

% test IR on audio data
% [sig,FS] = audioread('Voice_Test.wav');
% out = conv(IRdata_truncated,sig);
% soundsc(out,FS)

%% save the ir simulated

save('panel_IR_simulated_48k.mat', ...
     'IRdata', 'IRdata_truncated', 'IRdata_pressure', 'IRdata_pressure_truncated',...
     'velocity_sensor', 'pressure_far_pt', 'f', ...
     'Lx','Ly','H','D_x','D_y','S_x','S_y', ...
     'x_pos','y_pos','z_pos');
