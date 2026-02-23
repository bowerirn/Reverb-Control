clear;
% Load ir data
ir_data = load('measured_ir_data.mat');
%% 
panel2mic = ir_data.measurementData(1,:);
panel_Reflection = ir_data.measurementData(2,:);

%% plot the impulse response
figure;
plot(panel2mic.ImpulseResponse.Time, panel2mic.ImpulseResponse.Amplitude, 'b', 'DisplayName', 'Panel to microphone');
hold on;
plot(panel_Reflection.ImpulseResponse.Time, panel_Reflection.ImpulseResponse.Amplitude, 'r', 'DisplayName', 'Source hit Panel then Reflect');
xlabel('Time');
ylabel('Amplitude');
legend show;
grid on;
