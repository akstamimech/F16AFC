sys = ss(SS_lo.A, SS_lo.B, SS_lo.C, SS_lo.D);

G = tf(sys); 

G_19_2 = G(19,2); 


figure
t = 0:0.01:10;  % Time vector from 0 to 0.5 seconds
[y, t] = step(G_19_2, t);  % Get the step response within this time range
%[y, t] = step(G(18,2),t);
plot(t, -y);  % Plot the negative of the step response
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Negative Step Response (First 10 Seconds)');