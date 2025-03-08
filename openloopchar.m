%eigenvalue finder
%run after FindF16Dynamics

FindF16Dynamics

%20000, 600 

A_long = SS_long_lo.A;
B_long = SS_long_lo.B;
C_long= SS_long_lo.C; 
D_long = SS_long_lo.D;

disp(A_long)
%disp(B_long)
disp(C_long)
%disp(D_long)

A_lat = SS_lat_lo.A;
B_lat = SS_lat_lo.B;
C_lat = SS_lat_lo.C;
D_lat = SS_lat_lo.D;

A_long_red = A_long(1:end-2, 1:end-2);
B_long_red = A_long(:, end); 
C_long_red = C_long(1:end, 1:end);
D_long_red = D_long(:, 1:end-2);



function print_matrix_latex(name, M)
    fprintf('\\[\n');
    fprintf('%s = \\begin{bmatrix}\n', name);
    for i = 1:size(M,1)
        fprintf('%8.4f', M(i,1));
        for j = 2:size(M,2)
            fprintf(' & %8.4f', M(i,j));
        end
        fprintf(' \\\\\n');
    end
    fprintf('\\end{bmatrix}\n\\]\n\n');
end




print_matrix_latex('A', A_long_red);
print_matrix_latex('B', B_long_red);
print_matrix_latex('C', C_long_red);


disp(A_long_red)
disp(B_long_red)
disp(C_long_red)
disp(D_long_red)

A_lat_red = A_lat(1:end-2, 1:end-2);
B_lat_red = A_lat(1:4, end-1:end); 
C_lat_red = C_lat(1:end-2, :);
D_lat_red = D_lat(:, 1:end-2);



eig_long = eig(A_long_red);

disp(eig_long)

eig_lat = eig(A_lat_red);

disp(eig_lat)

%for short-period/Phugoid


eigenvalues = eig_lat;  %switch between eig_long and eig_lat

% Initialize storage
motion_data = [];

% Loop through eigenvalues (step by 2 to handle conjugate pairs)
for i = 1:length(eigenvalues)
    lambda = eigenvalues(i);
    
    % Extract real and imaginary parts
    sigma = real(lambda);
    omega_d = abs(imag(lambda));
    
    % Compute natural frequency
    omega_n = sqrt(sigma^2 + omega_d^2);
    
    % Compute damping ratio
    zeta = -sigma / omega_n;
    
    % Compute period P (if oscillatory)
    if omega_d ~= 0
        P = 2 * pi / omega_d;
    else
        P = NaN;
    end
    
    % Compute time constant tau
    if sigma ~= 0
        tau = -1 / sigma;
    else
        tau = NaN;  % If sigma is 0, tau is undefined
    end
    
    % Compute time to damp to half amplitude T1/2
    if zeta > 0
        T1_2 = log(2) / (zeta * omega_n);  % Only valid for zeta > 0
    else
        T1_2 = NaN;  % If zeta is 0 or less, T1/2 is undefined
    end

    % Store results
    motion_data = [motion_data; omega_n, zeta, P, tau, T1_2, lambda];
end

% Convert results to table
motion_table = array2table(motion_data, 'VariableNames', ...
    {'Omega_n', 'Damping_Zeta', 'Period_P', 'Tau', 'T1_2', 'Eigenvalue'});

% Display results
disp('Longitudinal Motion Characteristics:');
disp(motion_table);


% Define time vector
t = linspace(0, 100, 1000);

wn_short = 1.5459;  zeta_short = 0.4991;
wn_phugoid = 0.0677;  zeta_phugoid = 0.0588;
wn_dutch = 2.9452;  zeta_dutch = 0.10467;

% now we compute time responses so after definitiosn
short_period = exp(-zeta_short * wn_short * t) .* cos(wn_short * sqrt(1 - zeta_short^2) * t);
phugoid = exp(-zeta_phugoid * wn_phugoid * t) .* cos(wn_phugoid * sqrt(1 - zeta_phugoid^2) * t);
dutch_roll = exp(-zeta_dutch * wn_dutch * t) .* cos(wn_dutch * sqrt(1 - zeta_dutch^2) * t);
tau_roll = 0.45155;
tau_spiral = 99.189;

aperiodic_roll = exp(-t/tau_roll);
spiral_mode = exp(-t/tau_spiral);

% Plot responses
figure;

% Longitudinal motions
subplot(3,2,1);
plot(t, short_period, 'b', 'LineWidth', 1.5);
title('Short Period (Longitudinal)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,2,2);
plot(t, phugoid, 'b', 'LineWidth', 1.5);
title('Phugoid (Longitudinal)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Lateral motions
subplot(3,2,3);
plot(t, dutch_roll, 'r', 'LineWidth', 1.5);
title('Dutch Roll (Lateral)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,2,4);
plot(t, aperiodic_roll, 'r', 'LineWidth', 1.5);
title('Aperiodic Roll (Lateral)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,2,5);
plot(t, spiral_mode, 'r', 'LineWidth', 1.5);
title('Spiral Mode (Lateral)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Improve layout
sgtitle('F-16 Eigenmotion Time Responses');



%Computed Natural Frequencies (omega_n):
    %longitudinal
    %0.0677 corresponds to Phugoid with damping ratio 0.0588 (very
    %underdamped) with period 93.1239
    %1.5459 corresponds to Short Period with damping ratio 0.4991
    %(underdamped, but more damped) with period 4.6918

    %lateral
    %2.9452 corresponds to Roll Mode (Fast damping).
    %0.0101 corresponds to Spiral Mode (Very slow response).
    %0.0000 corresponds to Neutral Spiral Stability.
    %1.0000 corresponds to Dutch Roll with critical damping 1 with Inf Period??




%Computed damping ratios
    %0.0588
    %0.4991


    %lateral 
    %NaN
    %0.1047
    %1.0000
    %1.0000


%Computed periods

    %93.1239
    %4.6918
















