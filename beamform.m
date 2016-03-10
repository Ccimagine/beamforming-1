clear all
M = 8; % number of sensors
spacing = 1; % corresponds to lambda/2
numjam = 5; % number of jammers
numsigs = numjam + 1; % and one signal of interest
a = sqrt(0.5);
deg = pi/180;   % degrees to radians conversions
angled = 10*deg; % desired signal angle-of-arrival (AoA); 0 assumes signal comes from broadside
anglej = [-69 -30 -17 14 34]*deg; % jammer angles converted to radians
angle = [angled anglej]; % vector of all AoAs
sigpow = -10;   % signal power [in dB]
jampow = [30 32 27 30 29]; % jammer powers [in dB]

% Generate steering vectors for each jammer
for p = 1:numjam
    S(:,p) = exp(j*([0:M-1]'-(M-1)/2)*pi*spacing*sin(anglej(p)));
end

% Compute interference covariance matrix, i.e. jammers only
Rint = S*a*diag(10.^(jampow/10))*S';
% Compute total noise + interference covariance matrix
Rnoise = Rint + eye(M);

% Compute steering vector for signal-of-interest (SoI)
s = exp(j*([0:M-1]'-(M-1)/2)*pi*spacing*sin(angled));
s = s/sqrt(s'*s);

% Generate steering vectors
for m = 0:M-1
    steer(m+1,:) = exp(j*(m-(M-1)/2)*pi*spacing*sin(([1:0.5:360]-180)*deg));
end

% Compute and plot frequency response of conventional beamformer
for p = 1:(360/0.5-1)
    gainc(p) = s'*steer(:,p)*steer(:,p)'*s;
end

% Compute optimum Wiener solution and plot its frequency response
w = inv(Rnoise)*s/(s'*inv(Rnoise)*s);
for p = 1:(360/0.5-1)
    gain(p) = w'*steer(:,p)*steer(:,p)'*w;
end

% Plot resulting frequency responses to compare
figure;
plot([1:0.5:360]-180,10*log10(abs(gainc)),'m');
hold;
plot([1:0.5:360]-180,10*log10(abs(gain)));
plot(angled/deg*ones(1,101),-60:40,'gx');
for k = 1:numjam
    plot(anglej(k)/deg*ones(1,101),-60:40,'rx');
end
legend('Conventional Beamformer', 'Optimal Wiener Solution','True SoI DoA','True Jammer DoAs');
axis([-90 90 -60 40]);
xlabel('DoA [Degrees]');
ylabel('Gain [dB]');
%boldify