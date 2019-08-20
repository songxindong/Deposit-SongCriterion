function Mat = HumanCochleogram(waveform, SR, cf)
% PLOT HUMAN COCHLEOGRAM
% Given a waveform and SR, plot ERB gammatonegram for a series of cf
% based on Patterson's ear model (1987/1992) and Slaney's implementation (1993)
%   SR & cf in Hz

%% Determine the ERB (Equivalent Rectangular Bandwidths) for each cf the cf
% in the filterbank

cf =    reshape(cf, 1, [])';
ERB =   24.7*(4.37*cf/1000 +1); 

%% Filterbank Parameters
% Gammatone parameters: T, B, cf needed, 
% T     is a NUMBER indicates the sampling inteval
% B     is a VECTOR of bandwidths related to ERB
% cf    is a VECTOR of CFs of filters

T = 1/SR;
B = 1.019*2*pi*ERB;

% cascade filters implementation parameters
A0 = T *    ones(length(cf),1);
A2 = 0 *    ones(length(cf),1);
B0 = 1 *    ones(length(cf),1);
B1 = -2*cos(2*cf*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;

gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                         (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                          sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
           2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);

%% Determine filtering 
% filterring on waveform
Mat = single( zeros( size(gain,1), length(waveform) ) );
for chan = 1: size(gain,1)
%     chan
	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		   A2(chan)/gain(chan)], ...
				[B0(chan) B1(chan) B2(chan)], waveform);
	y2=filter([A0(chan) A12(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y1);
	y3=filter([A0(chan) A13(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y2);
	y4=filter([A0(chan) A14(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y3);
    Mat(chan, :) = y4;
end
    % Mat is a matrix of (filter# * waveform sample#)

























