function Spec_dB_Out = HumanExcitationPattern(Freq_In, Spec_Amp_In, Freq_Out)
% Human Excitation Pattern caclulates the steady excitation pattern of a
% sound witha a stationary spectrum in th auditory periphery
% For the sound:
%   Freq_In,        the discrete INPUT frequency series on linear scale
%   Spec_Amp_In,    the Amp(rms) at each INPUT frequency sampling bin,
% For the excitation pattern
%   Freq_Out,       the discrete OUTPUT frequency series,
%                   on any (linear/log/ERB) scale
%   Spec_dB_Out,    the Amp(dB relative) for each OUTPUT frequency sampling bin
%                   assume Amp = 1(rms) is X dB SPL
%                   the number in SpecdBOut is dB number relative to X in each bin
%
%   Technical NOTES:
%       (1) The model is a simplifed model in Moore & Glasberg 1983,
%           from which an auditory filter can be described as 
%               W(g) =              (1 + pg) * exp(-pg),        where
%               g =     (f - cf) / cf
%       (2) A full rounded-exponential filter is
%               W(g) =  (1 - r) *   (1 + pg) * exp(-pg) + r,    where
%               g =     (f - cf) / cf
%       (3) Is was said in Moore & Glasberg 1983: 
%           " For calculating excitation patterns fromfilter shapes, it is 
%           convenient to use a simple one-parameter approximation to the 
%           auditory filter shape, by dropping the parameter r from ...",
%       (4) and:
%           " For this simplifed filter shape, the value of ERB/f is 
%           exactly 4/p" 
%       (2) ERBs are calculated based on the equation in Glasberg & Moore
%           1990, as ERB = 24.7 * (4.37 * f/1000 +1)
%           as also appeared in Moore BCJ 2012 pp76
%
%   Code reorganized on 20190124


ERB_Out =       24.7 * (4.37 * Freq_Out/1000 + 1);  
ERB_p_Out =     4 * Freq_Out./ERB_Out; 
ERB_r_Out =     0*ERB_p_Out;            % in case to be used later
Spec_dB_Out =   0 * Freq_Out;

% figure
% loglog(FreqOut, ERB_Out);

for i = 1:length(Freq_Out)
    % for each output auditory filter,
        
    % it recieves frequencies as
    g_In =  abs( (Freq_In - Freq_Out(i)) / Freq_Out(i) );
    
    % each with the weight
    W_In =  (1 - ERB_r_Out(i)) * (1 + ERB_p_Out(i)*g_In) .* exp(-ERB_p_Out(i)*g_In) ...
            + ERB_r_Out(i);
    
    % so the output power at this filter together is
    Spec_dB_Out(i) = 10*log10( sum(Spec_Amp_In.^2 .*W_In) );
end;

% maybe for marmosets later
% ERB_r_raw       =[	0       0.03735 0.00848 0.00193 0.00095];
% ERB_p_Out   = interp1(ERB_freq_Raw,	ERB_p_raw, FreqOut,'spline');
% ERB_r_Out   = interp1(ERB_freq_Raw,	ERB_r_raw, FreqOut,'spline');
