function [MI,d] = mmi_psych(S1R1,S1R2,S2R1,S2R2)
% function [MI,d] = mmi_psych(S1R1,S1R2,S2R1,S2R2)
% 
% computes mutual information (MI) and sensitivity (d') for behavioral data in a two-stimulus (or two-category), two-choice task
% 
% INPUT
% hit     number of correct S1 trials (i.e., stimulus 1, response 1)
% miss    number of false S1 trials (i.e., stimulus 1, response 2)
% fa      number of false S2 trials (i.e., stimulus 2, response 1)
% cr      number of correct S2 trials (i.e., stimulus 2, response 2)
% 
% OUTPUT
% MI      behavioral categorization performance expressed as MI
% d       signal detection theory index d' (sensitivity)
% 
% MI and d' are correlated:
% p_hit = 0.6;p_fa = 0.3;
% for i = 1:100
%   ntrials1 = round(normrnd(100,40));
%   ntrials2 = round(normrnd(100,40));
%   nhits = binornd(ntrials1,p_hit);
%   nfa   = binornd(ntrials2,p_fa);
%   [MI(i,1),d(i,1)] = mmi_psych(nhits,ntrials1-nhits,nfa,ntrials2-nfa);
% end
% figure,scatter(d,MI)
% corr(d,MI,'rows','complete')
% 
% September 2021, Maik C. Stüttgen, University Medical Center Mainz, Germany
%% preparations
ntrials = S1R1+S1R2+S2R1+S2R2;
Pr1 = (S1R1+S2R1)/ntrials;
Pr2 = (S1R2+S2R2)/ntrials;
Ps1 = (S1R1+S1R2)/ntrials;
Ps2 = (S2R1+S2R2)/ntrials;
Ps1r1 = S1R1/(S1R1+S2R1);
Ps2r1 = S2R1/(S1R1+S2R1);
Ps1r2 = S1R2/(S1R1+S2R1);
Ps2r2 = S2R2/(S1R2+S2R2);

%% information
% I = sum over r,s: P(r)*P(s|r)*log2(P(s|r)/P(s))
I(1) = Pr1*Ps1r1*log2(Ps1r1/Ps1);
I(2) = Pr1*Ps2r1*log2(Ps2r1/Ps2);
I(3) = Pr2*Ps1r2*log2(Ps1r2/Ps1);
I(4) = Pr2*Ps2r2*log2(Ps2r2/Ps2);
MI = sum(I);

d = norminv(S1R1/(S1R1+S1R2)) - norminv(S2R1/(S2R1+S2R2));    % norminv(hit_rate_S1) - norminv (false_alarm_rate_S2)
end