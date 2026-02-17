function [MI,cMI1995,cMI1996] = mmi(s1,s2)
% function [MI,cMI1995,cMI1996] = mmi(s1,s2)
% 
% computes Shannon's mutual information uncorrected (MI) and corrected (cMI) between two spike count distributions from a single neuron.
% s1 and s2 must be vectors of integers (spikes per trial).
% 
% The formula is (according to Arabzadeh et al., J Neurosci 2007):
% MI = sum over r,s: P(r)*P(s|r)*log2(P(s|r)/P(s)
% 
% This is the formula according to Reyes-Puerta et al., PLoS Comp Biol 2014 (also see Panzeri & Treves, Network 1996):
% MI = H(R) - H(R|S) = sum over s: P(s) * sum over r: P(r|s)*log2(P(r|s)/P(r))
% H(R) is the response entropy (overall variability of the responses) and H(R|S) is the noise entropy (response variability to a given stimulus).
% 
% This is the formula in Treves and Panzeri, Neural Computation 1995:
% MI = sum over s: P(s) * sum over r: P(r|s)*log2(P(r|s)) - P(r)*log2(P(r))
% 
% Let's assume the following spike count distributions for a single neuron and two stimuli:
% s1 = [0,0,1,2,2,3,5]
% s2 = [0,0,1,2,3,4,4]
% P(r=3) = 2/14
% P(s=s1) = P(s=s2) = 7/14;
% P(s=s1|r=3) = 1/2
% 
% An apparent problem is that log2(0) yields negative infinity; this is unavoidable because either P(r|s) or P(s|r) will be 0 for
% some value of r in almost any data set! But I think one can ignore this - responses which do not occur are simply not used.
% Another issue is that the maximum information is determined by the probabilities of the two stimuli; if they are equal, max(I)=1.
% If they are unequal, max(I) is smaller: max(I) = H(S) (Borst & Theunissen, Nat Neurosci 1999).
% H(S) = -sum over s: P(s)*log2(Ps)
% For two equiprobable stimuli Ps1=0.5,Ps2=1-Ps1, H(S) = -sum([Ps1*log2(Ps1),Ps2*log2(Ps2)]) = 1.
% Ps1=1/6 and Ps2=1-Ps1, H(S) = 0.650
% Ps1=1/4 and Ps2=1-Ps1, H(S) = 0.811 
% Ps1=1/3 and Ps2=1-Ps1, H(S) = 0.918
% Ps1=6/10 and Ps2=1-Ps1, H(S) = 0.971
% Ps1=18/38 and Ps2=1-Ps1, H(S) = 0.998
% Thus, as long as we stay reasonable close to 50%, we end up in a regime where max(I) is close 1.
% We could temper this even further if we subsample and bootstrap when we stray too far from stimulus equiprobability.
% 
% The first output (MI) was confirmed through "information.m" from the Information Breakdown Toolbox.
% The second and third outputs match up with what they write in the 1996 and 2007 papers, they are however not identical to the output of the ITB.
% This is most probably because the effective number of responses Rs is estimated, not taken naively (see Panzeri 2007, p. 5). Since cMI1996 is higher
% than cMI1995, I advise to use the latter.
% 
% September 2021, Maik C. Stüttgen, University Medical Center Mainz, Germany
%% preparations
if ~any(int32(s1)==s1) || ~any(int32(s2)==s2),error('Please exclusively provide vectors with integers.'),end
s1 = s1(:);s2 = s2(:);  s1s2 = [s1;s2]; % ensure column vectors

minval = min(s1s2);maxval = max(s1s2);vals = minval:maxval;bins = minval-0.5:maxval+0.5;
ns1 = histcounts(s1,bins);
ns2 = histcounts(s2,bins);
ns1s2 = histcounts(s1s2,bins);
%% collect all probabilities
% get P(s=s1) and P(s=s2)
Ps1 = numel(s1)/numel(s1s2);Ps2 = 1-Ps1;
disp(['Stimulus probabilities are ',num2str(Ps1,'%0.3f'),' and ',num2str(Ps2,'%0.3f'),'.'])
disp(['Maximum information = H(S) = ',num2str(-sum([Ps1*log2(Ps1),Ps2*log2(Ps2)])),'.'])

% get P(r=x); sum(Pr)=1
Pr = (ns1s2/numel(s1s2))';

disp(['For S1, we have ',num2str(numel(unique(s1))),' different responses.'])
disp(['For S2, we have ',num2str(numel(unique(s2))),' different responses.'])
disp(['Across both S1 and S2, there are ',num2str(numel(unique(s1s2))),' different responses.'])

% get P(s=y|r=x) for all x; Ps1r+Ps2r=1 for each r
for i = 1:numel(vals)
  Ps1r(i,1) = sum(s1==vals(i)) / sum(s1s2==vals(i));
  Ps2r(i,1) = sum(s2==vals(i)) / sum(s1s2==vals(i));  % should be 1-Ps1r
end

% I = sum over r,s: P(r)*P(s|r)*log2(P(s|r)/P(s)
MI = nansum([Pr.*Ps1r.*log2(Ps1r/Ps1),Pr.*Ps2r.*log2(Ps2r/Ps2)],'all');
% this yields the same result
% for i = 1:numel(vals)
%   preMI(i,1) = Pr(i)*Ps1r(i)*log2(Ps1r(i)/Ps1);
%   preMI(i,2) = Pr(i)*Ps2r(i)*log2(Ps2r(i)/Ps2);
% end
% MI2 = nansum(preMI(:));

% formula in Panzeri et al., J Neurophysiol 2007:
% I = H(R) - H(R|S) = sum over r,s: P(s)*P(r|s)*log2(P(r|s)/P(r))
%% bias correction
% I computed three different types of bias correction according to formulas presented in Panzeri's 1995, 1996, and 2007 papers.
% The 1996 and 2007 factors yielded equal results.
% The 1995 factor is always higher but is an overestimation when Rs<R. See below for details
% 
% Treves and Panzeri, Neural Computation 1995: C = (1/(2N*ln(2))) * (S-1) * (R-1)
C1995 = (1/(2*numel(s1s2)*log(2))) * (2-1) * (sum(ns1s2~=0)-1);
% The first part is getting smaller with increasing N: for N=20, it is 0.036, for N=100, it is 0.007 (i.e., a fifth).
% The second part increases with S and R (S is the number of stimuli, R is the number of response bins).
% I am quite sure that I implemented this correctly.
% They however go on to say that this correction factor is too large, and a better one is:
% C = (1/(2*numel(s1s2)*log(2))) * (S-1)*(R'-1), where R' is an "effective" number of response bins.
% This however seems to apply for data which has been "regularized" in some way, and I would not know how to arrive at this for my data.
% 
% In the 1996 paper, they say that the C1995 measure is appropriate only when each response bin has a non-zero probability of being occupied 
% for every stimulus - this is rarely the case for us.
% So one should perhaps use the following expression from Panzeri and Treves, Network 1996:
% C = 1/(2N*log(2)) * ((sum over s: Rs) - R - (S-1))
C1996 = (1/(2*numel(s1s2)*log(2))) * (sum(ns1~=0) + sum(ns2~=0) - sum(ns1s2~=0)-(2-1));
% For Rs, we can simply take the number of relevant responses for the neuron, or one can estimate it with a Bayesian procedure (see the paper).
% If we do the former (the "naive count" - Ince et al., 2014), we are likely to underestimate Rs because of limited sampling.
% Not sure whether this is relevant for us anyway, because they mention it is important after regularization.
% Strange - the second term will always be negative!!!
% I am unsure whether I correctly grasped the meaning of the curly brackets in this context.
% In any case, C1996 gives weird results for small (but not very large) sample sizes.
% 
% The formula again looks different in Panzeri et al., J Neurophysiol 2007, equation 4:
% C = 1/(2N*log(2)) * (sum over s:(Rs-1)) - (R-1))
% C2007 = (1/(2*numel(s1s2)*log(2))) * ((sum(ns1~=0)-1)+(sum(ns2~=0)-1)-(sum(ns1s2~=0)-1));
% C2007 gives similarly weird results for small sample sizes.

disp(['Bias correction factor 1995: ', num2str(C1995)])
disp(['Bias correction factor 1996: ', num2str(C1996)])
% disp(['Bias correction factor 2007: ', num2str(C2007)])

% best bias correction according to Panzeri et al. 2007, equation 7
% H(R)      = -sum over r: P(r)*log2(P(r))
% H(R|S)    = -sum over r,s: P(s)*P(r|s)*log2(Pr|s)
% I(S;R)    = H(R)-H(R|S) = sum over r,s: P(s)*P(r|s)*log2(P(r|s)/P(r))
% Hind(R|S) = -sum over r,s: P(s)*Pind(r|s)*log2(Pind(r|s)
% Ish(S;R) = H(R) - Hind(R|S) + Hsh(R|S) - H(R|S)
% Is not straightforward because probability distributions have to be estimated.
%% get corrected MI
cMI1995 = MI-C1995;
cMI1996 = MI-C1996;
%% exploring things
% Borst & Theunissen, Nat Neurosci 1998:
% Information about stimulus x: I(R,Sx) = sum over r: P(r|sx) * log2(P(r|sx)/P(r))
% Prs1 = (ns1/sum(ns1))';
% Prs2 = (ns2/sum(ns2))';
% MI3 = nansum([Ps1.*Prs1.*log2(Prs1./Pr),Ps2.*Prs2.*log2(Prs2./Pr)],'all');
% % I(i,1) = nansum(Prs1.*log2(Prs1./Pr));
% % I(i,2) = nansum(Prs2.*log2(Prs2./Pr));
% 
% % MI = I(R,S) = H(S) - H(S|R) = H(R) - H(R|S)
% HR  = -nansum(Pr.*log2(Pr));                                        % validated with entropyPanzeri.m
% HRS = -sum(nansum([Ps1*Prs1.*log2(Prs1) , Ps2*Prs2.*log2(Prs2)]));  % validated with entropyPanzeri.m
% HS  = -(Ps1*log2(Ps1) + Ps2*log2(Ps2));
% HSR = -nansum(Pr.*nansum([Ps1r.*log2(Ps1r),Ps2r.*log2(Ps2r)],2));
% MI4 = HS-HSR;
% 
% thres = 1e-10;
% if MI-MI3>thres || MI-MI4>thres,error(['Calculations do not match, MI-MI3=',num2str(MI-MI3),' and MI-MI4=',num2str(MI-MI4),'.']),end 
end