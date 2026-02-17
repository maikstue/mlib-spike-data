function [vs,pvs] = mvecstrength(cycle_tspx,T)
% function [vs,pvs] = mvecstrength(cycle_tspx,T)
% 
% This function computes the so-called vector strength and performs a significance test on it.
% The measure was introduced first by Goldberg and Brown (1969, J Neurophysiol).
% 
% Vector strength (vs) is computed as follows:
% vs = sqrt ( sum(sin(phi))^2 + sum(cos(phi))^2 ) / n
% where n is the total number of spikes and
% phi = 2*pi*(t/T)
% where t is the time between cycle onset and an evoked spike, and T is the period of the stimulus frequency
% (e.g. the duration of a single sine wave). phi thus represents the phase of the spike in radians.
% 
% The measure has been used many times before, e.g. by Gao and Wehr (Neuron, 2015).
% Note that the definition of vector strength given in Khatri et al. (2004, J Neurophysiol) is incorrect.
% Further information can be found in Grün & Rotter (2010): Analysis of parallel spike trains, Chapter 4,
% written by Ashida, Wagner, and Carr: Processing of phase-locked spikes and periodic signals.
% Also see  https://www.cns.nyu.edu/~bijan/courses/sm10/Lectures/semple/NS2_07_Audition2.pdf.
% 
% INPUT
% cycle_tspx  a vector of spike timestamps relative to cycle onset in ms
% T           period length of the sine wave in ms
% 
% OUTPUT
% vs    vector strength, a number between 0 and 1
% pvs   the p-value from the Rayleigh test
% 
% HISTORY
% 2025 June     renamed input 'tspx' to 'cycle_tspx'
% 
% one could add optional bootstrapping the p-value for small sample sizes, see book by Grün & Rotter
% 
% October 2021
% Maik C. Stüttgen, University Medical Center Mainz
%% sanity checks
if any(cycle_tspx>T)
  disp('Some spike times exceed the period length, I will eliminate these.')
  cycle_tspx(cycle_tspx>T) = [];
end
if numel(cycle_tspx)<50
  disp('Less than 50 spikes provided, significance test may not be robust.')
end
%% the works
% compute phi
phi = nan(numel(cycle_tspx),1);
for i = 1:numel(phi)
  phi(i,1) = 2*pi*cycle_tspx(i)/T;
end

% compute vector strength
vs = sqrt(sum(cos(phi))^2 + sum(sin(phi))^2) / numel(phi);

% compute significance: exp(-n*vs^2)
pvs = exp(-numel(phi)*vs^2);
end