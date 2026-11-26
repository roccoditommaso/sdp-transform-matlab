%% SDPT_EXAMPLE  Worked example for the SDP time-frequency representation.
%
%   Loads the companion file sdpt_example_data.mat and feeds it to sdpt.m.
%   Run it as a script: it needs nothing else besides sdpt.m and the .mat file.
%
%   The data are a synthetic three-storey response whose fundamental frequency
%   drops from 2.50 to about 1.47 Hz during a transient and partially recovers
%   to 1.90 Hz, with a second mode and measurement noise. Because the true
%   instantaneous frequency is stored in the file, every step can be checked.
%
%   See SDPT_guide.pdf for an illustrated walk-through of these steps.
%
%   Funded by the Italian Civil Protection Department, project DPC-ReLUIS 2024-26
%   - WP6 "Monitoring and satellite data".
%
%   Reference
%     Ditommaso R., Ponzo F.C. (2026). Time-frequency representation of
%     nonlinear structural response through rationally sampled Stockwell
%     analysis. Mechanical Systems and Signal Processing.
%     DOI: https://doi.org/10.1016/j.ymssp.2026.114920

close all      % this is a script: it adds its variables to your workspace
               % and does not clear it

%% 1 -- load the example record ------------------------------------------------
S = load('sdpt_example_data.mat');
fprintf('%s\n\n',S.description);
fprintf('record: %g s at %g Hz, %d channels\n',S.t(end),S.fs,size(S.X,2));

fmin = 0.8;      % lowest frequency of interest [Hz]
fmax = 12;       % highest frequency of interest [Hz]

%% 2 -- default call: map and ridge, everything automatic ----------------------
% With no output argument sdpt plots the record on top and the map with its
% ridge below, on aligned time axes.
sdpt(S.x,fmin,fmax,S.fs);

%% 3 -- read the ridge back and check it against the truth ---------------------
[M1,f1,t1,info1] = sdpt(S.x,fmin,fmax,S.fs,'Plot',false);
fTrue = interp1(S.t,S.fTrue,t1(:));
rmse  = 1000*sqrt(mean((info1.ridge - fTrue).^2));
fprintf('\nsupport: D = %d, %d nodes between %.3f and %.3f Hz\n',info1.D,info1.nodes,f1(1),f1(end));
fprintf('ridge:   RMSE against the true instantaneous frequency = %.0f mHz\n',rmse);
fprintf('         initial %.2f Hz, minimum %.2f Hz, final %.2f Hz\n', ...
        median(info1.ridge(t1<20)),min(info1.ridge),median(info1.ridge(t1>100)));

figure('Color','w','Name','ridge vs truth');
plot(t1,info1.ridge,'LineWidth',1.3); hold on
plot(S.t,S.fTrue,'k--','LineWidth',1.1);
xlabel('Time [s]'); ylabel('Frequency [Hz]'); grid on; axis tight
legend('SDP ridge','true instantaneous frequency','Location','best');
title('Example 2  --  extracted ridge against the truth');

%% 4 -- the property that motivates the method ---------------------------------
% The support is fixed by (fs, D, band): two records of different length share
% the same nodes, so their coefficients can be compared position by position.
% Declare the conductor explicitly whenever records must be commensurable.
D = 120;
[Ma,fa] = sdpt(S.x     ,fmin,fmax,S.fs,'Conductor',D,'Ridge',false,'Plot',false);
[Mb,fb] = sdpt(S.xShort,fmin,fmax,S.fs,'Conductor',D,'Ridge',false,'Plot',false);
nb = size(Mb,2);
fprintf('\nfull record  : %d nodes\n',numel(fa));
fprintf('short record : %d nodes  (%.0f%% of the samples)\n',numel(fb),100*nb/size(Ma,2));
fprintf('identical support: %d\n',isequal(fa,fb));
w = t1(1:nb) > 10 & t1(1:nb) < 70;      % window well inside both records
rel = 100*norm(Ma(:,w)-Mb(:,w),'fro')/norm(Ma(:,w),'fro');
fprintf('coefficient-wise difference over a common window: %.2e %%\n',rel);

%% 5 -- three channels at once -------------------------------------------------
[Mc,fc,tc,infoC] = sdpt(S.X,fmin,fmax,S.fs,'Conductor',D,'TimeStep',2,'Plot',false);
fprintf('\nthree channels: map %s, ridge %s\n',mat2str(size(Mc)),mat2str(size(infoC.ridge)));
figure('Color','w','Name','channels');
plot(tc,infoC.ridge,'LineWidth',1.2); grid on; axis tight
xlabel('Time [s]'); ylabel('Ridge [Hz]'); legend(S.channels{:},'Location','best');
title('Example 3  --  ridge of each floor');

%% 6 -- band energy through the Voronoi weights --------------------------------
% info.weights are the cell widths df/f: they integrate the map over the band
% with the measure the Stockwell transform is built on.
E = sum((abs(M1).^2).*info1.weights,1);
figure('Color','w','Name','band energy');
plot(t1,E,'LineWidth',1.2); grid on; axis tight
xlabel('Time [s]'); ylabel('\Sigma |S|^2 df/f');
title('Example 4  --  band energy from the Voronoi weights');
fprintf('\nsum of the weights = %.4f   (log(fmax/fmin) = %.4f)\n',sum(info1.weights),log(fmax/fmin));
fprintf('\nDone. See the four figures.\n');