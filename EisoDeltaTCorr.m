% Correlation Coefficient With Measurement Error

clear;

%% Data

total = readtable('total1021.csv','CommentStyle','#')

x = horzcat(log10(total.energy),log10(total.deltaT));
sigmaerror = horzcat(log10(total.energy+total.energyErr)-log10(total.energy),log10(total.deltaT+total.deltaTErr)-log10(total.deltaT));


%% Sampling
% MCMC Parameters
WinBUGS = 'C:/Program Files/WinBUGS14'; % WinBUGS installation Path
nchains = 2; % How Many Chains?
nburnin = 5e4; % How Many Burn-in Samples?
nsamples = 5e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Constants
[n,~] = size(x);
lambdaerror = 1./sigmaerror.^2; % Precision of Measurements

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n,'lambdaerror',lambdaerror);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.r = 0;
    S.mu = zeros(1,2);
    S.lambda = ones(1,2);
    init0(i) = S;
end

% Use WinBUGS to Sample
tic
[samples, stats] = matbugs(datastruct, ...
fullfile(pwd, 'Correlation.txt'), ...
'init', init0, ...
'nChains', nchains, ...
'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
'monitorParams', {'r','mu','sigma'}, ...
'Bugdir', WinBUGS);
toc

%% Analysis
figure();clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
eps = .02; binsc = -1+eps/2:eps:1-eps/2; binse = -1:eps:1;
subplot(121);hold on;
ph = plot(x(:,1),x(:,2),'ko');
set(ph,'markersize',4,'markerfacecolor','k');
for i = 1:n
    ph = plot([x(i,1)-sigmaerror(i,1) x(i,1)+sigmaerror(i,1)],ones(1,2)*x(i,2),'k-');
    ph = plot(ones(1,2)*x(i,1),[x(i,2)-sigmaerror(i,2) x(i,2)+sigmaerror(i,2)],'k-');
end;
%axis([min(x(:,1))*0.98 max(x(:,1))*1.02 min(x(:,2))*0.99 max(x(:,2))*1.02]);
set(gca,'box','on','fontsize',14);
xlabel('Log(Eiso)','fontsize',16);
ylabel('Log(Duration)','fontsize',16);
subplot(122);hold on;
count = histc(reshape(samples.r,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
set(gca,'box','on','fontsize',14,'xtick',[-1:.5:1],'ytick',[]);
tmp = corrcoef(x);
ph = plot(ones(1,2)*tmp(1,2),[0 max(get(gca,'ylim'))],'k--');
xlabel('Correlation','fontsize',16);
ylabel('Posterior Density','fontsize',16);



