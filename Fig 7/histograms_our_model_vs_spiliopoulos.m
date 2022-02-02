close all

%% load data
load([cd() '\saved_escapes\bt_et.mat'])
load([cd() '\saved_escapes\bt_et_s.mat'])

bt = [bt_et.bt];
bts = [bt_et_s.bt];

tauA = bt_et(1).et;
tauB = bt_et(2).et;
tauC = bt_et(3).et;
tauD = bt_et(4).et;

tausA = bt_et_s(1).et;
tausB = bt_et_s(2).et;
tausC = bt_et_s(3).et;
tausD = bt_et_s(4).et;


%% plot 8 histograms
close all

ds = {'DisplayStyle','stairs'};
ps = {'PaperSize',[6.8 5.5]};
fs = {'FontSize',10};
fig_size = [250 200];
lwd = {'LineWidth',1.3};
ls = {'LineStyle','--'};
lbl = {'time','number of samples'};

pos1 = {'Position',[50 500 fig_size]};     %adaptation 1
pos2 = {'Position',[350 500 fig_size]};     %adaptation 2
pos3 = {'Position',[650 500 fig_size]};     %adaptation 3
pos4 = {'Position',[950 500 fig_size]};     %adaptation 4
pos1s = {'Position',[50 100 fig_size]};    %spilopoulos 1
pos2s = {'Position',[350 100 fig_size]};    %spiliopoulos 2
pos3s = {'Position',[650 100 fig_size]};    %spiliopoulos 3
pos4s = {'Position',[950 100 fig_size]};    %spiliopoulos 4

axs = {'XLim',[0 500],'YLim',[0 500]};
axss = {'XLim',[0 15],'YLim',[0 500]};

nbins = 20;
% firstly four with adaptation
fg1 = figure(pos1{:},ps{:});hold on
histogram(tauA(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tauA(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axs{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg1,[cd() '\histograms\hist1'],'-dpdf')


fg2 = figure(pos2{:},ps{:});hold on
histogram(tauB(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tauB(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axs{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg2,[cd() '\histograms\hist2'],'-dpdf')


fg3 = figure(pos3{:},ps{:});hold on
histogram(tauC(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tauC(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axs{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg3,[cd() '\histograms\hist3'],'-dpdf')


fg4 = figure(pos4{:},ps{:});hold on
histogram(tauD(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tauD(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axs{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg4,[cd() '\histograms\hist4'],'-dpdf')


% then four spiliopoulos
fg1s = figure(pos1s{:},ps{:});hold on
histogram(tausA(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tausA(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axss{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg1s,[cd() '\histograms\hist1s'],'-dpdf')


fg2s = figure(pos2s{:},ps{:});hold on
histogram(tausB(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tausB(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axss{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg2s,[cd() '\histograms\hist2s'],'-dpdf')


fg3s = figure(pos3s{:},ps{:});hold on
histogram(tausC(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tausC(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axss{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg3s,[cd() '\histograms\hist3s'],'-dpdf')


fg4s = figure(pos4s{:},ps{:});hold on
histogram(tausD(1,:),nbins,'EdgeColor','b',ds{:},lwd{:})
histogram(tausD(2,:),nbins,'EdgeColor','r',ds{:},lwd{:})
legend('\tau_1','\tau_2','Location','NorthEast')
set(gca,fs{:},axss{:});xlabel(lbl{1});ylabel(lbl{2});
set(gca, 'YScale', 'log')
box on
print(fg4s,[cd() '\histograms\hist4s'],'-dpdf')

