%% DefaultData
%
%
clear all

startdate='01-Jan-2000';
enddate='01-Jan-2016';

tab=readtable('poss_defaults_cleaned.csv');


hs=tab.health;
tick=tab.TICKER;
conm=tab.conm;
pd=tab.public_date;
dsel=logical((pd>=startdate).*(pd<enddate));

ticksel=tick(dsel);
conmsel=conm(dsel);
pdsel=pd(dsel);
hsel=hs(dsel);

%% Plot raw data
figure(1)
clf;
plot(pdsel,hsel,'*');

%% Sort into tickers
firms=unique(conmsel);

%% firm plot
f2=figure(2);
f2.Units='centimeters';
f2.OuterPosition=[2,2,14,9];
f2.PaperSize=[14 9];
clf;
hold on;

%firms=firms(2:end);
for i=1:size(firms,1)
    fsel=logical(strcmp(conmsel,firms{i}));
    pdf=pdsel(fsel);
    hsf=hsel(fsel);
    pl(i)=plot(pdf,hsf);
    pdend=pdf(end);
    c=get(pl(i),'Color');
    if pdend<enddate
        plot(pdf(end),hsf(end),'*','Color',c);
    end
end
yline(0);

xlabel('date');
ylabel('x = log inverse debt ratio');
ylim([-0.07,0.14]);
xlim([datetime('01-Jan-2000'), datetime('01-Jan-2014')]);

lgd=legend(pl,firms,'Location','southwest');
lgd.FontSize=5;

print('figdefaultdata.pdf','-dpdf');

keyboard;
    
    

