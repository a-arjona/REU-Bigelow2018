%% Plot results from cluster analysis in Sub-Arctic Atlantic
%
%
%

clear all
close all
clc

% open metadata
T = readtable('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/metadata.txt');

maxo        = table2array(T(:,1));
xo          = table2array(T(:,2));
yo          = table2array(T(:,3));
good_points = table2array(T(:,4));

T = readtable('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/clusters.txt');
cluster     = table2array(T);
T = readtable('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/si.txt');
si          = table2array(T);


% open timeseries
T = readtable('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/timeseries.txt');

timeseries  = table2array(T);

%% plot map
nclu = 6;
nclu = nclu - 1;
lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);

map_ccc = nan(size(LON));
for ii = 1:length(xo)
    map_ccc(yo(ii),xo(ii)) = cluster(ii,nclu);
%     map_ccc(yo(ii),xo(ii)) = cluster3(ii);
end

figure(2)
set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites

col = cbrewer('qual','Set1',max(cluster(:,nclu)));
ncol = [col(1,:); col(2,:); col(4,:); col(3,:); col(5,:); col(6,:)];
colormap(ncol)
[X,Y] = m_ll2xy(LON,LAT);
h = pcolor(X,Y,map_ccc);
set(h,'edgecolor','none')

tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])

m_grid()
m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
title("20-year Mean")

%% timeseries
figure(20)

M = [3:10];
D = ones(1,length(M));
Y = repmat(2000,1,length(M));
datem = datenum('01-01-2000')+ [3:8:365];
xtick = datenum(Y,M,D);

cpt = 0;
sub = [1 2 3 4 5 6];
for ccc = [1 2 4 3 5 6]
    f = find(cluster(:,nclu) == ccc);
    cpt = cpt + 1;
    subplot(2,3,sub(cpt))
    plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f)),'-','color',ncol(ccc,:),'LineWidth',2)
    hold on
    plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f))+std(timeseries(f,:))*mean(maxo(f)),'k-')
    plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f))-std(timeseries(f,:))*mean(maxo(f)),'k-')
    set(gca,'Ylim',[0 3],'Xtick',xtick,'Xticklabel',datestr(xtick,'mmm'),'Xlim',[xtick(1) xtick(end)+10],...
        'Xgrid','on','Ygrid','on','fontsize',9)
    ylabel('[Chl-a] (mg m^-^3)','fontsize',9)
    title("Mean Time Series")
    
end

%% Mapping all number of clusters
close all

lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;
[LON, LAT] = meshgrid(lon,lat);


for nclu = 2:9
    nclu = nclu - 1;
    
    subplot(2,5,nclu)
    
    map_ccc = nan(size(LON));
    for ii = 1:length(xo)
        map_ccc(yo(ii),xo(ii)) = cluster(ii,nclu);
    end

    set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites

    colormap(gca,ncol(1:(nclu+1),:))
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')

%     tb = m_tbase('contour',[-5000:1000:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])

    m_grid('fontsize',6)
    [X,Y] = m_ll2xy(-36,78.5);
    txt = ['Grp = ',num2str(nclu+1)];
    text(X,Y,txt,'BackgroundColor','w','edgecolor','k')
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    
end

%% All timeseries with all number of cluster
close all

M = [3:10];
D = ones(1,length(M));
Y = repmat(2000,1,length(M));
datem = datenum('01-01-2000')+ [3:8:365];
xtick = datenum(Y,M,D);
col = cbrewer('seq','Greys',10);
col = flipud(col(2:end,:));
col2 = cbrewer('qual','Set1',9);
H = [];

for nclu = 1:8
    cpt = 0;
    for ccc = 1:max(cluster(:,nclu))
        f = find(cluster(:,nclu) == ccc);
        cpt = cpt + 1;
        subplot(3,3,cpt)
        if nclu ~= 8
            h = plot(datem(9:35), mean(timeseries(f,:)),'-','color',col(nclu,:),'LineWidth',2);
        else
            h = plot(datem(9:35), mean(timeseries(f,:)),'-','color',col2(ccc,:),'LineWidth',2);
        end
        hold on
%         plot(datem(14:31), mean(timeseries(f,:))*mean(maxo(f))+std(timeseries(f,:))*mean(maxo(f)),'k-')
%         plot(datem(14:31), mean(timeseries(f,:))*mean(maxo(f))-std(timeseries(f,:))*mean(maxo(f)),'k-')
        set(gca,'Ylim',[0 1],'Xtick',xtick,'Xticklabel',datestr(xtick,'mmm'),'Xlim',[xtick(1) xtick(end)+10],...
            'Xgrid','on','Ygrid','on','fontsize',9)
        ylabel('[Chl-a]_n_o_r_m','fontsize',9)
        H = [H, h];
%         if ccc == 7
%             set(gca,'Ylim',[0 6])
%         elseif ccc == 2 ||ccc == 4 || ccc == 6
%             set(gca,'Ylim',[0 2])
%         end
    end
end
legend([H(1),H(22),H(36)],{'grp = 2','grp = 6','grp = 9'},'Location','NorthEast')


%% All Silhouette analysis
 close all

%col = cbrewer('qual','Set1',9);
cpt = 0;
sub = [1 3 5 7 2 4 6 8];
    
for nclu = 1:8
    x = 0;
    cpt = cpt + 1;
    for ccc = 1:max(cluster(:,nclu))
        subplot(4,2,sub(cpt))
        
        f = find(cluster(:,nclu) == ccc);
        Y = sort(si(f,nclu),'descend');
        X = [x:(x+length(Y)-1)];
        area(X,Y,'FaceColor',ncol(ccc,:),'LineStyle','none')
        hold on
        x = X(end)+1;
        set(gca,'Ylim',[-0.2 0.55],'Xlim',[-100 length(cluster)+100], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),25),[num2str(round(sum(Y < mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),75),[num2str(round(sum(Y >= mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
    
    end
    plot([0 length(cluster)],[mean(si(:,nclu)) mean(si(:,nclu))],'k-')
    ylabel('Silhouette value', 'fontsize',9)
end

legend([H(1),H(22),H(36)],{'grp = 2','grp = 6','grp = 9'},'Location','southeast')
%% save

set(gcf,'PaperPosition',[0 0 10 8])
txt_save = 'C:\Users\Ade\Documents\Bigelow\DINEOF-2018\silhouette-value-DINEOF.png';
print(txt_save,'-dpng','-r600')
