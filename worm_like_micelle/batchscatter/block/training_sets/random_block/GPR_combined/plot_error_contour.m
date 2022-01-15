clear
close all

%% load
filenames = {'model_GPR_b1.mat','model_GPR_rb.mat','model_GPR_f.mat'};
p_maxAE = zeros(3,3);
for j = 1:3
    data = load(filenames{j});
    Y_GT = data.Y_GT;
    Y_infer = data.Y_infer;
    Y_infer_std = data.Y_infer_std;
    
    AE = abs(Y_infer-Y_GT);
    [maxAE, i_maxAE] = max(AE);
    
    ymin = min(Y_GT);
    ymax = max(Y_GT);
    ypmin = min(Y_infer);
    ypmax = max(Y_infer);
    
    n_grid = 100;
    Dy = (ymax-ymin);
    dy = Dy/n_grid;
    X = ymin-Dy/2:dy:ymax+Dy/2;
    Y = ymin-Dy/2:dy:ymax+Dy/2;
    
    [XX,YY] = meshgrid(X,Y);

    D_all = zeros(size(XX,1),size(XX,2));
    
    for i = 1:length(Y_GT)
        index_X = ceil((Y_GT(i)-(min(X)))/dy);
        index_Y = ceil((Y_infer(i)-(min(Y)))/dy);
        if index_X<size(D_all,1)&&index_Y<size(D_all,2)
            D_all(index_X,index_Y) = D_all(index_X,index_Y)+1;
        end
    end
    
    D_all = D_all';
    D_all_sm = imgaussfilt(D_all,1);
    
    figure;
    axis equal
    hold on
    box on
    
    pcolor(XX((X>ymin-Dy/16)&(X<ymax+Dy/16),(X>ymin-Dy/16)&(X<ymax+Dy/16)),...
        YY((X>ymin-Dy/16)&(X<ymax+Dy/16),(X>ymin-Dy/16)&(X<ymax+Dy/16)),...
        D_all((X>ymin-Dy/16)&(X<ymax+Dy/16),(X>ymin-Dy/16)&(X<ymax+Dy/16)))
    plot([ymin-Dy/8,ymax+Dy/8],[ymin-Dy/8,ymax+Dy/8],'-','Color','#A0A0A0')
    plot([Y_GT(i_maxAE),Y_GT(i_maxAE)],...
        [Y_infer(i_maxAE),Y_GT(i_maxAE)],'--','LineWidth',2,'MarkerSize',20,'Color','r')
%     plot(Y_GT(i_maxAE),Y_infer(i_maxAE),'xr','LineWidth',2,'MarkerSize',20)
%     plot(Y_GT(i_maxAE),Y_GT(i_maxAE),'sk','LineWidth',8,'MarkerSize',24)
    plot(Y_GT(i_maxAE),Y_GT(i_maxAE),'ow','LineWidth',2,'MarkerSize',16,...
    'MarkerFaceColor','w','MarkerEdgeColor','k')
    plot(Y_GT(i_maxAE),Y_infer(i_maxAE),'or','LineWidth',2,'MarkerSize',16,...
    'MarkerFaceColor','r','MarkerEdgeColor','k')

%     errorbar(Y_GT(i_maxAE),Y_infer(i_maxAE),Y_infer_std(i_maxAE),...
%         'r','LineWidth',2,'MarkerSize',20,'CapSize',12)
    
    caxis([0 4])
    colormap(flipud(gray))
    shading flat
    xlim([ymin-Dy/8,ymax+Dy/8])
    ylim([ymin-Dy/8,ymax+Dy/8])
    
    xticks((10^floor(log10(abs(ymax+Dy/8))))*[-10:10]);
    yticks((10^floor(log10(abs(ymax+Dy/8))))*[-10:10]);
    
    if j==1
        xticks((10^floor(log10(abs(ymax+Dy/8))))*[-10:0.5:10]);
        yticks((10^floor(log10(abs(ymax+Dy/8))))*[-10:0.5:10]);
    end
    
    xlabel('','FontSize',24,'Interpreter','tex')
    ylabel('','FontSize',24,'Interpreter','tex')
    set(gcf,'Position',[100,100,600,600])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.22    0.22   0.72    0.72])
    set(gca,'FontSize',28,'FontName','Arial')
    saveas(gcf,[filenames{j},'.png'])
end

Y_GT_all = zeros(length(Y_GT),3);
Y_infer_all = zeros(length(Y_GT),3);
for j = 1:3
    data = load(filenames{j});
    Y_GT = data.Y_GT;
    Y_infer = data.Y_infer;
    
    AE = abs(Y_infer-Y_GT);
    [maxAE(j), i_maxAE(j)] = max(AE);
    
    Y_GT_all(:,j) = Y_GT;
    Y_infer_all(:,j) = Y_infer;
end
disp(['i_maxAE = ',num2str(i_maxAE)])
p_GT_maxAE = Y_GT_all(i_maxAE,:);
p_infer_maxAE = Y_infer_all(i_maxAE,:);

% %% plot prediction cloud
% figure(4);
% axis equal
% E_all = abs(Y_infer_all-Y_GT_all);
% RMSE_all = std(E_all);
% hold on
% plot3(E_all(:,1)/RMSE_all(1),E_all(:,2)/RMSE_all(2),E_all(:,3)/RMSE_all(3),...
%     'k.','MarkerSize',1)
% grid on