clear
close all

%% load
filenames = {'model_GPR_b1.mat','model_GPR_rb.mat','model_GPR_f.mat'};
data_p = load('../parameters_block.mat');
for j = 1:3
    data = load(filenames{j});
    Y_GT = data.Y_GT;
    Y_infer = data.Y_infer;
    
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
    plot([ymin-Dy/8,ymax+Dy/8],[ymin-Dy/8,ymax+Dy/8],'-r')
    
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

    AE = abs(Y_infer-Y_GT);
    [maxAE, i_maxAE] = max(AE);
    disp(['maxAE = ',num2str(maxAE)])
    disp(['i_maxAE = ',num2str(i_maxAE)])
    
end