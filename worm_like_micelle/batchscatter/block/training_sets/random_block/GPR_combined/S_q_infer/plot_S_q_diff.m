clear
close all

%% Load
data_rod = load('scatter_rod.mat');
data_GT = load('S_q_GT.mat');
data_infer = load('S_q_infer.mat');

q = data_rod.qq;
qL = q*1000;
S_q_rod = data_rod.S_q;
delta_S_q = data_rod.delta_S_q;

S_q_GT = data_GT.S_q_maxAE;
S_q_infer = data_infer.S_q;

S_q_rod = S_q_rod.*delta_S_q;
S_q_GT = S_q_GT'.*delta_S_q;
S_q_infer = S_q_infer'.*delta_S_q;

%% color
% x_color = (p_color-min(p_color))/(max(p_color)-min(p_color));
% color_parula = parula(101);
% color = arrayfun(@(x) color_parula(1+round(x*100),:),x_color,'UniformOutput',false);

%% plot S_q
filenames = {'S_q_b1.mat','S_q_rb.mat','S_q_f.mat'};
for i = 1:3
    % plot setup
    figure(i)
    hold on
    box on
    set(gcf,'Position',[2000,100,600,600])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.22    0.22   0.72    0.72])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
    ylabel('{\it S}({\itQL})','FontSize',24,'Interpreter','tex')
    set(gca,'FontSize',28,'FontName','Arial')
    
    % background
    for j = 1:9
        q_i = [1e-1 1e3];
        S_q_i_1 = 10^(j-5)./q_i;
        S_q_i_2 = 10^(2*j-10)./q_i.^2;
        plot(q_i,S_q_i_1,'--','LineWidth',1,'Color','#C0C0C0')
        plot(q_i,S_q_i_2,':','LineWidth',2,'Color','#C0C0C0')
    end
    
%     % S_q_rod
%     plot(qL,S_q_rod,'k-','LineWidth',2)
    
    % S_q
    plot(qL,S_q_GT(i,:),'k-','LineWidth',2)
    plot(qL,S_q_infer(i,:),'r-','LineWidth',2)
    
    
    xlim([1e-1,1e3])
    ylim([1e-3,2e0])
    xticks([1e-1,1e0,1e1,1e2,1e3])
    yticks([1e-3,1e-2,1e-1,1e0])
    saveas(gcf,[filenames{i},'.png'])
end

%% plot delta S_q
filenames = {'diff_b1.mat','diff_rb.mat','diff_f.mat'};
for i = 1:3
    % plot setup
    figure(i+3)
    hold on
    box on
    set(gcf, 'color', 'none')
    set(gcf,'Position',[2000,100,400,400])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.25    0.25   0.6    0.6])
    set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
    ylabel('diff. (%)','FontSize',24,'Interpreter','tex')
    set(gca,'FontSize',24,'FontName','Arial')

    
%     % S_q_rod
%     plot(qL,S_q_rod./S_q_rod,'k-','LineWidth',2)
    
    % S_q
    plot(qL,qL*0,'k-','LineWidth',2)
    plot(qL,(S_q_infer(i,:)-S_q_GT(i,:))./S_q_GT(i,:)*100,'r-','LineWidth',2)
    
    
    xlim([1e-1,1e3])
    ylim([-6e0,1e0])
    xticks([1e-1,1e0,1e1,1e2,1e3])
    % yticks([1e0,1e1])
    saveas(gcf,[filenames{i},'.svg'])
end