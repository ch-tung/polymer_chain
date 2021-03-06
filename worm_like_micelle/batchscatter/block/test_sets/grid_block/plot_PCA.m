clear
close all

%% Load
data = load('scatter_chain_block_0_g.mat');
data_rod = load('scatter_rod.mat');
data_homo = load('scatter_homo.mat');

q = data_rod.qq;
qL = q*1000;
S_q_rod = data_rod.S_q;
delta_S_q = data_rod.delta_S_q;

S_q_0 = data.S_q;

S_q = zeros(length(q),11*size(S_q_0,2));
p = zeros(3,11*size(S_q_0,2));
for i = 0:10
    filename = ['scatter_chain_block_', num2str(i) ,'_g.mat'];
    data = load(filename);
    S_q(:,1+i*size(S_q_0,2):(i+1)*size(S_q_0,2)) = data.S_q.*delta_S_q';
    p(:,1+i*size(S_q_0,2):(i+1)*size(S_q_0,2)) = data.p;
end

S_q_rod = S_q_rod.*delta_S_q;
S_q_homo = data_homo.S_q_homo;

F = S_q'./S_q_rod;
F = log(F);

F_homo = S_q_homo'./S_q_rod;
F_homo = log(F_homo);

X = F;
X_homo = F_homo;
%% indexing
b1 = p(2,:)/1000;
rb = p(1,:);
f  = p(3,:);

index_b1 = b1==0.01;
index_rb = rb==10;
index_f = f==0.5;

%% SVD
[U,S,V] = svd(X');
score_X = V*S';
score_X_homo = X_homo*U;

[sort_X,I] = sort(score_X(:,1));

N = size(X,1);
Dim = size(X,2);

Lambda_C = (S*S')/(N-1);

percentage_C = diag(Lambda_C)/Dim;

%% plot scree
figure;
hold on

Lambda_C = (S*S')/(N-1);

percentage_C = diag(Lambda_C);
singular_C = diag(S);

plot(1,singular_C(1),'o','MarkerSize',12,'LineWidth',2)
plot(2,singular_C(2),'o','MarkerSize',12,'LineWidth',2)
plot(3,singular_C(3),'o','MarkerSize',12,'LineWidth',2)
plot(4:Dim,singular_C(4:Dim),'o','Color','#808080','MarkerSize',12)
hold on

% figure settings
box on
xlabel('SVR','FontSize',24)
ylabel('\Sigma','FontSize',24)
% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')
set(gca,'LineWidth',2)
set(gca,'position',[0.22    0.22   0.72    0.72])

%% likelihood
list_q = 1:Dim;

% v = var(singular_C);

for i = 1:length(list_q)
    q = list_q(i);
    
    mu_1(i) = sum(singular_C(1:q))/q;
    mu_2(i) = sum(singular_C(q+1:Dim))/(Dim-q);
    
    v1(i) = var(singular_C(1:q));
    v2(i) = var(singular_C(q+1:Dim));
    
    if isnan(v2(i))
        v2(i)=0;
    end
    
    vq(i) = ((q-1)*v1(i) + (Dim-q-1)*v2(i))/(Dim-2);
    
    syms d
    syms mu
    syms v
    fun = @(d,mu,v) 1/sqrt(2*pi*v)*exp(-(d-mu)^2/2/v);
    
    lq1 = 0;
    lq2 = 0;
    
    for j = 1:q
        lq1 = lq1 + log(fun(singular_C(j),mu_1(i),vq(i)));
    end
    for k = q+1:Dim
        lq2 = lq2 + log(fun(singular_C(k),mu_2(i),vq(i)));
    end
    
    lq(i) = lq1 + lq2;
    
    
end

figure;
hold on
plot(list_q(1:20),lq(1:20),'-','LineWidth',2,'Color','#808080')
plot(list_q(4:20),lq(4:20),'o','LineWidth',2,'Color','#808080')
plot(list_q(1),lq(1),'o','LineWidth',2)
plot(list_q(2),lq(2),'o','LineWidth',2)
plot(list_q(3),lq(3),'o','LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')
set(gca,'LineWidth',2)
set(gca,'position',[0.22    0.22   0.72    0.72])

%% plot basis
QU = qL;
figure;
plot(repmat(QU,3,1)',U(:,1:3),'LineWidth',2)

xlim([1e-1,1e3])
xticks([1e-1,1e0,1e1,1e2,1e3])
ylim([-3e-1,3e-1])
yticks([-3e-1:1e-1:3e-1])
xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
ylabel('Singular Vector','FontSize',24,'Interpreter','tex')
set(gca, 'XScale', 'log')
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')
set(gca,'LineWidth',2)
set(gca,'position',[0.22    0.22   0.72    0.72])

%% plot SVD
index = 1:size(S_q,2);

for i = 1:3
    switch(i)
    case 1
%         index_rb = rb==1;
%         index = index_rb&index_f; % b1
        p_color = log10(b1); % b1
    case 2
%         index = index_b1&index_f; % rb
        p_color = log10(rb); % rb
    case 3
%         index = index_b1&index_rb; % f
        p_color = f; % f
    end
    
    % plot3D
    figure;
    hold on
    scatter3(score_X(index,1),score_X(index,2),score_X(index,3),200,p_color,'.')
    plot3(score_X_homo(:,1),score_X_homo(:,2),score_X_homo(:,3),'-r','LineWidth',2)
    plot3(0,0,0,'x','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k')
    
    % figure settings
    daspect([1 1 1])
    view(50,30)
%     view(90,0)
%     view(-360,90)
%     view(-360,0)
    
    box on
    grid on
    % caxis([0 1])
    % load tw_flag.mat

    xticks(-10:10)
    yticks(-10:10)
    zticks(-10:10)
    % xticklabels({'',''})
    % yticklabels({'',''})
    xlabel('SVD0','FontSize',24)
    ylabel('SVD1','FontSize',24)
    zlabel('SVD2','FontSize',24)
    axis equal
    axesLabelsAlign3D
%     xlim([-6 0])
%     ylim([-3 3])
%     zlim([-2 2])
    % shading flat
    % set(gca,'Position', get(gca, 'OuterPosition') - ...
    % get(gca,'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Arial')
    % legend('boxoff')
    set(gca,'LineWidth',2)
    
end
