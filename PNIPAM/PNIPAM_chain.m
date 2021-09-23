% clear
close all

%% load unimer
unimer_C = load('u_C.dat')';
unimer_H = load('u_H.dat')';
unimer_N = load('u_N.dat')';
unimer_O = load('u_O.dat')';

unimer_Ct = load('u_Ct.dat')';
unimer_Ht = load('u_Ht.dat')';
unimer_St = load('u_St.dat')';

nC = size(unimer_C,2);
nH = size(unimer_H,2);
nN = size(unimer_N,2);
nO = size(unimer_O,2);

nCt = size(unimer_Ct,2);
nHt = size(unimer_Ht,2);
nSt = size(unimer_St,2);


lambda = 7.523026017;
DP = round(250/6);
a = 4;
b = 2;
Cr = 16;
Nagg = 20; %20 or 30

%dodecahedron vertice
phi = (1+sqrt(5))/2;
D_C1 = [1 1 1;-1 1 1;1 -1 1;1 1 -1;1 -1 -1;-1 1 -1;-1 -1 1; -1 -1 -1];
D_C2 = [0 phi 1/phi;0 -phi 1/phi;0 phi -1/phi;0 -phi -1/phi];
D_C3 = [1/phi 0 phi;-1/phi 0 phi;1/phi 0 -phi;-1/phi 0 -phi];
D_C4 = [phi 1/phi 0;-phi 1/phi 0;phi -1/phi 0;-phi -1/phi 0];
D_C = [D_C1' D_C2' D_C3' D_C4'];

% %icosidodecahedron vertice
% D_Ci1 = phi*[1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
% D_Ci2 = 0.5*[1 phi phi^2;-1 phi phi^2;1 -phi phi^2;1 phi -phi^2;1 -phi -phi^2;-1 phi -phi^2;-1 -phi phi^2;-1 -phi -phi^2];
% D_Ci3 = 0.5*[phi phi^2 1;-phi phi^2 1;phi -phi^2 1;phi phi^2 -1;phi -phi^2 -1;-phi phi^2 -1;-phi -phi^2 1;-phi -phi^2 -1];
% D_Ci4 = 0.5*[phi^2 1 phi;-phi^2 1 phi;phi^2 -1 phi;phi^2 1 -phi;phi^2 -1 -phi;-phi^2 1 -phi;-phi^2 -1 phi;-phi^2 -1 -phi];
% D_Ci = [D_Ci1' D_Ci2' D_Ci3' D_Ci4'];
% plot3(0,0,0,'o','Color','k','MarkerSize',180,'MarkerFaceColor','#666666')
% hold on
% plot3(D_Ci(1,:),D_Ci(2,:),D_Ci(3,:),'o','Color','k','MarkerSize',60,'MarkerFaceColor','#FF0000')
% view([1 1 1])
% set (gcf,'Position',[0,0,1000,1000])
% axis equal

%orientation
GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
              norm(cross(A,B)) dot(A,B)  0;
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);

v0 = [1;0.1;0.1];

for i = 1:Nagg
vi = D_C(:,i);
U(:,:,i) = UU(FFi(v0,vi), GG(v0,vi));
end

for d = 1:Nagg
[Cc,Hc,Nc,Oc,CCN,NCN] = chain(DP,unimer_C,unimer_H,unimer_N,unimer_O,nC,nH,nN,nO,lambda,a,Cr);
[u_Ct,u_Ht] = chain_tail(unimer_Ct,unimer_Ht,b);
Cd(:,:,d) = U(:,:,d)*Cc;
Hd(:,:,d) = U(:,:,d)*Hc;
Nd(:,:,d) = U(:,:,d)*Nc;
Od(:,:,d) = U(:,:,d)*Oc;
Ct(:,:,d) = U(:,:,d)*u_Ct;
Ht(:,:,d) = U(:,:,d)*u_Ht;
St(:,:,d) = U(:,:,d)*unimer_St;
C_CN(:,d) = U(:,:,d)*CCN;
N_CN(:,d) = U(:,:,d)*NCN;
end

CC = reshape(Cd,[3,DP*nC*Nagg]);
HH = reshape(Hd,[3,DP*nH*Nagg]);
NN = reshape(Nd,[3,DP*nN*Nagg]);
OO = reshape(Od,[3,DP*nO*Nagg]);

CT = reshape(Ct,[3,nCt*Nagg]);
HT = reshape(Ht,[3,nHt*Nagg]);
ST = reshape(St,[3,nSt*Nagg]);

%% plot
% %head
% figure;
% box on
% plot3(0,0,0);
% hold on
% pC1=scatter3(CC(1,:),CC(2,:),CC(3,:),48,...  
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% pH1=scatter3(HH(1,:),HH(2,:),HH(3,:),12,...
%     'filled','o','MarkerFaceColor','#FFFFFF','MarkerFaceAlpha',0.99);
% pN1=scatter3(NN(1,:),NN(2,:),NN(3,:),48,...
%     'filled','o','MarkerFaceColor','#0000FF','MarkerFaceAlpha',0.99);
% pO1=scatter3(OO(1,:),OO(2,:),OO(3,:),48,...
%     'filled','o','MarkerFaceColor','#FF0000','MarkerFaceAlpha',0.99);
% 
% %tail
% pCT1=scatter3(CT(1,:),CT(2,:),CT(3,:),48,...  
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% pHT1=scatter3(HT(1,:),HT(2,:),HT(3,:),12,...
%     'filled','o','MarkerFaceColor','#FFFFFF','MarkerFaceAlpha',0.99);
% pST1=scatter3(ST(1,:),ST(2,:),ST(3,:),96,...
%     'filled','o','MarkerFaceColor','#EDB120','MarkerFaceAlpha',0.99);
% 
% %CN
% pCCN1=scatter3(C_CN(1,:),C_CN(2,:),C_CN(3,:),48,...
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% pNCN1=scatter3(N_CN(1,:),N_CN(2,:),N_CN(3,:),192,...
%     'filled','o','MarkerFaceColor','#7E2F8E','MarkerFaceAlpha',0.99);
% 
% axis off
% set (gcf,'Position',[0,0,1000,1000])
% axis equal
% % xlim([-600 600])
% % ylim([-600 600])
% % zlim([-600 600])
% view([1 0 0])
% % axis normal;

figure;
%head
% box on
plot3(0,0,0);
hold on
% pC2=plot3(CC(1,:),CC(2,:),CC(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% pH2=plot3(HH(1,:),HH(2,:),HH(3,:),  'o','Color','k','MarkerSize',4,'MarkerFaceColor','#FFFFFF');
% pN2=plot3(NN(1,:),NN(2,:),NN(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#0000FF');
% pO2=plot3(OO(1,:),OO(2,:),OO(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#FF0000');

pC2=plot3(CC(1,:),CC(2,:),CC(3,:),  'o','Color','#004060','MarkerSize',6,'MarkerFaceColor','#0080C0');
pH2=plot3(HH(1,:),HH(2,:),HH(3,:),  'o','Color','#004060','MarkerSize',4,'MarkerFaceColor','#0080C0');
pN2=plot3(NN(1,:),NN(2,:),NN(3,:),  'o','Color','#004060','MarkerSize',6,'MarkerFaceColor','#0080C0');
pO2=plot3(OO(1,:),OO(2,:),OO(3,:),  'o','Color','#004060','MarkerSize',6,'MarkerFaceColor','#0080C0');

%tail
% pCT2=plot3(CT(1,:),CT(2,:),CT(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% pHT2=plot3(HT(1,:),HT(2,:),HT(3,:),  'o','Color','k','MarkerSize',4,'MarkerFaceColor','#FFFFFF');
% pST2=plot3(ST(1,:),ST(2,:),ST(3,:),  'o','Color','k','MarkerSize',12,'MarkerFaceColor','#EDB120');

pCT2=plot3(CT(1,:),CT(2,:),CT(3,:),  'o','Color','#806000','MarkerSize',6,'MarkerFaceColor','#F08000');
pHT2=plot3(HT(1,:),HT(2,:),HT(3,:),  'o','Color','#806000','MarkerSize',4,'MarkerFaceColor','#F08000');
pST2=plot3(ST(1,:),ST(2,:),ST(3,:),  'o','Color','#806000','MarkerSize',6,'MarkerFaceColor','#F08000');

%CN
% pCCN2=plot3(C_CN(1),C_CN(2),C_CN(3),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% pNCN2=plot3(N_CN(1),N_CN(2),N_CN(3),  'o','Color','k','MarkerSize',12,'MarkerFaceColor','#AA00AA');

pCCN2=plot3(C_CN(1,:),C_CN(2,:),C_CN(3,:),  'o','Color','#004060','MarkerSize',6,'MarkerFaceColor','#0080C0');
pNCN2=plot3(N_CN(1,:),N_CN(2,:),N_CN(3,:),  'o','Color','#004060','MarkerSize',6,'MarkerFaceColor','#0080C0');

%CN
% pCCN2=plot3(C_CN(1,:),C_CN(2,:),C_CN(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% pNCN2=plot3(N_CN(1,:),N_CN(2,:),N_CN(3,:),  'o','Color','k','MarkerSize',12,'MarkerFaceColor','#AA00AA');

axis off
set (gcf,'Position',[0,0,800,800])
axis equal
% xlim([-600 600])
% ylim([-600 600])
% zlim([-600 600])
view([1 0.01 0.01])
% axis normal;

%%
% set(gcf,'color','none');
% set(gca,'color','none');
% axes.SortMethod = 'ChildOrder';
% addpath(genpath('D:\Documents\Project\Matlab\export_fig-master\export_fig-master'));
% export_fig PNIPAM.png -r300 -a2 -opengl -transparent

%%
% %coarse grained
% %head
% figure;
% box on
% plot3(0,0,0);
% hold on
% pC1=scatter3(CC(1,1:6:end),CC(2,1:6:end),CC(3,1:6:end),96,...  
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% % pH1=scatter3(HH(1,:),HH(2,:),HH(3,:),12,...
% %     'filled','o','MarkerFaceColor','#FFFFFF','MarkerFaceAlpha',0.99);
% % pN1=scatter3(NN(1,:),NN(2,:),NN(3,:),48,...
% %     'filled','o','MarkerFaceColor','#0000FF','MarkerFaceAlpha',0.99);
% % pO1=scatter3(OO(1,:),OO(2,:),OO(3,:),48,...
% %     'filled','o','MarkerFaceColor','#FF0000','MarkerFaceAlpha',0.99);
% 
% %tail
% pCT1=scatter3(CT(1,:),CT(2,:),CT(3,:),48,...  
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% pHT1=scatter3(HT(1,:),HT(2,:),HT(3,:),12,...
%     'filled','o','MarkerFaceColor','#FFFFFF','MarkerFaceAlpha',0.99);
% pST1=scatter3(ST(1,:),ST(2,:),ST(3,:),96,...
%     'filled','o','MarkerFaceColor','#EDB120','MarkerFaceAlpha',0.99);
% 
% %CN
% pCCN1=scatter3(C_CN(1,:),C_CN(2,:),C_CN(3,:),48,...
%     'filled','o','MarkerFaceColor','#666666','MarkerFaceAlpha',0.99);
% pNCN1=scatter3(N_CN(1,:),N_CN(2,:),N_CN(3,:),192,...
%     'filled','o','MarkerFaceColor','#7E2F8E','MarkerFaceAlpha',0.99);
% 
% axis off
% set (gcf,'Position',[0,0,500,500])
% axis equal
% % xlim([-600 600])
% % ylim([-600 600])
% % zlim([-600 600])
% view([1 0 0])
% % axis normal;
% 
% % figure;
% % %head
% % % box on
% % plot3(0,0,0);
% % hold on
% % pC2=plot3(CC(1,:),CC(2,:),CC(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% % pH2=plot3(HH(1,:),HH(2,:),HH(3,:),  'o','Color','k','MarkerSize',4,'MarkerFaceColor','#FFFFFF');
% % pN2=plot3(NN(1,:),NN(2,:),NN(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#0000FF');
% % pO2=plot3(OO(1,:),OO(2,:),OO(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#FF0000');
% % 
% % %tail
% % pCT2=plot3(CT(1,:),CT(2,:),CT(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% % pHT2=plot3(HT(1,:),HT(2,:),HT(3,:),  'o','Color','k','MarkerSize',4,'MarkerFaceColor','#FFFFFF');
% % pST2=plot3(ST(1,:),ST(2,:),ST(3,:),  'o','Color','k','MarkerSize',12,'MarkerFaceColor','#EDB120');
% % 
% % %CN
% % pCCN2=plot3(C_CN(1,:),C_CN(2,:),C_CN(3,:),  'o','Color','k','MarkerSize',6,'MarkerFaceColor','#888888');
% % pNCN2=plot3(N_CN(1,:),N_CN(2,:),N_CN(3,:),  'o','Color','k','MarkerSize',12,'MarkerFaceColor','#AA00AA');
% % 
% % axis off
% % set (gcf,'Position',[0,0,1000,1000])
% % axis equal
% % % xlim([-600 600])
% % % ylim([-600 600])
% % % zlim([-600 600])
% % view([-1 1 -1])
% % % axis normal;