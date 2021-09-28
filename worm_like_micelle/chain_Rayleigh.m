function [lc,Cc,O,n] = chain(DP,a,lambda,unit_C)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
n=zeros(3,DP);
l=zeros(3,DP);
lc=zeros(3,DP);
B=zeros(3,3);
C=zeros(3,3);
D=zeros(3,3);
R=zeros(3,3);
O=zeros(3,3,DP);
for i=1:DP
    if i==1
        n(:,i) = [1;0;0];
        l(:,i) = n(:,i);
        B(:,:,i) = eye(3);
        C(:,:,i) = eye(3);
        D(:,:,i) = eye(3);
        R(:,:,i) = eye(3);
        O(:,:,i) = R(:,:,i);
    else
        
%         % xyz convention
%         psi = 0; % Roll angle
%         theta = 2*(rand-0.5)*pi/a/2; % Pitch angle 
%         phi = 2*(rand-0.5)*pi/a/2; % Yaw angle
% 
%         D(:,:) = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]; % Yaw
%         C(:,:) = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]; % Pitch
%         B(:,:) = [1 0 0; 0 cos(psi) sin(psi); 0 -sin(psi) cos(psi)]; % Roll
%         
%         R = B*C*D;
        
        % quaternion
        phi_q = 2*(rand-0.5)*pi/2;
        theta_q = sqrt(-log(1-rand)/a);
        
        vq = O(:,2,i-1)*cos(phi_q) + O(:,3,i-1)*sin(phi_q);
        qr = cos(theta_q/2);
        qi = vq(1)*sin(theta_q/2);
        qj = vq(2)*sin(theta_q/2);
        qk = vq(3)*sin(theta_q/2);
        
        Rq = [1-2*(qj^2+qk^2) 2*(qi*qj-qk*qr) 2*(qi*qk+qj*qr);...
              2*(qi*qj+qk*qr) 1-2*(qi^2+qk^2) 2*(qj*qk-qi*qr);...
              2*(qi*qk-qj*qr) 2*(qj*qk+qi*qr) 1-2*(qj^2+qk^2)];
        
        R = Rq;
                
        O(:,:,i) = R*O(:,:,i-1);
        n(:,i) = R*n(:,i-1);
        l(:,i) = l(:,i-1) + n(:,i);
        
    end
end
lc = l*lambda;
%% map unimer
%C
nC = size(unit_C,2);
for j=1:DP
        for k=1:nC
            m_backbone_C(:,k,j) = O(:,:,j)*unit_C(:,k)+lc(:,j)+[0;0;0];
        end
end
Cc = reshape(m_backbone_C,[3,DP*nC]);
end