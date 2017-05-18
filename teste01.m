clear
close all

%====================================================================================
%====================================================================================
plotnodetext = 1;
plot_model = 1;

%Modelo
modelo05   % modelo xx - trelica a ser analisada

model_name = 'modelo05';

% materail properties (this case for alumninum)
E   = 70e9;  % young's modulus (modulo de elasticidade)
rho = 2710.0  % density

% cross section area
d = 0.01;
A = pi*d^2/4;
%====================================================================================
%====================================================================================

[ne,nc] = size(elem); % get ne = number of elements = number of lines in elements
[nn,ncn] = size(node); % get nn = number of nodes = number of lines in nodes

KS = zeros(nn*3);  % system stiffness matrix
MS = zeros(nn*3);  % system mass matrix

for e=1:ne
        E = elem(e,4);
        rho = elem(e,5);
        A = elem(e,6);
        x1 = node(elem(e,2),2);
        y1 = node(elem(e,2),3);
        z1 = node(elem(e,2),4);
        x2 = node(elem(e,3),2);
        y2 = node(elem(e,3),3);
        z2 = node(elem(e,3),4);
        L = sqrt((x2-x1)^2+(y2-y1)^2+(z1-z2)^2);
        Ke = (E*A/L)*[1 0 0 -1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; -1 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
        Me = (rho*A*L/6)*[2 0 0 1 0 0; 0 2 0 0 1 0; 0 0 2 0 0 1; 1 0 0 2 0 0; 0 1 0 0 2 0; 0 0 1 0 0 2];
        T = trans_matrix_beam3d(x1,y1,z1,x2,y2,z2);
        TT = [T zeros(3);zeros(3) T];
        Keg = transpose(TT)*Ke*TT;
        Meg = transpose(TT)*Me*TT;
        
        %Stiffnes of the local system
        Keg11 = Keg(1:3,1:3);
        Keg12 = Keg(1:3,4:6);
        Keg21 = Keg(4:6,1:3);
        Keg22 = Keg(4:6,4:6);
        
        %Mass of the local system
        Meg11 = Meg(1:3,1:3);
        Meg12 = Meg(1:3,4:6);
        Meg21 = Meg(4:6,1:3);
        Meg22 = Meg(4:6,4:6);
        
        no1 = node(elem(e,2),1);
        no2 = node(elem(e,3),1);
        
        %disp(sprintf('Elemento %d, no 1: %d, no 2: %d',e,no1,no2))
        
        %Add Stiffness of the element in the global matrix
        KS(no1*3-2:no1*3,no1*3-2:no1*3) = KS(no1*3-2:no1*3,no1*3-2:no1*3) + Keg11;
        KS(no1*3-2:no1*3,no2*3-2:no2*3) = KS(no1*3-2:no1*3,no2*3-2:no2*3) + Keg12;
        KS(no2*3-2:no2*3,no1*3-2:no1*3) = KS(no2*3-2:no2*3,no1*3-2:no1*3) + Keg21;
        KS(no2*3-2:no2*3,no2*3-2:no2*3) = KS(no2*3-2:no2*3,no2*3-2:no2*3) + Keg22;
        
        %Add Mass of the element in the global matrix
        MS(no1*3-2:no1*3,no1*3-2:no1*3) = MS(no1*3-2:no1*3,no1*3-2:no1*3) + Meg11;
        MS(no1*3-2:no1*3,no2*3-2:no2*3) = MS(no1*3-2:no1*3,no2*3-2:no2*3) + Meg12;
        MS(no2*3-2:no2*3,no1*3-2:no1*3) = MS(no2*3-2:no2*3,no1*3-2:no1*3) + Meg21;
        MS(no2*3-2:no2*3,no2*3-2:no2*3) = MS(no2*3-2:no2*3,no2*3-2:no2*3) + Meg22;
        
            
end
%disp (MS)

[ave,ava] = eig(KS,MS);

[AVA, i] = sort(diag(ava));

frequencies = sqrt(AVA);

AVE = ave(:,i);

mode = 12; %mode of vibration

factor = 0.4; 

for e=1:ne
    no1 = node(elem(e,2),1);
    no2 = node(elem(e,3),1);
    
    x1 = node(elem(e,2),2);
    y1 = node(elem(e,2),3);
    z1 = node(elem(e,2),4);
    x2 = node(elem(e,3),2);
    y2 = node(elem(e,3),3);
    z2 = node(elem(e,3),4);
    
    %new coordinates
    nx1 =x1 + factor*AVE(3*no1-2,mode); %3 = degrees of freedom
    ny1 =y1 + factor*AVE(3*no1-1,mode);
    nz1 =z1 + factor*AVE(3*no1,mode);
    nx2 =x2 + factor*AVE(3*no2-2,mode);
    ny2 =y2 + factor*AVE(3*no2-1,mode);
    nz2=z2 +factor*AVE(3*no2,mode);
    
    set(gcf,'color','w');
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    plot3([nx1 nx2], [nz1,nz2],[ny1,ny2],'k','LineWidth', 1.0);
    
    hold on
    
end

axis equal

% %%
% % nodex = node(:,2);
% % nodey = node(:,3);
% % nodez = node(:,4);
% % avex = AVE(1:3:end,mode);
% % avey = AVE(2:3:end,mode);
% % avez = AVE(3:3:end,mode);
% % newnodex = (nodex+avex)*factor;
% % newnodey = (nodey+avey)*factor;
% % newnodez = (nodez+avez)*factor;
% % 
% % t = linspace(-1,1,100);
% % 
% % for c1 = 1:1:length(t)
% %    newnodex2(:,c1) = nodex+avex*t(c1)*factor;
% %    newnodey2(:,c1) = nodey+avey*t(c1)*factor;
% %    newnodez2(:,c1) = nodez+avez*t(c1)*factor;
% % end
% t = 1:-0.2:-1;
% t = 0:.1:2*pi;
% figure(2)
% for c2 = 1:1:length(t)
%    
%    for e=1:ne
%         no1 = node(elem(e,2),1);
%         no2 = node(elem(e,3),1);
% 
%         x1 = node(elem(e,2),2);
%         y1 = node(elem(e,2),3);
%         z1 = node(elem(e,2),4);
%         x2 = node(elem(e,3),2);
%         y2 = node(elem(e,3),3);
%         z2 = node(elem(e,3),4);
% 
%         %new coordinates
%         nx1 = x1 + sin(t(c2))*factor*AVE(3*no1-2,mode); %3 = degrees of freedom
%         ny1 = y1 + sin(t(c2))*factor*AVE(3*no1-1,mode);
%         nz1 = z1 + sin(t(c2))*factor*AVE(3*no1,mode);
%         nx2 = x2 + sin(t(c2))*factor*AVE(3*no2-2,mode);
%         ny2 = y2 + sin(t(c2))*factor*AVE(3*no2-1,mode);
%         nz2 = z2 + sin(t(c2))*factor*AVE(3*no2,mode);
% 
%         set(gcf,'color','w');
%         axis off;
%         %set(gca,'xtick',[],'ytick',[],'ztick',[]);
%         plot3([nx1 nx2], [nz1,nz2],[ny1,ny2],'k','LineWidth', 1.0);
%         hold on
% 
%    end
%    xlim([-.2,10.2]);
%    ylim([-.2,1.2]);
%    zlim([-.2,1.2]);
%    axis equal
%    drawnow
%    mov(c2) = getframe(gcf);
%    hold off
%end

% video_filename = sprintf('%s_mode_%d.avi',model_name,mode)
% 
% movie2avi(mov, video_filename, 'compression', 'None');
