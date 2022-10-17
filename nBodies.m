%N-bodies problem 
%Given inputs: M,P,dt,tscale,duration
%M is a vector of masses in kilograms 
%P is a matrix of positions in meters 
%dt is in seconds, duration in dt's 

clear
figure('WindowState','Fullscreen','Color',[0,0,0])
BX = gca;
box; BX.BoxStyle = 'full'; BX.Color = 'k';
axis vis3d equal
BX.XTickLabel = [];BX.YTickLabel = [];BX.ZTickLabel = [];
BX.XTick = [];BX.YTick = [];BX.ZTick = [];
BX.XColor = [1,1,1];BX.YColor = [1,1,1];BX.ZColor = [1,1,1];
ROT = 20;
view(ROT,15)
load set
load masses
load positions
load colors4
load velocities
N = length(M);
if size(P,1) ~= length(M)
    disp('xwx')
    return
end
%%Initial velocities
%V = [
%    0, 0, 0
%    0,sqrt(G*M(1)/P(2,1)),0
%    sqrt(G*M(2)/(sqrt(sum((P(3,:)-P(2,:)).^2)))),sqrt(G*M(1)/P(3,1)) ,0
%    ];
%%We create the animatedlines
L = cell(N,1);
for j = 1:N
    L{j} = animatedline('LineStyle','none','MaximumNumPoint',tail,'MarkerFaceColor',C(j),'Marker','o','MarkerSize',5,'MarkerEdgeColor','none'); 
end
%%Text stuff
lgd = legend(string(M'),'Location','southeastoutside','EdgeColor',[1,1,1],'TextColor',[1,1,1]);
title(lgd,'Masses')
dimTE = [.83 .5 .4 .3];
TE = annotation('textbox',dimTE,'FitBoxToText','on','EdgeColor',[1,1,1],'Color',[1,1,1]');
dtHrs = sprintf('dt = %.1f hours',string(dt/3600));
E = 0;
tic
for k = 1:duration
    %%Animate
    for kk = 1:N
       addpoints(L{kk,1},P(kk,1),P(kk,2),P(kk,3));
    end
    if mod(k,speed)==0
        ROT = ROT+rotSpeed*.1; view(ROT,15);
        CAM = campos.*1.2;campos(CAM)
        drawnow
        tScale = sprintf('%.1f hours/second',dt*speed/(3600*toc));
        tElapsed = {sprintf('Time elapsed: %.1f days', string(k*dt/3600/24)),dtHrs,tScale};
        tic
        TE.String = tElapsed;
    end
    %%Potential Energy
    U = zeros(size(P,1),1);
    %%Velocity changes
    V2 = V;
    for kk = 1:N
        CM = 0;jjM = 0;
        for jj = 1:N
            if jj==kk
                continue
            else
                %%Resulting velocity of kk due to aceleration due to jj
                V(kk,:) = V(kk,:) + (G*M(jj)*(P(jj,:)-P(kk,:))/(sqrt(sum((P(jj,:)-P(kk,:)).^2))^3))*dt;
                U(kk,:) = U(kk,:) -G*M(jj)*M(kk)/sqrt(sum((P(jj,:)-P(kk,:)).^2));
                %CM = CM + M(jj)*P(jj,:);
                %jjM = jjM + M(jj);
            end
        end
        %U(kk) = -G*jjM*M(kk)/sqrt(sum(((CM/jjM)-P(kk,:)).^2));
    end
    %%Position changes
    for ii = 1:N
        P(ii,:) = P(ii,:) + V(ii,:)*dt;
    end
    LinearMomentum = sum(M.*V,'all')
    E = sum(U + (M.*(sum(V.^2,2)))./2)
end