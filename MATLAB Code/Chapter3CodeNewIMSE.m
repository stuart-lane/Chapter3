%%% STUART LANE
%%% 07/06/2023 
%%% PRODUCE TABLES AND GRAPHS FROM CHAPTER 3.

%% TAKES ABOUT 8 HOURS TO RUN ALL 

clear

tic;

rng(1029384);

nvals = [200 500];
nsimvals = [2000 2000];

fignum = 1;

printvals = zeros(1,10);

alphavals = [2 3];
betavals = [2 3];
deltavals = [0 0.5 1 2];

xgrid = [0.02:0.02:0.98];
nxgrid = length(xgrid);


for betaj = 1:length(betavals)
    beta = betavals(betaj);
    
    for alphaj = 1:length(alphavals)
        alpha = alphavals(alphaj);
    
        for deltaj = 1:length(deltavals)    
            deltamult = deltavals(deltaj);
            if alpha + 1/2 > beta
                deltabar = (2*beta-1)/(2*alpha+2*beta+1);
                delta = deltamult*deltabar;
            else
                deltabar = (alpha)/(2*alpha+1);
                delta = deltamult*deltabar;
            end

            for nj = 1:length(nvals)
                n = nvals(nj);
                nsim = nsimvals(nj); 

                J = 100;
                Jvec = [1:1:J];
                Jvalsall = [2:1:J];
                Jvalsodd = [1:2:J];

                warning('off')
                dfn = length(alphavals)*length(betavals)*length(deltavals)+1;
                fxzg = @(x,z) fun_fxz(J,x,z,alpha,n,delta);
                figure(dfn);
                S = fsurf(fxzg,[0 1 0 1]);
                ylabel("$z$",'Interpreter','Latex');
                xlabel("$x$",'Interpreter','Latex');
                zlabel("Density");
                zlim([0 5]);
                xticks([0:0.2:1]);
                yticks([0:0.2:1]);
                hold on;

                zdata = S.ZData ; 
                maxfxz = max(zdata(:));

                sb1 = zeros(nsim,nxgrid);
                sb2 = zeros(nsim,nxgrid);
                sb3 = zeros(nsim,nxgrid);
                sb4 = zeros(nsim,nxgrid);

                for si=1:nsim
                    jj = 1;
                    x = zeros(n,1);
                    z = zeros(n,1);

                    while jj < n+1
                        target = rand(2,1);
                        zval = maxfxz*rand;
                        val = fun_fxz(J,target(1),target(2),alpha,n,delta);
                        istrue = val > zval;
                        if istrue == 0
                           x(jj) = target(1);
                           z(jj) = target(2);
                        else 
                           x(jj) = target(1);
                           z(jj) = target(2);
                           jj = jj+1;
                        end
                    end

                    data = [x z];
                    data = sortrows(data,1);
                    x = data(:,1);
                    z = data(:,2);
                    
                    %%% endogenous errors
                    zregJ = 10;
                    zreg = zeros(n,zregJ);
                    for kk = 1:zregJ
                        zreg(:,kk) = sqrt(2)*cos((kk-1)*pi*z);
                    end
                    u = x - zreg*(zreg\x) + normrnd(0,0.3,n,1);                 
                    
                    phix = fun_phix(J,x,beta);
                    y = phix + u;
                    
                    eigs = zeros(J,1);
                    eigs(1) = 1;
                    for j = 2:J
                        eigs(j) = 0.2/(((j-1)^(alpha))*(n^(delta)));
                    end                    
                    
                    rjvec = zeros(J,1);
                    rjvec(1) = mean(y);
                    for jj = 2:J
                        rjvec(jj) = sqrt(2)*sqrt(eigs(jj))*(y'*cos((jj-1)*pi*z)/n);
                    end
                    
                    estcoeffs = rjvec./eigs;
                    
                    phihatsumvals = zeros(n,J);
                    phihatsumvals(:,1) = repmat(estcoeffs(1),n,1);
                    for kk = 2:J
                        phihatsumvals(:,kk) = estcoeffs(kk)*sqrt(2)*cos((kk-1)*pi*x);
                    end  

                    phiJ = 10;
                    phihatscmsevec = zeros(phiJ,1);
                    for kk=1:phiJ   
                        phihatscsumk = phihatsumvals(:,1:kk);
                        phihatsck = sum(phihatscsumk,2);
                        errorphihatsc1 = (phix-phihatsck).^2;
                        phihatscmsevec(kk) = sum(errorphihatsc1)/n;
                    end
                   [msescmin,I1sc] = min(phihatscmsevec);
                   phihatscsum = phihatsumvals(:,1:I1sc);
                   phihatsc = sum(phihatscsum,2);
                   

                   %%% TIKHONOV ESTIMATOR
                    anvec = logspace(-3,0,200)'; 
                    estcoeffstk = zeros(length(anvec),J);
                    estcoeffstk(:,1) = repmat(rjvec(1),length(anvec),1)./(1+anvec);
                    for jj = 2:J
                        estcoeffstk(:,jj) = repmat(rjvec(jj),length(anvec),1)./(eigs(jj)+anvec);
                    end 
                    
                    basisfunvals = zeros(n,J);
                    basisfunvals(:,1) = 1;
                    for jj = 2:J
                        basisfunvals(:,jj) = sqrt(2)*cos((jj-1)*pi*x);
                    end
                 
                    phihattkmsevec = zeros(length(anvec),1);
                    for kk=1:length(anvec)
                        phihattk = basisfunvals*estcoeffstk(kk,:)';
                        errorphihattk1 = (phix-phihattk).^2;
                        phihattkmsevec(kk) = sum(errorphihattk1)/n;
                    end
                   [msetkmin,I1tk] = min(phihattkmsevec);
                   phihattk = basisfunvals*estcoeffstk(I1tk,:)';
                   show = phihattkmsevec;
                   
                   % IMSE vectors
                   
                   % SPECTRAL CUTOFF

                    phihatsumvals = zeros(nxgrid,I1sc);
                    phihatsumvals(:,1) = repmat(estcoeffs(1),nxgrid,1);
                    for kk = 2:J
                        phihatsumvals(:,kk) = estcoeffs(kk)*sqrt(2)*cos((kk-1)*pi*xgrid);
                    end  

                   phihatscsum = phihatsumvals(:,1:I1sc);
                   phihatsc = sum(phihatscsum,2);
                    
                    basisfunvals = zeros(nxgrid,J);
                    basisfunvals(:,1) = 1;
                    for jj = 2:J
                        basisfunvals(:,jj) = sqrt(2)*cos((jj-1)*pi*xgrid);
                    end
                 
                   phihattk = basisfunvals*estcoeffstk(I1tk,:)';
                                      
                   sb3(si,:) = phihatsc';
                   sb4(si,:) = phihattk';

                end
                
            phigrid = fun_phix(J,xgrid,beta);      
            tkvals = sb4;
            scvals = sb3;
            
            isqbiastk = mean((mean(tkvals) - phigrid).^2);
            isqbiassc = mean((mean(scvals) - phigrid).^2);
            
            ivartk = mean(var(tkvals));
            ivarsc = mean(var(scvals));
            
            imsetk = isqbiastk + ivartk;
            imsesc = isqbiassc + ivarsc;

            printvals(end+1,:) = [n/100 deltamult alpha beta isqbiassc ivarsc imsesc isqbiastk ivartk imsetk]

            meanx = mean(sb1);
            meanphix = mean(sb2);
            meanphihatsca = mean(sb3);
            meanphihattka = mean(sb4);
            
            l_phihatsca = quantile(sb3,0.025);
            u_phihatsca = quantile(sb3,0.975);
            l_phihattka = quantile(sb4,0.025);
            u_phihattka = quantile(sb4,0.975);
            
            n1 = nvals(1);
            n2 = nvals(2);
%            n3 = nvals(3);

%             if n == n1;
%                 
%                 figure(fignum);
%                 p1 = plot(meanx,meanphix,'k','LineWidth',2);
%                 p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p2 = plot(meanx,meanphihatsca,'r','LineWidth',1.5);
%                 hold on
%                 p3 = plot(meanx,u_phihatsca,'r--','LineWidth',0.5);
%                 p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p4 = plot(meanx,l_phihatsca,'r--','LineWidth',0.5);
%                 p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p5 = plot(meanx,meanphihattka,'b','LineWidth',1.5);
%                 hold on
%                 p6 = plot(meanx,u_phihattka,'b--','LineWidth',0.5);
%                 p6.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p7 = plot(meanx,l_phihattka,'b--','LineWidth',0.5);
%                 p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 xlim([0 1]);
%                 ylim([-0.2 1.8]);
%                 ylabel("$\varphi(x)$",'Interpreter','Latex');
%                 xlabel("$x$",'Interpreter','Latex');
%                 hold on;
%                 legend('$(n = 200)$','Location','northwest','Interpreter','Latex');
% 
%             elseif n == n2;
% 
%                 figure(fignum);
%                 p8 = plot(meanx,meanphihatsca,'b','LineWidth',1.5);
%                 hold on
%                 p9 = plot(meanx,u_phihatsca,'b--','LineWidth',0.5);
%                 p9.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p10 = plot(meanx,l_phihatsca,'b--','LineWidth',0.5);
%                 p10.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p11 = plot(meanx,meanphihattka,'b','LineWidth',1.5);
%                 p11.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p12 = plot(meanx,u_phihattka,'b--','LineWidth',0.5);
%                 p12.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p13 = plot(meanx,l_phihattka,'b--','LineWidth',0.5);
%                 p13.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 xlim([0 1]);
%                 ylim([-0.2 1.8]);
%                 ylabel("$\varphi(x)$",'Interpreter','Latex');
%                 xlabel("$x$",'Interpreter','Latex');
%                 hold on;
%                 legend('$(n = 200)$','$(n = 500)$','Location','northwest','Interpreter','Latex');
% 
%             else
%  
%                 figure(fignum);
%                 p14 = plot(meanx,meanphihatsca,'Color',[0.09,0.77,0.09],'LineWidth',1.5);
%                 hold on
%                 p15 = plot(meanx,u_phihatsca,'Color',[0.09,0.77,0.09],'LineStyle','--','LineWidth',0.5);
%                 p15.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p16 = plot(meanx,l_phihatsca,'Color',[0.09,0.77,0.09],'LineStyle','--','LineWidth',0.5);
%                 p16.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p17 = plot(meanx,meanphihattka,'b','LineWidth',1.5);
%                 p17.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p18 = plot(meanx,u_phihattka,'b--','LineWidth',0.5);
%                 p18.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 hold on
%                 p19 = plot(meanx,l_phihattka,'b--','LineWidth',0.5);
%                 p19.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                 xlim([0 1]);
%                 ylim([-0.2 1.8]);
%                 ylabel("$\varphi(x)$",'Interpreter','Latex');
%                 xlabel("$x$",'Interpreter','Latex');
%                 hold on;
%                 legend('$(n = 200)$','$(n = 500)$','$(n = 100)$','Location','northwest','Interpreter','Latex');
% 
%             end
% % 
% %             if n == n1;
% %                 
% %                 figure(fignum);
% %                 p1 = plot(meanx,meanphix,'k','LineWidth',2);
% %                 p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p2 = plot(meanx,meanphihatsca,'r','LineWidth',1.5);
% %                 hold on
% %                 p3 = plot(meanx,u_phihatsca,'r--','LineWidth',0.5);
% %                 p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p4 = plot(meanx,l_phihatsca,'r--','LineWidth',0.5);
% %                 p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p5 = plot(meanx,meanphihattka,'b','LineWidth',1.5);
% %                 p5.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p6 = plot(meanx,u_phihattka,'b--','LineWidth',0.5);
% %                 p6.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p7 = plot(meanx,l_phihattka,'b--','LineWidth',0.5);
% %                 p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 xlim([0 1]);
% %                 ylim([-0.2 1.8]);
% %                 ylabel("$\varphi(x)$",'Interpreter','Latex');
% %                 xlabel("$x$",'Interpreter','Latex');
% %                 hold on;
% %                 legend('$(n = 200)$','Location','northwest','Interpreter','Latex');
% % 
% %             elseif n == n2;
% % 
% %                 figure(fignum);
% %                 p8 = plot(meanx,meanphihatsca,'r-+','LineWidth',1.5,'MarkerIndices',50:50:450);
% %                 hold on
% %                 p9 = plot(meanx,u_phihatsca,'r--+','LineWidth',0.5,'MarkerIndices',50:50:450);
% %                 p9.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p10 = plot(meanx,l_phihatsca,'r--+','LineWidth',0.5,'MarkerIndices',50:50:450);
% %                 p10.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p11 = plot(meanx,meanphihattka,'b--+','LineWidth',1.5,'MarkerIndices',50:50:450);
% %                 p11.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p12 = plot(meanx,u_phihattka,'b--+','LineWidth',0.5,'MarkerIndices',50:50:450);
% %                 p12.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p13 = plot(meanx,l_phihattka,'b--+','LineWidth',0.5,'MarkerIndices',50:50:450);
% %                 p13.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 xlim([0 1]);
% %                 ylim([-0.2 1.8]);
% %                 ylabel("$\varphi(x)$",'Interpreter','Latex');
% %                 xlabel("$x$",'Interpreter','Latex');
% %                 hold on;
% %                 legend('$(n = 200)$','$(n = 500)$','Location','northwest','Interpreter','Latex');
% % 
% %             else
% %  
% %                 figure(fignum);
% %                 p14 = plot(meanx,meanphihatsca,'r-o','LineWidth',1.5,'MarkerIndices',100:100:900);
% %                 hold on
% %                 p15 = plot(meanx,u_phihatsca,'r--o','LineWidth',0.5,'MarkerIndices',100:100:900);
% %                 p15.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p16 = plot(meanx,l_phihatsca,'r--o','LineWidth',0.5,'MarkerIndices',100:100:900);
% %                 p16.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p17 = plot(meanx,meanphihattka,'b-o','LineWidth',1.5,'MarkerIndices',100:100:900);
% %                 p17.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p18 = plot(meanx,u_phihattka,'b--o','LineWidth',0.5,'MarkerIndices',100:100:900);
% %                 p18.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 hold on
% %                 p19 = plot(meanx,l_phihattka,'b--o','LineWidth',0.5,'MarkerIndices',100:100:900);
% %                 p19.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                 xlim([0 1]);
% %                 ylim([-0.2 1.8]);
% %                 ylabel("$\varphi(x)$",'Interpreter','Latex');
% %                 xlabel("$x$",'Interpreter','Latex');
% %                 hold on;
% %                 legend('$(n = 200)$','$(n = 500)$','$(n = 1000)$','Location','northwest','Interpreter','Latex');
% % 
% %             end
% % 
% %             if n == nvals(end);
% %                 fignum = fignum + 1;
% %             else
% %                 fignum = fignum;
% %             end
% % 
             end

        end

    end

end

toc

header = {'n','delta','alpha','beta','ISQ','IV','IMSE'}

FINALVALUES = printvals(2:end,:)

% unblock line below to save table of results
csvwrite('finalsimulationsvalues.csv',FINALVALUES,0,1)

T1 = array2table(FINALVALUES,'VariableNames',{'n','delta','alpha','beta','ISQBS','IVARS','IMSES','ISQBT','IVART','IMSET'})
table2latex(T1)

toc

alpha = 3;
beta = 4;
delta = 0;
anrate = -(alpha + 2*delta*(alpha+beta))/(alpha+2*beta)

data = [x y z u phix]
file_name = 'my_data.csv';

% Use csvwrite to save the matrix to a CSV file
csvwrite(file_name, data);
