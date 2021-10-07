close all;
do_rixs = 0;

B = dlmread('./Data/chi_Nk64_T170K.dat');
I = find(B(:,1) == 0); nw = length(I); w = B(I,2);
I = find(B(:,2) == 0); nq = length(I); q = B(I,1);
Sqw = reshape(B(:,4),nw,nq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the RIXS intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(do_rixs == 1)
    A = dlmread('./Data/rixs_OKedge_ARPESbands_Zeq1_Nk64_T170K.dat');
    
    Ibr = reshape(A(:,4),nw,nq);
    Ib1g = reshape(A(:,5),nw,nq);
    Ia1g = reshape(A(:,6),nw,nq);
    Iapex = reshape(A(:,7),nw,nq);
    Iac = reshape(A(:,8),nw,nq);
    
    set(0,'defaulttextInterpreter','latex')
    figure(1); hold on;
    imagesc(q/2,1000*w,6*Ibr+10*Ib1g+5*Ia1g+1*Iapex+Iac);
    A = load('24K_data_dispersion.txt');
    errorbar(A(:,1),1000*A(:,2),1000*A(:,3),'sk','MarkerFaceColor','k')
    errorbar(A(:,1),1000*A(:,4),1000*A(:,5),'db','MarkerFaceColor','b')
    errorbar(A(:,1),1000*A(:,6),1000*A(:,7),'vr','MarkerFaceColor','r')
    axis([0.1,0.4,0,100]);
    set(gca,'FontSize',25);
    xlabel('$q_\parallel~(\pi/2a)$','FontSize',25,'FontName','Times');
    ylabel('Energy (meV)','FontSize',25);
    caxis([0,2])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the phonon dispersions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = dlmread('ImDChi.dat');
Ibr = reshape(C(:,9),nw,nq);
Ib1g = reshape(C(:,7),nw,nq);
Ia1g = reshape(C(:,5),nw,nq);
Iapex = reshape(C(:,3),nw,nq);
Iac = reshape(C(:,11),nw,nq);

I = Ibr+Ib1g+Ia1g+Iapex+Iac;

fh = figure('Renderer', 'painters', 'Position', [10 10 1000 500]); 
subplot(1,2,1); hold on;
imagesc(q/2,1000*w,I); hold on;
axis([0,0.5,0,100]);
set(gca,'FontSize',25);
xlabel('$q_\parallel~(\pi/2a)$','FontSize',25,'FontName','Times');
ylabel('Energy (meV)','FontSize',25);
A = load('24K_data_dispersion.txt');
errorbar(A(:,1),1000*A(:,2),1000*A(:,3),'sk','MarkerFaceColor','k')
errorbar(A(:,1),1000*A(:,4),1000*A(:,5),'db','MarkerFaceColor','b')
errorbar(A(:,1),1000*A(:,6),1000*A(:,7),'vr','MarkerFaceColor','r')
title('$-\mathrm{Im}D(q,\omega)$')

subplot(1,2,2); hold on;
imagesc(q/2,1000*w,Sqw); hold on;
axis([0,0.5,0,100]);
set(gca,'FontSize',25);
xlabel('$q_\parallel~(\pi/2a)$','FontSize',25,'FontName','Times');
ylabel('Energy (meV)','FontSize',25);
title('$\mathrm{Im}\chi(q,\omega)$')

%saveas(fh,'phonons_Nk64_T170K.eps','epsc');
