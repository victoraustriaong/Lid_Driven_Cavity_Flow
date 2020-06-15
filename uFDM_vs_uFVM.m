%This code and the data included graphs the Lid-Driven Cavity Ansys Fluent u-velocity results (FVM) versus the
%results from FDM approach using Matlab.
%By: Victor Ong

load 'uFDM_vs_uFVM.mat'

plot(uu,Y)
hold on
plot(uuu,YY,'LineWidth',3)

title('FDM u-velocity vs FVM u-velocity ')
xlabel('U velocity (m/s)');
ylabel('Y (m)');



% load 'uFDM_vs_uFVM_RE1000.mat'
% plot(UU,YY)
% hold on
% plot(uu,yy,'LineWidth',3)
% yy=linspace (0,1,50);
% title('FDM u-velocity vs FVM u-velocity ')
% xlabel('U velocity (m/s)');
% ylabel('Y (m)');