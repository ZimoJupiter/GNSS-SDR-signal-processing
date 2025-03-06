function TrackingPlot(TckResult,prn)

figure(prn);

subplot(4,3,1);
plot(TckResult.P_i(:),TckResult.P_q(:),'b*');
title('Discrete-Time Scatter Plot');
xlabel('I prompt');
ylabel('Q prompt');
grid on;
subplot(4,3,2:3);
plot(TckResult.P_i(:),'r-');
hold on;
plot(TckResult.P_q(:),'b-');
title('Ip and Qp of Tracking Result');
legend('I_P','Q_P');
grid on;

subplot(4,3,4);
plot(TckResult.PLLdiscri(:),'b-');
title('Raw PLL Discriminator');
xlabel('Time(ms)');
ylabel('Amplitude');
grid on;

subplot(4,3,5:6)
plot(sqrt(TckResult.E_i(:).^2+TckResult.E_q(:).^2),'b*-');hold on
plot(sqrt(TckResult.L_i(:).^2+TckResult.L_q(:).^2),'g*-');hold on
plot(sqrt(TckResult.P_i(:).^2+TckResult.P_q(:).^2),'r*-');hold on
legend('Early','Late','Prompt');
title('Correlation');
xlabel('Time(ms)');
grid on;

subplot(4,3,7);
plot(TckResult.doppler,'b');
title('Doppler Shift');
xlabel('Time(ms)');
ylabel('Hz');
grid on;

subplot(4,3,8);
plot(TckResult.DLLdiscri,'r-');
title('Raw DLL Discriminator');
xlabel('Time(ms)');
ylabel('Amplitude');
grid on;

subplot(4,3,9);
plot(TckResult.codedelay,'r*-');
title('Code Delay Sample');
xlabel('Time(ms)');
ylabel('Sample');
grid on;

RawNavigationData = TckResult(1).P_i(:);
NaviData(find(RawNavigationData>=0)) = 1;
NaviData(find(RawNavigationData<0)) = -1;
% NaviData = kron(NaviData,ones(1,4));
datalength = length(TckResult(1).P_i);
subplot(4,3,10:12);
stairs(NaviData);
axis([1 datalength -1.2 1.2])
title('Navigation Data');
xlabel('Time(ms)');
grid on;

% figure(100+prn);
% subplot(311)
% plot(TckResult.doppler,'b');
% title('Doppler Shift');
% ylabel('Hz');
% grid on;
% subplot(312)
% plot(TckResult.doppler2,'g');
% ylabel('Hz');
% grid on;
% subplot(313)
% plot(TckResult.doppler3,'r');
% xlabel('Time(ms)');
% ylabel('Hz');
% grid on;

end