clc; clear;

%% Correlation analysis
for i = 1:2
    if i == 1
        load('TckResult_EphOpensky_90.mat');
        ncase = 'Opensky';
    else
        load('TckResult_EphUrban_90.mat');
        ncase = 'Urban';
    end
    TckResultCT= TckResult_Eph;


    prn = 3;
    Spacing = -0.6:0.05:0.6;

    E = sqrt(TckResultCT(prn).E_i.^2 + TckResultCT(prn).E_q.^2);
    P = sqrt(TckResultCT(prn).P_i.^ 2 + TckResultCT(prn).P_q.^2);
    L = sqrt(TckResultCT(prn).L_i.^2 + TckResultCT(prn).L_q.^2);
    E_5 = sqrt(TckResultCT(prn).E_5_i.^2 + TckResultCT(prn).E_5_q.^2);
    L_5 = sqrt(TckResultCT(prn).L_5_i.^2 + TckResultCT(prn).L_5_q.^2);


    figure;
    plot([-0.5, -0.25, 0, 0.25, 0.5], [E(25000), E_5(25000), P(25000), L_5(25000), L(25000)], 'Color', "#0072BD", 'LineWidth', 1); hold on;
    plot([-0.5, -0.25, 0, 0.25, 0.5], [E(50000), E_5(25000), P(50000), L_5(50000), L(25000)], 'Color', "#D95319", 'LineWidth', 1); 
    plot([-0.5, -0.25, 0, 0.25, 0.5], [E(75000), E_5(25000), P(75000), L_5(75000), L(25000)], 'Color', "#EDB120", 'LineWidth', 1); 
    plot([-0.5, -0.25, 0, 0.25, 0.5], [E(90000), E_5(25000), P(90000), L_5(90000), L(25000)], 'Color', "#7E2F8E", 'LineWidth', 1); 
    xlabel('Time delay (chips)');
    ylabel('Correlation');
    title(['Correlation analysis for ', ncase, ' case']);
    grid on;

    legend('t = 25s', 't = 50s', 't = 75s', 't = 90s', 'location', 'best');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print(['Figures/correlation_analysis_plot', ncase], '-dpng', '-r600');
end





%% WLS plots
    OpenskyGT = [22.328444770087565, 114.1713630049711];
    UrbanGT = [22.3198722, 114.209101777778];
    load('navSolCT_1ms_Opensky.mat')
    OpenskySolutionCT = navSolutionsCT;
    load('navSolCT_1ms_Urban.mat')
    UrbanSolutionCT = navSolutionsCT;

    figure;
    plot(OpenskySolutionCT.localTime(1:5000), OpenskySolutionCT.usrPosLLH(1:5000,1), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Open sky latitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Open sky latitude', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.localTime(1:5000), UrbanSolutionCT.usrPosLLH(1:5000,1), 'Color', "#D95319", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Urban latitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban latitude', '-dpng', '-r600');

    figure;
    plot(OpenskySolutionCT.localTime(1:5000), OpenskySolutionCT.usrPosLLH(1:5000,2), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Opensky longitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Opensky longitude', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.localTime(1:5000), UrbanSolutionCT.usrPosLLH(1:5000,2), 'Color', "#D95319", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Urban longitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban longitude', '-dpng', '-r600');

    figure;
    plot(OpenskySolutionCT.localTime(1:5000), OpenskySolutionCT.usrVelENU(1:5000,3), 'Color', "#EDB120", 'LineWidth', 1); hold on;
    plot(OpenskySolutionCT.localTime(1:5000), OpenskySolutionCT.usrVelENU(1:5000,2), 'Color', "#D95319", 'LineWidth', 1);
    plot(OpenskySolutionCT.localTime(1:5000), OpenskySolutionCT.usrVelENU(1:5000,1), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('m/s');
    title('Opensky velocity');
    grid on;
    legend('U', 'N', 'E');
    h = legend;
    legend_strings = h.String;
    reversed_legend_strings = flip(legend_strings);
    legend(reversed_legend_strings, 'location', 'best');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Opensky velocity', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.localTime(1:5000), UrbanSolutionCT.usrVelENU(1:5000,3), 'Color', "#EDB120", 'LineWidth', 1);hold on;
    plot(UrbanSolutionCT.localTime(1:5000), UrbanSolutionCT.usrVelENU(1:5000,2), 'Color', "#D95319", 'LineWidth', 1);
    plot(UrbanSolutionCT.localTime(1:5000), UrbanSolutionCT.usrVelENU(1:5000,1), 'Color', "#0072BD", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('m/s');
    title('Urban velocity');
    grid on;
    legend('U', 'N', 'E');
    h = legend;
    legend_strings = h.String;
    reversed_legend_strings = flip(legend_strings);
    legend(reversed_legend_strings, 'location', 'best');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban velocity', '-dpng', '-r600');

    figure;
    geoscatter(UrbanSolutionCT.usrPosLLH(1:5000,1), UrbanSolutionCT.usrPosLLH(1:5000,2), 1,'g','filled'); hold on
    geoscatter(UrbanGT(1), UrbanGT(2), 5,'r','filled');
    geoscatter(OpenskySolutionCT.usrPosLLH(1:5000,1), OpenskySolutionCT.usrPosLLH(1:5000,2), 1,'c','filled');
    geoscatter(OpenskyGT(1), OpenskyGT(2), 5,'b','filled');
    title('Positioning Results');
    legend('Urban','Urban GT','Open sky','Open sky GT', 'Location','northwest');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Map', '-dpng', '-r600');





%% KF plots
    OpenskyGT = [22.328444770087565, 114.1713630049711];
    UrbanGT = [22.3198722, 114.209101777778];
    load('navSolCT_KF_Opensky.mat')
    OpenskySolutionCT = navSolutionsCT;
    load('navSolCT_KF_Urban.mat')
    UrbanSolutionCT = navSolutionsCT;

    figure;
    plot(OpenskySolutionCT.usrTime, OpenskySolutionCT.usrPosLLH(:,1), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Open sky latitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Open sky latitude_KF', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.usrTime, UrbanSolutionCT.usrPosLLH(:,1), 'Color', "#D95319", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Urban latitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban latitude_KF', '-dpng', '-r600');

    figure;
    plot(OpenskySolutionCT.usrTime, OpenskySolutionCT.usrPosLLH(:,2), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Opensky longitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Opensky longitude_KF', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.usrTime, UrbanSolutionCT.usrPosLLH(:,2), 'Color', "#D95319", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('^\circ');
    title('Urban longitude');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban longitude_KF', '-dpng', '-r600');

    figure;
    plot(OpenskySolutionCT.usrTime, OpenskySolutionCT.usrVelENU(:,3), 'Color', "#EDB120", 'LineWidth', 1); hold on;
    plot(OpenskySolutionCT.usrTime, OpenskySolutionCT.usrVelENU(:,2), 'Color', "#D95319", 'LineWidth', 1);
    plot(OpenskySolutionCT.usrTime, OpenskySolutionCT.usrVelENU(:,1), 'Color', "#0072BD", 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('m/s');
    title('Opensky velocity');
    grid on;
    legend('U', 'N', 'E');
    h = legend;
    legend_strings = h.String;
    reversed_legend_strings = flip(legend_strings);
    legend(reversed_legend_strings, 'location', 'best');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Opensky velocity_KF', '-dpng', '-r600');

    figure;
    plot(UrbanSolutionCT.usrTime, UrbanSolutionCT.usrVelENU(:,3), 'Color', "#EDB120", 'LineWidth', 1);hold on;
    plot(UrbanSolutionCT.usrTime, UrbanSolutionCT.usrVelENU(:,2), 'Color', "#D95319", 'LineWidth', 1);
    plot(UrbanSolutionCT.usrTime, UrbanSolutionCT.usrVelENU(:,1), 'Color', "#0072BD", 'LineWidth', 1); 
    xlabel('Time (ms)');
    ylabel('m/s');
    title('Urban velocity');
    grid on;
    legend('U', 'N', 'E');
    h = legend;
    legend_strings = h.String;
    reversed_legend_strings = flip(legend_strings);
    legend(reversed_legend_strings, 'location', 'best');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Urban velocity_KF', '-dpng', '-r600');

    figure;
    geoscatter(UrbanSolutionCT.usrPosLLH(:,1), UrbanSolutionCT.usrPosLLH(:,2), 1,'g','filled'); hold on
    geoscatter(UrbanGT(1), UrbanGT(2), 5,'r','filled');
    geoscatter(OpenskySolutionCT.usrPosLLH(:,1), OpenskySolutionCT.usrPosLLH(:,2), 1,'c','filled');
    geoscatter(OpenskyGT(1), OpenskyGT(2), 5,'b','filled');
    title('Positioning Results');
    legend('Urban','Urban GT','Open sky','Open sky GT', 'Location','northwest');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 3.25, 3.25*3/4]);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman'); 
    print('Figures/Map_KF', '-dpng', '-r600');
