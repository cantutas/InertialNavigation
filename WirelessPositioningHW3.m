clc; clear; close all;

% 参数设置
fs = 25e6;          % 采样频率 25 MHz
fIF = 6.25e6;       % 中频频率 6.25 MHz
f_clk = 1.023e6;    % 码片速率 1.023 MHz
codeLength = 1023;  % C/A 码长度
fd_step = -10e3:500:10e3;    % 多普勒频移搜索范围 (-10 kHz 到 +10 kHz，步长 500 Hz)
code_phase_step = 0:1:1022;  % 码相位搜索范围 (0 到 1022 码片)

% 加载实际的 GPS 信号数据
filepath = "D:\python\AI_hw\dsp\2. 附件1：无线定位大作业题目3-数据.xlsx";
dataTable = readtable(filepath);
gpsSignal = dataTable{:, 1}';
t=(0:length(gpsSignal)-1)/fs;

% 生成本地 PRN8 伪码
localCode = CAcode_generation([2 9],f_clk,[t(1),t(end)],fs);

% 初始化捕获结果矩阵
capture_result = zeros(length(fd_step), length(code_phase_step));

% 串行捕获算法
for i1 = 1:length(fd_step)
    fd = fd_step(i1);    % 当前多普勒频移
    wave_carrier_receive_I = cos(2 * pi * (fIF + fd) * t);
    wave_carrier_receive_Q = sin(2 * pi * (fIF + fd) * t);
    
    for i2 = 1:length(code_phase_step)
        code_phase = code_phase_step(i2);    % 当前码相位
        shiftedCode = circshift(localCode, round(-(1+(code_phase_step(i2)/f_clk-t(1))*fs)));
        
        % 相关运算
        R_I = gpsSignal .* wave_carrier_receive_I .* shiftedCode;
        R_Q = gpsSignal .* wave_carrier_receive_Q .* shiftedCode;
        
        % 积分并计算相关值
        NB = 1;    % 判决视域点数
        N = round(length(R_I) / NB);    % 相干积分点数
        R_I_sum = zeros(1, NB);
        R_Q_sum = zeros(1, NB);
        
        for i3 = 1:NB
            index_1 = (i3 - 1) * N + 1;
            index_2 = i3 * N;
            R_I_sum(i3) = sum(R_I(index_1:index_2));
            R_Q_sum(i3) = sum(R_Q(index_1:index_2));
        end
        
        M_I = sum(R_I_sum.^2);
        M_Q = sum(R_Q_sum.^2);
        capture_result(i1, i2) = (M_I + M_Q) / N;
    end
end

% 绘制捕获结果
[FD, PHASE_code] = meshgrid(fd_step, code_phase_step);
figure;
surf(FD, PHASE_code, capture_result'); shading interp;
xlabel("多普勒频移 (Hz)");
ylabel("码相位 (码片)");
zlabel("相关值");
title("GPS 信号捕获结果 (PRN8)");
colorbar;

% 找到最大相关值及其对应的多普勒频移和码相位
[max_val, max_idx] = max(capture_result(:));
[max_fd_idx, max_code_phase_idx] = ind2sub(size(capture_result), max_idx);
estimated_fd = fd_step(max_fd_idx);
estimated_code_phase = code_phase_step(max_code_phase_idx);

fprintf('Estimated Doppler Frequency: %.2f Hz\n', estimated_fd);
fprintf('Estimated Code Phase: %d chips\n', estimated_code_phase);

function [wave,code]=CAcode_generation(phase_tag,F_CLK,T,FS)
%phase_tag为G2码相位选择参数，F_CLK为码片速率设置
% T(时间范围)，FS为采样频率)
% 实际输出的伪码波形对应的时间横坐标为：t=[T(1):(1/fs):(T(2)-1/fs)];
%输出：wave(伪码时域波形)，code(伪码二进制序列)
    reg_G1=-ones(1,10);    %-1代表1，1代表0
    reg_G2=-ones(1,10);    
    index_code=0;
    t_chip=1/F_CLK;
    wave=zeros(1,(T(2)+1/FS)*FS);
    code=zeros(1,(T(2)+1/FS)*F_CLK);
    for i0=1:(T(2)+1/FS)*FS
        if ((T(1)+(i0-1)/FS)>=T(1)+index_code*t_chip)
            sum_G1=reg_G1(3)*reg_G1(10);
            reg_G1(2:10)=reg_G1(1:9);
            reg_G1(1)=sum_G1;
            sum_G2=reg_G2(2)*reg_G2(3)*reg_G2(6)*reg_G2(8)*reg_G2(9)*reg_G2(10);
            reg_G2(2:10)=reg_G2(1:9);
            reg_G2(1)=sum_G2;
            index_code=index_code+1;
            code(index_code)=reg_G1(10)*reg_G2(phase_tag(1))*reg_G2(phase_tag(2));
        end
        wave(i0)=code(index_code);
    end
end
