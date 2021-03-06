clear all;
close all;
clc;
[S,F]=rdsamp('f:\Database\CUDB\cu02\cu02',[],(60*0+5)*250);
% load('f:\Database\CUDB\cu02\cu02m.mat');
B = 1/5*ones(5,1);
out = filter(B,1,F);                                                       %滑动均值滤波
% t=S((60*2+10)*250:end-1,1);
f=out((60*0+0)*250+1:end,1);
plot(f)
hold on
Hd=untitled();
output=filter(Hd,f);                                                       %滤波器使用
plot(output)
hold on
Hd=untitled1();
output1=filter(Hd,output);
plot(output1)
hold on
Hd=untitled2();
output2=filter(Hd,output1);
plot(output2)
% f=F((60*8+10)*250:end-1,1);
% p=F((60*2+10)*250:end-1,1);
% B = 1/5*ones(5,1);
% out = filter(B,1,f);
% plot(p)
% hold on
% plot(f)
% grid on
% legend('Input Data','Filtered Data','Location','NorthWest')              %5阶滑动均值滤波
% title('Plot of Input and Filtered Data')
% [c1,l1]=wavedec(w1,50,'db3');
% decmp1=wrcoef('a',c1,l1,'db3',7);                                        %提取趋势项
% corrcoef(a(1:len),b(1:len))                                              %比较两个波形相关性                                        

%--------------- Preparation初步建议不要使用各种滤波器造成模态频率偏移
alpha = 5000;        % moderate bandwidth constraint
tau = 1;            % noise-tolerance (no strict fidelity enforcement)
K = 4;              % 3 modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;
%--------------- Visualization

% Time Domain 0 to T
T = length(f);
fs = 1/T;
t = (1:T)/T;
freqs = 2*pi*(t-0.5-1/T)/(fs);


% for visualization purposes
% fsub = {};
% wsub = {};
% fsub{1} = v_1;
% fsub{2} = v_2;
% fsub{3} = v_3;
% wsub{1} = 2*pi*f_1;
% wsub{2} = 2*pi*f_2;
% wsub{3} = 2*pi*f_3;
%--------------- Run actual VMD code

[u, u_hat, omega] = VMD(f, alpha, tau, K, DC, init, tol);
f_hat = fftshift((fft(f)));

%--------------- Visualization

% For convenience here: Order omegas increasingly and reindex u/u_hat
[~, sortIndex] = sort(omega(end,:));
omega = omega(:,sortIndex);
u_hat = u_hat(:,sortIndex);
u = u(sortIndex,:);
linestyles = {'b', 'g', 'm', 'c', 'c', 'r', 'k'};

figure('Name', 'Composite input signal' );%输出合成波
plot(t,f, 'k');
% set(gca, 'XLim', [0 2]);  %截取横坐标0到2范围内的图形

% for sub = 1:length(fsub)
%     figure('Name', ['Input signal component ' num2str(sub)] );
%     plot(t,fsub{sub}, 'k');
%     set(gca, 'XLim', [0 1]);
% end

% figure('Name', 'Input signal spectrum' );
% loglog(freqs(T/2+1:end), abs(f_hat(T/2+1:end)), 'k');
% set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% ylims = get(gca, 'YLim');
% hold on;
% for sub = 1:length(wsub)
%     loglog([wsub{sub} wsub{sub}], ylims, 'k--');
% end
% set(gca, 'YLim', ylims);

figure('Name', 'Evolution of center frequencies omega');
for k=1:K
    semilogx(2*pi/fs*omega(:,k), 1:size(omega,1), linestyles{k});          %横坐标用10对数表示，纵坐标直接表示画图
    hold on;
end
set(gca, 'YLim', [1,size(omega,1)]);
set(gca, 'XLim', [2*pi,0.5*2*pi/fs], 'XGrid', 'on', 'XMinorGrid', 'on');
% corrcoef(a(1:len),b(1:len))相关系数求解

% figure('Name', 'Spectral decomposition');
% loglog(freqs(T/2+1:end), abs(f_hat(T/2+1:end)), 'k:');                     %横坐标和纵坐标都用10对数表示画图
% % semilogx(freqs(T/2+1:end), abs(f_hat(T/2+1:end)),'k:');
% set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% hold on;
% for k = 1:K
%     loglog(freqs(T/2+1:end), abs(u_hat(T/2+1:end,k)), linestyles{k});
% end
% set(gca, 'YLim', ylims);


for k = 1:K
    figure('Name', ['Reconstructed mode ' num2str(k)]);
    plot(t,u(k,:), linestyles{k});   hold on;
%     plot(t,u(2,:)+u(3,:)+u(4,:));
    plot(t,f);
%     corrcoef(u(k,:),u(2,:)+u(3,:)+u(4,:)) 
    corrcoef(f,u(k,:)) 
%     if ~isempty(fsub)
%         plot(t, fsub{min(k,length(fsub))}, 'k:');
%     end
    set(gca, 'XLim', [0 1]);
end

