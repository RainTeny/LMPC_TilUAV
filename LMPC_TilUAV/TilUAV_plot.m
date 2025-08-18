%% 画图
% 修改saveas函数后的路径用于保存图片和路径
p=0:N_sim-1;
% 输出与参考值
figure
for i=1:4 
    subplot(4,1,i)
    plot(p,y_vector(i,:),'LineWidth', 1);
    hold on
    grid on
    plot(p,r(i,:));
end
subplot(4,1,1); xlabel('x (m)');legend('x','xd');
subplot(4,1,2); xlabel('y (m)');legend('y','yd');
subplot(4,1,3); xlabel('z (m)');legend('z','zd');
subplot(4,1,4); xlabel('yaw (rad)');legend('\psi','\psi d');
% saveas(gcf,['C:\Users\LJY\Desktop\confer\try3\fig1_output.png'])
%  alpha1 与 alpha2 曲线
% figure
% subplot(2,1,1)
% plot(p, alpha1_vector, 'LineWidth', 1.5)
% xlabel('Time step')
% ylabel('\alpha_1 (rad)')
% title('Tilt Angle \alpha_1 vs Time')
% grid on
% 
% subplot(2,1,2)
% plot(p, alpha2_vector, 'LineWidth', 1.5)
% xlabel('Time step')
% ylabel('\alpha_2 (rad)')
% title('Tilt Angle \alpha_2 vs Time')
% grid on

% u
figure
for i=1:6
    subplot(6,1,i)
    plot(p,u_vector(i,:),'LineWidth', 1);
    xlabel('u'+string(i));
    grid on
end
% saveas(gcf,['C:\Users\LJY\Desktop\confer\try3\fig2_u.jpg'])
% 
% % deltau
% figure
% for i=1:4
%     subplot(4,1,i)
%     plot(p,deltau_vector(i,:),'LineWidth', 1);
%     xlabel('deltau'+string(i));
%     grid on
% end
% saveas(gcf,['C:\Users\LJY\Desktop\confer\try3\fig3_deltau.jpg'])

% %状态量 x  画图
%    figure  
%     plot(p,xm_vector(1,:),'r','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(m/s)')
%     legend('x')
%     
%     figure
%    plot(p,xm_vector(2,:),'g','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(m/s)')
%     legend('y')
%     
%      figure
%     plot(p,xm_vector(3,:),'r','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(m/s)')
%     legend('z')
%     
%      figure
%   plot(p,xm_vector(4,:),'b','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('m')
%     legend('v_x')
%     
%      figure
%     plot(p,xm_vector(5,:),'c','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(rad/s)')
%     legend('v_y')
%     
%      figure
%     plot(p,xm_vector(6,:),'m','linewidth',1.0);
%   xlabel('Time(s)')
%     ylabel('rad/s')
%     legend('v_z')
%     
     figure
    plot(p,xm_vector(7,:),'b','linewidth',1.0);
    xlabel('Time(s)')
    ylabel('(m)')
    legend('\psi')
    
    
     figure
   plot(p,xm_vector(8,:),'r','linewidth',1.0);
    xlabel('Time(s)')
    ylabel('(rad/s)')
   legend('\phi')
    
     figure
   plot(p,xm_vector(9,:),'g','linewidth',1.0);
    xlabel('Time(s)')
    ylabel('(rad/s)')
    legend('\Theta')
%      
%     
%      figure
%    plot(p,xm_vector(10,:),'b','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(rad/s^2)')
%     legend('dot\psi')
%     
%      figure
%   plot(p,xm_vector(11,:),'c','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(rad/s^2)')
%     legend('dot\phi')
%     
%      figure
%     plot(p,xm_vector(12,:),'m','linewidth',1.0);
%     xlabel('Time(s)')
%     ylabel('(rad/s^2)')
%     legend('dot\Theta')
%  
% figure
% for i=1:12
%     subplot(4,3,i)
%     plot(p,xm_vector(i,:),'LineWidth', 1); 
%     hold on
%     grid on
%     xlabel('x'+string(i));    
% end
% % saveas(gcf,['C:\Users\LJY\Desktop\confer\try3\fig4_state12.jpg'])

%%图
% figure 
% plot3(y_vector(1,:),y_vector(2,:),y_vector(4,:),'b','LineWidth', 1.5)
% hold on
% plot3(r(1,:),r(2,:),r(4,:),'r--','LineWidth', 1.0)   
% xlabel('Position x(m)')
% ylabel('Position y(m)')
% zlabel('Position z(m)')
% legend('实际轨迹','参考轨迹');


 %%三维图 
figure
plot3(y_vector(1,:),y_vector(2,:),y_vector(3,:),'b','LineWidth', 1.0)
hold on
plot3(r(1,:),r(2,:),r(3,:),'r:','LineWidth', 1.5)
xlabel('Position x(m)');ylabel('Position y(m)');zlabel('Position z(m)');
xlabel('Position x(m)')
ylabel('Position y(m)')
zlabel('Position z(m)')
legend('实际轨迹','参考轨迹');
% saveas(gcf,['C:\Users\LJY\Desktop\confer\try3\fig5_tracking.jpg'])

errorsum=sum(sqrt((y_vector(1,:)-r(1,:)).^2)+sqrt((y_vector(2,:)-r(2,:)).^2)+sqrt((y_vector(3,:)-r(3,:)).^2)+sqrt((y_vector(4,:)-r(4,:)).^2))
% save(['C:\Users\LJY\Desktop\confer\try3\matlab.mat'])



%% Nc=10; Np=50; 1.0213e+03
%% Nc=20; Np=50;  938.4251