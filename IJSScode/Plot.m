function Plot(cellnumber, Micsig, inc, Maceps_save,...
    Macsig_save, T_save,  Macporosity, Maceps_v,...
    nonconnect_porosity, cellnumberA, connect_porosity)

Micsig_distri(1:6) = 0;
for i = 1 : cellnumber
    micsig_m = Micsig(:,:,i);
    [V, Prin] = eig(micsig_m);
    eigenvalue = [Prin(1,1) Prin(2,2) Prin(3,3)];
    prinsig = sort(eigenvalue);
    if prinsig(1) == prinsig(2)
        majorsigone = prinsig(3);
        mediansigone = prinsig(2) + 1e-10;
        minorsigone = prinsig(1);
    elseif prinsig(2) == prinsig(3)
        majorsigone = prinsig(3) + 1e-10;
        mediansigone = prinsig(2) ;
        minorsigone = prinsig(1);
    else
        majorsigone = prinsig(3);
        mediansigone = prinsig(2) ;
        minorsigone = prinsig(1);
    end
    for k = 1:3
        if eigenvalue(1,k) == minorsigone
            Vabs = abs(V(:,k));
            r = sqrt(Vabs(2)^2+Vabs(3)^2);
            z = Vabs(1);
            sinw_minor = r;
            cosw_minor = z;
        elseif eigenvalue(1,k) == majorsigone
            Vabs = abs(V(:,k));
            r = sqrt(Vabs(2)^2+Vabs(3)^2);
            z = Vabs(1);
            sinw_major = r;
            cosw_major = z;
        elseif eigenvalue(1,k) == mediansigone
            Vabs = abs(V(:,k));
            r = sqrt(Vabs(2)^2+Vabs(3)^2);
            z = Vabs(1);
            sinw_median = r;
            cosw_median = z;
        end
    end
    
    if majorsigone >=0
        majorsigr = abs(majorsigone*sinw_major);
        majorsigz = abs(majorsigone*cosw_major);
    else
        majorsigr = -abs(majorsigone*sinw_major);
        majorsigz = -abs(majorsigone*cosw_major);
    end 
    if mediansigone >=0
        mediansigr = abs(mediansigone*sinw_median);
        mediansigz = abs(mediansigone*cosw_median);
    else
        mediansigr = -abs(mediansigone*sinw_median);
        mediansigz = -abs(mediansigone*cosw_median);
    end
    if minorsigone >=0
        minorsigr = abs(minorsigone*sinw_minor);
        minorsigz = abs(minorsigone*cosw_minor);
    else
        minorsigr = -abs(minorsigone*sinw_minor);
        minorsigz = -abs(minorsigone*cosw_minor);
    end
    
    Micsig_distri(i,:) = [majorsigr majorsigz mediansigr mediansigz minorsigr minorsigz];
    
end
figure(1)
plot(Micsig_distri(:, 1),Micsig_distri(:, 2),'xb');hold on;
plot(Micsig_distri(:, 3),Micsig_distri(:, 4),'or');hold on;
plot(Micsig_distri(:, 5),Micsig_distri(:, 6),'+g');hold on;
x = max([max(abs(Micsig_distri(:, 1))); max(abs(Micsig_distri(:, 3))); max(abs(Micsig_distri(:, 5)))]);
y = max([max(abs(Micsig_distri(:, 2))); max(abs(Micsig_distri(:, 4))); max(abs(Micsig_distri(:, 6)))]);
x1 = [-3*x 3*x]; y1 = [0 0]; plot(x1,y1,'--k','lineWidth',1);
x2 = [0 0]; y2 = [-y y]; plot(x2,y2,'--k','lineWidth',2);
set(gca,'fontsize',12);
ylim([-y, y])
legend('Major Principal Micro-stress','Median Principal Micro-stress','Minor Principal Micro-stress','Location','Best','Orientation','Vertical');
axis equal;

saveas(gcf, 'result1.jpg');


figure(2)
set(gcf, 'Units', 'centimeter', 'Position', [5 10 30 10] )
subplot(1, 2, 1)
if cellnumberA > 0
plot(Micsig_distri(1:cellnumberA, 1),Micsig_distri(1:cellnumberA, 2),'xb');hold on;
plot(Micsig_distri(1:cellnumberA, 3),Micsig_distri(1:cellnumberA, 4),'or');hold on;
plot(Micsig_distri(1:cellnumberA, 5),Micsig_distri(1:cellnumberA, 6),'+g');hold on;
x = max([max(abs(Micsig_distri(:, 1))); max(abs(Micsig_distri(:, 3))); max(abs(Micsig_distri(:, 5)))]);
y = max([max(abs(Micsig_distri(:, 2))); max(abs(Micsig_distri(:, 4))); max(abs(Micsig_distri(:, 6)))]);
x1 = [-3*x 3*x]; y1 = [0 0]; plot(x1,y1,'--k','lineWidth',1);
x2 = [0 0]; y2 = [-y y]; plot(x2,y2,'--k','lineWidth',2);
set(gca,'fontsize',12);
ylim([-y, y])
legend('Major Principal Micro-stress','Median Principal Micro-stress','Minor Principal Micro-stress','Location','Best','Orientation','Vertical');
axis equal;
end

subplot(1, 2, 2)
if cellnumberA ~= cellnumber
plot(Micsig_distri((cellnumberA + 1):cellnumber, 1),Micsig_distri((cellnumberA + 1):cellnumber, 2),'xb');hold on;
plot(Micsig_distri((cellnumberA + 1):cellnumber, 3),Micsig_distri((cellnumberA + 1):cellnumber, 4),'or');hold on;
plot(Micsig_distri((cellnumberA + 1):cellnumber, 5),Micsig_distri((cellnumberA + 1):cellnumber, 6),'+g');hold on;
x = max([max(abs(Micsig_distri(:, 1))); max(abs(Micsig_distri(:, 3))); max(abs(Micsig_distri(:, 5)))]);
y = max([max(abs(Micsig_distri(:, 2))); max(abs(Micsig_distri(:, 4))); max(abs(Micsig_distri(:, 6)))]);
x1 = [-3*x 3*x]; y1 = [0 0]; plot(x1,y1,'--k','lineWidth',1);
x2 = [0 0]; y2 = [-y y]; plot(x2,y2,'--k','lineWidth',2);
set(gca,'fontsize',12);
ylim([-y, y])
legend('Major Principal Micro-stress','Median Principal Micro-stress','Minor Principal Micro-stress','Location','Best','Orientation','Vertical');
axis equal;
end
saveas(gcf, 'result2.jpg');  


sz = inc-1;
figure(3)
set(gcf, 'Units', 'centimeter', 'Position', [5 10 30 10]    )
subplot(1, 2, 1)
plot(reshape(Maceps_save(1, 1, :), sz, 1), reshape(Macsig_save(1, 1, :), sz, 1),'-b','LineWidth',3);
hold on;
plot(reshape(Maceps_save(2, 2, :), sz, 1), reshape(Macsig_save(1, 1, :), sz, 1), '--r','LineWidth',4);
hold on;
plot(reshape(Maceps_save(3, 3, :), sz, 1), reshape(Macsig_save(1, 1, :), sz, 1),'-k','LineWidth',1.5);
legend('\epsilon_{11}','\epsilon_{22}','\epsilon_{33}','Location','Best','Orientation','Horizontal');
xlabel('Strain, \epsilon','FontSize',12);
ylabel('Axial stress, \sigma','FontSize',12);
grid on;
set(gca,'FontSize',12);

subplot(1, 2, 2)
plot(T_save, reshape(Maceps_save(1, 1, :), sz, 1),'-b','LineWidth',3);
hold on;
plot(T_save, reshape(Maceps_save(2, 2, :), sz, 1),'--r','LineWidth',4);
hold on
plot(T_save, reshape(Maceps_save(3, 3, :), sz, 1),'-k','LineWidth',1.5);
legend('\epsilon_{11}','\epsilon_{22}','\epsilon_{33}','Location','Best','Orientation','horizontal')
xlabel('Time, t(s)','FontSize',12)
ylabel('Strain, \epsilon','FontSize',12)
grid on
set(gca,'FontSize',12)

saveas(gcf, 'result3.jpg');

Spiers = xlsread('Spiers1990');
figure(4)
set(gcf, 'Units', 'centimeter', 'Position', [5 10 45 10])
subplot(1, 3, 1)
plot(T_save, Maceps_v(1:sz),'r','Linewidth',3);
hold on
plot(Spiers(:, 5)*24*3600, -(Spiers(:, 6))/100, 'or', 'LineWidth', 1);
legend('Volumetric strain','Location','Best','Orientation','horizontal')
xlabel('Time, t(s)','FontSize',12)
ylabel('\epsilon_v, \Phi','FontSize',12)
grid on
set(gca,'FontSize',12)

subplot(1, 3, 2)
plot(T_save, Macporosity(1:sz),'b-','Linewidth',3);
hold on
legend('\Phi','Location','Best','Orientation','horizontal')
xlabel('Time, t(s)','FontSize',12)
ylabel('Porosity, \Phi','FontSize',12)
grid on
set(gca,'FontSize',12)

subplot(1, 3, 3)
plot(T_save, nonconnect_porosity(1:sz),'b-','Linewidth',3)
hold on 
plot(T_save, connect_porosity(1:sz), 'b-','Linewidth',3)
legend('non- \Phi','connect \Phi','Location','Best','Orientation','horizontal')
xlabel('Time, t(s)','FontSize',12)
ylabel('Porosity, \Phi','FontSize',12)
grid on
set(gca,'FontSize',12)


saveas(gcf, 'result4.jpg');

Ter = xlsread('Ter2005');
figure(5)
set(gcf, 'Units', 'centimeter', 'Position', [5 10 30 10]    )
subplot(1, 2, 1)
plot(reshape(Maceps_save(1, 1, :), sz, 1), reshape(Macsig_save(1, 1, :), sz, 1),'-b','LineWidth',3);
hold on;
plot(reshape(Maceps_save(1, 1, :), sz, 1), reshape(Macsig_save(2, 2, :), sz, 1),'--r','LineWidth',4);
hold on;
plot(reshape(Maceps_save(1, 1, :), sz, 1), reshape(Macsig_save(3, 3, :), sz, 1),'-k','LineWidth',1.5);
hold on;
plot(-Ter(:, 1), -Ter(:, 2)-50, 'or', 'LineWidth', 1);
legend('\sigma_{11}','\sigma_{22}','\sigma_{33}','Location','Best','Orientation','Horizontal');
xlabel('Axial strain, \epsilon','FontSize',12);
ylabel('Stress, \sigma','FontSize',12);
grid on;
set(gca,'FontSize',12);

subplot(1, 2, 2)
plot(T_save, reshape(Macsig_save(1, 1, :), sz, 1),'-b','LineWidth',3);
hold on;
plot(T_save, reshape(Macsig_save(2, 2, :), sz, 1),'--r','LineWidth',4);
hold on;
plot(T_save, reshape(Macsig_save(3, 3, :), sz, 1),'-k','LineWidth',1.5);
legend('\sigma_{11}','\sigma_{22}','\sigma_{33}','Location','Best','Orientation','Horizontal');
xlabel('Time, t(s)','FontSize',12)
ylabel('Stress, \sigma','FontSize',12)
grid on
set(gca,'FontSize',12)

saveas(gcf, 'result5.jpg');
end
