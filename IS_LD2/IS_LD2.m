%% PAGRINDINE
clc; clear all; close all; 
% Pasirinkti DNT struktura:
% vienas iejimas x
% vienas pasleptasis sluosnis su 5 neuronais
% vienas isejimas y1

% 1. Surinkti/paruosti mokymo duomenu rinkini
x = 0.1:1/22:1;
y1 = (1 + 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))/2;

% 2. Sugeneruoti pradines svoriu reiksmes
% 1 sluoksnio parametrai
w11_1 = randn(1); b1_1 = randn(1);
w21_1 = randn(1); b2_1 = randn(1);
w31_1 = randn(1); b3_1 = randn(1);
w41_1 = randn(1); b4_1 = randn(1);
w51_1 = randn(1); b5_1 = randn(1);
% 2 sluoksnio parametrai
w11_2 = randn(1); b1_2 = randn(1);
w12_2 = randn(1);
w13_2 = randn(1); 
w14_2 = randn(1);
w15_2 = randn(1); 

mok_zing = 0.15;

for indx = 1:10000
    for n = 1:20
        % 3. Apskaiciuoti tinklo atsaka (momentini)
        % pasverta suma 1 sluoksnyje
        v1_1 = x(n)*w11_1 + b1_1;
        v2_1 = x(n)*w21_1 + b2_1;
        v3_1 = x(n)*w31_1 + b3_1;
        v4_1 = x(n)*w41_1 + b4_1;
        v5_1 = x(n)*w51_1 + b5_1;
        % aktyvavimo f-ja 1 sluoksnyje
        y1_1 = 1/(1+exp(-v1_1));
        y2_1 = 1/(1+exp(-v2_1));
        y3_1 = 1/(1+exp(-v3_1));
        y4_1 = 1/(1+exp(-v4_1));
        y5_1 = 1/(1+exp(-v5_1));
        % pasverta suma isejimo sluoksnyje
        v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2 + y4_1*w14_2 + y5_1*w15_2 + b1_2;
       
        % skaiciuojame isejima/tinklo atsaka pritaikydami aktyvavimo f-jas
        y1_apskaiciuota = v1_2;
        
        % 4. Palyginti su norimu atsaku ir paskaiciuoti klaida
        e1 = y1(n) - y1_apskaiciuota;
        
        % 5. Atnaujinti rysiu svorius taip, kad klaida mazetu (atlikti tinklo mokyma)
        % formule: w = w + mok_zing*delta*input
        delta_out_1 = e1;
        
        delta_1_1 = y1_1*(1-y1_1)*(delta_out_1*w11_2);
        delta_2_1 = y2_1*(1-y2_1)*(delta_out_1*w12_2);
        delta_3_1 = y3_1*(1-y3_1)*(delta_out_1*w13_2);
        delta_4_1 = y4_1*(1-y4_1)*(delta_out_1*w14_2);
        delta_5_1 = y5_1*(1-y5_1)*(delta_out_1*w15_2);
        
        % atnaujinti svorius isejimo sluoksnyje
        w11_2 = w11_2 + mok_zing*delta_out_1*y1_1;
        w12_2 = w12_2 + mok_zing*delta_out_1*y2_1;
        w13_2 = w13_2 + mok_zing*delta_out_1*y3_1;
        w14_2 = w14_2 + mok_zing*delta_out_1*y4_1;
        w15_2 = w15_2 + mok_zing*delta_out_1*y5_1;
        b1_2 = b1_2 + mok_zing*delta_out_1*1;
        
        % atnaujinti svorius pasleptajame sluoksnyje
        w11_1 = w11_1 + mok_zing*delta_1_1*x(n);
        w21_1 = w21_1 + mok_zing*delta_2_1*x(n);
        w31_1 = w31_1 + mok_zing*delta_3_1*x(n);
        w41_1 = w41_1 + mok_zing*delta_4_1*x(n);
        w51_1 = w51_1 + mok_zing*delta_5_1*x(n);
        b1_1 = b1_1 + mok_zing*delta_1_1;
        b2_1 = b2_1 + mok_zing*delta_2_1;
        b3_1 = b3_1 + mok_zing*delta_3_1;
        b4_1 = b4_1 + mok_zing*delta_4_1;
        b5_1 = b5_1 + mok_zing*delta_5_1;
    end
end

x2 = linspace(0.1,1,100);
for m = 1:100
    % 3. Apskaiciuoti tinklo atsaka (momentini)
    % pasverta suma 1 sluoksnyje
    v1_1 = x2(m)*w11_1 + b1_1;
    v2_1 = x2(m)*w21_1 + b2_1;
    v3_1 = x2(m)*w31_1 + b3_1;
    v4_1 = x2(m)*w41_1 + b4_1;
    v5_1 = x2(m)*w51_1 + b5_1;
    % aktyvavimo f-ja 1 sluoksnyje
    y1_1 = 1/(1+exp(-v1_1));
    y2_1 = 1/(1+exp(-v2_1));
    y3_1 = 1/(1+exp(-v3_1));
    y4_1 = 1/(1+exp(-v4_1));
    y5_1 = 1/(1+exp(-v5_1));
    % pasverta suma isejimo sluoksnyje
    v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2 + y4_1*w14_2 + y5_1*w15_2 + b1_2;
    % skaiciuojame isejima/tinklo atsaka pritaikydami aktyvavimo f-jas
    y1_apskaiciuota(m) = v1_2;
end

figure,

plot(x,y1,'kx');
hold on
plot(x2,y1_apskaiciuota,'rx');
%% PAPILDOMA
clc; clear all; close all;
[X,Y,Z] = peaks(25);
C0(:,:,1) = zeros(25);
C0(:,:,2) = ones(25).*linspace(0.5,0.6,25);
C0(:,:,3) = ones(25).*linspace(0,1,25);
mesh(X,Y,Z,C0);
x1 = linspace(-3,3,25);
x2 = linspace(-3,3,25);


% Centru reiksmes ir spinduliu
c1_x = 0;
c1_y = 1.5;
r1 = 0.50;

c2_x = 0.5;
c2_y = 0.75;
r2 = 0.50;

c3_x = 1.25;
c3_y = 0;
r3 = 0.50;

c4_x = -1.25;
c4_y = 0.5;
r4 = 0.50;

c5_x = 0.25;
c5_y = -1.5;
r5 = 0.50;

c6_x = 0.25;
c6_y = 0.25;
r6 = 0.50;

% Rysiu svoriai
w1 = randn(1);
w2 = randn(1);
w3 = randn(1);
w4 = randn(1);
w5 = randn(1);
w6 = randn(1);
w0 = randn(1);
mok_zing = 0.1;
% Skaiciuojame SB f-ju atsakus/isejimus
for m = 1:25
    for n = 1:25
        Phi_1 = exp(-((x1(n) - c1_x)^2 + (x2(m) - c1_y)^2)/(2*r1^2));
        Phi_2 = exp(-((x1(n) - c2_x)^2 + (x2(m) - c2_y)^2)/(2*r2^2));
        Phi_3 = exp(-((x1(n) - c3_x)^2 + (x2(m) - c3_y)^2)/(2*r3^2));
        Phi_4 = exp(-((x1(n) - c4_x)^2 + (x2(m) - c4_y)^2)/(2*r4^2));
        Phi_5 = exp(-((x1(n) - c5_x)^2 + (x2(m) - c5_y)^2)/(2*r5^2));
        Phi_6 = exp(-((x1(n) - c6_x)^2 + (x2(m) - c6_y)^2)/(2*r6^2));
        % Skaiciuojama pasvertaja suma
        v = Phi_1*w1 + Phi_2*w2 + Phi_3*w3 + Phi_4*w4 + Phi_5*w5 + Phi_6*w6 + w0;
        y = v;
        e = Z(n,m) - y;
        % Rysiu svoriu atnaujinimas
        w1 = w1 + mok_zing*e*Phi_1;
        w2 = w2 + mok_zing*e*Phi_2;
        w3 = w3 + mok_zing*e*Phi_3;
        w4 = w4 + mok_zing*e*Phi_4;
        w5 = w5 + mok_zing*e*Phi_5;
        w6 = w6 + mok_zing*e*Phi_6;
        w0 = w0 + mok_zing*e;
    end
end

for n = 1:25
    for m = 1:25
        Phi_1 = exp(-((x1(n) - c1_x)^2 + (x2(m) - c1_y)^2)/(2*r1^2));
        Phi_2 = exp(-((x1(n) - c2_x)^2 + (x2(m) - c2_y)^2)/(2*r2^2));
        Phi_3 = exp(-((x1(n) - c3_x)^2 + (x2(m) - c3_y)^2)/(2*r3^2));
        Phi_4 = exp(-((x1(n) - c4_x)^2 + (x2(m) - c4_y)^2)/(2*r4^2));
        Phi_5 = exp(-((x1(n) - c5_x)^2 + (x2(m) - c5_y)^2)/(2*r5^2));
        Phi_6 = exp(-((x1(n) - c6_x)^2 + (x2(m) - c6_y)^2)/(2*r6^2));
        % Skaiciuojama pasvertaja suma
        v = Phi_1*w1 + Phi_2*w2 + Phi_3*w3 + Phi_4*w4 + Phi_5*w5 + Phi_6*w6 + w0;
        Z_est(n,m) = v;
    end
end
figure
mesh(X,Y,Z_est,C0);
