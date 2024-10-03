clear,clc,close;

u = 10^-6;
m=10^-3;
n=10^-9;
p=10^-12;
M=10^6;

core_radius = 5*u;
ac = pi*(core_radius^2);
gamma_p = 0.4;
gamma_s = 0.8;
nt = 1*(10^19)*M;   %m^3
tau = 10*m;
ase_bw = 0.1*n;
c=300*M;
h=6.62606957*(10^-34);

delta_z = 1;       %steps or delta in fiber
fiber_length = 0:delta_z:50;

sigma_a = csvread('absorption.csv'); %step in lambda is 0.1nm
sigma_a = sigma_a(:,2)/10000;
sigma_abs=@(x) sigma_a((round(round(x-1450,1)/0.1))+1); %takes input without nm

sigma_e = csvread('emission.csv'); %step in lambda is 0.1nm
sigma_e = sigma_e(:,2)/10000;
sigma_ems=@(x) sigma_e((round(round(x-1450,1)/0.1))+1); %takes input without nm

delta_pump = 0.01;
p_in_pump = 0.1:delta_pump:1; %Watts 41 is 0.5W
lambda_pump = 1480*n;

delta_signal = 1;
p_in_sig = -40:delta_signal:10; %Watts 41 is 1mW
p_in_sig = 10.^((p_in_sig-30)/10);

delta_lambda = 0.1*n;
lambda = 1525*n:delta_lambda:1565*n; %1550 at 251

figure(1)
plot(lambda,sigma_e(751:1151));
hold on;
plot(lambda,sigma_a(751:1151));
grid on; legend('Sigma Emission', 'Sigma Absorbtiom'); title('Emission and Absorbtion cross-sections');
xlabel('Wavelength (nm)');ylabel('Cross Section (m^2)');ax=gca; ax.XAxis.Exponent = -9;

for s = 1:1:5
    switch s
        
        case 1
            p_pump_z = zeros(length(lambda),length(fiber_length));
            p_sig_z = zeros(length(lambda),length(fiber_length));
            n2 = zeros(length(lambda),length(fiber_length));
            n1 = zeros(length(lambda),length(fiber_length));
            gain = zeros(length(lambda),length(fiber_length));
            ase_noise = zeros(length(lambda),length(fiber_length));
            p_pump_z(:,1) = p_in_pump(41);
            p_sig_z(:,1) = p_in_sig(41);
            
            for j =1:1:length(lambda)
                
                up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,1)/(h*c)) + (lambda(j)*gamma_s*sigma_abs(lambda(j)/n)*p_sig_z(j,1)/(h*c));
                down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,1)/(h*c)) + (lambda(j)*gamma_s*(sigma_abs(lambda(j)/n)+sigma_ems(lambda(j)/n))*p_sig_z(j,1)/(h*c));
                n2(j,1) = (up/(down+ac/tau))*nt;
                n1(j,1) = nt - n2(j,1);
                
                gain(j,1) = p_sig_z(j,1)/p_sig_z(j,1);
                
                ase_noise(j,1) = gamma_s*sigma_ems(lambda(j)/n)*n2(j,1)*(2*h*c^2*ase_bw/(lambda(j)^3));
                
                for i = 1:1:length(fiber_length)-1
                    
                    [zp,pp] = ode45(@(z,p) p*gamma_p*(n2(j,i)*sigma_ems(lambda_pump/n) - n1(j,i)*sigma_abs(lambda_pump/n)), [fiber_length(i) fiber_length(i+1)],p_pump_z(j,i));
                    p_pump_z(j,i+1) = (pp(end)-pp(1)) + p_pump_z(j,i);
                    
                    [zs,ps] = ode45(@(z,p) p*gamma_s*(n2(j,i)*sigma_ems(lambda(j)/n) - n1(j,i)*sigma_abs(lambda(j)/n)) + gamma_s*sigma_ems(lambda(j)/n)*n2(j,i)*(2*h*c^2*ase_bw/(lambda(j)^3)), [fiber_length(i) fiber_length(i+1)],p_sig_z(j,i));
                    p_sig_z(j,i+1) = (ps(end)-ps(1)) + p_sig_z(j,i);
                    
                    up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,i+1)/(h*c)) + (lambda(j)*gamma_s*sigma_abs(lambda(j)/n)*p_sig_z(j,i+1)/(h*c));
                    down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,i+1)/(h*c)) + (lambda(j)*gamma_s*(sigma_abs(lambda(j)/n)+sigma_ems(lambda(j)/n))*p_sig_z(j,i+1)/(h*c));
                    n2(j,i+1) = (up/(down+ac/tau))*nt;
                    n1(j,i+1) = nt - n2(j,i+1);
                    
                    gain(j,i+1) = p_sig_z(j,i+1)/p_sig_z(j,1);
                    
                    ase_noise(j,i+1) = gamma_s*sigma_ems(lambda(j)/n)*n2(j,i+1)*(2*h*c^2*ase_bw/(lambda(j)^3));
                    
                end
            end
            
            figure(2)
            plot(fiber_length,10*log10(1000*p_pump_z(251,:))); hold on; plot(fiber_length,10*log10(1000*p_sig_z(251,:)));
            grid on; legend('Pump Power', 'Signal Power'); title('Pump and Signal Power at Pump Power 0.5W or 26dBm, Signal Power 1mW or 0dBm');
            ylabel('Power (dBm)');xlabel('Fiber length (m)');
            
            figure(7)
            plot(lambda,10*log10(gain(:,20)));
            grid on; legend('Output Signal Gain'); title('Output Signal Gain at Pump Power 0.5W or 26dBm, Signal Power 1mW or 0dBm and Fiber Length 19m');
            ylabel('Gain (dB)');xlabel('Wavelength (nm)'); ax=gca; ax.XAxis.Exponent = -9;
            
        case 2   %%optimum is at 19 meters index = 20 at 0.5W
            delta_pump = 0.2;
            p_in_pump = 0.1:delta_pump:1; %Watts 3 is 0.5W
            lambda_pump = 1480*n;
            
            p_pump_z = zeros(length(p_in_pump),length(fiber_length));
            p_sig_z = zeros(length(p_in_pump),length(fiber_length));
            n2 = zeros(length(p_in_pump),length(fiber_length));
            n1 = zeros(length(p_in_pump),length(fiber_length));
            gain = zeros(length(p_in_pump),length(fiber_length));
            ase_noise = zeros(length(p_in_pump),length(fiber_length));
            p_sig_z(:,1) = p_in_sig(41);
            
            for j =1:1:length(p_in_pump)
                
                p_pump_z(j,1) = p_in_pump(j);
                
                up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,1)/(h*c));
                down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,1)/(h*c));
                n2(j,1) = (up/(down+ac/tau))*nt;
                n1(j,1) = nt - n2(j,1);
                
                gain(j,1) = p_sig_z(j,1)/p_sig_z(j,1);
                
                ase_noise(j,1) = gamma_s*sigma_ems(lambda(251)/n)*n2(j,1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                
                for i = 1:1:length(fiber_length)-1
                    
                    [zp,pp] = ode45(@(z,p) p*gamma_p*(n2(j,i)*sigma_ems(lambda_pump/n) - n1(j,i)*sigma_abs(lambda_pump/n)), [fiber_length(i) fiber_length(i+1)],p_pump_z(j,i));
                    p_pump_z(j,i+1) = p_pump_z(j,i)+(pp(end)-pp(1));
                    
                    [zs,ps] = ode45(@(z,p) p*gamma_s*(n2(j,i)*sigma_ems(lambda(251)/n) - n1(j,i)*sigma_abs(lambda(251)/n)) + gamma_s*sigma_ems(lambda(251)/n)*n2(j,i)*(2*h*c^2*ase_bw/(lambda(251)^3)), [fiber_length(i) fiber_length(i+1)],p_sig_z(j,i));
                    p_sig_z(j,i+1) = p_sig_z(j,i)+(ps(end)-ps(1));
                    
                    up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,i+1)/(h*c));
                    down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,i+1)/(h*c));
                    n2(j,i+1) = (up/(down+ac/tau))*nt;
                    n1(j,i+1) = nt - n2(j,i+1);
                    
                    gain(j,i+1) = p_sig_z(j,i+1)/p_sig_z(j,1);
                    
                    ase_noise(j,i+1) = gamma_s*sigma_ems(lambda(251)/n)*n2(i+1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                end
            end
            
            figure(3)
            subplot(2,1,1)
            plot(fiber_length,10*log10(gain));
            grid on; legend('Pump Power = 0.1', 'Pump Power = 0.3', 'Pump Power = 0.5', 'Pump Power = 0.7', 'Pump Power = 0.9'); title('Gain at Diffrent Pump Powers Across Fiber Length');
            ylabel('Gain (dB)');xlabel('Fiber length (m)');
            subplot(2,1,2)
            plot(fiber_length,10*log10(gain(3,:)));
            grid on; legend('Pump Power = 0.5'); title('Gain at 0.5W Pump Powers Across Fiber Length');
            ylabel('Gain (dB)');xlabel('Fiber length (m)');
            
            
        case 3 %%keep z(20) = 19m, lambda const, pin pump --- pin signal is variable
            delta_pump = 0.01;
            p_in_pump = 0.1:delta_pump:1; %Watts 41 is 0.5W
            lambda_pump = 1480*n;
            
            p_pump_z = zeros(length(p_in_sig),length(fiber_length));
            p_sig_z = zeros(length(p_in_sig),length(fiber_length));
            n2 = zeros(length(p_in_sig),length(fiber_length));
            n1 = zeros(length(p_in_sig),length(fiber_length));
            gain = zeros(length(p_in_sig),length(fiber_length));
            ase_noise = zeros(length(p_in_sig),length(fiber_length));
            p_pump_z(:,1) = p_in_pump(41);
            
            
            for j =1:1:length(p_in_sig)
                
                p_sig_z(j,1) = p_in_sig(j);
                
                up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,1)/(h*c));
                down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,1)/(h*c));
                n2(j,1) = (up/(down+ac/tau))*nt;
                n1(j,1) = nt - n2(j,1);
                
                gain(j,1) = p_sig_z(j,1)/p_sig_z(j,1);
                
                ase_noise(j,1) = gamma_s*sigma_ems(lambda(251)/n)*n2(j,1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                
                for i = 1:1:length(fiber_length)-1
                    
                    [zp,pp] = ode45(@(z,p) p*gamma_p*(n2(j,i)*sigma_ems(lambda_pump/n) - n1(j,i)*sigma_abs(lambda_pump/n)), [fiber_length(i) fiber_length(i+1)],p_pump_z(j,i));
                    p_pump_z(j,i+1) = p_pump_z(j,i)+(pp(end)-pp(1));
                    
                    [zs,ps] = ode45(@(z,p) p*gamma_s*(n2(j,i)*sigma_ems(lambda(251)/n) - n1(j,i)*sigma_abs(lambda(251)/n)) + gamma_s*sigma_ems(lambda(251)/n)*n2(j,i)*(2*h*c^2*ase_bw/(lambda(251)^3)), [fiber_length(i) fiber_length(i+1)],p_sig_z(j,i));
                    p_sig_z(j,i+1) = p_sig_z(j,i)+(ps(end)-ps(1));
                    
                    up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,i+1)/(h*c));
                    down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,i+1)/(h*c));
                    n2(j,i+1) = (up/(down+ac/tau))*nt;
                    n1(j,i+1) = nt - n2(j,i+1);
                    
                    gain(j,i+1) = p_sig_z(j,i+1)/p_sig_z(j,1);
                    
                    ase_noise(j,i+1) = gamma_s*sigma_ems(lambda(251)/n)*n2(i+1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                end
            end
            
            figure(4)
            plot(10*log10(1000*p_in_sig),10*log10(1000*p_sig_z(:,20)))
            grid on; legend('Output Signal Power'); title('Output Signal Power at 0.5W or 26dBm pump power at Fiber legth = 19m');
            ylabel('Output Signal Power (dBm)');xlabel('Input Signal Power (dBm)');
            
            figure(5)
            plot(10*log10(1000*p_in_sig),10*log10(gain(:,20)))
            grid on; legend('Output Signal Gain'); title('Output Signal Gain at 0.5W or 26dBm pump power at Fiber legth = 19m');
            ylabel('Output Signal Gain (dB)');xlabel('Input Signal Power (dBm)');
            
        case 4
            p_pump_z = zeros(length(lambda),length(fiber_length));
            p_sig_z = zeros(length(lambda),length(fiber_length));
            n2 = zeros(length(lambda),length(fiber_length));
            n1 = zeros(length(lambda),length(fiber_length));
            gain = zeros(length(lambda),length(fiber_length));
            ase_noise = zeros(length(lambda),length(fiber_length));
            p_pump_z(:,1) = p_in_pump(41);
            p_sig_z(251,1) = p_in_sig(41);
            
            for j =1:1:length(lambda)
                
                up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,1)/(h*c)) + (lambda(j)*gamma_s*sigma_abs(lambda(j)/n)*p_sig_z(j,1)/(h*c));
                down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,1)/(h*c)) + (lambda(j)*gamma_s*(sigma_abs(lambda(j)/n)+sigma_ems(lambda(j)/n))*p_sig_z(j,1)/(h*c));
                n2(j,1) = (up/(down+ac/tau))*nt;
                n1(j,1) = nt - n2(j,1);
                
                gain(j,1) = p_sig_z(j,1)/p_sig_z(j,1);
                
                ase_noise(j,1) = gamma_s*sigma_ems(lambda(j)/n)*n2(j,1)*(2*h*c^2*ase_bw/(lambda(j)^3));
                
                for i = 1:1:length(fiber_length)-1
                    
                    [zp,pp] = ode45(@(z,p) p*gamma_p*(n2(j,i)*sigma_ems(lambda_pump/n) - n1(j,i)*sigma_abs(lambda_pump/n)), [fiber_length(i) fiber_length(i+1)],p_pump_z(j,i));
                    p_pump_z(j,i+1) = (pp(end)-pp(1)) + p_pump_z(j,i);
                    
                    [zs,ps] = ode45(@(z,p) p*gamma_s*(n2(j,i)*sigma_ems(lambda(j)/n) - n1(j,i)*sigma_abs(lambda(j)/n)) + gamma_s*sigma_ems(lambda(j)/n)*n2(j,i)*(2*h*c^2*ase_bw/(lambda(j)^3)), [fiber_length(i) fiber_length(i+1)],p_sig_z(j,i));
                    p_sig_z(j,i+1) = (ps(end)-ps(1)) + p_sig_z(j,i);
                    
                    up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,i+1)/(h*c)) + (lambda(j)*gamma_s*sigma_abs(lambda(j)/n)*p_sig_z(j,i+1)/(h*c));
                    down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,i+1)/(h*c)) + (lambda(j)*gamma_s*(sigma_abs(lambda(j)/n)+sigma_ems(lambda(j)/n))*p_sig_z(j,i+1)/(h*c));
                    n2(j,i+1) = (up/(down+ac/tau))*nt;
                    n1(j,i+1) = nt - n2(j,i+1);
                    
                    gain(j,i+1) = p_sig_z(j,i+1)/p_sig_z(j,1);
                    
                    ase_noise(j,i+1) = gamma_s*sigma_ems(lambda(j)/n)*n2(j,i+1)*(2*h*c^2*ase_bw/(lambda(j)^3));
                    
                end
            end
            
            figure(8)
            plot(lambda,10*log10(1000*p_sig_z(:,21)));
            grid on; legend('Output ASE Cloud'); title('Output ASE Cloud at Pump Power 0.5W, Signal Power 1mW or 0dBm and Fiber Length 19m at 1550nm Wavelength');
            ylabel('Output Power (dBm)');xlabel('Wavelength (nm)'); ax=gca; ax.XAxis.Exponent = -9;
            
        case 5 %optimum length for pump power 0.1W is 10m, index z(11)
            delta_pump = 0.01;
            p_in_pump = 0.1:delta_pump:1; %Watts 41 is 0.5W
            lambda_pump = 1480*n;
            
            p_pump_z = zeros(length(p_in_sig),length(fiber_length));
            p_sig_z = zeros(length(p_in_sig),length(fiber_length));
            n2 = zeros(length(p_in_sig),length(fiber_length));
            n1 = zeros(length(p_in_sig),length(fiber_length));
            gain = zeros(length(p_in_sig),length(fiber_length));
            ase_noise = zeros(length(p_in_sig),length(fiber_length));
            p_pump_z(:,1) = p_in_pump(1);
            
            
            for j =1:1:length(p_in_sig)
                
                p_sig_z(j,1) = p_in_sig(j);
                
                up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,1)/(h*c));
                down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,1)/(h*c));
                n2(j,1) = (up/(down+ac/tau))*nt;
                n1(j,1) = nt - n2(j,1);
                
                gain(j,1) = p_sig_z(j,1)/p_sig_z(j,1);
                
                ase_noise(j,1) = gamma_s*sigma_ems(lambda(251)/n)*n2(j,1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                
                for i = 1:1:length(fiber_length)-1
                    
                    [zp,pp] = ode45(@(z,p) p*gamma_p*(n2(j,i)*sigma_ems(lambda_pump/n) - n1(j,i)*sigma_abs(lambda_pump/n)), [fiber_length(i) fiber_length(i+1)],p_pump_z(j,i));
                    p_pump_z(j,i+1) = p_pump_z(j,i)+(pp(end)-pp(1));
                    
                    [zs,ps] = ode45(@(z,p) p*gamma_s*(n2(j,i)*sigma_ems(lambda(251)/n) - n1(j,i)*sigma_abs(lambda(251)/n)) + gamma_s*sigma_ems(lambda(251)/n)*n2(j,i)*(2*h*c^2*ase_bw/(lambda(251)^3)), [fiber_length(i) fiber_length(i+1)],p_sig_z(j,i));
                    p_sig_z(j,i+1) = p_sig_z(j,i)+(ps(end)-ps(1));
                    
                    up = (lambda_pump*gamma_p*sigma_abs(lambda_pump/n)*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*sigma_abs(lambda(251)/n)*p_sig_z(j,i+1)/(h*c));
                    down = (lambda_pump*gamma_p*(sigma_abs(lambda_pump/n)+sigma_ems(lambda_pump/n))*p_pump_z(j,i+1)/(h*c)) + (lambda(251)*gamma_s*(sigma_abs(lambda(251)/n)+sigma_ems(lambda(251)/n))*p_sig_z(j,i+1)/(h*c));
                    n2(j,i+1) = (up/(down+ac/tau))*nt;
                    n1(j,i+1) = nt - n2(j,i+1);
                    
                    gain(j,i+1) = p_sig_z(j,i+1)/p_sig_z(j,1);
                    
                    ase_noise(j,i+1) = gamma_s*sigma_ems(lambda(251)/n)*n2(i+1)*(2*h*c^2*ase_bw/(lambda(251)^3));
                end
            end
            
            figure(6)
            plot(10*log10(1000*p_in_sig),10*log10(gain(:,11)))
            grid on; legend('Output Signal Gain'); title('Output Signal Gain at 0.1W or 20dBm pump power at Fiber legth = 10m');
            ylabel('Output Signal Gain (dB)');xlabel('Input Signal Power (dBm)');
            
    end
    
end





