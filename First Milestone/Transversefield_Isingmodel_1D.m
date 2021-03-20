clc 
clear
tic
%=============================================================
%            1D Transverse-field Ising model             
%=============================================================
%Parameters:
N = 12; %number of sites
J = 1;
bnd = 1; %periodic bc = 1 & %open bc = 0

%pauli matrices:
sigma_x = [0 1;1 0];sigma_y = [0 -1i;1i 0];sigma_z = [1 0;0 -1];I = eye(2);
X = sparse(sigma_x);Y = sparse(sigma_y);Z = sparse(sigma_z);I = sparse(I);

first_h = 0; last_h = 7; steps = 0.1;
Ns = ((last_h - first_h)/steps)+1; %Number of samples
Nz = (N*(N-1))/2; %Number of correlation functions
Nf = Nz + 2*N + 2^N +1;   %Number of features
y = linspace(0,last_h,Ns); %'h' which is y in our data
y = reshape(y,Ns,1);

%% Construction of S Matrices
Sx = cell(1,N);Sz = cell(1,N);
I_others = 1;
for i=1:N
    Sx{i} = kron(I_others,X);
    Sz{i} = kron(I_others,Z);
    I_others = kron(I_others,I);
    
    for j=1:i-1
        Sx{j} = kron(Sx{j},I);
        Sz{j} = kron(Sz{j},I);
    end
    
    for j=i+1:N
        Sx{j} = I_others;
        Sz{j} = I_others;
    end
end

%% Construction of mio Matrices
mio = cell(1,N);
X_cell = cell(1,N); X_cell{1} = X; I_cell=cell(1,N);
for i=2:N
    X_cell{i} = kron(X,X_cell{i-1});
end
for i=1:N
    I_cell{i} = eye(2^(i));
end 
for i=1:N
    if (i~=N)
        mio{i} = kron(X_cell{i},I_cell{N-i});
    end
    if (i==N)
        mio{i} = X_cell{i};
    end
end
%% One dimentional Quantum Ising model & Data generation:

%Desired Quantities:
G_Energy = zeros(1,Ns);    %Ground state energy 
Energy_gap = zeros(Ns,1);  %Energy gap
G_State = zeros(2^N,Ns);   %Ground state wave function
G_Data = zeros(Ns,2^N);    %Ground state wave function 
exp_Sx = zeros(Ns,1);      %Expectation values for Sx
corr_SzSz = zeros(Ns,Nz);  %Correlation functions for Sz
exp_mio = zeros(Ns,N);     %Expectation values for mio
RDM = cell(1,N-1);         %Reduced Density matrix for each subsystem 
Entanglemennt_Entropy = zeros(Ns,N-1); %Entanglemennt_Entropy
Datas = zeros(Ns+1,Nf);    %Datas in form of (X,Y)
k=0; M =0;
%Main Loop
for h=first_h:steps:last_h
    k=k+1;
    
    %=======================================
    %   1D Quantum ising model Hamiltonian             
    %=======================================
    H = -J*bnd*(Sz{N}*Sz{1}) -h*(Sx{N}); %boundary conditions
    for i=1:N-1
        H = H -J*(Sz{i}*Sz{i+1}) -h*(Sx{i});
    end  
    H = (H+H')/2; %making sure that it's hermitian
   
    %======================
    %   Diagonalization           
    %======================
    [States,Energys] = eigs(H,3,'sa'); 
    energy = diag(Energys); %Energy spectrum
    G_Energy(:,k) = energy(1); %Ground state energy
    Energy_gap(k,1) = abs(energy(3) - G_Energy(:,k)); %Energy gap
    G_State(:,k) = States(:,1);  %Ground state wavefunction
    norm = abs(G_State(:,k)'*G_State(:,k));
    G_State(:,k)= G_State(:,k)/norm;%normalization
    
    %===========================
    %   Entanglemennt Entropy         
    %===========================
    for i=1:N-1
        N_left = i;
        N_right = N - i;
        psi = reshape(G_State(:,k),2^(N_right),2^(N_left));
        RDM{i} = psi'*psi ;
        [U,D]=eig(RDM{i});
        p = diag(D);
        p (p< 10^-15) = 10 ^ -15;
        Entanglemennt_Entropy(k,i) = -dot(p,log(p));
    end
        
    %==================================================
    %   Expectation values & Correlation functions         
    %==================================================
    exp_Sx(k,1) = (G_State(:,k)')*Sx{N/2}*(G_State(:,k));
    c=1;
    for i=1:N
        for j=i+1:N
            corr_SzSz(k,c) = ((G_State(:,k)'*(Sz{i}*Sz{j})*G_State(:,k))- (G_State(:,k)'*(Sz{i})*G_State(:,k))*(G_State(:,k)'*(Sz{j})*G_State(:,k)));
            c = c + 1;
        end
        exp_mio(k,i) = (G_State(:,k)')*mio{i}*(G_State(:,k));
    end 
    G_Data(k,:) = reshape(G_State(:,k),1,2^N);
    Datas(k+1,1:Nz) = corr_SzSz(k,:);%66
    Datas(k+1,Nz+1) = exp_Sx(k,1);%67 
    Datas(k+1,Nz+2:Nz+12) = exp_mio(k,1:N-1); %78
    Datas(k+1,Nz+13) = Energy_gap(k,1); %79
    Datas(k+1,Nz+14:Nz+24) = Entanglemennt_Entropy(k,:);%90
    Datas(k+1,Nz+25:Nf-1) = G_Data(k,:); %4186
    
end
Datas(2:Ns+1,Nf) = y;%4187

%% plots:
%%(All SzSz correlation functions)
figure(1)
x_plot=linspace(1,Nz,Nz);
plot(x_plot,corr_SzSz(6,:),'-*');
legend('h=0.5')
hold on
plot(x_plot,corr_SzSz(11,:),'-*','DisplayName','h=1');
plot(x_plot,corr_SzSz(16,:),'-*','DisplayName','h=1.5');plot(x_plot,corr_SzSz(21,:),'-*','DisplayName','h=2');
plot(x_plot,corr_SzSz(31,:),'-*','DisplayName','h=3');plot(x_plot,corr_SzSz(41,:),'-*','DisplayName','h=4');
plot(x_plot,corr_SzSz(51,:),'-*','DisplayName','h=5');hold all
t=title('Correlation Functions(S_zS_z)');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([1 Nz])
ylabel({'$<\sigma_z^i \sigma_z^{i+1}> - <\sigma_z^i><\sigma_z^{i+1}>$'},'Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('Number of sites');

%% (SzSz correlation functions for first site)
figure(2)
x_plot=linspace(1,N-1,N-1);
plot(x_plot,corr_SzSz(6,1:N-1),'-s','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
legend('h=0.5')
hold on
plot(x_plot,corr_SzSz(9,1:N-1),'-s','MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'DisplayName','h=0.8');
plot(x_plot,corr_SzSz(11,1:N-1),'-s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'DisplayName','h=1');
plot(x_plot,corr_SzSz(13,1:N-1),'-s','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880],'DisplayName','h=1.2');
plot(x_plot,corr_SzSz(16,1:N-1),'-s','MarkerEdgeColor',	[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'DisplayName','h=1.5');
plot(x_plot,corr_SzSz(21,1:N-1),'-s','MarkerEdgeColor',[0.75, 0, 0.75],'MarkerFaceColor',[0.75, 0, 0.75],'DisplayName','h=2');
plot(x_plot,corr_SzSz(31,1:N-1),'-s','MarkerEdgeColor',	[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'DisplayName','h=3');
plot(x_plot,corr_SzSz(41,1:N-1),'-s','MarkerEdgeColor',[0, 0.5, 0],'MarkerFaceColor',[0, 0.5, 0],'DisplayName','h=4');
plot(x_plot,corr_SzSz(51,1:N-1),'-s','MarkerEdgeColor',[0, 0.75, 0.75],'MarkerFaceColor',[0, 0.75, 0.75],'DisplayName','h=5');hold all
t=title('Correlation Functions');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([1 N-1])
ylabel({'$<\sigma_z^1 \sigma_z^{i+1}> - <\sigma_z^1><\sigma_z^{i+1}>$'},'Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('Number of sites');

%% (Expectation values for mio_x vs sites)
figure(3)
x_plot=linspace(1,N-1,N-1);
plot(x_plot,exp_mio(6,1:N-1),'-s','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
legend('h=0.5')
hold on
plot(x_plot,exp_mio(9,1:N-1),'-s','MarkerEdgeColor',[0.75, 0.75, 0],'MarkerFaceColor',[0.75, 0.75, 0],'DisplayName','h=0.8');
plot(x_plot,exp_mio(11,1:N-1),'-s','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'DisplayName','h=1');
plot(x_plot,exp_mio(16,1:N-1),'-s','MarkerEdgeColor',[0.75, 0, 0.75],'MarkerFaceColor',	[0.75, 0, 0.75],'DisplayName','h=1.5');
plot(x_plot,exp_mio(21,1:N-1),'-s','MarkerEdgeColor',[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'DisplayName','h=2');
plot(x_plot,exp_mio(31,1:N-1),'-s','MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'DisplayName','h=3');
plot(x_plot,exp_mio(41,1:N-1),'-s','MarkerEdgeColor',[0, 0.75, 0.75],'MarkerFaceColor',[0, 0.75, 0.75],'DisplayName','h=4');
plot(x_plot,exp_mio(51,1:N-1),'-s','MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'DisplayName','h=5');
hold all
xlim([1 N-1])
t=title('Expectation values for \mu_l');
set(0,'DefaultAxesTitleFontWeight','normal');
ylabel('$<\mu_l>$','Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('Number of sites');

%% (Expectation values for mio_x vs h)
figure(4)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,exp_mio(:,1),'-*');
legend('l=1');
hold on
plot(h_plot,exp_mio(:,2),'-*','DisplayName','l=2');plot(h_plot,exp_mio(:,3),'-*','DisplayName','l=3');
plot(h_plot,exp_mio(:,4),'-*','DisplayName','l=4');plot(h_plot,exp_mio(:,5),'-*','DisplayName','l=5');
plot(h_plot,exp_mio(:,6),'-*','DisplayName','l=6');plot(h_plot,exp_mio(:,7),'-*','DisplayName','l=7');
plot(h_plot,exp_mio(:,8),'-*','DisplayName','l=8');plot(h_plot,exp_mio(:,9),'-*','DisplayName','l=9');
plot(h_plot,exp_mio(:,10),'-*','DisplayName','l=10');plot(h_plot,exp_mio(:,11),'-*','DisplayName','l=11');
xline(1,'--k','DisplayName','h = 1');
hold all
xlim([0 last_h])
t=title('Expectation values for \mu_l');
set(0,'DefaultAxesTitleFontWeight','normal');
ylabel('$<\mu_l>$','Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('h');

%% (SzSz correlation functions vs h)
figure(5)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,corr_SzSz(:,1),'-s','MarkerEdgeColor','b','MarkerFaceColor',[0, 0.75, 0.75]);
xline(1,'--k');
t=title('Correlation Functions');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([0 last_h])
ylabel({'$<\sigma_z^1 \sigma_z^{2}> - <\sigma_z^1><\sigma_z^{2}>$'},'Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('h');

%% (Entanglemennt_Entropy vs size of subsystem)
figure(6)
x_plot = linspace(1,N-1,N-1);
plot(x_plot,Entanglemennt_Entropy(11,:),'-s','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
legend('h=1')
hold on
plot(x_plot,Entanglemennt_Entropy(16,:),'-s','MarkerEdgeColor',	[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'DisplayName','h=1.5');
plot(x_plot,Entanglemennt_Entropy(21,:),'-s','MarkerEdgeColor',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'DisplayName','h=2');
plot(x_plot,Entanglemennt_Entropy(31,:),'-s','MarkerEdgeColor',	[0.75, 0, 0.75],'MarkerFaceColor',	[0.75, 0, 0.75],'DisplayName','h=3');
plot(x_plot,Entanglemennt_Entropy(41,:),'-s','MarkerEdgeColor',[0, 0.75, 0.75],'MarkerFaceColor',[0, 0.75, 0.75],'DisplayName','h=4');
plot(x_plot,Entanglemennt_Entropy(51,:),'-s','MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'DisplayName','h=5');hold all
t=title('Entanglement Entropy');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([1 N-1])
ylabel('S');
xlabel('size of subsystem');

%% (Energy difference)
figure(7)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,Energy_gap(:,1),'-s','MarkerEdgeColor','b','MarkerFaceColor',[0, 0.75, 0.75]);
xline(1,'--k');
t=title('Energy difference');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([0 last_h])
ylabel('|E_2 - E_0|');
xlabel('h');

%% (Entanglemennt_Entropy vs h)
figure(8)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,Entanglemennt_Entropy(:,1),'-*');
legend('size=1')
hold on
plot(h_plot,Entanglemennt_Entropy(:,2),'-*','DisplayName','size=2');plot(h_plot,Entanglemennt_Entropy(:,3),'-*','DisplayName','size=3');
plot(h_plot,Entanglemennt_Entropy(:,4),'-*','DisplayName','size=4');plot(h_plot,Entanglemennt_Entropy(:,5),'-*','DisplayName','size=5');
plot(h_plot,Entanglemennt_Entropy(:,6),'-*','DisplayName','size=6');plot(h_plot,Entanglemennt_Entropy(:,7),'-*','DisplayName','size=7');
plot(h_plot,Entanglemennt_Entropy(:,8),'-*','DisplayName','size=8');plot(h_plot,Entanglemennt_Entropy(:,9),'-*','DisplayName','size=9');
plot(h_plot,Entanglemennt_Entropy(:,10),'-*','DisplayName','size=10');plot(h_plot,Entanglemennt_Entropy(:,11),'-*','DisplayName','size=11');
xline(1,'--k','DisplayName','h = 1');
hold all
t=title('Entanglement Entropy');
set(0,'DefaultAxesTitleFontWeight','normal');
xlim([0 last_h])
ylabel('S');
xlabel('h');
%%
%%(Expectation values for S_x vs h)
figure(9)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,exp_Sx(:,1),'-s','MarkerEdgeColor','b','MarkerFaceColor',[0, 0.75, 0.75]);
xline(1,'--k');
xlim([0 last_h])
t=title('Expectation values for \sigma_x');
set(0,'DefaultAxesTitleFontWeight','normal');
ylabel('$<\sigma_x>$','Interpreter','latex','FontSize',15,'FontWeight','bold');
xlabel('h');
toc
csvwrite('Data_XY.csv',Datas);

