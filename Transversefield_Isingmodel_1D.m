clc 
clear
tic
%=============================================================
%           Transverse-field Ising model              
%=============================================================
%Parameters:
N = 14; %number of sites
J = 1;
bnd = 1; %periodic bc = 1 & %open bc = 0

%pauli matrices:
sigma_x = [0 1;1 0];sigma_y = [0 -1i;1i 0];sigma_z = [1 0;0 -1];I = eye(2);
X = sparse(sigma_x);Y = sparse(sigma_y);Z = sparse(sigma_z);I = sparse(I);

first_h = 0; last_h = 5; steps = 0.1;
Ns = ((last_h - first_h)/steps)+1; %Number of samples
Nz = (N*(N-1))/2; %Number of correlation functions
Nf = Nz + N +1;   %Number of features

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
G_State = zeros(2^N,Ns);   %Ground state wave function 
exp_Sx = zeros(Ns,1);      %Expectation values for Sx
corr_SzSz = zeros(Ns,Nz);  %Correlation functions for Sz
exp_mio = zeros(Ns,N);     %Expectation values for mio
Datas = zeros(Ns+1,Nf);    %Datas in form of (X,Y)
k=0;
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
    [U,D] = eigs(H,5,'sa'); 
    energy = diag(D); %Energy spectrum
    G_Energy(:,k) = energy(1); %Ground state energy
    G_State(:,k) = U(:,1); %Ground state wavefunction
    
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
    Datas(k+1,1:Nz) = corr_SzSz(k,:);
    Datas(k+1,Nz+1:Nf-2) = exp_mio(k,1:N-1);
    Datas(k+1,Nf-1) = exp_Sx(k,1);
end
Datas(2:Ns+1,Nf) = y;

%% plots:
%% (SzSz correlation functions)
figure(1)
x_plot=linspace(1,Nz,Nz);
plot(x_plot,corr_SzSz(2,:),'-*');
legend('h=0.1')
hold on
plot(x_plot,corr_SzSz(6,:),'-*','DisplayName','h=0.5');plot(x_plot,corr_SzSz(9,:),'-*','DisplayName','h=0.8');
plot(x_plot,corr_SzSz(11,:),'-*','DisplayName','h=1');plot(x_plot,corr_SzSz(13,:),'-*','DisplayName','h=1.2');
plot(x_plot,corr_SzSz(16,:),'-*','DisplayName','h=1.5');plot(x_plot,corr_SzSz(21,:),'-*','DisplayName','h=2');
plot(x_plot,corr_SzSz(31,:),'-*','DisplayName','h=3');plot(x_plot,corr_SzSz(41,:),'-*','DisplayName','h=4');
plot(x_plot,corr_SzSz(51,:),'-*','DisplayName','h=5');hold all
t=title('Correlation Functions(S_zS_z)');
xlim([1 Nz])
ylabel('Correlation functions');
xlabel('Number of sites');
grid on

%% (SzSz correlation functions for l = 1)
figure(2)
x_plot=linspace(1,N-1,N-1);
plot(x_plot,corr_SzSz(2,1:N-1),'-*');
legend('h=0.1')
hold on
plot(x_plot,corr_SzSz(6,1:N-1),'-*','DisplayName','h=0.5');plot(x_plot,corr_SzSz(9,1:N-1),'-*','DisplayName','h=0.8');
plot(x_plot,corr_SzSz(11,1:N-1),'-*','DisplayName','h=1');plot(x_plot,corr_SzSz(13,1:N-1),'-*','DisplayName','h=1.2');
plot(x_plot,corr_SzSz(16,1:N-1),'-*','DisplayName','h=1.5');plot(x_plot,corr_SzSz(21,1:N-1),'-*','DisplayName','h=2');
plot(x_plot,corr_SzSz(31,1:N-1),'-*','DisplayName','h=3');plot(x_plot,corr_SzSz(41,1:N-1),'-*','DisplayName','h=4');
plot(x_plot,corr_SzSz(51,1:N-1),'-*','DisplayName','h=5');hold all
t=title('Correlation Functions(S_zS_z) for site 1');
xlim([1 N-1])
ylabel('Correlation functions');
xlabel('Number of sites');
grid on

%% (expectation values for Sx)
figure(3)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,exp_Sx(:,1),'-*');
t=title('Expectation Value(S_x)');
xlim([0 last_h])
ylabel('<S_x>');
xlabel('h');
grid on

%% (expectation values for mio_x)
figure(4)
x_plot=linspace(1,N-1,N-1);
plot(x_plot,exp_mio(2,1:N-1),'-*');
legend('h=0.1')
hold on
plot(x_plot,exp_mio(9,1:N-1),'-*','DisplayName','h=0.8');plot(x_plot,exp_mio(11,1:N-1),'-*','DisplayName','h=1');
plot(x_plot,exp_mio(13,1:N-1),'-*','DisplayName','h=1.2');plot(x_plot,exp_mio(16,1:N-1),'-*','DisplayName','h=1.5');
plot(x_plot,exp_mio(21,1:N-1),'-*','DisplayName','h=2');plot(x_plot,exp_mio(31,1:N-1),'-*','DisplayName','h=3');
plot(x_plot,exp_mio(41,1:N-1),'-*','DisplayName','h=4');plot(x_plot,exp_mio(51,1:N-1),'-*','DisplayName','h=5');
hold all
t=title('Expectation Value(mio_l)');
xlim([1 N-1])
ylabel('<mio_l>');
xlabel('Number of sites');
grid on

%% (expectation values for mio_x)
figure(5)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,exp_mio(:,1),'-*');
legend('l=1');
hold on
plot(h_plot,exp_mio(:,2),'-*','DisplayName','l=2');plot(h_plot,exp_mio(:,3),'-*','DisplayName','l=3');
plot(h_plot,exp_mio(:,4),'-*','DisplayName','l=4');plot(h_plot,exp_mio(:,5),'-*','DisplayName','l=5');
plot(h_plot,exp_mio(:,6),'-*','DisplayName','l=6');plot(h_plot,exp_mio(:,7),'-*','DisplayName','l=7');
plot(h_plot,exp_mio(:,8),'-*','DisplayName','l=8');plot(h_plot,exp_mio(:,9),'-*','DisplayName','l=9');
plot(h_plot,exp_mio(:,10),'-*','DisplayName','l=10');plot(h_plot,exp_mio(:,11),'-*','DisplayName','l=11');
plot(h_plot,exp_mio(:,12),'-*','DisplayName','l=12');plot(h_plot,exp_mio(:,13),'-*','DisplayName','l=13');
hold all
t=title('Expectation Value(mio)');
xlim([0 last_h])
ylabel('<mio_l>');
xlabel('h');
grid on
%% (SzSz correlation functions vs h)
figure(6)
h_plot=linspace(0,last_h,Ns);
plot(h_plot,corr_SzSz(:,1),'-*');
t=title('Correlation Functions(S_zS_z)');
xlim([0 last_h])
ylabel('Correlation functions');
xlabel('h');
grid on
toc
%% 
%csvwrite('Data_XY.csv',Datas);