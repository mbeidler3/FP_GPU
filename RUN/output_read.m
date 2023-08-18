format longE

%% LOAD OUTPUT DATA

filename1='./output.korc';

fileID = fopen(filename1,'r');

tmp = fgetl(fileID);

tmp = textscan(fileID,'%s %s %s %u',1);
nRE=tmp{4};

tmp = fgetl(fileID);
tmp = fgetl(fileID);
tmp = fgetl(fileID);
tmp = fgetl(fileID);

tmp = textscan(fileID,'%s %s %f',1);
sim_time=tmp{3};

tmp = textscan(fileID,'%s %s %s %d',1);
ndump=tmp{4};

fclose(fileID);

%% LOAD OUTPUT DATA

KE=zeros(1,ndump);

filename2='./data.korc';

fileID = fopen(filename2,'r');


for ii=1:ndump    
    tmp = textscan(fileID,'%s %s %s %f %f',1);
    KE(ii)=tmp{4};
end

fclose(fileID);

%% Analysis

c=299792458;
me=9.109382E-31;
qe=1.602176E-19;
mu0=4.0*pi*1E-7;
ep0= 1.0/(mu0*c^2);

time=linspace(0,sim_time,ndump);

KE0=10.e6*qe; 
E0=KE0+(me*c^2);
ne=1.E20; 
Te=1.E0*qe; 

vth=sqrt(2*Te/me);
Clog0=14.9-log(ne/1E20)/2+log(Te/qe/1E3);

gammac=qe^4/(4*pi*ep0^2);

K_NB_model1=(KE0-gammac*ne*Clog0/(me*c)*time)/qe;

%% Plotting

fig=figure;
plot(time,KE)
hold on
plot(time,K_NB_model1)

legend({'FP','model'})
