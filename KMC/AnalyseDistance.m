

dist=load('distance.dat');
tck=load('tcheck.dat');

step=10
T=283
n=floor(size(dist,1)/step);

matrices=zeros(3*(n),3);
matrices_V=zeros(3*(n),3)
avedir=zeros(3);
aveeig=zeros(3);

for i=1:n
    
    sta=(i-1)*step+1;
    
    fini=i*step;
    
    mobmat=cov(dist(sta:fini,:))*(1e-8)^2/(2*mean(tck(sta:fini))*8.617343e-5*T);
    
     [V,D]=eig(mobmat);
    
    matrices(((3*i)-2):(3*i),:) = D; %[mobmat(1,1),mobmat(2,2),mobmat(3,3), abs(onmin/offmax)];
    matrices_V(((3*i)-2):(3*i),:) = V;
    

    avedir=avedir+V;
    aveeig=aveeig+D;
    
end    
    mean(matrices)
    
    
    %
    
    dlmwrite('g.txt',matrices,'delimiter', '\t');
    dlmwrite('g_V.txt',matrices_V,'delimiter', '\t')

    
    
    mx=matrices(1:3:end,1);
    my=matrices(2:3:end,2);
    mz=matrices(3:3:end,3);
   m=[mean(mx),mean(my),mean(mz)]
