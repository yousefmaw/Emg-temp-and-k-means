function [timestamps,muaps] = K_means( datadir, nmeans)
emg=fopen(datadir);
data = fscanf(emg,'%f');
nsamples=size(data,1);
absdata=zeros(nsamples,1);
for i=1:1:nsamples
    absdata(i)=abs(data(i));
end
avrdata=zeros(nsamples,1);
for i = 1:1:nsamples
    
    if(i<20+1)
        for j = 20 : -1: 20-i+1
        avrdata(i)= avrdata(i)+ absdata(i-(20-j));
        end
         avrdata(i)=(1/i)*avrdata(i);
    else
        for j = 1:1:20
         avrdata(i)= avrdata(i)+ absdata(i-(20-j));    
        end
         avrdata(i)=(1/20)*avrdata(i);
    end
        
end
%%%%% threshold value calculations
threshold=3*std(absdata(500:800));
%%% range from 500:800 determined visually 
%%%% peak detection

aps= zeros(nsamples,1);nap=1; apstart=0;apend=0;
timestamps= zeros(nsamples,1);flag2=0;flag1=0;i=1;
%%action potentials detection
while (i<=nsamples-20)
    for j = 0:1:20
        if(avrdata(i+j)>threshold)
            timestamps(nap)= i+j;
            aps(nap)=data(i+j);
            nap=nap+1;
            i=i+j+20;
            break;
        end  
    end
    if (flag1==0)
       if i > 30000
           apstart=nap;
           flag1=1;
       end
    end
    if (flag2==0)
       if i > 35000
           apend=nap;
           flag2=1;
       end
    end
   i=i+20;
end
%%%%centring
for i=1:1:nap-1
    [M,I]=max(transpose(data(timestamps(i):timestamps(i)+20)));
    timestamps(i)=I+timestamps(i)-1;
    aps(i)=M;
end
muaps=zeros(20,1);
reaps=zeros(20,nap-1);
for i=1:1:nap-1
    reaps(:,i)=data(timestamps(i)-10:timestamps(i)+9);
end
%%%%%%%%%%%% K-means
knew=rand(20,nmeans);kold=knew;
r=zeros(nmeans,nap-1);summ=zeros(1,3);error=5;
while(error >(0.000000000001))
    kold=knew;
    for i=1:1:nmeans
       r(i,:) =sqrt(sum((reaps- repmat(knew(:,i), [1,nap-1])).^2));
    end
    for i=1:1:nap-1
      [M I]= min(r(:,i));
      r(:,i)=0;
      r(I,i)=1;
    end
     knew=reaps*transpose(r);
    summ=sum(transpose(r));
    for i=1:1:nmeans
    knew(:,i)=knew(:,i)./summ(1,i);
    end
   error=sum(sum((knew-kold).^2));
end
timestamps=transpose(timestamps(1:nap-1,1));
timestamps=repmat(timestamps,[nmeans,1]);
timestamps = r.*timestamps;
%%%%%%%%stats
colors=['b' 'r' 'k' 'g' 'y' 'm' 'c' 'w'];
figure();
plot(knew);
title('detected ap in the data');
xlabel('# sample');
ylabel('voltage');
%%%
figure();
plot((30000:35000),data(30000:35000),'g');
hold on
for i=1:1:nmeans
plot(timestamps(i,apstart:apend),1000,'color',colors(i),'marker','*');
hold on
end
title('muap');
xlabel('#samples');
ylabel('emg');
xlim([30000 35000])
hold off
fclose(emg);
end

