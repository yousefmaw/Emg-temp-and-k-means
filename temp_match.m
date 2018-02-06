function [tempedmuap,muaps] = temp_match( datadir )
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
%%%%%%%%%%%% temp matching
muaps(:,1)=reaps(:,1);D=zeros(20,nap-1);tempmuap=zeros(1,nap-1);
D1=zeros(1,nap-1);D2=zeros(1,nap-1);j=1;D3=zeros(1,nap-1);
while(sum(muaps(:,j))~=0)
    
for i=1:1:nap-1
D(:,i)=reaps(:,i)-muaps(:,j);
end
D = D.^2;
D1=sum(D);
D1 =(D1<=(12.65^5));
%D1used for thersholding
tempmuap=tempmuap+j*D1; %to retrieve the temp for each muap
[M,I]=min(tempmuap); %to detect the next unmatched muap
muaps=[muaps reaps(:,I)];
%zeroing muaps that belong to the detected template
D1=~D1;
for k=1:1:20
reaps(k,:)=reaps(k,:).*D1;
end
j=j+1;
end
%%%allocating a temp for each muap
ntemps= size(muaps,2)-1;tempedmuap=zeros(ntemps,nap-1);
muaps=muaps(:,1:ntemps);
timestamps=transpose(timestamps(1:nap-1,1));
tempmuap=tempmuap+1;
for i= ntemps+1:-1:2
tempedmuap(i,:)=timestamps.*(~mod(tempmuap,i));
end
tempedmuap=tempedmuap(2:ntemps+1,:);
%%%%%%%%stats
colors=['b' 'r' 'k' 'g' 'y' 'm' 'c' 'w'];
figure();
plot(muaps);
title('detected ap in the data');
xlabel('# sample');
ylabel('voltage');
%%%
figure();
plot((30000:35000),data(30000:35000),'g');
hold on
for i=1:1:ntemps
plot(tempedmuap(i,apstart:apend),1000,'color',colors(i),'marker','*');
hold on
end
title('muap');
xlabel('#samples');
ylabel('emg');
xlim([30000 35000])
hold off
fclose(emg);
end

