clear all
clc

disp("Sprain detection simulation");

Datanum = 490;
data = ReadIMUFunction("imu_test1.txt",Datanum);
% data{1} = Gyro data{2} = Linear data{3} = etc
Gyro = data{1};
Linear = data{2};
etc = data{3};


Xlocal=[1,0,0]'; %sagittal vector
Zffhat=zeros(1,3);
i=1;
j=1;
k=100;
Datanum = 490;
oldRmatrix=eye(3);
pitchDATA=0;
rollDATA=0;

diff = 0;
oldtime= etc(i,1)-20;
oldnum = 0;
SwingNum = -1;
SwingCheck = 0;
IsSwing = -1;
IsSwingData =0;



delta=0
c=0;
d=0;

while(1)
 %vertical axis 구하기
    if (j<k)
        while(1)
            Zffhat=Zffhat+Linear(j,1:3);
            j=j+1;
            if (Linear(j,3)>=Linear(1,3)+0.015 || Linear(1,3)<=Linear(1,3)-0.015) %z 좌표축으로 변경하는게 좋을 듯
                k=j;
                break;
            end
        end
    end
Zff=Zffhat./sqrt(sum(Zffhat.^2));


%mediolateral axis 구하기
Yfhat=cross(Zff',Xlocal);
Yfleft=-Yfhat./sqrt(sum(Yfhat.^2));%left foot
Yright=+Yfhat./sqrt(sum(Yfhat.^2));%right foot

%anterposterior axis 구하기
Xhat=cross(cross(Zff',Xlocal),Zff');
Xf=Xhat./sqrt(sum(Xhat.^2));

%rotation matrix 구하기
    S=[0 -Gyro(j,3) -Gyro(j,2);Gyro(j,3) 0 -Gyro(j,1);-Gyro(j,2) Gyro(j,1) 0];
    timeinterval=(etc(j,1)-etc(j-1,1))*0.001;%msec
    A=S.*timeinterval;
    sigma=sqrt(Gyro(j,1)^2+Gyro(j,2)^2+Gyro(j,3)^2);%rad/sec
    Rmatrixhat=oldRmatrix*(eye(3)+A.*(sin(sigma)/sigma)+(A*A).*(((1-cos(sigma))/sigma^2)));%Rotationmatrix
    
    Rmatrix=normalize(Rmatrixhat);%normalization
    %Rmatrix=Rmatrixhat;

    pitch=asind(Zff*Rmatrix*Xf)
    roll=90-asind(Zff*Rmatrix*Yfleft);

    %swing check
    time=(etc(i,1)-oldtime);
    if(SwingCheck>0)
        SwingCheck=SwingCheck-time;
    end
    pause(time*0.001);

    diff(i) =  (Gyro(i,2)-oldnum)/ time;

    if(diff(i)>0.06 && SwingCheck<=0)
        IsSwing=IsSwing*(-1);
        SwingCheck=200;
    end
    SwingData(i)=IsSwing;

    %rotating matrix 누적 오류 보정
    if i>1
    if IsSwing==1 && SwingData(i-1)==-1
        c=etc(i,1);
    end

    if IsSwing==-1 && SwingData(i-1)==1
        d=etc(i,1);
    end
    end

    if d>c
        oldRmatrix=eye(3);
    else
        oldRmatrix=Rmatrix;
    end
    d=0;

    %DATA값
    pitchDATA(j)=pitch;
    rollDATA(j)=roll;
    

    if i<k
        figure(5),plot((etc(1:i,1)-etc(1,1)),[Gyro(1:i,2),zeros(i,1),zeros(i,1)] );
    else
        figure(5),plot((etc(1:i,1)-etc(1,1)),[Gyro(1:i,2),rollDATA',pitchDATA'] );%roll 빨간색,pitch 노란색
    end
    oldnum=Gyro(i,2);
    oldtime= etc(i,1);
    
    if i>k-1
        j=j+1;
    end    

    i=i+1;

    if(j>Datanum)
        break;
    end
end
