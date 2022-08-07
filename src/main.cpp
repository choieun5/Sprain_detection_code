#include <Vector.h>
#include <MatrixMath.h>
#include <math.h>
#define N 3
float LOOP_DELAY=50;

void setup(){
 Serial.begin(9600);
}

void loop(){
  sprainDetect();
  delay(2000);
}

//sprain detection function
  int icount=1;
  int kcount=100;
  int jcount;
  int lcount;
  float firstLinear_z=0;
  float res = 0;
  float magnitude=0;
  float oldswing=0;//이거 위로 빼기
  mtx_type Zffhat[N];
  mtx_type Zff[N];
  mtx_type Yfhat[N];
  mtx_type Yfleft[N];
  mtx_type Yfright[N];
  mtx_type Xfhat[N];
  mtx_type Xf[N];
  mtx_type oldRmatrix[N][N]={1,0,0,0,1,0,0,0,1};
  mtx_type Rmatrix[N][N]{0,0,0,0,0,0,0,0,0};



void sprainDetect(){

  sensors_event_t angVelocityData, linearAccelData;
  bno.getEvent(&angVelocityData, Adafruit_BNO055::VECTOR_GYROSCOPE);
  bno.getEvent(&linearAccelData, Adafruit_BNO055::VECTOR_LINEARACCEL);

  float Gyro_x= angVelocityData.gyro.x;
  float Gyro_y= angVelocityData.gyro.y;
  float Gyro_z= angVelocityData.gyro.z;

  float Linear_x= linearAccelData.acceleration.x;
  float Linear_y= linearAccelData.acceleration.y;
  float Linear_z= linearAccelData.acceleration.z;
  
/*
  float Gyro_x=1.5;
  float Gyro_y=1.2;
  float Gyro_z=-2.3;
  float Linear_x=3;
  float Linear_y=5;
  float Linear_z=9;
  */
  
   //vertical axis
  if (icount<kcount && icount>1){
    
    if (Linear_z >= (firstLinear_z + 0.015) || Linear_z <= (firstLinear_z - 0.015) ){
      res = 0;
      for (jcount = 0; jcount < 3; jcount++){
        res += Zffhat[jcount] * Zffhat[jcount];
      }
      magnitude = sqrt(res);
      Zff[0]=Zffhat[0]/magnitude;
      Zff[1]=Zffhat[1]/magnitude;
      Zff[2]=Zffhat[2]/magnitude;
      kcount=icount;
  }
    else{
      Zffhat[0]+=Linear_x;
      Zffhat[1]+=Linear_y;
      Zffhat[2]+=Linear_z;
      icount+=1;
  }
  
 }

   if(icount=1){
    firstLinear_z= Linear_z;
    Zffhat[0]=Linear_x;
    Zffhat[1]=Linear_y;
    Zffhat[2]=Linear_z;
    icount+=1;
    }
    
   //mediolateral axis & anterposterior axis

  if (icount = kcount){
      Yfhat[0]= Zff[1] * 0 - Zff[2] * 0;
      Yfhat[1]= Zff[2] * 1 - Zff[0] * 0;
      Yfhat[2]= Zff[0] * 0 - Zff[1] * 1;
      
      res = 0;
      for (jcount = 0; jcount < 3; jcount++){
        res += Yfhat[jcount] * Yfhat[jcount];
      }
      magnitude = sqrt(res);
        
      Yfleft[0]=-Yfhat[0]/magnitude;
      Yfleft[1]=-Yfhat[1]/magnitude;
      Yfleft[2]=-Yfhat[2]/magnitude;
      Yfright[0]=Yfhat[0]/magnitude;
      Yfright[1]=Yfhat[1]/magnitude;
      Yfright[2]=Yfhat[2]/magnitude;
     
      Xfhat[0]= Yfhat[1] * Zff[2] - Yfhat[2] * Zff[1];
      Xfhat[1]= Yfhat[2] * Zff[0] - Yfhat[0] * Zff[2];
      Xfhat[2]= Yfhat[0] * Zff[1] - Yfhat[1] * Zff[0];

      res = 0;
      for (jcount = 0; jcount < 3; jcount++){
        res += Xfhat[jcount] * Xfhat[jcount];
      }
      magnitude = sqrt(res);
      
      Xf[0]=Xfhat[0]/magnitude;
      Xf[1]=Xfhat[1]/magnitude;
      Xf[2]=Xfhat[2]/magnitude;

      icount+=1;
    }

    //Rotation matrix

    
      if(icount>kcount){
        mtx_type S[N][N]={0,-Gyro_z,-Gyro_y,Gyro_z,0,-Gyro_x,-Gyro_y,Gyro_x,0};
        mtx_type A[N][N];
        mtx_type B[N][N];
        mtx_type C[N][N];
        mtx_type D[N][N];
        mtx_type Rmatrixhat[N][N]{0,0,0,0,0,0,0,0,0};
        mtx_type eye[N][N]={1,0,0,0,1,0,0,0,1};
        mtx_type AA[N][N]={0,0,0,0,0,0,0,0,0};
        float sigma= sqrt(Gyro_x*Gyro_x+Gyro_y*Gyro_y+Gyro_z*Gyro_z);

      for (jcount = 0; jcount < N; jcount++){
        for (lcount = 0; lcount <N; lcount++){
          A[jcount][lcount]=(S[jcount][lcount]*50*0.001);//50=loop_delay
        }
      }
       
      Matrix.Multiply((mtx_type*)A, (mtx_type*)A, N, N, N, (mtx_type*)AA);


      for (jcount = 0; jcount < N; jcount++){
        for (lcount = 0; lcount <N; lcount++){
          B[jcount][lcount]=A[jcount][lcount]*sin(sigma)/sigma;
          C[jcount][lcount]=AA[jcount][lcount]*(1-cos(sigma))/sigma/sigma;
        }
      }

      Matrix.Add((mtx_type*) B, (mtx_type*) C, N, N, (mtx_type*) C);
      Matrix.Add((mtx_type*) eye, (mtx_type*) C, N, N, (mtx_type*) D); //숫자가 작아서 0으로 나와요~
      Matrix.Multiply((mtx_type*)oldRmatrix, (mtx_type*)D, N, N, N, (mtx_type*)Rmatrixhat);
      
//normalize
       for (jcount = 0; jcount < N; jcount++){
        float mid=(Rmatrixhat[0][jcount]+Rmatrixhat[1][jcount]+Rmatrixhat[2][jcount])/3;
        float norm=sqrt((Rmatrixhat[0][jcount]*Rmatrixhat[0][jcount]+Rmatrixhat[1][jcount]*Rmatrixhat[1][jcount]+Rmatrixhat[2][jcount]*Rmatrixhat[2][jcount])/2);
        for (lcount = 0; lcount <N; lcount++){
          Rmatrix[lcount][jcount]=(Rmatrixhat[lcount][jcount]-mid)/norm;
        }
      }
      
      mtx_type W[N];
      mtx_type X[N];
      Matrix.Multiply((mtx_type*)Rmatrix, (mtx_type*)Xf, N, N, 1, (mtx_type*)W);
      Matrix.Multiply((mtx_type*)Rmatrix, (mtx_type*)Yfleft, N, N, 1, (mtx_type*)X);
      float x =0, y=0;
      for (jcount = 0; jcount < N; jcount++){
        x+=Zff[jcount]*W[jcount];
        y+=Zff[jcount]*X[jcount];
      }
      
      float pitch= asin(x)*180/PI;
      float roll= 90-asin(y)*180/PI;
      
      
float IsSwing=0; // 지울거임     
//누적오류보정
    if (icount>1){
      if ((IsSwing-oldswing)==1){
        Matrix.Copy((mtx_type*)eye, N, N, (mtx_type*)oldRmatrix);
      }
      else{
        Matrix.Copy((mtx_type*)Rmatrix, N, N, (mtx_type*)oldRmatrix);
      }
      oldswing=IsSwing;
    }

      }
      
      else{
        Serial.println("looking for a vertical axis");
      }
      

  
}