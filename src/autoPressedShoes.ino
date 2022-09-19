#include <Wire.h>
extern TwoWire Wire1;
#include "pinHeader.h"
#include "encoder.h"
#include <MatrixMath.h>
#include <math.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>

#define DEVELOPMENT 1

//차압센서 함수
void DataFetch_ISEN_P10K(void);
void RawToDecimal_ISEN_P10k(void);
void Calculate_ISEN_P10k(void);

//차압센서 변수
unsigned long press_decimal, temp_decimal;
float pressure, temperature;
unsigned char dat[4];

float predictValue[1000];
int CheckingStep = 0;
int stepstart = 0;
int step=0;

//matrix
mtx_type Zffhat[3]={0,0,0};
mtx_type Zff[3];

//BNO55 센서
Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);

//Step 분석
void readStep() {
    sensors_event_t angVelocityData, linearAccelData;
    bno.getEvent(&angVelocityData, Adafruit_BNO055::VECTOR_GYROSCOPE);
    bno.getEvent(&linearAccelData, Adafruit_BNO055::VECTOR_LINEARACCEL);

    float y = angVelocityData.gyro.y;

    float diff = (y - oldgyroY) / LOOP_DELAY;

//Swing 분석
    if (diff > 0.06 && swingCheck <= 0) {
        IsSwing = !IsSwing;
        swingCheck = 200;
        //paze 변경
        if (!IsSwing) {
            //조여주기
            target_pressure = highPressure;
            Serial.println("target_pressure is highPressure");
        }
        else {
            //풀어주기
            target_pressure = lowPressure;
            Serial.println("target_pressure is lowPressure");
        }
    }else if((y - oldgyroY)<0.005){
        Serial.print(y);
        Serial.println(" : 값이 변하지 않음.");
    }

    if(swingCheck>0)
        swingCheck-= LOOP_DELAY;

    oldgyroY = y;
//step 예측
/*  초기값
    float s_timeinterval={-1,-1,-1,-1,-1,-1};
    float s_time={-1,-1,-1,-1,-1};
    float v_high={1,1};
    float v_low=-2;
*/
}
//Sprain detect
void sprain(){
  sensors_event_t orientationData, accelerometerData,angVelocityData;
  bno.getEvent(&orientationData, Adafruit_BNO055::VECTOR_EULER);
  bno.getEvent(&angVelocityData, Adafruit_BNO055::VECTOR_GYROSCOPE);
  bno.getEvent(&accelerometerData, Adafruit_BNO055::VECTOR_ACCELEROMETER);


  float Linear_x= accelerometerData.acceleration.x;
  float Linear_y= accelerometerData.acceleration.y;
  float Linear_z= accelerometerData.acceleration.z;
  float Orient_z= orientationData.orientation.z;
  Gyro_integ+=(angVelocityData.gyro.y*50);

//find the z-axis of a flat foot
  if (ic<kc && ic>1){

    if (Orient_z >= (firstOrient_z + 500) || Orient_z <= (firstOrient_z - 500) ){
      res = 0;
      for (jc = 0; jc < 3; jc++){
        res += (Zffhat[jc] * Zffhat[jc]);
      }
      magnitude = sqrt(res);
      Zff[0]=Zffhat[0]/magnitude;
      Zff[1]=Zffhat[1]/magnitude;
      Zff[2]=Zffhat[2]/magnitude;
      ic=kc;

      ic+=1;
  }
    else{
      Zffhat[0]+=Linear_x;
      Zffhat[1]+=Linear_y;
      Zffhat[2]+=Linear_z;
      ic+=1;

  }

 }
   if(ic==1)
   {
    firstOrient_z= Orient_z;
    Zffhat[0]=Linear_x;
    Zffhat[1]=Linear_y;
    Zffhat[2]=Linear_z;
    ic+=1;
    }

    if (ic == kc){
      res = 0;
      for (jc = 0; jc < 3; jc++){
        res += (Zffhat[jc] * Zffhat[jc]);
      }
      magnitude = sqrt(res);
      Zff[0]=Zffhat[0]/magnitude;
      Zff[1]=Zffhat[1]/magnitude;
      Zff[2]=Zffhat[2]/magnitude;
      ic+=1;
    }

// find an angle between Zff and present gravity vector
    if(ic>kc){
      float pitch;
      float roll;
      pitch=acos((Linear_y*Zff[1]+Linear_z*Zff[2])/sqrt(Linear_y * Linear_y+Linear_z * Linear_z)/sqrt(Zff[1]*Zff[1]+Zff[2]*Zff[2]))*180/3.141592;
      roll=acos((Linear_x*Zff[0]+Linear_z*Zff[2])/sqrt(Linear_x * Linear_x+Linear_z * Linear_z)/sqrt(Zff[0]*Zff[0]+Zff[2]*Zff[2]))*180/3.141592;

      Serial.print("pitch:");
      Serial.print(pitch);
      Serial.print("  ");
      Serial.print("roll:");
      Serial.println(roll);


      if(pitch>inversion || (roll>dorsi&&Gyro_integ<0) || (roll>plantar&&Gyro_integ>0)){
        target_pressure = highPressure;
      }
      else{
        target_pressure = lowPressure;

      }


      ic+=1;

    }
    else{
        Serial.println("looking for the z-axis of a flat foot");
        target_pressure = lowPressure;
        Serial.println(ic);
      }

}
//루프마다 모터
void motorWork() {
    if (target_pressure - current_pressure != 0 && target_pressure != 0) {
        //모터의 방향 설정
        if (target_pressure - current_pressure > 0) {
            digitalWrite(MOTOR_DIR, HIGH);
#ifdef DEBUG
            Serial.print(" motor HIGH ");
#endif
        }else {
            digitalWrite(MOTOR_DIR, LOW);
#ifdef DEBUG
            Serial.print(" motor LOW ");
#endif
        }
        //모터가 반대방향으로 변경하는 것 방지
        if (motor_current_degree < 10 && target_pressure - current_pressure <= 0) {
            float motorspeed = 0;
            analogWrite(MOTOR_PWM, motorspeed);
            return;
        }
        //모터 PID코드
       float motorspeed = pControl + iControl + dControl;
       //Serial.println(motorspeed);
       //float motorspeed = 250;

        analogWrite(MOTOR_PWM, motorspeed);
        //디버그 모드 프린트
        Serial.print("pressure : ");
        Serial.print(target_pressure);
        Serial.print(" ");
        Serial.print(current_pressure);
        Serial.print(" ");
        Serial.print(motor_current_degree);
        Serial.print(" ");
        Serial.print(motorspeed);
        Serial.println(" ");

    }

}

//모터 리셋
void motorReset() {
  Serial.print("motor Reset");
  digitalWrite(MOTOR_DIR, LOW);
  while(motor_current_degree > 10){
    analogWrite(MOTOR_PWM, 100);
    Serial.print("motor_current_degree : ");
    Serial.println(motor_current_degree);
  }
  analogWrite(MOTOR_PWM, 0);
}

//셋업
void setup() {
    Wire.begin(); // I2C 초기화(Master Mode)
    Wire1.begin(); // I2C 초기화(Master Mode)
    Serial.begin(115200); // UART 초기화 9600bps

    pinMode(13, OUTPUT);
    pinMode(MOTOR_DIR, OUTPUT);
    pinMode(MOTOR_PWM, OUTPUT);
    pinMode(MOTOR_ENCODER_A, INPUT_PULLUP);
    pinMode(MOTOR_ENCODER_B, INPUT_PULLUP);
    attachInterrupt(digitalPinToInterrupt(MOTOR_ENCODER_A), doEncoderA, CHANGE);
    attachInterrupt(digitalPinToInterrupt(MOTOR_ENCODER_B), doEncoderB, CHANGE);
    Serial.println("startup the autoPressing");

    if (!bno.begin()) {
        /* There was a problem detecting the BNO055 ... check your connections */
        Serial.println("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    }

    if(DEVELOPMENT){
        lowPressure = 2900;
        highPressure = 3200;

    }else{
        //최저압 측정
        Serial.println("checking lowPressure");
        delay(1000);
        lowPressure = 0;
        for (int i = 0; i < 10; i++) {
            DataFetch_ISEN_P10K();
            RawToDecimal_ISEN_P10k();
            Calculate_ISEN_P10k();
            lowPressure += press_decimal;
            delay(100);
        }
        lowPressure = lowPressure / 10;
        Serial.print("lowPressure ");
        Serial.println(lowPressure);

        //최고압 측정
        Serial.println("checking highPressure");
        digitalWrite(MOTOR_DIR, HIGH);
        float motorspeed = 200;
        analogWrite(MOTOR_PWM, motorspeed);
        int tmp = 0;
        while (1) {
            unsigned long mfirst = millis();
            DataFetch_ISEN_P10K();
            RawToDecimal_ISEN_P10k();
            Calculate_ISEN_P10k();
            if (press_decimal == highPressure && tmp > 100) break;
            highPressure = press_decimal;
            unsigned long msecond = millis();
            unsigned long code_time = msecond - mfirst;
            tmp++;
            delay(LOOP_DELAY - code_time);
        }
        Serial.print("highPressure ");
        Serial.println(highPressure);

        //다시 최저압으로 변경
        Serial.println("checking complete");
        tmp = 0;
        target_pressure = lowPressure;
        while (1) {
            unsigned long mfirst = millis();
            DataFetch_ISEN_P10K();
            RawToDecimal_ISEN_P10k();
            Calculate_ISEN_P10k();
            current_pressure = press_decimal;
            motorWork();
            if (abs(target_pressure - current_pressure) < 5 && tmp > 20) break;
            unsigned long msecond = millis();
            unsigned long code_time = msecond - mfirst;
            tmp++;
            delay(LOOP_DELAY - code_time);
        }
    }
    Serial.println("start autopressing");
}

//압력센서 값을 받아 걸음걸이 분석
void checkingPressureSensor(int val) {
    if (val >= SENSOR_CRITICAL_VALUE) {//압력센서가 일정한 값 이상일때 작동
        if (preMillis != 0) {
            unsigned long millidelay = millis() - preMillis;
            if (millidelay > 10000) {
                Serial.println("걸음걸이에 오류가 있습니다.");
            }
            else {
                distance_walking = millis() - preMillis;
                Serial.print("걸음걸이 간격 시간은 ");
                Serial.println(distance_walking);
                walking_percentage = distance_walking;
            }
        }
        sensorDelay = 20;
        preMillis = millis();
    }
}

//루프
void loop() {
    unsigned long mfirst = millis();

    //motor reset
   /*
    while (Serial.available() > 0) {
      int type = Serial.parseInt();
      Serial.println(type);
      if(type==1){
        motorReset();
        while (true) {
          delay(1000);
          Serial.println("delaying");
        }
      }
    } */

    /*압력센서 쓰면 쓰기
    //현재 압력값 받아오기
    int val = 1023- analogRead(A0);
    Serial.print(" Sensor : ");
    Serial.println(val);
    if(sensorDelay<=0){
      checkingPressureSensor(val);
    }else{
      sensorDelay--;
      if(sensorDelay<0) sensorDelay=0;
    }
    */

    //readStep();
    sprain();

    //모터 함수
    //motorWork();


    //현재 압력값 받아오기

    DataFetch_ISEN_P10K();
    RawToDecimal_ISEN_P10k();
    Calculate_ISEN_P10k();

    current_pressure = press_decimal;
    Serial.println(current_pressure);


    //PID 값 알아내기 위한 코드

    errorGap = target_pressure - current_pressure - realError;
    realError = target_pressure - current_pressure;
    accError = accError + realError;

    pControl = mP * realError;
    iControl = mI * (accError * LOOP_DELAY);  // 전에 에러의 합과 지금 에러를 더해서 딜레이를 곱해서 더함
    dControl = mD * (errorGap/ LOOP_DELAY);  // 전에 error와 지금 error의 차이를 딜레이로 나눔

    unsigned long msecond = millis();
    unsigned long code_time = msecond - mfirst;
    delay(LOOP_DELAY - code_time);
}

//--------------------------------------------------------------------------------------------------------------------
void printEvent(sensors_event_t* event) {
    double x = -1000000, y = -1000000, z = -1000000; //dumb values, easy to spot problem
    if (event->type == SENSOR_TYPE_ACCELEROMETER) {
        Serial.print("Accl:");
        x = event->acceleration.x;
        y = event->acceleration.y;
        z = event->acceleration.z;
    }
    else if (event->type == SENSOR_TYPE_ORIENTATION) {
        Serial.print("Orient:");
        x = event->orientation.x;
        y = event->orientation.y;
        z = event->orientation.z;
    }
    else if (event->type == SENSOR_TYPE_MAGNETIC_FIELD) {
        Serial.print("Mag:");
        x = event->magnetic.x;
        y = event->magnetic.y;
        z = event->magnetic.z;
    }
    else if (event->type == SENSOR_TYPE_GYROSCOPE) {
        Serial.print("Gyro:");
        x = event->gyro.x;
        y = event->gyro.y;
        z = event->gyro.z;
    }
    else if (event->type == SENSOR_TYPE_ROTATION_VECTOR) {
        Serial.print("Rot:");
        x = event->gyro.x;
        y = event->gyro.y;
        z = event->gyro.z;
    }
    else if (event->type == SENSOR_TYPE_LINEAR_ACCELERATION) {
        Serial.print("Linear:");
        x = event->acceleration.x;
        y = event->acceleration.y;
        z = event->acceleration.z;
    }
    else if (event->type == SENSOR_TYPE_GRAVITY) {
        Serial.print("Gravity:");
        x = event->acceleration.x;
        y = event->acceleration.y;
        z = event->acceleration.z;
    }
    else {
        Serial.print("Unk:");
    }

    Serial.print("\tx= ");
    Serial.print(x);
    Serial.print(" |\ty= ");
    Serial.print(y);
    Serial.print(" |\tz= ");
    Serial.println(z);
}
//--------------------------------------------------------------------------------------------------------------------
// 압력센서 모듈에서 압력 및 온도 Raw Data를 Fetch
void DataFetch_ISEN_P10K(void) {
    int i = 0;
    /****************************************************************************
    * 압력센서는 특별한 데이타 요청이 없을 때 소모전력을 최소화 하기 위해서 Sleep 상태로 들어감.
    * 센서를 Wake up 하기 위해서는 두 가지 명령이 존재.
    * * 1. Read MR -> slave address 0x28 + I2C read command bit를 적용하여 Wake up
    * --> 압력 및 온도 데이터를 Read할때 마다 update 됨.
    * 2. Write MR -> slave address 0x28 + I2C write command bit를 적용햐여 Wake up
    * --> 압력 값만 Update 되고 온도 값은 이전 값을 유지함.
    * 압력센서 Application 에 따라 Wake up 명령을 적절해 사용할 필요 있음.
    ****************************************************************************/
    Wire1.write(P10K_i2c_address + 0x01); // I2C Read_MR
    //Wire.write(P10K_i2c_address); // I2C Write_MR
    Wire1.endTransmission();
    delay(10);
    Wire1.requestFrom(0x28, 4);
    while (Wire1.available()) {
        if (i >= 4) i = 0;
        dat[i] = Wire1.read();
        i++;
    }
}
/******************************************************************
// -.I2C라인을 통해 읽어들인 압력센서의 Raw 값을 decimal 값으로 Bit Arrange
// -.pressure : 2^14(16384),temperature : 2^11(2048)
*******************************************************************/
void RawToDecimal_ISEN_P10k(void) {
    press_decimal = (((dat[0] & 0x3F) << 8) | dat[1]);
    temp_decimal = ((dat[2] << 3) | (dat[3] >> 5));
}
/*******************************************************************
// Decimal 값을 10kPa 범위의 압력값으로 환산
// - 압력[Pa]=10000*(decimal값-0Pa일때 Offset)(16383.0-Offset), Offset=3000
// - 온도[C]=(200*decimal값/2047)-50, 온도범위: -50~+150[C]
********************************************************************/
void Calculate_ISEN_P10k(void) {
    pressure = (float)((press_decimal - 3000.0) / (16383.0 - 3000)) * 10000;
    if (pressure < 0) pressure = 0.0; // 0Pa에서 음수 값이면 ->0
    if (pressure > 10000.0) pressure = 10000.0; // 10kPa 이상 값이면 -> 10,000 Truncation
    temperature = 200.0 * (float)(temp_decimal / 2047.0) - 50.0;
}
