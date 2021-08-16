#include "time.h"
#include "ACMSim.h"
#include "controller.h"
#define VVVF_CONTROL 0

struct PMSMSimulated PM;
struct ControllerForExperiment CTRL;
struct InductionMachine im;

//DTC
double HysPhi;
double HysTem;
double DeltaPhi;
double DeltaTem;

int HysPhicheck;
int HysTemcheck;
int HysPhicheckback;
int HysTemcheckback;

int checkx;
int checky;

double freq; // 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）

//定子电压查表
  int vvect[6][6] = {{2,3,4,5,6,1},
                     {7,0,7,0,7,0},
                     {6,1,2,3,4,5},
                     {3,4,5,6,1,2},
                     {0,7,0,7,0,7},
                     {5,6,1,2,3,4}};//电压向量

 //  int vvect[6][6] = {{6,2,3,1,5,4},
 //                     {7,0,7,0,7,0},
 //                     {5,4,6,2,3,1},
 //                     {2,3,1,5,4,6},
 //                     {0,7,0,7,0,7},
 //                     {1,5,4,6,2,3}};//电压向量

double volt;
double theta;

void CTRL_DTC();
void test();


void PM_init(){
    int i;
    for(i=0;i<5;++i){
        PM.x[i] = 0.0;
    }
    PM.rpm = 0.0;

    PM.iq = 0.0;
    PM.id = 0.0;

    PM.Tload = 0.0;
    PM.rpm_cmd = 0.0;

    PM.R  = 0.45;
    PM.Ld = 4.15*1e-3;
    PM.Lq = 16.74*1e-3;
    PM.KE = 0.504; // Vs/rad
    PM.L0 = 0.5*(PM.Ld + PM.Lq);
    PM.L1 = 0.5*(PM.Ld - PM.Lq);

    PM.Js = 0.06; 
    PM.npp = 2;
    PM.mu_m = PM.npp/PM.Js; //惯量系数

    PM.Ts  = PM_TS;

    PM.ual = 0.0;
    PM.ube = 0.0;

}

void rK5_dynamics(double t, double *x, double *fx){
        // electromagnetic model
        fx[0] = (PM.ud - PM.R * x[0] + x[2]*PM.Lq*x[1]) / PM.Ld;   //id电流状态方程
        fx[1] = (PM.uq - PM.R * x[1] - x[2]*PM.Ld*x[0] - x[2]*PM.KE) / PM.Lq;  //iq电流状态方程

        // mechanical model
        PM.Tem = PM.npp*(x[1]*PM.KE + (PM.Ld - PM.Lq)*x[0]*x[1]);    //转矩
        fx[2] = (PM.Tem - PM.Tload)*PM.mu_m; // elec. angular rotor speed
        fx[3] = x[2];                           // elec. angular rotor position
}
void rK555_Lin(double t, double *x, double hs){
    double k1[4], k2[4], k3[4], k4[4], xk[4];
    double fx[4];
    int i;

    rK5_dynamics(t, x, fx); // timer.t,
    for(i=0;i<4;++i){        
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<4;++i){        
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<4;++i){        
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs, 
    for(i=0;i<4;++i){        
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
    }
}
int machine_simulation(){
    rK555_Lin(CTRL.timebase, PM.x, PM.Ts);

    PM.id = PM.x[0];
    PM.iq = PM.x[1];
    PM.rpm = PM.x[2] * 60 / (2 * M_PI * PM.npp);

    PM.ial = MT2A(PM.id, PM.iq, cos(PM.theta_d), sin(PM.theta_d)); //DQ反变化
    PM.ibe = MT2B(PM.id, PM.iq, cos(PM.theta_d), sin(PM.theta_d)); //DQ反变化

    if(isNumber(PM.rpm)) //过速方跑飞
        return false;
    else
        return true;
}
void measurement(){
    US_C(0) = CTRL.ual;
    US_C(1) = CTRL.ube;
    US_P(0) = US_C(0);
    US_P(1) = US_C(1);

    IS_C(0) = PM.ial;
    IS_C(1) = PM.ibe;

    PM.omg = PM.x[2];
    PM.theta_d = PM.x[3];     //转子磁链角度
    PM.theta_r = PM.theta_d; //转子位置=定子磁链 磁场定向

}
void inverter_model(){
    PM.ual = CTRL.ual;
    PM.ube = CTRL.ube;

    PM.ud = AB2M(PM.ual, PM.ube, cos(PM.theta_d), sin(PM.theta_d)); //DQ变换，印加电压
    PM.uq = AB2T(PM.ual, PM.ube, cos(PM.theta_d), sin(PM.theta_d)); //DQ变换，印加电压
}
int main(){
    printf("NUMBER_OF_LINES: %d\n\n", NUMBER_OF_LINES);

    /* Initialization */
    PM_init();
    CTRL_init();

    FILE *fw;
    fw = fopen("algorithm.dat", "w");

    /* MAIN LOOP */
    clock_t  begin, end;
    begin = clock();   
    int _; // _ for the outer iteration
    int dfe=0; // dfe for down frequency execution
    for(_=0;_<NUMBER_OF_LINES;++_){

        /* Command and Load Torque */
        PM.rpm_cmd = 1500;
        PM.Tload = 30;//CTRL.Tem_cmd;//0.1;

        /* SPMulated PM */
        if(machine_simulation()){ 
            printf("Break the loop.\n");
            break;
        } 

        if(++dfe==DOWN_FREQ_EXE){
            dfe = 0;

            /* time */
            CTRL.timebase += TS;

            measurement();

            // observation();

            write_data_to_file(fw);

            #if VVVF_CONTROL == TRUE
                #define VF_RATIO 18 //18.0 // 8 ~ 18 shows saturated phenomenon
                double freq = 2; // 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
                double volt = VF_RATIO*freq;
                CTRL.ual = volt*cos(2*M_PI*freq*CTRL.timebase);
                CTRL.ube = volt*sin(2*M_PI*freq*CTRL.timebase);
            #else
                //freq = 50;
                volt = CTRL.phi_cmd; //励磁电压为磁链的5倍，PM磁石影响
                //CTRL.vecter_cmd = 6;000
                CTRL_DTC();
            #endif
        }

        inverter_model();
    }
    end = clock(); printf("The sPMulation in C costs %g sec.\n", (double)(end - begin)/CLOCKS_PER_SEC);
    fclose(fw);

    /* Fade out */
    system("python ./ACMPlot.py"); 
    // getch();
    // system("pause");
    // system("exit");
    return 0; 
}


/* Utility */
void write_data_to_file(FILE *fw){
    static int j=0,jj=0; // j,jj for down sampling

    // if(CTRL.timebase>20)
    {
    if(++j == 10)
    {
        j=0;
        fprintf(fw, "%g,%g,%g,%g,%g,%d,%d\n",
                PM.x[0], PM.x[1],CTRL.phial,CTRL.phibe,PM.Tem,HysPhicheck,HysTemcheck
                );
    }
    }
}

bool isNumber(double x){
    // This looks like it should always be true, 
    // but it's false if x is a NaN (1.#QNAN0).
    return (x == x); 
    // see https://www.johndcook.com/blog/IEEE_exceptions_in_cpp/ cb: https://stackoverflow.com/questions/347920/what-do-1-inf00-1-ind00-and-1-ind-mean
}

/* Initialization */
void CTRL_init(){
    int i=0,j=0;

    CTRL.timebase = 0.0;

        /* Parameter (including speed) Adaptation */ 
        CTRL.R     = PM.R;
        CTRL.Ld = PM.Ld;
        CTRL.Lq  = PM.Lq;
        CTRL.Js     = PM.Js;
        CTRL.Js_inv = 1.0/PM.Js;    

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    CTRL.rpm_cmd = 0.0;

    CTRL.phial = 0.0;
    CTRL.phibe = 0.0;

    CTRL.Tem = 0.0;
    CTRL.phi = 0.0;


    CTRL.theta_M = 0.0;

    CTRL.vecter_cmd = 0.0;

    checky = 0;
    checkx = 0;

    HysPhi = 0.0;
    HysTem = 0.0;

    HysPhicheck = 0;
    HysTemcheck = 0;

    DeltaPhi = 0.0;
    DeltaTem = 0.0;

}

void CTRL_DTC()
{
//指令
    CTRL.phi_cmd = 0.1;
    CTRL.Tem_cmd = 30; //直接转矩控制，无负载空转，反馈的转矩始终小于指令转矩
    //CTRL.Tem_cmd = PM.rpm_cmd - PM.rpm; //速度控制正常

//反馈电压电流获取

   CTRL.ual_fb = PM.ual; 
   CTRL.ube_fb = PM.ube; 
   CTRL.ial_fb = PM.x[0]; //d轴电流，FOC准确的情况
   CTRL.ibe_fb = PM.x[1]; //q轴电流，FOC准确的情况


//定子磁链以及电磁转矩估算
   CTRL.phial += (PM.ual - CTRL.R*PM.ial)*TS; //定子α磁链估算 （电压-电流损耗）*时间常数积分
   CTRL.phibe += (PM.ube - CTRL.R*PM.ibe)*TS; //定子β磁链估算 （电压-电流损耗）*时间常数积分

   CTRL.theta_M += PM.x[4]*TS;  //theta估算 测试用

   CTRL.phi = sqrt(CTRL.phial*CTRL.phial + CTRL.phibe*CTRL.phibe);  //定子总磁链估算
   CTRL.Tem = PM.npp*(PM.x[1]*PM.KE + (PM.Ld - PM.Lq)*PM.x[0]*PM.x[1]);  //转矩估算

   if (atan2(CTRL.phibe,CTRL.phial) >= 0)
   {
        CTRL.theta = atan2(CTRL.phibe,CTRL.phial); //获得磁链角度
   //     theta=0;
   }
   else 
   {
        CTRL.theta = 2*M_PI + atan2(CTRL.phibe,CTRL.phial); //获得磁链角度
    //    theta=1;
   }
   

//Hysteresis控制
   HysPhi = CTRL.phi_cmd - CTRL.phi;
   HysTem = CTRL.Tem_cmd - CTRL.Tem;

   
   DeltaPhi = CTRL.phi_cmd * TS;  //滞环比较器磁链调整分辨率，实际加载的电压向量微分
   DeltaTem = CTRL.Tem_cmd * TS;  //滞环比较器转矩调整分辨率，转矩指令的微分

//转码 滞环域为一个指令所能控制的范围
//电流滞环
    if(HysPhi>=(DeltaPhi))  
    {
        HysPhicheck = 0;
    }
    if(HysPhi<=(-DeltaPhi))
    {
        HysPhicheck = 3;   
    }
    if(((-DeltaPhi)<HysPhi)&&(HysPhi<(DeltaPhi)))
    {
        HysPhicheck = HysPhicheck;
    }
    

//转矩滞环 先进行状态判断
    if(HysTem > (DeltaTem))
    {
        HysTemcheck = 0;
    }
    if (HysTem < (-DeltaTem))
    {
        HysTemcheck = 2;
    }
    if( HysTem == 0 )
    {
        HysTemcheck = 1;
    }
    if( (DeltaTem) < HysTem < (-DeltaTem))
    {
        HysTemcheck = 1;
    }


//转矩滞环内 维持状态
    if(((-DeltaTem) <= HysTem) &&(HysTem < 0)&&(HysTemcheckback == 0))
    {
        HysTemcheck = 1;
    }  
    if(((-DeltaTem) <= HysTem) &&(HysTem < 0)&&(HysTemcheckback == 2))
    {
        HysTemcheck = 2;
    } 

    if((0 < HysTem) && (HysTem <= DeltaTem)&&(HysTemcheckback == 0))
    {
        HysTemcheck = 0;
    }
    if((0 < HysTem) && (HysTem <= DeltaTem)&&(HysTemcheckback == 2))
    {
        HysTemcheck = 1;
    }


//保存前回值
    HysTemcheckback = HysTemcheck;

//寻找扇区
    if((0<=CTRL.theta)&&(CTRL.theta<(60*M_PI/180))) //第一扇区
    {
        checkx = 0;
    }
    if(((60*M_PI/180)<=CTRL.theta)&&(CTRL.theta<(120*M_PI/180)))
    {
        checkx = 1;
    }
    if(((120*M_PI/180)<=CTRL.theta)&&(CTRL.theta<(180*M_PI/180)))
    {
        checkx = 2;
    }
    if(((180*M_PI/180)<=CTRL.theta)&&(CTRL.theta<(240*M_PI/180))) 
    {
        checkx = 3;
    }
    if(((240*M_PI/180)<=CTRL.theta)&&(CTRL.theta<(300*M_PI/180)))
    {
        checkx = 4;
    }
    if(((300*M_PI/180)<=CTRL.theta)&&(CTRL.theta<=(360*M_PI/180)))
    {
        checkx = 5;
    }

//寻表
    checky = HysPhicheck + HysTemcheck;
    CTRL.vecter_cmd = vvect[checky][checkx];



    
//电压矢量V0~V7
    if(CTRL.vecter_cmd == 0)
    {
        CTRL.ual = 0;
        CTRL.ube = 0;
    }
    if(CTRL.vecter_cmd == 1)
    {
        CTRL.ual = volt*cos(0);
        CTRL.ube = volt*sin(0);
    }
    if(CTRL.vecter_cmd == 2)
    {
        CTRL.ual = volt*cos(60*M_PI/180);
        CTRL.ube = volt*sin(60*M_PI/180);
    }
    if(CTRL.vecter_cmd == 3)
    {
        CTRL.ual = volt*cos(120*M_PI/180);
        CTRL.ube = volt*sin(120*M_PI/180);
    }
    if(CTRL.vecter_cmd == 4)
    {
        CTRL.ual = volt*cos(180*M_PI/180);
        CTRL.ube = volt*sin(180*M_PI/180);
    }
    if(CTRL.vecter_cmd == 5)
    {
        CTRL.ual = volt*cos(240*M_PI/180);
        CTRL.ube = volt*sin(240*M_PI/180);
    }
    if(CTRL.vecter_cmd == 6)
    {
        CTRL.ual = volt*cos(300*M_PI/180);
        CTRL.ube = volt*sin(300*M_PI/180);
    }
    if(CTRL.vecter_cmd == 7)
    {
        CTRL.ual = 0;
        CTRL.ube = 0;
    }

    //CTRL.ual = volt*cos(2*M_PI*freq*CTRL.timebase);
    //CTRL.ube = volt*sin(2*M_PI*freq*CTRL.timebase);

    //CTRL.vecter_cmd += 1;
    //test
    //if(CTRL.vecter_cmd >= 8)
    //{
    //   CTRL.vecter_cmd = 0;
    //}
//
//印加在电机上vecter_cmd
    //CTRL.ual = MT2A(CTRL., CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);  //指令帕克反变换
   // CTRL.ube = MT2B(CTRL.vecter_cmd, CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);  //指令帕克反变换
}

