
#ifndef CONTROLLER_H
#define CONTROLLER_H
//PI定义
struct PI_Reg{
   double   Kp;
   double   Ti;
   double   Ki; // Ki = Kp/Ti*TS
   double   i_state;
   double   i_limit;
};
double PI(struct PI_Reg *r, double err);


struct ControllerForExperiment{

    double timebase;

    double ual;
    double ube;

    double rs;
    double rreq;
    double Lsigma;
    double alpha;
    double Lmu;
    double Lmu_inv;

    double Tload;
    double rpm_cmd;

    double Js;
    double Js_inv;
//增加VCC控制
    double ud_cmd;
    double uq_cmd;

    double id_cmd;
    double iq_cmd;

    double iq_fb;

    double ial_fb;
    double ibe_fb;

    double sinT;
    double cosT;

    double omega_fb;

    double theta_M;

    struct PI_Reg pi_speed;
    struct PI_Reg pi_iMs;
    struct PI_Reg pi_iTs;
//DTC
    double phi_cmd;
    double Tem_cmd;

    double phi;
    double Tem;

    double theta;

    double phial;
    double phibe;

    double ual_fb;
    double ube_fb;

    double vecter_cmd;

    



};
extern struct ControllerForExperiment CTRL;



void CTRL_init();


#endif
