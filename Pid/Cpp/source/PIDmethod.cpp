#include "../include/PIDmethod.h"

SystemTick_Fun PIDtimer::Get_SystemTick = NULL;//静态变量必须实现

/**
 * @brief Update time stamp
 * @return none
 */
uint8_t PIDtimer::UpdataTimeStamp(void)
{
    uint32_t now_time;

    /*Check `Get_SystemTick` */
    if (PIDtimer::Get_SystemTick != NULL)
    {
        /*Convert to system time*/
        if (last_time == 0)
        {
            last_time = PIDtimer::Get_SystemTick();
            return 1;//第一步先不跑
        }
        now_time = PIDtimer::Get_SystemTick();

        /*Overflow*/
        if (now_time < last_time)
            dt = (double)(now_time + (0xFFFFFFFF - last_time));
        else
            dt = (double)(now_time - last_time);

        last_time = now_time;

        dt *= (double)0.000001;

        return 0;
    }
    else {
        dt = 0.001f;
        return 1;
    }
}

/**
 * @brief  Regist get time function(1Tick = 1us)
 * @param  realTime_fun: Pointer of function to get system real time
 * @retval 1: success
           0: error input param
 * @author
 */
uint8_t PIDtimer::getMicroTick_regist(uint32_t(*getTick_fun)(void))
{
    if (getTick_fun != NULL)
    {
        PIDtimer::Get_SystemTick = getTick_fun;
        return 1;
    }
    else
        return 0;
}


/*设置PID控制模式*/
PIDmethod::PIDmethod(Params_Mode mode, double _timeStep)
{
    params_mode = mode;
    timeStep = _timeStep;
}

void PIDmethod::PID_Init(Params_Mode mode, double _timeStep)
{
    params_mode = mode;
    timeStep = _timeStep;
}

/****************************************** PID参数录入 ******************************************/
void PIDmethod::Params_Config(Fit_Params _fun_p, Fit_Params _fun_i, Fit_Params _fun_d, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    fun_p = _fun_p;
    fun_i = _fun_i;
    fun_d = _fun_d;
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

void PIDmethod::Params_Config(Fit_Params _fun_p, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    fun_p = _fun_p;
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

void PIDmethod::Params_Config(PID_Mode mode, Fit_Params _fun_p, Fit_Params _fun_id, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    fun_p = _fun_p;
    if (mode == PID_Mode::IS_PI)
    {
        fun_i = _fun_id;
    }
    else if (mode == PID_Mode::IS_PD)
    {
        fun_d = _fun_id;
    }
    else {}
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

void PIDmethod::Params_Config(double _kp, double _ki, double _kd, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    kp = _kp;
    ki = _ki;
    kd = _kd;
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

void PIDmethod::Params_Config(double _kp, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    kp = _kp;
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

void PIDmethod::Params_Config(PID_Mode mode, double _kp, double _kid, double _I_Term_Max, double _Output_Max, double _Output_Min)
{
    kp = _kp;
    if (mode == PID_Mode::IS_PI)
    {
        ki = _kid;
    }
    else if (mode == PID_Mode::IS_PD)
    {
        kd = _kid;
    }
    else {}
    I_Term_Max = _I_Term_Max;
    Output_Max = _Output_Max;
    Output_Min = _Output_Min;
}

/**************************************** 拟合三次函数 *********************************************/
double PIDmethod::fit_function(Fit_Params param, double x)
{
    return param.a * pow(x, 3) + param.b * pow(x, 2) + x * param.c + param.d;
}
/****************************************** PID运算 **********************************************/
double PIDmethod::Adjust(double _x)
{
    if (timeStep > 0)
    {
        this->dt = timeStep;
    }
    else
    {
        //if (this->UpdataTimeStamp())
        //    return 0;//如果时间栈出错则不执行pid
        this->dt = this->UpdataTimeStamp();
    }

    if (params_mode == Fit)
    {
        fact_kp = fit_function(fun_p, _x);
        fact_ki = fit_function(fun_i, _x);
        fact_kd = fit_function(fun_d, _x);
    }
    else if (params_mode == Common)
    {
        fact_kp = kp;
        fact_ki = ki;
        fact_kd = kd;
    }
    else
    {
        return 0;//参数出错不执行
    }

    error = upper::constrain(target - current, Error_Max);
    /*error = target - current;*/
    d_error = (error - last_error) / this->dt;
    d_current = (current - last_current) / this->dt;

    P_Term = error * fact_kp;

    integral += error * this->dt;
    integral = upper::constrain(integral, I_Term_Max / fact_ki);
    I_Term = integral * fact_ki;
    if (abs(I_Term) > I_SeparThresh)
    {
        I_Term = 0;
    }
    else {}

    if (d_of_current)
    {
        D_Term = d_current * fact_kd;
    }
    else
    {
        D_Term = d_error * fact_kd;
    }

    out = P_Term + I_Term + D_Term;

    if (Output_Min >= Output_Max)
    {
        out = upper::constrain(out, Output_Max);
    }
    else
    {
        out = upper::constrain(out, Output_Min, Output_Max);
    }


    last_current = current;
    last_error = error;
    return out;
}

double PIDmethod::Adjust(double _x, double extern_d)
{

    if (timeStep > 0)
    {
        this->dt = timeStep;
    }
    else
    {
        //if (this->UpdataTimeStamp())
        //    return 0;//如果时间栈出错则不执行pid
        this->dt = this->UpdataTimeStamp();
    }

    if (params_mode == Fit)
    {
        fact_kp = fit_function(fun_p, _x);
        fact_ki = fit_function(fun_i, _x);
        fact_kd = fit_function(fun_d, _x);
    }
    else if (params_mode == Common)
    {
        fact_kp = kp;
        fact_ki = ki;
        fact_kd = kd;
    }
    else
    {
        return 0;//参数出错不执行
    }

    error = upper::constrain(target - current, Error_Max);
    /*error = target - current;*/
    d_error = (error - last_error) / this->dt;
    d_current = (current - last_current) / this->dt;

    P_Term = error * fact_kp;

    integral += error * this->dt;
    integral = upper::constrain(integral, I_Term_Max / fact_ki);
    I_Term = integral * fact_ki;
    if (abs(I_Term) > I_SeparThresh)
    {
        I_Term = 0;
    }
    else {}

    D_Term = extern_d * fact_kd;

    out = P_Term + I_Term + D_Term;
    if (Output_Min >= Output_Max)
    {
        out = upper::constrain(out, Output_Max);
    }
    else
    {
        out = upper::constrain(out, Output_Min, Output_Max);
    }

    last_current = current;
    last_error = error;
    return out;
}

