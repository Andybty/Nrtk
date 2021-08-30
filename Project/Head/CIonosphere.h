#ifndef CIONOSPHERE_H
#define CIONOSPHERE_H
/*
@@brief 电离层类，里面包含了不同的电离层改正模型，
通过不同的模型计算天顶方向或者倾斜方向上的电离层值（米）
*/
class CIonosphere
{
public:
	CIonosphere();
	~CIonosphere();
public:
	/**
	@brief  设置使用的电离层模型类型
	@param
	@返回值
	*/
	void SetIonoModel(int nionmodel);  
	/**
	@brief  Klobuchar电离层模型
	@param
	@返回值
	*/
	bool KlobIonoModel(); 
	/**
	@brief  GIM电离层模型
	@param
	@返回值
	*/
	bool GimIonoModel();   
	/**
	@brief  电离层投影系数
	@param
	@返回值
	*/
	void IonoMapFunc();  
	/**
	@brief  估计电离层值
	@param
	@返回值
	*/
	bool IonoEstModel();  
	/**
	@brief  消电离层电离层模型
	@param
	@返回值
	*/
	bool IonoFreeModel();   
	/**
	@brief  得到倾斜方向的电离层延迟和对应的方差
	@param
	@返回值
	*/
	void GetSlantIono(double dSlantIono, double dVarSIono);    
	/**
	@brief  得到天顶方向的电离层延迟和对应的方差
	@param
	@返回值
	*/
	void GetZeniIono(double  dZenIono, double dVarZIono);   

	double DDIonoGen();
private:
	int     m_nIonoModel   ;                  /*电离层模型类型*/
	double  m_dSlantIono   ;                  /*L1倾斜方向上的电离层延迟*/
	double  m_dZenIono     ;                  /*L1天顶方向上的电离层延迟*/
	double  m_dVarSIono    ;                  /*L1倾斜方向上的电离层延迟方差*/
	double  m_dVarZIono    ;                  /*L1天顶方向上的电离层延迟方差*/
	double  m_dMapCoffe    ;                  /*L1倾斜方向上的电离层投影函数*/
};




#endif // !CIONOSPHERE_H


