#ifndef CANTENNA_H
#define CANTENNA_H

/*
@@brief 天线类，主要对卫星天线误差进行改正；
*/
class CAntenna
{
public:
	CAntenna();
	~CAntenna();
public:
	/**
	@brief  对卫星天线相位中心变化进行改正
	@param
	@返回值
	*/
	double* SatPcvCorr();
	double* SatPcoCorr();
	/**
	@brief  对接收机天线相位中心进行改正
	@param
	@返回值
	*/
	double* RecPcvCorr();
	double* RecPcoCorr();
private:
	//AntData;

};
#endif // !CANTENNA_H
