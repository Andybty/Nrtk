#ifndef CTIDE_H
#define CTIDE_H

/*
@@brief 潮汐类，采用潮汐模型对海潮、固体潮对基站位置的影响进行改正；
*/
class CTide
{
public:
	CTide();
	~CTide();
public:
	/**
	@brief  固体潮改正
	@param
	@返回值 
	*/
	double* SolidTideCorr();
	/**
	@brief  海洋潮汐改正
	@param
	@返回值 
	*/
	double* OceanTideCorr();
};


#endif // !CTIDE_H


