#ifndef CTIDE_H
#define CTIDE_H

/*
@@brief ��ϫ�࣬���ó�ϫģ�ͶԺ��������峱�Ի�վλ�õ�Ӱ����и�����
*/
class CTide
{
public:
	CTide();
	~CTide();
public:
	/**
	@brief  ���峱����
	@param
	@����ֵ 
	*/
	double* SolidTideCorr();
	/**
	@brief  ����ϫ����
	@param
	@����ֵ 
	*/
	double* OceanTideCorr();
};


#endif // !CTIDE_H


