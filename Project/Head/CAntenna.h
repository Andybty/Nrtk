#ifndef CANTENNA_H
#define CANTENNA_H

/*
@@brief �����࣬��Ҫ���������������и�����
*/
class CAntenna
{
public:
	CAntenna();
	~CAntenna();
public:
	/**
	@brief  ������������λ���ı仯���и���
	@param
	@����ֵ
	*/
	double* SatPcvCorr();
	double* SatPcoCorr();
	/**
	@brief  �Խ��ջ�������λ���Ľ��и���
	@param
	@����ֵ
	*/
	double* RecPcvCorr();
	double* RecPcoCorr();
private:
	//AntData;

};
#endif // !CANTENNA_H
