#ifndef CIONOSPHERE_H
#define CIONOSPHERE_H
/*
@@brief ������࣬��������˲�ͬ�ĵ�������ģ�ͣ�
ͨ����ͬ��ģ�ͼ����춥���������б�����ϵĵ����ֵ���ף�
*/
class CIonosphere
{
public:
	CIonosphere();
	~CIonosphere();
public:
	/**
	@brief  ����ʹ�õĵ����ģ������
	@param
	@����ֵ
	*/
	void SetIonoModel(int nionmodel);  
	/**
	@brief  Klobuchar�����ģ��
	@param
	@����ֵ
	*/
	bool KlobIonoModel(); 
	/**
	@brief  GIM�����ģ��
	@param
	@����ֵ
	*/
	bool GimIonoModel();   
	/**
	@brief  �����ͶӰϵ��
	@param
	@����ֵ
	*/
	void IonoMapFunc();  
	/**
	@brief  ���Ƶ����ֵ
	@param
	@����ֵ
	*/
	bool IonoEstModel();  
	/**
	@brief  �����������ģ��
	@param
	@����ֵ
	*/
	bool IonoFreeModel();   
	/**
	@brief  �õ���б����ĵ�����ӳٺͶ�Ӧ�ķ���
	@param
	@����ֵ
	*/
	void GetSlantIono(double dSlantIono, double dVarSIono);    
	/**
	@brief  �õ��춥����ĵ�����ӳٺͶ�Ӧ�ķ���
	@param
	@����ֵ
	*/
	void GetZeniIono(double  dZenIono, double dVarZIono);   

	double DDIonoGen();
private:
	int     m_nIonoModel   ;                  /*�����ģ������*/
	double  m_dSlantIono   ;                  /*L1��б�����ϵĵ�����ӳ�*/
	double  m_dZenIono     ;                  /*L1�춥�����ϵĵ�����ӳ�*/
	double  m_dVarSIono    ;                  /*L1��б�����ϵĵ�����ӳٷ���*/
	double  m_dVarZIono    ;                  /*L1�춥�����ϵĵ�����ӳٷ���*/
	double  m_dMapCoffe    ;                  /*L1��б�����ϵĵ����ͶӰ����*/
};




#endif // !CIONOSPHERE_H


