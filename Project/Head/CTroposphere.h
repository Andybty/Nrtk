#ifndef CTROPOSPHERE_H
#define CTROPOSPHERE_H

/*
@@brief �������࣬��������˲�ͬ�Ķ��������ģ�ͣ�

*/
class CTroposphere
{
public:
	CTroposphere();
	~CTroposphere();
public:
	/**
	@brief  ���ö�����ģ������
	@param
	@����ֵ
	*/
	void SetTropModel(int nTropModel);
	/**
	@brief  ��˹��Ī��ģ��
	@param
	@����ֵ
	*/
	double SaasTropModel();
	/**
	@brief  ������ӳ�亯��
	@param
	@����ֵ ӳ�亯��ֵ
	*/
	double TropMapNeil();

	/**
	@brief  ����˫�������ֵ
	@param
	@����ֵ ӳ�亯��ֵ
	*/
	double DDTropGen();
};

#endif // !CTROPOSPHERE_H

