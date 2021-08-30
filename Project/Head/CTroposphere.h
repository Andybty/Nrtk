#ifndef CTROPOSPHERE_H
#define CTROPOSPHERE_H

/*
@@brief 对流层类，里面包含了不同的对流层改正模型，

*/
class CTroposphere
{
public:
	CTroposphere();
	~CTroposphere();
public:
	/**
	@brief  设置对流层模型类型
	@param
	@返回值
	*/
	void SetTropModel(int nTropModel);
	/**
	@brief  萨斯塔莫宁模型
	@param
	@返回值
	*/
	double SaasTropModel();
	/**
	@brief  对流程映射函数
	@param
	@返回值 映射函数值
	*/
	double TropMapNeil();

	/**
	@brief  生成双差对流层值
	@param
	@返回值 映射函数值
	*/
	double DDTropGen();
};

#endif // !CTROPOSPHERE_H

