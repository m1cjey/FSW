#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"MMM_CONFIG.h"	//class CON定義

MMM_config::MMM_config()
{
	///////解析条件
	Hf_type=0;		//外部磁場　0:一様磁場 1:磁石
	Hf_H=30000;		//外部磁場Hの強さ
}