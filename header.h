#ifndef HEADER
#define HEADER


#include<stdio.h>

#define _CRTDBG_MAP_ALLOC//メモリリーク検出用
#include<stdlib.h>
#include <crtdbg.h>//メモリリーク検出用
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)//メモリリーク検出用
#include<string.h>
#include<tchar.h>
#include<time.h>
#include<assert.h>
#include<conio.h>
#include<math.h>
#include<windows.h>
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include "tetgen.h"		//TetGenヘッダー
#include <string>

using namespace std;

#include"define.h"		//#define 格納
#include"CONFIG.h"		//設定用変数 定義
#include"PART.h"	//PART class 定義
#include"FEM3Dclass.h"	//FEM class 定義

#include "tetclass.h"	//TetGen用自作クラス
#include "tetfunc.h"	//TetGen用自作関数


#endif