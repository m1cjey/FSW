#ifndef HEADER
#define HEADER


#include<stdio.h>

#define _CRTDBG_MAP_ALLOC//���������[�N���o�p
#include<stdlib.h>
#include <crtdbg.h>//���������[�N���o�p
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)//���������[�N���o�p
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
#include "tetgen.h"		//TetGen�w�b�_�[
#include <string>

using namespace std;

#include"define.h"		//#define �i�[
#include"CONFIG.h"		//�ݒ�p�ϐ� ��`
#include"PART.h"	//PART class ��`
#include"FEM3Dclass.h"	//FEM class ��`

#include "tetclass.h"	//TetGen�p����N���X
#include "tetfunc.h"	//TetGen�p����֐�


#endif