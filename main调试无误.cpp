//
//			static analysis of plane trusses
//
//					main program
//
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include<iomanip>
/* 头文件为#include <iomanip>
其中io代表输入输出，manip是manipulator（操纵器）的缩写
iomanip的作用:
主要是对cin,cout之类的一些操纵运算子，比如setfill,setw,setbase,setprecision等等。它是I/O流控制头文件,就像C里面的格式化输出一样。 */

using namespace std;

// Functions declaration
void INPUT(double X[], double Y[], int NCO[], double PROP[],
	double AL[], int IB[], double REAC[]);
void ASSEM(double X[], double Y[], int NCO[], double PROP[],
	double TK[][20], double ELST[][5], double AL[]);
void STIFF(int NEL, double X[], double Y[], int NCO[], double PROP[],
	double ELST[][5], double AL[]);
void ELASS(int NEL, int NCO[], double TM[][20], double ELMAT[][5]);
void BOUND(double TK[][20], double AL[], double REAC[], int IB[]);
void SLBSI(double A[][20], double B[], double D[], int N, int MS,
	int NRMX, int NCMX);
void FORCE(int NCO[], double PROP[], double FORC[], double REAC[],
	double X[], double Y[], double AL[]);
void OUTPT(double AL[], double FORC[], double REAC[]);

//initialization of global variables
int NN, NE, NLN, NBN, N, MS;
double E, G;

// assign data set numbers to in,for input,and io for output 
ifstream READ_IN;
ofstream WRITE_IO;
// Initialization of program parameters
int NRMX = 200;
int NCMX = 20;
int NDF = 2;
int NNE = 2;
int NDFEL = NDF * NNE;

int main()
{

	double X[100], Y[100], PROP[100], TK[200][20], AL[200], FORC[100],
		REAC[200], ELST[5][5], V[20];
	int NCO[200], IB[60];
	//char file1[20], file2[20];
	//
	// open all files
	//cout << "请输入读入文件名" << endl;
	//cin >> file1;

	/*cout << "请输入数据文件名:" << endl;
	cin >> file2;*/
	//WRITE_IO.open(file1);
	WRITE_IO.open("output.dat");
	READ_IN.open("input2.dat");
	//
	// data input 
	INPUT(X, Y, NCO, PROP, AL, IB, REAC);

	//assembling of total stiffness matrix
	ASSEM(X, Y, NCO, PROP, TK, ELST, AL);
	//
	// introduction of boundary conditions
	BOUND(TK, AL, REAC, IB);
	// solution of the system of equations
	SLBSI(TK, AL, V, N, MS, NRMX, NCMX);
	//
	// computation of member forces
	FORCE(NCO, PROP, FORC, REAC, X, Y, AL);
	//
	// output
	OUTPT(AL, FORC, REAC);

	//close all files
	READ_IN.close();
	WRITE_IO.close();
	return 0;

}

// 

void INPUT(double X[], double Y[], int NCO[], double PROP[],
	double AL[], int IB[], double REAC[])
{
	// input program 
	//
	int I, NUM, N1, IC[2], K, L, L1, L2, N2;
	double W[3];
	WRITE_IO.setf(ios::fixed);
	WRITE_IO.setf(ios::showpoint);
	WRITE_IO << "" << "****************************************************" << endl;

	//
	// read basic parameters
	READ_IN >> NN >> NE >> NLN >> NBN >> E;

	WRITE_IO << "\n\nINTERNAL DATA \n\n" << "Number of Nodes         :"
		<< setw(5) << NN << "\n"
		<< "Number of Elements      :" << setw(5) << NE << "\n"
		<< "Number of Loaded Nodes  :" << setw(5) << NLN << "\n"
		<< "Number of Support Nodes :" << setw(5) << NBN << "\n"
		<< "Modulus of Elasticity   :" << setw(6) << setprecision(0) << E << "\n\n"
		<< "Nodal coordinates \n" << setw(11) << "Node" << setw(7)
		<< "X" << setw(10) << "Y\n";
	//使用setprecision(n)可控制输出流显示浮点数的数字个数。
	//
	//read nodal coordinates in array X and Y 
	for (I = 1; I <= NN; I++)
	{
		NUM = I;
		READ_IN >> NUM >> X[NUM] >> Y[NUM];
		WRITE_IO.precision(2);
		WRITE_IO << setw(10) << NUM << setw(10) << X[I] << setw(10) << Y[I] << "\n";
	}

	//for (I = 1; I <= NN; I++)
	//{
	//	WRITE_IO.precision(2);
	//	WRITE_IO << setw(10) << NUM << setw(10) << X[I] << setw(10) << Y[I] << "\n";
	//}

	//read element connectivity in array NCO and element properties in array PROP
	WRITE_IO << "\n Element connectivity and properties\n" << setw(11)
		<< "Element" << setw(23) << "Start Node  End Node" << setw(9)
		<< "AREA" << endl;
	for (I = 1; I <= NE; I++)
	{
		//NUM = I;
		//READ_IN >> NUM >> IC[0] >> IC[1] >> PROP[NUM];
		READ_IN >> I >> IC[0] >> IC[1] >> PROP[I];
		WRITE_IO.precision(5);
		WRITE_IO << setw(10) << I << setw(10) << IC[0] << setw(10) << IC[1]
			<< setw(15) << PROP[I] << "\n";
		N1 = NNE * (I - 1);

		//调试输出系数***************
		/*cout << "调试输出系数N1+1=" << N1 + 1 << "\n";
		cout << "调试输出系数N1+2=" << N1 + 2 << "\n";*/
		//调试输出系数***************

		NCO[N1 + 1] = IC[0];
		NCO[N1 + 2] = IC[1];

		//调试输出系数对应的每个单元节点号************************************
		/*cout << "NCO=" << NCO[N1+1]<<"\n";
		cout << "NCO=" << NCO[N1 +2] << "\n";*/
		//调试输出系数对应的每个单元节点号************************************

		//
	}
	//
	// compute actual number of unknowns and clear the load vector
	N = NN * NDF;
	for (I = 1; I <= N; I++)
	{
		REAC[I] = 0.0;
		AL[I] = 0.0;
	}
	//
	// Read the nodal loads and store then in array AL 
	WRITE_IO << "\n Nodal loads \n" << setw(11) << "Node" << setw(7) << "PX"
		<< setw(10) << "PY" << endl;
	for (I = 1; I <= NLN; I++)
	{

		READ_IN >> NUM >> W[0] >> W[1];
		//cout << setw(10) << J << setw(10) << W[0] << setw(10) << W[1] << "\n";
		WRITE_IO.precision(2);
		WRITE_IO << setw(10) << NUM << setw(10) << W[0] << setw(10) << W[1] << "\n";
		for (K = 1; K <= NDF; K++)
		{
			L = NDF * (NUM - 1) + K;
			AL[L] = W[K - 1];
		}
	}
	// Read boundary nodes data.store unknown status indicators
	//in array IB,and prescribed unknown values in array REAC
	WRITE_IO << "\n Boundary conditon data\n" << setw(29)
		<< "Status" << setw(31) << "Prescribed values\n" << setw(37)
		<< "(0:prescribed,1:FREE) \n" << setw(11) << "Node" << setw(9)
		<< "U" << setw(10) << "V" << setw(17) << "U" << setw(10) << "V"
		<< endl;
	for (I = 1; I <= NBN; I++)
	{
		NUM = I;
		READ_IN >> NUM >> IC[0] >> IC[1] >> W[0] >> W[1];
		//cout << " NUM="<<NUM<<setw(10)<<"IC[0]="<<IC[0]<<setw(10)<<"IC[1]="<<IC[1]<<"\n";
		WRITE_IO.precision(4);
		WRITE_IO << setw(10) << NUM << setw(10) << IC[0] << setw(10) << IC[1]
			<< setw(20) << W[0] << setw(10) << W[1] << "\n";
		L1 = (NDF + 1)*(I - 1) + 1;
		L2 = NDF * (NUM - 1);
		IB[L1] = NUM;
		for (K = 1; K <= NDF; K++)
		{
			N1 = L1 + K;
			N2 = L2 + K;
			IB[N1] = IC[K - 1];
			REAC[N2] = W[K - 1];
		}
	}
	return;
}



void ASSEM(double X[], double Y[], int NCO[], double PROP[],
	double TK[][20], double ELST[][5], double AL[])
{
	// Assembling of the total matrix for the problem
	int N1, I, L1, J, L2, J1, K, L3, L, NEL;
	// Compute half band width and store in MS
	N1 = NNE - 1;
	MS = 0;
	for (I = 1; I <= NE; I++)
	{
		L1 = NNE * (I - 1);
		for (J = 1; J <= N1; J++)
		{
			L2 = L1 + J;
			J1 = J + 1;
			for (K = J1; K <= NNE; K++)
			{
				L3 = L1 + K;
				L = abs(NCO[L2] - NCO[L3]);
				if ((MS - L) <= 0)
				{
					MS = L;
				}
			}
		}
	}
	MS = NDF * (MS + 1);

	//调试输出*********************************
	//cout << "调试输出半带宽 MS=" << MS << "\n";
	//调试输出********************************

	// clear the total stiffness matrix
	for (I = 1; I <= N; I++)
	{
		for (J = 1; J <= MS; J++)
		{
			TK[I][J] = 0;
		}
	}


	//调试输出初始化的半带矩阵***************************
	//cout << "调试输出初始化的半带矩阵!" << "\n";
	//for (I = 1; I <= N; I++)
	//{
	//	for (J = 1; J <= MS; J++)
	//	{
	//		cout << TK[I][J] << setw(5);
	//	}
	//	cout << "\n";
	//}
	//*********************************************

	for (NEL = 1; NEL <= NE; NEL++)
	{
		// compute the stiffness matrix for element NEL
		STIFF(NEL, X, Y, NCO, PROP, ELST, AL);
		// Place the matrix in the total stiffness matrix
		ELASS(NEL, NCO, TK, ELST);
	}

	//调试输出总刚矩阵***************************
	/*cout << "调试输出总刚矩阵!" << "\n";
	for (I = 1; I <= N; I++)
	{
		for (J = 1; J <= MS; J++)
		{
			cout << TK[I][J] << setw(10);
		}
		cout << "\n";
	}*/
	//*********************************************
	return;
}

// 计算单元刚度矩阵
void STIFF(int NEL, double X[], double Y[], int NCO[],
	double PROP[], double ELST[][5], double AL[])
{
	// computationof element stiffness matrix for current element
	//
	int L, N1, N2, I, J, K1, K2;
	double D, CO, SI, COEF;
	//
	L = NNE * (NEL - 1);
	N1 = NCO[L + 1];
	N2 = NCO[L + 2];
	//
	// compute length of element, and sine and cosine of its local X axis
	D = sqrt(pow((X[N2] - X[N1]), 2) + pow((Y[N2] - Y[N1]), 2));//element length

	//调试输出单元长度***************************
	//cout << "D=" << D << "\n";
	// 调试输出单元长度 ***************************

	CO = (X[N2] - X[N1]) / D;
	SI = (Y[N2] - Y[N1]) / D;

	//调试输出正弦和余弦值***************************
	/*cout << "调试输出正弦和余弦值\n";
	cout << "CO=" << CO<<setw(5)<<"SI="<<SI<< "\n";*/
	// ********************************

	//
	// compute element stiffness matrix
	COEF = E * PROP[NEL] / D;   // COEF=E*A/L

	//调试输出弹性系数***************************
	//cout << "E*A/L=" << COEF << "\n";

	//初始化单元刚度矩阵
	for (I = 1; I <= NDFEL; I++)
	{
		for (J = 1; J <= NDFEL; J++)
		{
			ELST[I][J] = 0;
		}
	}
	ELST[1][1] = COEF * CO*CO;
	ELST[1][2] = COEF * CO*SI;
	ELST[2][2] = COEF * SI*SI;
	for (I = 1; I <= 2; I++)
	{
		for (J = 1; J <= 2; J++)
		{
			K1 = I + NDF;
			K2 = J + NDF;
			ELST[K1][K2] = ELST[I][J];
			ELST[I][K2] = -ELST[I][J];
		}
	}
	ELST[2][3] = -ELST[1][2];

	//调试输出单元刚度矩阵***************************
	/*cout << "调试输出单元刚度矩阵\n";
	for (I = 1; I <= 4; I++)
	{
		for (J = 1; J <= 4; J++)
		{
			cout <<setprecision(3)<< ELST[I][J] << "\t";
		}
		cout << "\n";
	}*/
	//*********************************************

	return;
}



// 将单元刚度矩阵生成总刚度矩阵
void ELASS(int NEL, int NCO[], double TK[][20], double ELMAT[][5])
{
	// store the element matrix for element NEL in the total matrix
	int I, L2, N1, I1, J1, J, N2, I2, J2, K, KI, KR, IC, K1, K2, L, KC;
	//L1 = NNE * (NEL - 1);

	//调试输出系数***************************
	//cout << "调试输出系数 L1=" << L1 << "\n";
	//调试输出系数***************************

	for (I = 1; I <= NNE; I++)
	{
		//L1 = NNE * (NEL - 1);

		L2 = NNE * (NEL - 1) + I;
		N1 = NCO[L2];

		//调试输出系数***************************
		//cout << "调试输出系数 N1=" << N1 << "\n";
		//调试输出系数***************************

		I1 = NDF * (I - 1);
		J1 = NDF * (N1 - 1);
		for (J = 1; J <= NNE; J++)
		{
			//L2 = L1 + J;
			L2 = NNE * (NEL - 1) + J;
			N2 = NCO[L2];
			I2 = NDF * (J - 1);
			J2 = NDF * (N2 - 1);

			for (K = 1; K <= NDF; K++)
			{
				KI = 1;
				if ((N1 - N2) == 0)
				{
					// store A diagonal submatrix
					KI = K;
				}
				if ((N1 - N2) <= 0)
				{
					// store an off diagonal submatrix
					KR = J1 + K;
					IC = J2 - KR + 1;
					K1 = I1 + K;
				}
				else
				{
					// store the transpose of an off diagonal matrix
					KR = J2 + K;
					IC = J1 - KR + 1;
					K2 = I2 + K;
				}
				for (L = KI; L <= NDF; L++)
				{
					KC = IC + L;
					if ((N1 - N2) <= 0)
					{
						K2 = I2 + L;
					}
					else
					{
						K1 = I1 + L;
					}
					TK[KR][KC] = TK[KR][KC] + ELMAT[K1][K2];

					//调试输出**************************************
					//cout << "调试输出 ELMAT["<<K1<<"]["<<K2<<"]=" << ELMAT[K1][K2] << "\n";
					//调试输出**************************************

				}
			}
		}
	}
	//调试输出刚度矩阵***************************
	//cout << "调试输出单元刚度矩阵第"<<I<<"次叠加"<<"\n";
	//for (I = 1; I <= 6; I++)
	//{
	//	for (J = 1; J <= 4; J++)
	//	{
	//		cout << TK[I][J] << setw(10);
	//	}
	//	cout << "\n";
	//}
	//*********************************************

	return;
}

//引入边界条件
void BOUND(double TK[][20], double AL[], double REAC[], int IB[])
{
	// introduction of the boundary conditions
	int L, L1, NO, K1, I, L2, KR, J, KV;
	//
	//调试输出划边界前总刚矩阵***************************
	/*cout << "调试输出划边界前总刚矩阵!" << "\n";
	for (I = 1; I <= N; I++)
	{
		for (J = 1; J <= MS; J++)
		{
			cout << TK[I][J] << setw(10);
		}
		cout << "\n";
	}*/
	//*********************************************
	for (L = 1; L <= NBN; L++)
	{
		L1 = (NDF + 1)*(L - 1) + 1;
		NO = IB[L1];
		K1 = NDF * (NO - 1);
		for (I = 1; I <= NDF; I++)
		{
			L2 = L1 + I;
			if (IB[L2] == 0)
			{
				// prescribed unknown to be considered
				KR = K1 + I;
				for (J = 2; J <= MS; J++)
				{
					KV = KR + J - 1;
					if ((N - KV) >= 0)
					{
						// modify row of TK and corresponding elements in AL
						AL[KV] = AL[KV] - TK[KR][J] * REAC[KR];
						TK[KR][J] = 0.0;
					}
					KV = KR - J + 1;
					if (KV > 0)
					{
						// modify column in TK and conrresponding elements in AL
						AL[KV] = AL[KV] - TK[KV][J] * REAC[KR];
						TK[KV][J] = 0.0;
					}
				}
				// set diagonal coefficient of TK equal to 1 place prescribed unknown
				// value in AL
				TK[KR][1] = 1.0;
				AL[KR] = REAC[KR];
			}
		}
	}
	//调试输出划边界后总刚矩阵***************************
	/*cout << "调试输出划边界后总刚矩阵!" << "\n";
	for (I = 1; I <= N; I++)
	{
		for (J = 1; J <= MS; J++)
		{
			cout << TK[I][J] << setw(10);
		}
		cout << "\n";
	}*/
	//*********************************************
	return;
}

//求解线性方程组
void SLBSI(double A[][20], double B[], double D[], int N, int MS, int NRMX, int NCMX)
{
	// solution of simutaneous systems of equations by the gauss elimunation
	// methond for symmetric banded matrices
	int N1, K, K1, NI, L, J, K2, I, K3;
	double C;
	// 
	N1 = N - 1;
	for (K = 1; K <= N1; K++)
	{
		C = A[K][1];
		K1 = K + 1;

		//cout << "C=" << C << "\n";

		if (C <= 0.000001 && C >= -0.000001)
		{
			WRITE_IO << "IN SLBSI***** singularity in row" << setw(5) << K;
			return;
		}
		else
		{
			//divide row by diagonal coefficient
			NI = K1 + MS - 2;
			if (N1 <= N)
			{
				L = NI;
			}
			else
			{
				L = N;
			}
			for (J = 2; J <= MS; J++)
			{
				D[J] = A[K][J];
			}
			for (J = K1; J <= L; J++)
			{
				K2 = J - K + 1;
				A[K][K2] = A[K][K2] / C;
			}
			B[K] = B[K] / C;

			// eliminate unknown X(K) from row I
			for (I = K1; I <= L; I++)
			{
				K2 = I - K1 + 2;
				C = D[K2];
				for (J = I; J <= L; J++)
				{
					K2 = J - I + 1;
					K3 = J - K + 1;
					A[I][K2] = A[I][K2] - C * A[K][K3];
				}
				B[I] = B[I] - C * B[K];
			}
		}
	}
	//compute last unknown
	if (A[N][1] <= 0.00001 && A[N][1] >= 0.000001)
	{
		WRITE_IO << " ***** singularity in row" << setw(5) << K;
		return;
	}
	else
	{
		B[N] = B[N] / A[N][1];
		//
		// apply backsubstitute process to compute remaining unknows
		for (I = 1; I <= N1; I++)
		{
			K = N - I;
			K1 = K + 1;
			NI = K1 + MS - 2;
			if (NI <= N)
			{
				L = NI;
			}
			else
			{
				L = N;
			}
			for (J = K1; J <= L; J++)
			{
				K2 = J - K + 1;
				B[K] = B[K] - A[K][K2] * B[J];
			}
		}
	}
	return;
}
//计算节点力
void FORCE(int NCO[], double PROP[], double FORC[], double REAC[], double X[], double Y[], double AL[])
{
	// computation of element force
	int I, NEL, L, N1, N2, K1, K2;
	double D, CO, SI, COEF;
	//
	//clear the reaction array
	for (I = 1; I <= N; I++)
	{
		REAC[I] = 0.0;
	}
	for (NEL = 1; NEL <= NE; NEL++)
	{
		L = NNE * (NEL - 1);
		N1 = NCO[L + 1];
		N2 = NCO[L + 2];
		K1 = NDF * (N1 - 1);
		K2 = NDF * (N2 - 1);
		//
		// compute length of element,and sine/cosine of its local x axis
		D = sqrt(pow((X[N2] - X[N1]), 2) + pow((Y[N2] - Y[N1]), 2));
		CO = (X[N2] - X[N1]) / D;
		SI = (Y[N2] - Y[N1]) / D;
		COEF = E * PROP[NEL] / D;

		//
		// compute member axial force and store in array force
		FORC[NEL] = COEF * ((AL[K2 + 1] - AL[K1 + 1])*CO + (AL[K2 + 2] - AL[K1 + 2])*SI);
		//
		// compute nodal resultants
		REAC[K1 + 1] = REAC[K1 + 1] - FORC[NEL] * CO;
		REAC[K1 + 2] = REAC[K1 + 2] - FORC[NEL] * SI;
		REAC[K2 + 1] = REAC[K2 + 1] + FORC[NEL] * CO;
		REAC[K2 + 2] = REAC[K2 + 2] + FORC[NEL] * SI;
	}
	return;
}
//输出结果
void OUTPT(double AL[], double FORC[], double REAC[])
{
	// output program
	int I, K1, K2, J;
	//
	//write nodal displaements
	WRITE_IO <<
		"\n\n ***************************************************\n\n"
		<< "Results\n\n" << "Nodal displacements\n" << setw(11) << "Node"
		<< setw(12) << "U" << setw(15) << "V" << endl;
	for (I = 1; I <= NN; I++)
	{
		K1 = NDF * (I - 1) + 1;
		K2 = K1 + NDF - 1;
		WRITE_IO << setw(10) << I;
		for (J = K1; J <= K2; J++)
		{
			WRITE_IO << setw(15) << AL[J];
		}
		WRITE_IO << endl;
	}
	//
	//write nodal reactions
	WRITE_IO << "\n Nodal Reactions\n" << setw(11) << "Node" << setw(12) << "PX"
		<< setw(15) << "PY\n";
	for (I = 1; I <= NN; I++)
	{
		K1 = NDF * (I - 1) + 1;
		K2 = K1 + NDF - 1;
		WRITE_IO << setw(10) << I;
		for (J = K1; J <= K2; J++)
		{
			WRITE_IO << setw(15) << REAC[J];
		}
		WRITE_IO << endl;
	}
	//
	// write member axial forces
	WRITE_IO << "\nMember forces" << setw(27) << "Member axial force\n";
	for (I = 1; I <= NE; I++)
	{
		WRITE_IO << setw(10) << I << setw(15) << FORC[I] << endl;
	}
	WRITE_IO <<
		"\n\n****************************************************\n";
	return;
}