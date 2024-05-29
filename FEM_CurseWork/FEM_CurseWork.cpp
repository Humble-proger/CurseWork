#include <iostream>
#include <Windows.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <functional>
#include <tuple>
#include <fstream>
#include <set>
#include <vector>

using namespace std;
// Реализация функции sign
template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}
template <typename T> int ElemInList(T* List, T Elem, int n) {
	for (int i = 0; i < n; i++)
		if (List[i] == Elem) return i;
	return -1;
}

#pragma region Глобальная матрица и вектор правой части
vector<vector<vector<int>>> edge;
vector<vector<double>> A_1, A_2, A_3;
vector<double> ggl, di, b, q, L_sq, di_sq;
vector<int> ig, jg;
//вектор весов, результат работы программы
int* choise;
#pragma endregion
#pragma region Конечноэлементная сетка
int Wn, Dn; // Размерность массивов W, D
function<double(double, double)> *W; //Базисные функции
function<double(double, double)>* D[2]; //Производные от базисных функций
double koef_r; //Коэффициент растяжения по r
int nr; // Количество разбиений по r и z  соответственно
int* numfe[9];	     // номера узлов конечных элементов
double* rz[2];	     // координаты узлов
int NodesCount;      //количество узлов
int FECount;	     //количество КЭ
int NodesCountR;
#pragma endregion
#pragma region Локальные матрицы Жёсткости и Массы и вестора правой части
//Локальные матрицы жёсткости и массы
int GLen, MLen; //Размерность массивов G, M
double alpha[3]; //Элементы альфа которые участвуют при счёте Якоби
double Myu[6]; //Коэффициенты Мю которые участвуют после решения Крамера
double* MkXi,
	  * MkSig,
      * Gk;
double* bk;
#pragma endregion
#pragma region Сетка по времени

int num_t; // Количество узлов в сетке
vector<double> tj; //Массив с значениями tj
double koef_t; // Коэффициент растяжения по оси времени

#pragma endregion
#pragma region Вектор решений или вектор qj и константы
vector<vector<double>> qj; // Размер данного вектора будет ограничен до 4 необходимых. Чтобы меньше тратить памяти.

#pragma endregion


#pragma region Прототипы
#pragma region Прототипы Гаусса
double Gauss2rod(function<double(double, double, tuple<int, int, int>, double)> Ifun, tuple<int, int, int> FuncIndex, double time);
double Gauss2rod(function<double(double, double, tuple<int, int, int>)> Ifun, tuple<int, int, int> FuncIndex);
#pragma endregion
#pragma region Прототипы параметров уравнения
double lambda(double z);
double RightFunction(double r, double z, double time);
#pragma endregion
#pragma region Прототипы краевых условий
double ug(int IndexFunc, double r, double z, double t);
double teta(int IndexFunc, double r, double z, double t);
double beta();
double sigma(double r, double z);
double xi_d(double r, double z);
double u_beta(int IndexFunc, double r, double z, double t);
double u0(double r, double z);
double u1(double r, double z);
double u2(double r, double z);
void EnterBoundCondit();
int GetStepEndEl(int i);
int GetStartPoint(int i);
bool FindInd(int i);
void ConsiderBoundConditFirstType(int node_num, double t);
void ConsiderBoundConditSecType(int num_edge, double t);
void AddLocalMartBound_3(int num_edge, vector<vector<double>> matr);
void AddLocalVecBound(int bound_condit, int num_edge, vector<double> vec);
void ConsiderBoundConditThirdType(int num_edge, double t);
void BuildBoundMatrices();
void BoundCondit();
void ConsiderBoundCondit(double t);
#pragma endregion
#pragma region Прототипы загрузки и обработки сетки
void LoadNet(bool debag = false); // Загрузка сетки
void LoadTNet();
void InitQJ();
void InitFirstQ();
tuple<double, double> GetCordIK(int i, int IndexFE);
#pragma endregion
#pragma region Прототипы загрузки базисных функций и их производных
double FiniteFunc(int i, double var);
double DerivFiniteFunc(int i, double var);
void InitBaseFunction();
void InitDРerivativeBaseFunction();
#pragma endregion
#pragma region Прототипы построений векторов Альфа и Мю
void CalcVecAlpha(unsigned int indexFE);
void CalcVecMyu(unsigned int indexFE);
#pragma endregion

#pragma region Построение глобальной матрицы и вектора
void LocalCalcGK(int IndexFE);
double RPsyEta(double Psy, double Eta, int IndexFe); // r с заменой меременной
double ZPsyEta(double Psy, double Eta, int IndexFE);
void LocalCalcVectorFEK(int IndexFE, double time);
void SumVector(double* &FirstVec, double *SecondVec, int n);
void MultVectorConstant(double*& FirstVec, double Constant, int n);
void InitVectors();

#pragma endregion
#pragma endregion
#pragma region Реализация функций
#pragma region Реализация Загрузок и обработок сеток
void LoadNet(bool debag) {
	/*
	int L; // Число подобластей
    int* numfe[4],	// номера узлов конечных элементов
       * numsubarea[2];		// номера подобластей конечных элементов
    double koef_r; //Коэффициент растяжения по r от 0 до 2 не включая 0 и 2
    double koef_z; //Коэффициент растяжения по z 
    int nr, nz; // Количество разбиений по r и z  соответственно
	double* rz[2];	// координаты узлов
    int NodesCount; //количество узлов
    int FECount;	//количество КЭ
	//2 род
    int* numbc2[2], //ребра, на которых заданы
    int* numbc2[2], //ребра, на которых заданы
         numbc2_n;  //размерность
    //3 род
    int* numbc3[2], //ребра, на которых заданы
         numbc3_n;  //размерность
	*/
	
	ifstream NodesFile("Nodes.txt");
	if (!NodesFile.is_open())
		throw "Файл Nodes.txt не удалось открыть";
	double* rz_area[2]; // Координаты области
	rz_area[0] = new double[4]; rz_area[1] = new double[4];
	
	for (int i = 0; i < 4; i++) {
		NodesFile >> rz_area[0][i] >> rz_area[1][i];
	}
	NodesFile.close();
	ifstream FEFile("fe.txt");
	if (!FEFile.is_open()) {
		throw "Файл fe.txt не удалось открыть";
	}
	FEFile >> koef_r;
	//Генерация сетки
	double* r0[2], * r1[2];
	FEFile >> nr;
	FEFile.close();
	NodesCountR = 2 * nr + 1;
	NodesCount = NodesCountR * NodesCountR;
	FECount = nr * nr;
	for (int i = 0; i < 2; i++) {
		rz[i] = new double[NodesCount];
		r0[i] = new double[NodesCountR];
		r1[i] = new double[NodesCountR];
	}
	double rh, zh, max_val; 
	if (abs(koef_r - 1) < 1e-15) {
		rh = (rz_area[0][2] - rz_area[0][0]) / (NodesCountR - 1);
		zh = (rz_area[1][2] - rz_area[1][0]) / (NodesCountR - 1);
		for (int i = 0; i < NodesCountR; i++) {
			r0[0][i] = rz_area[0][0] + rh * i;
			r0[1][i] = rz_area[1][0] + zh * i;
		}
		rh = (rz_area[0][3] - rz_area[0][1]) / (NodesCountR - 1);
		zh = (rz_area[1][3] - rz_area[1][1]) / (NodesCountR - 1);
		for (int i = 0; i < NodesCountR; i++) {
			r1[0][i] = rz_area[0][1] + rh * i;
			r1[1][i] = rz_area[1][1] + zh * i;
		}
		for (int i = 0, p = 0; i < NodesCountR; i++) {
			rh = (r1[0][i] - r0[0][i]) / (NodesCountR - 1);
			zh = (r1[1][i] - r0[1][i]) / (NodesCountR - 1);
			for (int j = 0; j < NodesCountR; j++, p++) {
				rz[0][p] = r0[0][i] + rh * j;
				rz[1][p] = r0[1][i] + zh * j;

			}
		}
	}
	else {
		rh = (rz_area[0][2] - rz_area[0][0]);
		zh = (rz_area[1][2] - rz_area[1][0]);
		if (koef_r > 1) {
			max_val = pow(koef_r, NodesCountR - 2);
			r0[0][0] = rz_area[0][0];
			r0[1][0] = rz_area[1][0];
			for (int i = 0; i < NodesCountR - 1; i++) {
				r0[0][i + 1] = rz_area[0][0] + rh * pow(koef_r, i) / max_val;
				r0[1][i + 1] = rz_area[1][0] + zh * pow(koef_r, i) / max_val;
			}
			rh = (rz_area[0][3] - rz_area[0][1]);
			zh = (rz_area[1][3] - rz_area[1][1]);
			r1[0][0] = rz_area[0][1];
			r1[1][0] = rz_area[1][1];
			for (int i = 0; i < NodesCountR - 1; i++) {
				r1[0][i + 1] = rz_area[0][1] + rh * pow(koef_r, i) / max_val;
				r1[1][i + 1] = rz_area[1][1] + zh * pow(koef_r, i) / max_val;
			}
			int p = 0;
			for (int i = 0; i < NodesCountR; i++) {
				rh = (r1[0][i] - r0[0][i]);
				zh = (r1[1][i] - r0[1][i]);
				rz[0][p] = r0[0][i];
				rz[1][p] = r1[1][i];
				p += 1;
				for (int j = 0; j < NodesCountR - 1; j++) {
					rz[0][p] = r0[0][i] + rh * pow(koef_r, j) / max_val;
					rz[1][p] = r0[1][i] + zh * pow(koef_r, j) / max_val;
					p += 1;
				}
			}

		}
		else {
			koef_r += 1;
			double max_val = pow(koef_r, NodesCountR - 2);
			r0[0][NodesCountR - 1] = rz_area[0][2];
			r0[1][NodesCountR - 1] = rz_area[1][2];
			rh *= -1; zh *= -1;
			for (int i = 0; i < NodesCountR - 1; i++) {
				r0[0][NodesCountR - i - 2] = rz_area[0][2] + rh * pow(koef_r, i) / max_val;
				r0[1][NodesCountR - i - 2] = rz_area[1][2] + zh * pow(koef_r, i) / max_val;
			}
			rh = -(rz_area[0][3] - rz_area[0][1]);
			zh = -(rz_area[1][3] - rz_area[1][1]);
			r1[0][NodesCountR - 1] = rz_area[0][3];
			r1[1][NodesCountR - 1] = rz_area[1][3];
			for (int i = 0; i < NodesCountR - 1; i++) {
				r1[0][NodesCountR - i - 2] = rz_area[0][3] + rh * pow(koef_r, i) / max_val;
				r1[1][NodesCountR - i - 2] = rz_area[1][3] + zh * pow(koef_r, i) / max_val;
			}
			int p = NodesCount - 1;
			for (int i = 0; i < NodesCountR; i++) {
				rh = -(r1[0][NodesCountR - i - 1] - r0[0][NodesCountR - i - 1]);
				zh = -(r1[1][NodesCountR - i - 1] - r0[1][NodesCountR - i - 1]);
				rz[0][p] = r0[0][NodesCountR - i - 1];
				rz[1][p] = r1[1][NodesCountR - i - 1];
				p--;
				for (int j = 0; j < NodesCountR - 1; j++, p--) {
					rz[0][p] = r1[0][NodesCountR - i - 1] + rh * pow(koef_r, j) / max_val;
					rz[1][p] = r1[1][NodesCountR - i - 1] + zh * pow(koef_r, j) / max_val;
				}
			}
		}
	}
	
	//Заполнение массива rz_fe
	for (int i = 0; i < 9; i++) numfe[i] = new int[FECount];
	int xn = 0;
	for (int i = 0; i < FECount; i++) {
		if (i != 0 && i % nr == 0) xn += NodesCountR + 1;
		for (int j = 0; j < 3; j++) {
			numfe[j][i] = xn; numfe[j + 3][i] = xn + NodesCountR; numfe[j + 6][i] = xn + 2 * NodesCountR;
			if (j == 2) continue;
			xn++;
		}
	}

	if (debag) {
		printf("Print Num Global Function\n");
		for (int k = 0; k < FECount; k++) {
			printf("%d %d %d\n", numfe[6][k], numfe[7][k], numfe[8][k]);
			printf("%d %d %d\n", numfe[3][k], numfe[4][k], numfe[5][k]);
			printf("%d %d %d\n", numfe[0][k], numfe[1][k], numfe[2][k]);
		}
	}
}
tuple<double, double> GetCordIK(int i, int IndexFE) {
	int NumGlobalFunc = numfe[i][IndexFE];
	return make_tuple(rz[0][NumGlobalFunc], rz[1][NumGlobalFunc]);
}
void LoadTNet() {
	ifstream FileT("time.txt");
	if (!FileT.is_open()) {
		throw "Файл time.txt не открылся";
	}
	double max_t, min_t;
	FileT >> min_t >> max_t >> koef_t >> num_t;
	FileT.close();
	if (koef_t <= 0 || num_t <= 0) {
		throw "Сетка по времени заданна не верно";
	}
	tj.resize(num_t);
	if (abs(koef_t - 1) < 1e-15) {
		double h0 = (max_t - min_t) / (num_t - 1);
		for (int i = 0; i < num_t; i++) {
			tj[i] = min_t + h0 * i;
		}
	}
	else {
		double h0 = (max_t - min_t) * (koef_t - 1) / (pow(koef_t, num_t - 1) - 1);
		for (int i = 0; i < num_t; i++) {
			tj[i] = min_t + h0 * (pow(koef_r, i) - 1) / (koef_t - 1);
		}
	}
}
void InitQJ() {
	qj.resize(3);
	for (int i = 0; i < 3; i++) {
		qj[i].resize(NodesCount, 0);
	}
}
void InitFirstQ() {
	for (int i = 0; i < NodesCount; i++) {
		qj[0][i] = u0(rz[0][i], rz[1][i]);
		qj[1][i] = u1(rz[0][i], rz[1][i]);
		//qj[2][i] = u2(rz[0][i], rz[1][i]);
	}
}
#pragma endregion
#pragma region Реализация загрузок базисных функций и и производных
double FiniteFunc(int i, double var) {
	if (i == 0) {
		return 2 * (var - 0.5) * (var - 1.0);
	}
	else if (i == 1) {
		return -4 * var * (var - 1.0);
	}
	else return 2 * var * (var - 0.5);
}

double DerivFiniteFunc(int i, double var) { // производная финитной функции
	if (i == 0) {
		return 4 * var - 3;
	}
	else if (i == 1) {
		return -8 * var + 4;
	}
	else return 4 * var - 1;
}
void InitBaseFunction() {
	/*int Wn, Dn; // Размерность массивов W, D
      function<double(double)> *W; //Базисные функции
      function<double(double)>* D; //Производные от базисных функций
    */
	Wn = 9;
	W = new function<double(double, double)> [Wn];
	for (int i = 0, p = 0; i < 3; i++)
		for (int j = 0; j < 3; j++, p++) {
			W[p] = [i, j](double Psy, double Eta) {return FiniteFunc(j, Psy) * FiniteFunc(i, Eta);  };
		}
}
void InitDРerivativeBaseFunction() {
	Dn = Wn * 2;
	D[0] = new function<double(double, double)>[9]; D[1] = new function<double(double, double)>[9];
	for (int i = 0, p = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++, p++) {
			D[0][p] = [i, j](double Psy, double Eta) {return DerivFiniteFunc(j, Psy) * FiniteFunc(i, Eta);  };
			D[1][p] = [i, j](double Psy, double Eta) {return FiniteFunc(j, Psy) * DerivFiniteFunc(i, Eta);  };
		}
	}
}

#pragma endregion

#pragma region Реализация функций краевых условий
double ug(int IndexFunc, double r, double z, double t) {
	return exp(3 * t);
}
double teta(int IndexFunc, double r, double z, double t) {
	if (IndexFunc == 0) {
		return (-2 * r) / sqrt(26);
		
	}
	else {
		return (-2 * r) / sqrt(26);
	}
}
double beta() {
	return 1;
}


double u0(double r, double z) {
	return exp(3 * tj[0]);
}

double u1(double r, double z) {
	return exp(3 * tj[1]);
}
double u2(double r, double z) {
	return exp(3 * tj[2]);
}

double u_beta(int IndexFunc, double r, double z, double t) {
	if (IndexFunc == 0) {
		return pow(r, 2) + t + 18 * r / sqrt(82);
	}
	else {
		return pow(r, 2) + t - 10 * r / sqrt(26);
	}
}
void EnterBoundCondit() { // ввод краевых условий
	cout << "Введите тип краевого условия на соответствующем ребре:" << endl;
	cout << "1 - первое краевое условие;" << endl;
	cout << "2 - второе краевое условие;" << endl;
	cout << "3 - третье краевое условие;" << endl << endl;
	choise = new int[4];
	cout << "Краевое условие нижнего ребра:" << endl;
	cin >> choise[0];
	cout << endl;
	cout << "Краевое условие правого ребра:" << endl;
	cin >> choise[1];
	cout << endl;
	cout << "Краевое условие верхнего ребра:" << endl;
	cin >> choise[2];
	cout << endl;
	cout << "Краевое условие левого ребра:" << endl;
	cin >> choise[3];
	cout << endl;
}
int GetStepEndEl(int i) { // задаем шаг до след конеч. эл.
	if (i == 0 || i == 2) {
		return 2;
	}
	else return 2 * NodesCountR;
}
int GetStartPoint(int i) { // задаем начальную точку ребра
	switch (i)
	{
	case 0:
		return 0;
		break;
	case 1:
		return NodesCountR - 1;
		break;
	case 2:
		return NodesCount - NodesCountR;
		break;
	case 3:
		return 0;
		break;
	}
}
bool FindInd(int i) {
	for (int j = 0; j < edge[0].size(); j++) {
		if (edge[0][j][0] == i) {
			return true;
		}
	}
	return false;
}
void ConsiderBoundConditFirstType(int node_num, double t) { // учет краевых условий первого типа
	double _r, _z;
	_r = rz[0][edge[0][node_num][0]], _z = rz[1][edge[0][node_num][0]];
	b[edge[0][node_num][0]] = ug(edge[0][node_num][1], _r, _z, t);
	di[edge[0][node_num][0]] = 1;
	for (int i = ig[edge[0][node_num][0]]; i < ig[edge[0][node_num][0] + 1]; i++) {
		int _i = jg[i];
		if (FindInd(_i)) {
			ggl[i] = 0;
			continue;
		}
		b[_i] -= b[edge[0][node_num][0]] * ggl[i];
		ggl[i] = 0;
	}
	for (int i = edge[0][node_num][0]; i < NodesCount; i++) {
		int k = 0;
		for (int j = ig[i]; j < ig[i + 1]; j++) {
			if (jg[j] == edge[0][node_num][0]) {
				if (FindInd(i)) {
					ggl[j] = 0;
					continue;
				}
				b[i] -= b[edge[0][node_num][0]] * ggl[j];
				ggl[j] = 0;
			}
		}
	}
}

void ConsiderBoundConditSecType(int num_edge, double t) { // учет краевых условий второго типа
	vector<double> _r(3), _z(3);
	_r[0] = rz[0][edge[1][num_edge][0]], _z[0] = rz[1][edge[1][num_edge][0]];
	_r[1] = rz[0][edge[1][num_edge][1]], _z[1] = rz[1][edge[1][num_edge][1]];
	_r[2] = rz[0][edge[1][num_edge][2]], _z[2] = rz[1][edge[1][num_edge][2]];
	double h = sqrt(pow(_r[0] - _r[2], 2) + pow(_z[0] - _z[2], 2));
	vector<double> b_s2(3, 0);
	vector<double> _theta(3);
	_theta[0] = teta(edge[1][num_edge][3], _r[0], _z[0], t);
	_theta[1] = teta(edge[1][num_edge][3], _r[1], _z[1], t);
	_theta[2] = teta(edge[1][num_edge][3], _r[2], _z[2], t);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			b_s2[i] += (A_1[i][j] * _r[0] + A_2[i][j] * _r[1] + A_3[i][j] * _r[2]) * _theta[j];
		}
		b_s2[i] *= h * beta() / double(420);
	}
	AddLocalVecBound(1, num_edge, b_s2);
}
void AddLocalMartBound_3(int num_edge, vector<vector<double>> matr) { // добавление матрицы из третьего краевого в глобальную
	for (int i = 0; i < 3; i++) {
		di[edge[2][num_edge][i]] += matr[i][i];
	}
	for (int i = 0; i < 3; i++) {
		int i_beg = ig[edge[2][num_edge][i]];
		for (int j = 0; j < i; j++) {
			int i_end = ig[edge[2][num_edge][i] + 1];
			while (jg[i_beg] != edge[2][num_edge][j]) {
				int ind = (i_beg + i_end) / 2;
				if (jg[ind] <= edge[2][num_edge][j]) i_beg = ind;
				else i_end = ind;
			}
			ggl[i_beg] += matr[i][j];
			i_beg++;
		}
	}
}

void AddLocalVecBound(int bound_condit, int num_edge, vector<double> vec) { // добавление вектора из третьего или второго краевого в глобальный
	for (int i = 0; i < 3; i++) {
		b[edge[bound_condit][num_edge][i]] += vec[i];
	}
}
void ConsiderBoundConditThirdType(int num_edge, double t) { // учет краевых условий третьего типа
	vector<double> _r(3), _z(3);
	_r[0] = rz[0][edge[2][num_edge][0]], _z[0] = rz[1][edge[2][num_edge][0]];
	_r[1] = rz[0][edge[2][num_edge][1]], _z[1] = rz[1][edge[2][num_edge][1]];
	_r[2] = rz[0][edge[2][num_edge][2]], _z[2] = rz[1][edge[2][num_edge][2]];
	double h = sqrt(pow(_r[0] - _r[2], 2) + pow(_z[0] - _z[2], 2));
	vector<vector<double>> _A(3);
	_A[0].resize(1);
	_A[1].resize(2);
	_A[2].resize(3);
	vector<double> b_s3(3, 0);
	vector<double> _u_beta(3);
	_u_beta[0] = u_beta(edge[2][num_edge][3], _r[0], _z[0], t);
	_u_beta[1] = u_beta(edge[2][num_edge][3], _r[1], _z[1], t);
	_u_beta[2] = u_beta(edge[2][num_edge][3], _r[2], _z[2], t);
	for (int i = 0; i < 3; i++) {
		for (int j = i; j < 3; j++) {
			_A[j][i] = h * beta() / double(420) * (A_1[j][i] * _r[0] + A_2[j][i] * _r[1] + A_3[j][i] * _r[2]);
		}
	}
	AddLocalMartBound_3(num_edge, _A);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < i; j++) {
			b_s3[i] += _A[i][j] * _u_beta[j];
		}
		for (int j = i; j < 3; j++) {
			b_s3[i] += _A[j][i] * _u_beta[j];
		}
	}
	AddLocalVecBound(2, num_edge, b_s3);
}
void BuildBoundMatrices() { // построение "краевых" матриц
	A_1.assign(3, vector<double>(3));
	A_2.assign(3, vector<double>(3));
	A_3.assign(3, vector<double>(3));
	A_1[0][0] = 39, A_2[0][0] = 20, A_3[0][0] = -3;
	A_1[0][1] = A_1[1][0] = 20, A_2[0][1] = A_2[1][0] = 16, A_3[0][1] = A_3[1][0] = -8;
	A_1[0][2] = A_1[2][0] = -3, A_2[0][2] = A_2[2][0] = -8, A_3[0][2] = A_3[2][0] = -3;
	A_1[1][1] = 16, A_2[1][1] = 192, A_3[1][1] = 16;
	A_1[1][2] = A_1[2][1] = -8, A_2[1][2] = A_2[2][1] = 16, A_3[1][2] = A_3[2][1] = 20;
	A_1[2][2] = -3, A_2[2][2] = 20, A_3[2][2] = 39;
}
void BoundCondit() { // краевые условия
	edge.resize(3);
	int num_func_1 = 0, num_func_2 = 0, num_func_3 = 0;
	for (int i = 0; i < 4; i++) {
		int start_point = GetStartPoint(i);
		int step = GetStepEndEl(i);
		switch (choise[i])
		{
		case 1:
			for (int j = 0; j < NodesCountR; j++) {
				edge[0].resize(NodesCountR * (num_func_1 + 1));
				edge[0][j + NodesCountR * num_func_1].push_back(start_point + j * step / 2);
				edge[0][j + NodesCountR * num_func_1].push_back(num_func_1);
			}
			num_func_1++;
			break;
		case 2:
			for (int j = 0; j < nr; j++) {
				edge[1].resize(nr * (num_func_2 + 1));
				edge[1][j + nr * num_func_2].push_back(start_point);
				edge[1][j + nr * num_func_2].push_back(start_point + step / 2);
				edge[1][j + nr * num_func_2].push_back(start_point + step);
				edge[1][j + nr * num_func_2].push_back(num_func_2);
				start_point += step;
			}
			num_func_2++;
			break;
		case 3:
			for (int j = 0; j < nr; j++) {
				edge[2].resize(nr * (num_func_3 + 1));
				edge[2][j + nr * num_func_3].push_back(start_point);
				edge[2][j + nr * num_func_3].push_back(start_point + step / 2);
				edge[2][j + nr * num_func_3].push_back(start_point + step);
				edge[2][j + nr * num_func_3].push_back(num_func_3);
				start_point += step;
			}
			num_func_3++;
			break;
		}
	}
}
void ConsiderBoundCondit(double t) { // учет всех краевых
	for (int i = 0; i < edge[2].size(); i++) { // учет третьих краевых
		ConsiderBoundConditThirdType(i, t);
	}
	for (int i = 0; i < edge[1].size(); i++) { // учет вторых краевых
		ConsiderBoundConditSecType(i, t);
	}
	for (int i = 0; i < edge[0].size(); i++) { // учет первых краевых(не в верхнем цикле, тк должен быть в самом конце)
		ConsiderBoundConditFirstType(i, t);
	}
}

#pragma endregion

#pragma region Реализация функций параметров уравнения
double lambda(double z) {
	return 1;
}
double RightFunction(double r, double z, double time) {
	return 21 * exp(3 * time);
}
double sigma(double r, double z) {
	return 1;
}
double xi_d(double r, double z) {
	return 2;
}
#pragma endregion
#pragma region Реализация построений Альфа и Мю
void CalcVecAlpha(unsigned int indexFE) {
	/*
	int* numfe[4],	     // номера узлов конечных элементов
       * numsubarea[2];  // номера подобластей конечных элементов
      double* rz[2];	     // координаты узлов
      int NodesCount;      //количество узлов
      int FECount;	     //количество КЭ
	*/
	if (FECount != NULL and NodesCount != NULL) {
		if (indexFE < FECount) {
			alpha[0] = (get<0>(GetCordIK(2, indexFE)) - get<0>(GetCordIK(0, indexFE))) * (get<1>(GetCordIK(6, indexFE)) - get<1>(GetCordIK(0, indexFE))) -
				(get<1>(GetCordIK(2, indexFE)) - get<1>(GetCordIK(0, indexFE))) * (get<0>(GetCordIK(6, indexFE)) - get<0>(GetCordIK(0, indexFE)));
			alpha[1] = (get<0>(GetCordIK(2, indexFE)) - get<0>(GetCordIK(0, indexFE))) * (get<1>(GetCordIK(8, indexFE)) - get<1>(GetCordIK(6, indexFE))) -
				(get<1>(GetCordIK(2, indexFE)) - get<1>(GetCordIK(0, indexFE))) * (get<0>(GetCordIK(8, indexFE)) - get<0>(GetCordIK(6, indexFE)));
			alpha[2] = (get<0>(GetCordIK(8, indexFE)) - get<0>(GetCordIK(2, indexFE))) * (get<1>(GetCordIK(6, indexFE)) - get<1>(GetCordIK(0, indexFE))) -
				(get<1>(GetCordIK(8, indexFE)) - get<1>(GetCordIK(2, indexFE))) * (get<0>(GetCordIK(6, indexFE)) - get<0>(GetCordIK(0, indexFE)));
		}
		else throw "Ошибка: указан индекс за диапазоном";
	}
	else throw "Ошибка: вычисление вектора альфа невозможно без инициализации области";
}
void CalcVecMyu(unsigned int indexFE) {
	/*
	int* numfe[4],	     // номера узлов конечных элементов
	   * numsubarea[2];  // номера подобластей конечных элементов
	  double* rz[2];	     // координаты узлов
	  int NodesCount;      //количество узлов
	  int FECount;	     //количество КЭ
	*/
	if (FECount != NULL and NodesCount != NULL) {
		if (indexFE < FECount) {
			Myu[0] = get<0>(GetCordIK(6, indexFE)) - get<0>(GetCordIK(0, indexFE));
			Myu[1] = get<0>(GetCordIK(2, indexFE)) - get<0>(GetCordIK(0, indexFE));
			Myu[2] = get<1>(GetCordIK(6, indexFE)) - get<1>(GetCordIK(0, indexFE));
			Myu[3] = get<1>(GetCordIK(2, indexFE)) - get<1>(GetCordIK(0, indexFE));
			Myu[4] = get<0>(GetCordIK(0, indexFE)) - get<0>(GetCordIK(2, indexFE)) - get<0>(GetCordIK(6, indexFE)) + get<0>(GetCordIK(8, indexFE));
			Myu[5] = get<1>(GetCordIK(0, indexFE)) - get<1>(GetCordIK(2, indexFE)) - get<1>(GetCordIK(6, indexFE)) + get<1>(GetCordIK(8, indexFE));
		}
		else throw "Ошибка: указан индекс за диапазоном";
	}
	else throw "Ошибка: вычисление вектора мю невозможно без инициализации области";
}
#pragma endregion

#pragma region Реализация Гаусса
double Gauss2rod(function<double(double, double, tuple<int, int, int>, double)> Ifun, tuple<int, int, int> FuncIndex, double time) {
	// Расчёт констант
	double Tau[3]; Tau[0] = 8.0 / 9; Tau[1] = Tau[2] = 5.0 / 9;
	double t[3] = {0, 0.77459666924148337, -0.77459666924148337 };
	double result = 0.0;
	//Расчёт метода Гаусса
	for (int i = 0; i < 3; i++) {
		double PsyI = (1 + t[i]) / 2;
		for (int j = 0; j < 3; j++) {
			double EtaJ = (1 + t[j]) / 2;
			result += Tau[i] * Tau[j] * Ifun(PsyI, EtaJ, FuncIndex, time);
		}
	}
	return result/double(4);
}
double Gauss2rod(function<double(double, double, tuple<int, int, int>)> Ifun, tuple<int, int, int> FuncIndex) {
	// Расчёт констант
	double Tau[3]; Tau[0] = 8.0 / 9; Tau[1] = Tau[2] = 5.0 / 9;
	double t[3] = { 0, 0.77459666924148337, -0.77459666924148337 };
	double result = 0.0;
	//Расчёт метода Гаусса
	for (int i = 0; i < 3; i++) {
		double PsyI = (1 + t[i]) / 2;
		for (int j = 0; j < 3; j++) {
			double EtaJ = (1 + t[j]) / 2;
			result += Tau[i] * Tau[j] * Ifun(PsyI, EtaJ, FuncIndex);
		}
	}
	return result / double(4);
}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора
void GeneratePortrait() { // генерация портрета матрицы
	vector<set<int>> list(NodesCount);
	for (int i = 0; i < FECount; i++) {
		for (int j = 8; j >= 0; j--) {
			for (int m = j - 1; m >= 0; m--) {
				list[numfe[j][i]].insert(numfe[m][i]);
			}
		}
	}
	ig.resize(NodesCount + 1, 0);
	for (int i = 2; i <= NodesCount; i++) {
		ig[i] = ig[i - 1] + list[i - 1].size();
	}
	jg.resize(ig[NodesCount]);
	auto iter = jg.begin();
	for (int i = 0; i < NodesCount; i++) {
		copy(list[i].begin(), list[i].end(), iter);
		iter += list[i].size();
	}
	ggl.resize(ig[NodesCount]);
	di.resize(NodesCount);
	b.resize(NodesCount);
	L_sq.resize(ig[NodesCount]);
	di_sq.resize(NodesCount);
}
void add_local_to_global(int IndexFE, int debag = false)//занесение локальной матрицы и локального вектора правой части в  глобальные
{
	for (int i = 0; i < 9; i++)
	{
		int ind = numfe[i][IndexFE];
		di[ind] += Gk[i];
		b[ind] += bk[i];
	}
	for (int i = 0, p = 9; i < 9; i++)
	{
		int d = numfe[i][IndexFE];//глобальный номер узла ii-го КЭ с локальным номером i
		int ibeg = ig[d];
		for (int j = 0; j < i; j++, p++)
		{
			int iend = ig[d + 1] - 1;
			while (jg[ibeg] != numfe[j][IndexFE])//дихотомия
			{
				int ind = (ibeg + iend) / 2 + 1;
				if (jg[ind] <= numfe[j][IndexFE]) ibeg = ind;
				else iend = (iend == ind) ? ind - 1 : ind;
			}
			ggl[ibeg] += Gk[p];
			ibeg++;
		}
	}
}
void LocalCalcMKXi(int IndexFE) {
	if (FECount == NULL || NodesCount == NULL)
		throw "Ошибка: Нельзя посчитать локальную матрицу не зная область";
	//Создаём функцию по которой интегрируем
	function<double(double, double, tuple<int, int, int> Index)> IntegratedFunction = [](double Psy, double Eta, tuple<int, int, int> IndexFunc) { return RPsyEta(Psy, Eta, get<2>(IndexFunc)) * xi_d(RPsyEta(Psy, Eta, get<2>(IndexFunc)), ZPsyEta(Psy, Eta, get<2>(IndexFunc))) * W[get<0>(IndexFunc)](Psy, Eta) * W[get<1>(IndexFunc)](Psy, Eta) * (alpha[0] + alpha[1] * Psy + alpha[2] * Eta); };
	//Вычисляем M для первого конечного элемента
	for (int i = 0; i < 9; i++)
		MkXi[i] = sign(alpha[0]) * Gauss2rod(IntegratedFunction, make_tuple(i, i, IndexFE));
	for (int i = 1, p = 9; i < 9; i++)
		for (int j = 0; j < i; j++, p++)
			MkXi[p] = sign(alpha[0]) * Gauss2rod(IntegratedFunction, make_tuple(i, j, IndexFE));
}
void LocalCalcMKSig(int IndexFE) {
	if (FECount == NULL || NodesCount == NULL)
		throw "Ошибка: Нельзя посчитать локальную матрицу не зная область";
	//Создаём функцию по которой интегрируем
	function<double(double, double, tuple<int, int, int> Index)> IntegratedFunction = [](double Psy, double Eta, tuple<int, int, int> IndexFunc) { return RPsyEta(Psy, Eta, get<2>(IndexFunc)) * sigma(RPsyEta(Psy, Eta, get<2>(IndexFunc)), ZPsyEta(Psy, Eta, get<2>(IndexFunc))) * W[get<0>(IndexFunc)](Psy, Eta) * W[get<1>(IndexFunc)](Psy, Eta) * (alpha[0] + alpha[1] * Psy + alpha[2] * Eta); };
	//Вычисляем M для первого конечного элемента
	for (int i = 0; i < 9; i++)
		MkSig[i] = sign(alpha[0]) * Gauss2rod(IntegratedFunction, make_tuple(i, i, IndexFE));
	for (int i = 1, p = 9; i < 9; i++)
		for (int j = 0; j < i; j++, p++)
			MkSig[p] = sign(alpha[0]) * Gauss2rod(IntegratedFunction, make_tuple(i, j, IndexFE));
}
void LocalCalcGK(int IndexFE) {
	if (FECount == NULL || NodesCount == NULL)
		throw "Ошибка: Нельзя посчитать локальную матрицу не зная область";
	//Создаём функцию по которой интегрируем
	function<double(double, double, tuple<int, int, int>)> IntegretedFunc = [](double Psy, double Eta, tuple<int, int, int> IndexFunc) { return lambda(ZPsyEta(Psy, Eta, get<2>(IndexFunc))) * RPsyEta(Psy, Eta, get<2>(IndexFunc)) * ((D[0][get<0>(IndexFunc)](Psy, Eta) * (Myu[5] * Psy + Myu[2]) - D[1][get<0>(IndexFunc)](Psy, Eta) * (Myu[5] * Eta + Myu[3])) * (D[0][get<1>(IndexFunc)](Psy, Eta) * (Myu[5] * Psy + Myu[2]) - D[1][get<1>(IndexFunc)](Psy, Eta) * (Myu[5] * Eta + Myu[3])) + (D[1][get<0>(IndexFunc)](Psy, Eta) * (Myu[4] * Eta + Myu[1]) - D[0][get<0>(IndexFunc)](Psy, Eta) * (Myu[4] * Psy + Myu[0])) * (D[1][get<1>(IndexFunc)](Psy, Eta) * (Myu[4] * Eta + Myu[1]) - D[0][get<1>(IndexFunc)](Psy, Eta) * (Myu[4] * Psy + Myu[0]))) / (alpha[0] + alpha[1]*Psy + alpha[2]*Eta); };
	//Выделяем память
	Gk = new double[45];
	//Вычисляем G для IndexFE конечного элемента
	for (int i = 0; i < 9; i++)
		Gk[i] = sign(alpha[0]) * Gauss2rod(IntegretedFunc, make_tuple(i, i, IndexFE));
	for (int i = 1, p = 9; i < 9; i++)
		for (int j = 0; j < i; j++, p++)
			Gk[p] = sign(alpha[0]) * Gauss2rod(IntegretedFunc, make_tuple(i, j, IndexFE));
}
double RPsyEta(double Psy, double Eta, int IndexFE) {
	return (1 - Psy) * (1 - Eta) * get<0>(GetCordIK(0, IndexFE)) + Psy * (1 - Eta) * get<0>(GetCordIK(2, IndexFE)) +
		(1 - Psy) * Eta * get<0>(GetCordIK(6, IndexFE)) + Psy * Eta * get<0>(GetCordIK(8, IndexFE));
}
double ZPsyEta(double Psy, double Eta, int IndexFE) {
	return (1 - Psy) * (1 - Eta) * get<1>(GetCordIK(0, IndexFE)) + Psy * (1 - Eta) * get<1>(GetCordIK(2, IndexFE)) +
		(1 - Psy) * Eta * get<1>(GetCordIK(6, IndexFE)) + Psy * Eta * get<1>(GetCordIK(8, IndexFE));
}
void SumVector(double*& FirstVec, double* SecondVec, int n) {
	for (int i = 0; i < n; i++) {
		FirstVec[i] += SecondVec[i];
	}
}
void SubVector(double*& FirstVec, double* SecondVec, int n) {
	for (int i = 0; i < n; i++) {
		FirstVec[i] -= SecondVec[i];
	}
}
void EqVector(double*& FirstVec, double* SecondVec, int n) {
	for (int i = 0; i < n; i++) {
		FirstVec[i] = SecondVec[i];
	}
}
void MultMKVector(double*& FirstVec, double* SecVector, double * Mk) {
	for (int i = 0; i < 9; i++) {
		FirstVec[i] = Mk[i] * SecVector[i];
	}
	for (int i = 0, p = 9; i < 9; i++)
		for (int j = 0; j < i; j++, p++) {
			FirstVec[i] += Mk[p] * SecVector[j];
			FirstVec[j] += Mk[p] * SecVector[i];
		}
}
void MultVectorConstant(double*& FirstVec, double Constant, int n) {
	for (int i = 0; i < n; i++) FirstVec[i] *= Constant;
}
void LocalCalcVectorFEK(int IndexFE, double time) { // реализовано
	if (FECount == NULL || NodesCount == NULL)
		throw "Ошибка: Нельзя посчитать вектор правой части без заданной области";
	//Создаём функцию по которой интегрируем
	function<double(double, double, tuple<int, int, int>, double)> IntegretedFunc = [] (double Psy, double Eta, tuple<int, int, int> IndexFunc, double time) { return RightFunction(RPsyEta(Psy, Eta, get<2>(IndexFunc)), ZPsyEta(Psy, Eta, get<2>(IndexFunc)), time) * RPsyEta(Psy, Eta, get<2>(IndexFunc)) * W[get<1>(IndexFunc)](Psy, Eta) * W[get<0>(IndexFunc)](Psy, Eta) * (alpha[0] + alpha[1] * Psy + alpha[2] * Eta); };
	for (int i = 0; i < 9; i++) bk[i] = 0.;
	CalcVecAlpha(IndexFE);
	for (int i = 0; i < 9; i++)
		for (int j = 0; j < 9; j++)
			bk[i] += sign(alpha[0]) * Gauss2rod(IntegretedFunc, make_tuple(i, j, IndexFE), time);
}

void InitVectors() {
	MkXi = new double[45];
	MkSig = new double[45];
	Gk = new double[45];
	bk = new double[9];
}
static void Calc_deltat_3SCHEME(double*& Vector, int j) {
	if (j < 2) {
		throw "j < 2";
	}
	if (num_t < 2 || j >= num_t) {
		throw "Wrong j";
	}
	Vector[0] = tj[j] - tj[j - 1];
	Vector[1] = tj[j - 1] - tj[j - 2];
	Vector[2] = tj[j] - tj[j - 2];
}

void GetQpj(double *& FirstVec, int IndexFe, int mj) {
	int index;
	for (int i = 0; i < 9; i++) {
		index = numfe[i][IndexFe];
		FirstVec[i] = qj[mj][index];
	}
}
void GetJLayer_3SCHEME(bool debag, int j) {
	if (j < 2) {
		throw "j < 2";
	}
	InitVectors();
	double* TempM = new double[45];
	double* deltat = new double[3];
	Calc_deltat_3SCHEME(deltat, j);
	double con1 = 2 / (deltat[2] * deltat[0]);
	double con2 = 1 / deltat[0] + 1 / deltat[2];
	double con3 = 2 / (deltat[1] * deltat[0]);
	double con4 = 2 / (deltat[2] * deltat[1]);
	double con5 = deltat[2] / (deltat[1] * deltat[0]);
	double con6 = deltat[0] / (deltat[1] * deltat[2]);
	double* qp1 = new double[9];
	double* qp2 = new double[9];
	for (int p = 0; p < FECount; p++)//проходим по всем КЭ
	{
		CalcVecAlpha(p);
		CalcVecMyu(p);
		LocalCalcGK(p);
		LocalCalcMKXi(p);
		LocalCalcMKSig(p);
		if (debag) {
			printf("Gk Mk до суммирования\nGk: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", Gk[i]);
			printf("\nMkXi: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", MkXi[i]);
			printf("\nMkSig: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", MkSig[i]);
			printf("\n");
		}
		EqVector(TempM, MkXi, 45);
		MultVectorConstant(TempM, con1, 45);
		SumVector(Gk, TempM, 45);
		EqVector(TempM, MkSig, 45);
		MultVectorConstant(TempM, con2, 45);
		SumVector(Gk, TempM, 45);

		if (debag) {
			printf("Gk после суммирования\nGk: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", Gk[i]);
			printf("\n");
		}
		LocalCalcVectorFEK(p, tj[j]);
		GetQpj(qp1, p, 1);
		MultVectorConstant(qp1, con3, 9);
		GetQpj(qp2, p, 0);
		MultVectorConstant(qp2, con4, 9);
		SubVector(qp1, qp2, 9);
		MultMKVector(qp2, qp1, MkXi);
		SumVector(bk, qp2, 9);
		GetQpj(qp1, p, 1);
		MultVectorConstant(qp1, con5, 9);
		GetQpj(qp2, p, 0);
		MultVectorConstant(qp2, con6, 9);
		SubVector(qp1, qp2, 9);
		MultMKVector(qp2, qp1, MkSig);
		SumVector(bk, qp2, 9);
		if (debag) {
			printf("bk: ");
			for (int i = 0; i < 9; i++)
				printf("%f ", bk[i]);
			printf("\n");
		}
		add_local_to_global(p);
	}
	delete[] MkXi; delete[] MkSig; delete[] Gk; delete[] bk;
}

void Calc_deltat_4SCHEME(double*& Vector, int j) {
	if (j < 3) {
		throw "j < 3";
	}
	if (num_t < 3 || j >= num_t) {
		throw "Wrong j";
	}
	Vector[0] = tj[j - 3] - tj[j - 2];
	Vector[1] = tj[j - 3] - tj[j - 1];
	Vector[2] = tj[j - 3] - tj[j];
	Vector[3] = tj[j - 2] - tj[j - 1];
	Vector[4] = tj[j - 2] - tj[j];
	Vector[5] = tj[j - 1] - tj[j];
}

void GetJLayer_4SCHEME(bool debag, int j) {
	if (j < 3) {
		throw "j < 3";
	}
	InitVectors();
	double* TempM = new double[45];
	double* deltat = new double[6];
	Calc_deltat_4SCHEME(deltat, j);
	double con1 = deltat[4] * deltat[5]/(deltat[0] * deltat[1] * deltat[2]);
	double con2 = -deltat[2] * deltat[5] / (deltat[0] * deltat[3] * deltat[4]);
	double con3 = deltat[2]*deltat[4]/(deltat[1]*deltat[3]*deltat[5]);
	double con4 = -(1 / deltat[2] + 1 / deltat[4] + 1 / deltat[5]);
	double con5 = -2 * (deltat[4] + deltat[5]) / (deltat[0] * deltat[1] * deltat[2]);
	double con6 = 2 * (deltat[2] + deltat[5]) / (deltat[0] * deltat[3] * deltat[4]);
	double con7 = -2 * (deltat[2] + deltat[4]) / (deltat[1] * deltat[3] * deltat[5]);
	double con8 = 2 * (deltat[2] + deltat[5] + deltat[4]) / (deltat[2] * deltat[4] * deltat[5]);
	double* qp1 = new double[9];
	double* qp2 = new double[9];
	for (int p = 0; p < FECount; p++)//проходим по всем КЭ
	{
		CalcVecAlpha(p);
		CalcVecMyu(p);
		LocalCalcGK(p);
		LocalCalcMKXi(p);
		LocalCalcMKSig(p);
		if (debag) {
			printf("Gk Mk до суммирования\nGk: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", Gk[i]);
			printf("\nMkXi: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", MkXi[i]);
			printf("\nMkSig: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", MkSig[i]);
			printf("\n");
		}
		EqVector(TempM, MkXi, 45);
		MultVectorConstant(TempM, con8, 45);
		SumVector(Gk, TempM, 45);
		EqVector(TempM, MkSig, 45);
		MultVectorConstant(TempM, con4, 45);
		SumVector(Gk, TempM, 45);

		if (debag) {
			printf("Gk после суммирования\nGk: ");
			for (int i = 0; i < 45; i++)
				printf("%f ", Gk[i]);
			printf("\n");
		}
		LocalCalcVectorFEK(p, tj[j]);
		GetQpj(qp1, p, 1);
		MultVectorConstant(qp1, con6, 9);
		GetQpj(qp2, p, 0);
		MultVectorConstant(qp2, con5, 9);
		SumVector(qp1, qp2, 9);
		GetQpj(qp2, p, 2);
		MultVectorConstant(qp2, con7, 9);
		SumVector(qp1, qp2, 9);
		MultMKVector(qp2, qp1, MkXi);
		SubVector(bk, qp2, 9);

		GetQpj(qp1, p, 1);
		MultVectorConstant(qp1, con2, 9);
		GetQpj(qp2, p, 0);
		MultVectorConstant(qp2, con1, 9);
		SumVector(qp1, qp2, 9);
		GetQpj(qp2, p, 2);
		MultVectorConstant(qp2, con3, 9);
		SumVector(qp1, qp2, 9);
		MultMKVector(qp2, qp1, MkSig);
		SubVector(bk, qp2, 9);
		if (debag) {
			printf("bk: ");
			for (int i = 0; i < 9; i++)
				printf("%f ", bk[i]);
			printf("\n");
		}
		add_local_to_global(p);
	}
	delete[] MkXi; delete[] MkSig; delete[] Gk; delete[] bk;
}

#pragma endregion
#pragma region Реалиация MSG function
double ScalarMult(vector<double>& first_vec, vector<double>& second_vec, int n) {
	double res = 0;
	for (int i = 0; i < n; i++) {
		res += first_vec[i] * second_vec[i];
	}
	return res;
}

void MultMartVec(vector<double>& input_vec, vector<double>& res_vec, int n) {
	for (int i = 0; i < n; i++) {
		res_vec[i] = di[i] * input_vec[i];
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			int j = jg[k];
			res_vec[i] += ggl[k] * input_vec[j];
			res_vec[j] += ggl[k] * input_vec[i];
		}
	}
}

void LU_sq(int n) {
	for (int i = 0; i < n; i++) {
		double s_d = 0;
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			double s_l = 0;
			int j = jg[k];
			int j0 = ig[j];
			int j1 = ig[j + 1];
			int ki = ig[i];
			int kj = j0;
			for (; ki < k && kj < j1;) {
				int jl = jg[ki];
				int ju = jg[kj];
				if (jl == ju) {
					s_l += L_sq[kj] * L_sq[ki];
					ki++, kj++;
				}
				else if (jl < ju) ki++;
				else kj++;
			}
			L_sq[k] = (ggl[k] - s_l) / di_sq[j];
			s_d += L_sq[k] * L_sq[k];
		}
		di_sq[i] = sqrt(di[i] - s_d);
	}
}


void Ly_f(vector<double>& f, int n) {
	for (int i = 0; i < n; i++) {
		double y = 0;
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			int j = jg[k];
			y += L_sq[k] * f[j];
		}
		f[i] = (f[i] - y) / di_sq[i];
	}
}

void Ly_f_transp(vector<double>& f, int n) {
	for (int i = n - 1; i >= 0; i--) {
		f[i] /= di_sq[i];
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			int j = jg[k];
			f[j] -= L_sq[k] * f[i];
		}
	}
}

void PrMinusR(vector<double>& r, int n) {
	for (int i = 0; i < n; i++) {
		r[i] = b[i] - r[i];
	}
}

void FirstEqualSecond(vector<double>& z, vector<double>& r, int n) {
	for (int i = 0; i < n; i++) {
		z[i] = r[i];
	}
}

double CulcResid(vector<double>& r, double norma_f, int n) {
	return sqrt(ScalarMult(r, r, n)) / norma_f;
}

void CulcX(vector<double>& x, vector<double>& z, double alpha, int n) {
	for (int i = 0; i < n; i++) {
		x[i] += alpha * z[i];
	}
}
void CulcR(vector<double>& r, vector<double>& Az, double alpha, int n) {
	for (int i = 0; i < n; i++) {
		r[i] -= alpha * Az[i];
	}
}

void CulcZ(vector<double>& Mult, vector<double>& z, double beta, int n) {
	for (int i = 0; i < n; i++) {
		z[i] = Mult[i] + beta * z[i];
	}
}

double CulcAlpha(vector<double>& Az, vector<double>& z, vector<double>& Mult, vector<double>& r, int n) {
	return ScalarMult(Mult, r, n) / ScalarMult(Az, z, n);
}

double CulcBeta(double numerator, double denominator) {
	return numerator / denominator;
}

void CulcMult(vector<double>& Mult, int n) {
	Ly_f(Mult, n);
	Ly_f_transp(Mult, n);
}

void LU_sq_MSG(vector<double>& x, vector<double>& r, vector<double>& z, vector<double>& Az, vector<double>& Mult, int n, double eps, int max_iter) {
	double resid, norma_f;
	int counter = 1;
	LU_sq(n);
	MultMartVec(x, r, n);
	PrMinusR(r, n);
	norma_f = sqrt(ScalarMult(b, b, n));
	resid = CulcResid(r, norma_f, n);
	cout << "Start resid: " << resid << endl;
	FirstEqualSecond(Mult, r, n);
	CulcMult(Mult, n);
	FirstEqualSecond(z, Mult, n);
	for (; counter < max_iter && resid > eps; counter++) {
		double alpha, beta, numerator, denominator;
		resid = 0;
		MultMartVec(z, Az, n);
		alpha = CulcAlpha(Az, z, Mult, r, n);
		denominator = ScalarMult(Mult, r, n);
		CulcX(x, z, alpha, n);
		CulcR(r, Az, alpha, n);
		FirstEqualSecond(Mult, r, n);
		CulcMult(Mult, n);
		numerator = ScalarMult(Mult, r, n);
		beta = CulcBeta(numerator, denominator);
		CulcZ(Mult, z, beta, n);
		resid = sqrt(ScalarMult(r, r, n)) / norma_f;
		cout << setprecision(15) << resid << " " << counter << endl;
	}
	MultMartVec(x, r, n);
	PrMinusR(r, n);
	resid = CulcResid(r, norma_f, n);
	cout << "End resid: " << resid << endl;
}

void ClearVectorsLU() {
	ggl.resize(ig[NodesCount]);
	di.resize(NodesCount);
	b.resize(NodesCount);
	L_sq.resize(ig[NodesCount]);
	di_sq.resize(NodesCount);
	for (int i = 0; i < ig[NodesCount]; i++) {
		ggl[i] = 0.0;
		L_sq[i] = 0.0;
	}
	for (int i = 0; i < NodesCount; i++) {
		di[i] = 0.0;
		b[i] = 0.0;
		di_sq[i] = 0.0;
	}
}
#pragma endregion
#pragma endregion
#pragma region Тестирование
void Test(vector<double> q, int layer) {
	vector<double> q_u(NodesCount, 0);
	for (int i = 0, p = 0; i < NodesCountR; i++) {
		for (int j = 0; j < NodesCountR; j++, p++) {
			q_u[p] = ug(0, rz[0][p], rz[1][p], tj[layer]);
		}
	}
	
	double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
	for (int i = 0; i < NodesCount; i++) {
		norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
		norm_vec_q_u += (q_u[i]) * (q_u[i]);
	}
	cout << endl;
	cout << "Относительная норма вектора погрешности полученного решения на слое " << layer + 1 << " при t = " << tj[layer] << " равна: " << endl;
	cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl;
	cout.precision(15);
	/*
	for (int i = 0; i < NodesCount; i++)
		cout << q[i] << " " << q_u[i] << endl;
	*/
}
#pragma endregion
int main()
{
	setlocale(LC_ALL, "Russian");
	EnterBoundCondit();
	InitBaseFunction();
	InitDРerivativeBaseFunction();
	LoadNet();
	LoadTNet();
	InitQJ();
	InitFirstQ();
	BoundCondit();
	BuildBoundMatrices();
	GeneratePortrait();
	GetJLayer_3SCHEME(false, 2);
	ConsiderBoundCondit(tj[2]);
	int max_iter = 1000;
	double eps = 1e-15;
	vector<double> r(NodesCount);
	vector<double> z(NodesCount);
	vector<double> Mult(NodesCount);
	vector<double> Az(NodesCount);
	vector<double> q(NodesCount, 0);
	LU_sq_MSG(q, r, z, Az, Mult, NodesCount, eps, max_iter);
	qj[2] = q;
	Test(q, 2);
	for (int j = 3; j < num_t; j++) {
		ClearVectorsLU();
		GetJLayer_4SCHEME(false, j);
		ConsiderBoundCondit(tj[j]);
		vector<double> r(NodesCount);
		vector<double> z(NodesCount);
		vector<double> Mult(NodesCount);
		vector<double> Az(NodesCount);
		vector<double> q(NodesCount, 0);
		LU_sq_MSG(q, r, z, Az, Mult, NodesCount, eps, max_iter);
		qj[0] = qj[1]; qj[1] = qj[2]; qj[2] = q;
		Test(q, j);
	}
	return 0;
}