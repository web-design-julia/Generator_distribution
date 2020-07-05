#pragma once

class Minstd {
public:
	// Конструктор
	Minstd();
	// Сброс генератора к начальному состоянию
	void At_the_start();
	// Новое число
	double Generate_new_number();
	//Распределение Бернулли
	int Bern(double entry_p);
	void DistBern(double entry_p);
	//Биноминальное распределение
	int Bin(int n_event, double entry_p);
	int Factorial(int n);
	//Формула сочетания из n по k
	double C(int k, int n);
	void DistBin(int n_event, double entry_p);
	//Геометрическое распределение
	int Geom(double entry_p);
	void DistGeom(double entry_p);
	//Гипергеометрическое распределение
	int HGeom(int N, int n_event, int K);
	void DistHgeom(int N, int n_event, int K);
	//Распределение Пуассона
	double Pois(double lambda);
	void DistPois(double lambda);
	//Равномерное распределение
	double  U(double a, double b);
	void DistU(double a, double b);
	//Нормальное распределение
	double N(double exp, double dis);
	double Gamma(double param_k);
	void DistGamma(double entry_p);
	//Аппроксимация нормального распределения
	double Approx(double x);
	void DistN(double exp, double dis);
	//Экспоненциальное распределение
	double Exp(double lambda);
	void DistExp(double lambda);
	double Beta(double alfa, double beta);
	void DistBeta(double alfa, double beta);
	double Gauss_beta(double alfa, double beta);
	double Gauss_beta_probability(double alfa, double beta, double a, double b);
	double Lagerr_gamma(double param_k);
	double Gauss_gamma_probability(double param_k, double a, double b);
	//Хи-квадрат
	double Chi_Test(const std::vector<int> count,
		const std::vector<double> p, const unsigned long long n);
	//Табличное значение хи-квадрата
	void Hi_tabl(const double n_event);
	//Вывод распределения
	void Distribution(const std::vector<int> count);


private:
	const int m{ static_cast<int>(pow(2, 31) - 1) };
	const int a{ 48271 }, c{ 0 }, seed_{ 1 };
	//Кол-во сгенерированных чисел
	int n{ 1000000 };
	// Текущее случайное значение
	unsigned long long seed_current{ 0 };

};