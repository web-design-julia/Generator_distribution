#include "Header.h"
#include <math.h>
#include <iostream>
#include <vector>

Minstd::Minstd()
  : seed_current(seed_) {
}

void Minstd::At_the_start() {
  seed_current = seed_;
}

double Minstd::Generate_new_number() {
  seed_current = (a * seed_current + c) % m;
  return static_cast<double>(seed_current) / static_cast<double>(m);

}
//Распределение Бернулли
int Minstd::Bern(double p) {
  
  double r = Generate_new_number();
  if (r >= p) {
    return 1;
  }
  else {
  return  0;
  }
}


void Minstd::DistBern(double entry_p) {
  At_the_start();
  std::vector<int> v(2,0);
  for (int i(0); i < n; i += 1) {
    if (Bern(entry_p) == 0) {
      v[0] += 1;
    }
    else {
      v[1] += 1;
    }
  }
  
  //Вероятность
  std::vector<double> p(2, 0);
  p[0] = entry_p;
  p[1] = 1 - entry_p;
  std::cout << "Bern_Chi_test( "<<entry_p<<" ):"<< Chi_Test(v, p, n) << std::endl;
  Hi_tabl(p.size()-1);
  Distribution(v);
}

//Биноминальное распределение
int Minstd::Bin(int n_event, double entry_p) {
  //Некруткин "Моделирование распределений" стр.22
  double c = entry_p / (1 - entry_p);
  double s = pow((1 - entry_p), n_event);
  double p_current = s;
  int k(0);
  double random = Generate_new_number();
  while (random > s) {
    k += 1;
    p_current = p_current*c*(n_event - k + 1) / k;
    s += p_current;
  }
  return k;
}

int Minstd::Factorial(int n) {
  int f(1);
  for (int i(2); i < n + 1; i +=1) {
    f *= i;
  }
  return f;
}
//Число сочетаний n по k
double Minstd::C(int k, int n)
{
  return (Factorial(n) / (Factorial(k) *Factorial(n - k)));
}


void Minstd::DistBin(int n_event, double entry_p) {
  At_the_start();
  std::vector<int> vec(n_event + 1, 0);
  for (int i(0); i < n; i += 1) {
    vec[Bin(n_event, entry_p)] += 1;
  }
  std::vector<double> p(n_event + 1, 0);
  p[0] = pow(1 - entry_p, n_event);
  for (int i(1); i <= n_event; i += 1) {
    p[i] = C(i, n_event)*pow(entry_p, i)*pow(1 - entry_p, n_event - i);
  }
  std::cout << "Bin_Chi_test(" << n_event << "," << entry_p << "): " << Chi_Test(vec,p,n) << std::endl;
  Hi_tabl(n_event);
  Distribution(vec);
}

//Геометрическое распределение
int Minstd::Geom(double entry_p) {
  double s = entry_p, r = entry_p, q = 1 - entry_p;
  double random = Generate_new_number();
  int k(0);
  //Возвращаем минимальное значение k, 
  //для которого сумма больше случайного числа
  while (random > s) {
    k += 1;
    r *= q;
    s += r;
  }
  return k;
}

void Minstd::DistGeom(double entry_p) {
  At_the_start();
  std::vector<double> dist(n, 0);
  double x_max(0);
  //Найдем максимальное число
  for (int i(0); i < n; i += 1) {
    dist[i] = Geom(entry_p);
    if (dist[i] > x_max) x_max = dist[i];
  }
  //Сколько раз получилось число= dist[i]
  std::vector<int> v(x_max + 1, 0);
  for (int i(0); i < n; i += 1) {
    v[dist[i]]++;
  }
  //Вероятность
  std::vector<double> p(x_max + 1, 0);
  for (int i(0); i <= x_max; i += 1) {
    p[i] = entry_p*pow(1 - entry_p, i);
  }
  
  std::cout << "Geom(" << entry_p << "): " << Chi_Test(v,p,n) << std::endl;
  Hi_tabl(v.size()-1);
  Distribution(v);
}
//Гипергеометрическое распределение
int Minstd::HGeom(int N, int n_event, int K) {
  int sum(0), p = static_cast<double>(K) / static_cast<double>(N);
  for (int i(1); i <= n_event; i += 1) {
    if (Bern(p) && ++sum == K)
      return sum;
    p = (K - sum) / (N - i);
  }
  return sum;
}

void Minstd::DistHgeom(int N, int n_event, int K) {
  At_the_start();
  n = 3000;
  std::vector<int> dist(n, 0);
  
  dist[0] = HGeom(N, n_event, K);
  int min(dist[0]);
  int max = min;
  //Находим минимальное и максимальное значение
  for (int i(1); i < n; i += 1) {
    dist[i] = HGeom(N, n_event, K);
    if (dist[i] > max) max = dist[i];
    if (dist[i] < min) min = dist[i];
  }
  //Диапазон значений
  int range = max - min + 1;
  std::vector<int> v(range, 0);
  //Считаем кол-во появлений dist[i]
  for (int i(0); i < n; i +=1) {
    v[dist[i] - min] += 1;
  }
  //Вероятность
  std::vector<double> p(range, 0);
  for (int i(0); i < range; i += 1) {
    p[i] = (C(i + min, K)*C(n_event - i - min, N - K)) / C(n_event, N);
  }
  
  std::cout << "HGeom(" << N << "," << n_event << "," << K << "): " << Chi_Test(v,p,n)<< std::endl;
  Hi_tabl(range);
  Distribution(v);
}
//Распределение Пуассона
double Minstd::Pois(double lambda) {
  double s = exp(-lambda), r = exp(-lambda);
  int k(0);
  double random = Generate_new_number();
  while (random > s) {
    k++;
    r = r*(lambda / k);
    s += r;
  }
  return k;
}

void Minstd::DistPois(double lambda) {
  At_the_start();
  std::vector<double> dist(n, 0);
  int x_max(0);
  for (int i(0); i < n; i += 1) {
    dist[i] = Pois(lambda);
    if (dist[i] > x_max) x_max = dist[i];
  }
  std::vector<int> v(x_max + 1, 0);
  for (int i(0); i < n; i += 1) {
    v[dist[i]]++;
  }
  std::vector<double> p(x_max + 1, 0);
  for (int i(0); i < x_max + 1; i += 1) {
    p[i] = (exp(-lambda)*pow(lambda, i)) / Factorial(i);
  }
  
  std::cout << "Pois(" << lambda << "): " << Chi_Test(v, p,n) << std::endl;
  Hi_tabl(x_max);
  Distribution(v);
}

//Равномерное распределение
double  Minstd::U(double a, double b) {
  //F(x)=1/(b-a)*(x-a)
  return (b - a)*Generate_new_number() + a;
}

void Minstd::DistU(double a, double b) {
  At_the_start();
  const int partition{100};
  std::vector<int> v(partition, 0);
  double u = U(a, b);
  //делим промежуток на 100 частей 
  double step = (b - a) / partition;

  double sum_step = step;
  int x(0);
  for (int i(0); i < n; i += 1) {
    //Если U меньше a
    if (u < step + a) v[0] += 1;
    else {
      //Если U больше b
      if (u >a + (step * partition)) v[partition - 1] += 1;
      else {
        while (u >a + sum_step) {
          sum_step += step;
          x += 1;
        }
        v[x] += 1;
        sum_step = step;
        x = 0;
      }
    }
    u = U(a, b);

  }

  std::vector<double> p(partition, 0);
  for (int i(0); i < partition; i += 1) {
    p[i] = (step*(i + 1)) / (b - a) - (step*(i)) / (b - a);
  }

  std::cout << "U(" << a << "," << b << "): " << Chi_Test(v,p,n) << std::endl;
  Hi_tabl(v.size()-1);
  Distribution(v);
}

//Нормировка
//Теорема: при сложении 6 с.ч. с разным распределением
//получим число с нормальным распределением
double Minstd::N(double exp, double dis) {
  double sum(0);
  for (int i(0); i < 12; i += 1) sum += Generate_new_number();

  return (sum - 6.0) * sqrt(dis) + exp;
}
//Кобзарь "Прикладная_математическая_статистика" стр.30
double Minstd::Approx(double x) {
  return 1 - 0.852*exp(-(pow((x + 1.5774) / (2.0637), 2.34)));
}

void Minstd::DistN(double exp, double dis) {
  At_the_start();
  const int partition{1000};
  n = 3000;
  std::vector<double> dist(n, 0);
  dist[0] = N(exp, dis);
  double min(dist[0]);
  double max = min;
  for (int i(1); i < n; i += 1) {
    dist[i] = N(exp, dis);
    if (dist[i] > max) max = dist[i];
    if (dist[i] < min) min = dist[i];
  }
  double step = (max - min) / partition;
  double sum_step = 0;
  int x(0);
  std::vector<int> v(partition, 0);
  for (int i(1); i < n; i += 1) {
    while (dist[i] > min + sum_step) {
      sum_step += step;
      x += 1;
    }
    if (x <= partition - 1) {
      v[x] += 1;
    }
    else {
      v[partition - 1] += 1;
    }
    sum_step = 0;
    x = 0;
  }
  std::vector<double> p(partition, 0);
  double right(0), left(0);
  for (int i = 0; i < partition; i++) {
    right = (min + (i + 1)*step)*(1 - exp) / sqrt(dis);
    left = (min + (i)*step)*(1 - exp) / sqrt(dis);
    p[i] = abs(Approx(abs(right)) - Approx(abs(left)));
  }
  std::cout << "N(" << exp << "," << dis << "): " << Chi_Test(v, p, n) << std::endl;
  Hi_tabl(v.size()-1);
}
//Экспоненциальное распределение
double Minstd::Exp(double lambda) {
  //f(x)=e^(-lambda*x)
  return -(1 / lambda)*log(Generate_new_number());
}

void Minstd::DistExp(double lambda) {
  At_the_start();
  const int partition{ 1000 };
  n = 3000;
  std::vector<double> dist(n, 0);
  dist[0] = Exp(lambda);
  double min(dist[0]);
  double max = min;
  for (int i(1); i < n; i += 1) {
    dist[i] = Exp(lambda);
    if (dist[i]  > max) max = dist[i];
    if (dist[i]  < min) min = dist[i];
  }

  double step = (max - min) / partition;
  double sum_step(0);
  int x(0);

  std::vector<int> v(partition, 0);
  for (int i(1); i < n; i += 1) {
    while (dist[i] > min + sum_step) {
      sum_step += step;
      x += 1;
    }

    if (x <= partition - 1)
      v[x]++;
    else v[partition - 1] += 1;
    sum_step = 0;
    x = 0;
  }

  std::vector<double> p(partition, 0);

  for (int i = 0; i < partition; i++) {
    //F(x)=1-e^(-lambda*x)
    p[i] = (1 - exp((-(i + 1)*step)*lambda)) - (1 - exp((-i*step) * lambda));
  }

  std::cout << "Exp(" << lambda << "): " << Chi_Test(v, p, n) << std::endl;
  Hi_tabl(v.size()-1);


}

//Гамма-распределение Gamma(k,1), где к-целое число
double Minstd::Gamma(double param_k) {
	
	double sum(0);
	for (int i(0); i < param_k; i++) {
		sum = Exp(1);
	}
	return sum;
}

void Minstd::DistGamma(double param_k) {
	At_the_start();
	const int partition{ 1000 };
	n = 3000;
	std::vector<double> dist(n, 0);
	dist[0] = Gamma(param_k);
	double min(dist[0]);
	double max = min;
	for (int i(1); i < n; i += 1) {
		dist[i] = Gamma(param_k);
		if (dist[i]  > max) max = dist[i];
		if (dist[i]  < min) min = dist[i];
	}

	double step = (max - min) / partition;
	double sum_step(0);
	int x(0);

	std::vector<int> v(partition, 0);
	for (int i(1); i < n; i += 1) {
		while (dist[i] > min + sum_step) {
			sum_step += step;
			x += 1;
		}

		if (x <= partition - 1)
			v[x]++;
		else v[partition - 1] += 1;
		sum_step = 0;
		x = 0;
	}

	std::vector<double> p(partition, 0);

	for (int i = 0; i < partition; i++) {
		
		p[i] = Gauss_gamma_probability(param_k,i,i+1);
	}

	std::cout << "Gamma(" << param_k << "): " << Chi_Test(v, p, n) << std::endl;
	Hi_tabl(v.size() - 1);
	
}

//метод для вычисления значения несобственного интеграла
//по формуле Лагерра
//Gamma(k)=от 0 до 1 x^(alfa-1)*(1-x)^(beta-1)
double Minstd::Lagerr_gamma(double param_k)
{
	int n(4);
	std::vector<double> A = {0.603154,0.357419,0.0388879,0.000539295};
	std::vector<double> x = {0.322548,1.745761,4.536620,9.395071};

	double sum(0);
	double y(0);
	for (int i(0); i < n; i++)
	{
		y = pow(x[i],param_k-1);
		sum += A[i]*y;
	}
	
	return sum;
}

//метод для вычисления значения интеграла по формуле Гаусса
//Интеграл Gamma для вероятности
double Minstd::Gauss_gamma_probability(double param_k, double a, double b)
{
	int n(4);
	std::vector<double> t = { -0.86113631,-0.33998104, 0.86113631, 0.33998104 };
	std::vector<double> A = { 0.34785484, 0.65214516 };
	double sum(0);
	double x(0);
	double c(0);
	double y(0);
	for (int i(0); i < n; i++)
	{
		x = (b + a) / 2 + (b - a)*t[i] / 2;
		if(i>1) c = (b - a)*A[1] / 2;
		else 
			c = (b - a)*A[0] / 2;
		y = pow(x,param_k -1)*exp(-x)/Lagerr_gamma(param_k);
		sum += c*y;
	}

	return sum;
}


//Бета-распределение Beta(alfa,beta)
double Minstd::Beta(double alfa, double beta) {

	double gamma1 = Gamma(alfa);
	double gamma2 = Gamma(beta);

	return gamma1 / (gamma1 + gamma2);
}

void Minstd::DistBeta(double alfa, double beta) {
	At_the_start();
	const int partition{ 1000 };
	n = 3000;
	std::vector<double> dist(n, 0);
	dist[0] = Beta(alfa,beta);
	double min(dist[0]);
	double max = min;
	for (int i(1); i < n; i += 1) {
		dist[i] = Beta(alfa,beta);
		if (dist[i]  > max) max = dist[i];
		if (dist[i]  < min) min = dist[i];
	}

	double step = (max - min) / partition;
	double sum_step(0);
	int x(0);

	std::vector<int> v(partition, 0);
	for (int i(1); i < n; i += 1) {
		while (dist[i] > min + sum_step) {
			sum_step += step;
			x += 1;
		}

		if (x <= partition - 1)
			v[x]++;
		else v[partition - 1] += 1;
		sum_step = 0;
		x = 0;
	}

	std::vector<double> p(partition, 0);
	
	for (int i = 0; i < partition; i++) {
		//f(x)=(1-x)^(beta-1)*x^(alfa-1)/B(alfa,beta)
	p[i] = Gauss_beta_probability(alfa,beta,i,i+1);
	}

	std::cout << "Beta(" << alfa <<"," << beta << "): " << Chi_Test(v, p, n) << std::endl;
	Hi_tabl(v.size() - 1);
	
}

//метод для вычисления значения интеграла по формуле Гаусса
//B(alfa,beta)=от 0 до 1 x^(alfa-1)*(1-x)^(beta-1)
double Minstd::Gauss_beta(double alfa, double beta)
{
	int n(4);
	int a(0), b(1);
	std::vector<double> t = {-0.86113631,-0.33998104, 0.86113631, 0.33998104 };
	std::vector<double> A = {0.34785484, 0.65214516};
	double sum(0);
	double x(0);
	double c(0);
	double y(0);
	for (int i(0); i< n; i++)
	{
		x = (b + a) / 2 + (b - a)*t[i] / 2;
		if(i>1) c = (b - a)*A[1] / 2;
		else 
			c = (b - a)*A[0] / 2;
		y = pow(x, alfa - 1)*pow(1 - x, beta - 1);
		sum += c*y;
	}

	return sum;
}

//метод для вычисления значения интеграла по формуле Гаусса
//Интеграл Beta для вероятности
double Minstd::Gauss_beta_probability(double alfa, double beta, double a, double b)
{
	int n(4);
	std::vector<double> t = { -0.86113631,-0.33998104, 0.86113631, 0.33998104 };
	std::vector<double> A = { 0.34785484, 0.65214516 };
	double sum(0);
	double x(0);
	double c(0);
	double y(0);
	for (int i(0); i < n; i++)
	{
		x = (b + a) / 2 + (b - a)*t[i] / 2;
		if(i>1) c = (b - a)*A[1] / 2;
		else 
			c = (b - a)*A[0] / 2;
		y = pow(1 - x, beta - 1)*pow(x, alfa - 1) / Gauss_beta(alfa, beta);
		sum += c*y;
	}

	return sum;
}


double Minstd::Chi_Test(const std::vector<int> count,
  const std::vector<double> p, const unsigned long long n) {

  double ans(0.0);
  for (int i(0); i < count.size(); i += 1) {
    ans += pow(count[i] - n * p[i], 2) / (n * p[i]);
  }

  return ans;
}

void Minstd::Hi_tabl(const double n_event) {
	std::cout << "Hi_tabl" << std::endl;
	std::cout << " 25%  " << n_event+sqrt(2 * n_event)*(-0.674) + (double)(2 / 3)*(-0.674)*(-0.674) - (double)(2 / 3)
    << "  50%  " << n_event+sqrt(2 * n_event)*(0) + (double)(2 / 3)*(0) - (double)(2 / 3) 
    << "  75%  " << n_event+sqrt(2 * n_event)*(0.674) + (double)(2 / 3)*(0.674)*(0.674) - (double)(2 / 3) << std::endl;
	std::cout << std::endl;
}

void Minstd::Distribution(const std::vector<int> count) {
  for (int i(0); i < count.size(); i += 1) {
	  std::cout << count[i] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

int main()
{
  Minstd d;
  
  d.DistGamma(10);
  d.DistBeta(2, 8);

  system("pause");
}

