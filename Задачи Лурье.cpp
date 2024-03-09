#include <iomanip>
#include <clocale>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <fixed/fixed.h>
#include <fixed/fixed_nonlinear_solver.h>
#include <pde_solvers/pde_solvers.h>  // подключение библиотеки pde_solvers (В.В. Южанин)
#include <gtest/gtest.h>
using namespace std;
using namespace pde_solvers;

/// @brief Данные по задаче 1
struct zadacha_1 {
	double L;		// длина трубопровода
	double D;		// внешний диаметр трубопровода
	double b;		// толщина стенки трубопровода
	double z_0;		// начальная высотная отметка профиля трубопровода 
	double z_l;		// конечная высотная отметка профиля трубопровода 
	double density; // плоность нефти
	double nu;		// кинематическая вязкость нефти 
	double Q;		// объемный расход [м^3/час]
	double p_k;		// давление в конце трубопровода 
	double d;		// внешний диаметр трубопровода [мм]
	double d_m;		// внешний диаметр трубопровода [м]
	double q;		// Объемный расход [м^3/с]

	// дальше данные для решения методом Эйлера

	double n;		//кол-во шагов стеки
	double h;		//шаг по координате расчетной сетки, в метрах


	zadacha_1() {

		L = 80e3;	// км
		D = 720;				// мм
		b = 10;					// мм
		z_0 = 50;				// м 
		z_l = 100;				// м
		density = 870;			// кг/м3
		nu = 15;				// сСт
		Q = 3500;				// м3/ч
		p_k = 0.6 * pow(10, 6); // Па

		d = D - 2 * b;			// мм
		d_m = d / 1000;			// м
		q = Q / 3600;			// м3/с

		// для решения методом Эйлера

		n = 100;
		

	}
};


/// @brief Данные по задаче 2
struct zadacha_2 {

	double Z; 
	double lamd_2;
	double epsilon;		// допустимая погрешность 
	double lamd; 
	double sher;		// относительная шероховатость 
	double pp_n;		// давление в начале 
	double pp_k;		// давление в конце
	

  zadacha_2(zadacha_1& iniz1) {

		Z = 0.15; // мм
		lamd_2 = 0.02; 
		epsilon = 5e-5; //1 * pow(10, (- 3.21));
		pp_n = 5 * pow(10, 6);
		pp_k = 0.8 * pow(10, 6); 

		sher = Z / iniz1.d;


	}
};


struct massiv {
	vector<double> parameter;
};


/// @brief Задача 1, QP
/// @param iniz1 ссылка на данные по задаче 1
/// @return 
/// 
TEST(Task_1, QP_Lurie) {

		const zadacha_1 iniz1;

		double V = 4 * iniz1.q / (M_PI * pow(iniz1.d_m, 2));    // скорость в трубопровде
		double Re = V * iniz1.d_m / (iniz1.nu * pow(10, -6));   // число Рейнольдса  
		double e = iniz1.nu * pow(10, -2) / iniz1.d;            // шероховатость 


		double resistance = hydraulic_resistance_isaev(Re, e);  // ф-ция расчета значения лямбды

		double p_n = (iniz1.p_k / (iniz1.density * M_G) + iniz1.z_0 - iniz1.z_l + resistance
			* (iniz1.L / iniz1.d_m * pow(V, 2) / 2 / M_G)) * (iniz1.density * M_G);	// значение распределения начльного давления, МПа

		double pogr = (6 * pow(10, 6) - p_n) / (6 * pow(10, 6)); // Вычисление расхождения в расчетах программы и М.В.Лурье

		cout << "Задача QP: " << p_n << "Па" << endl;

	double error = 0.1e6;
	EXPECT_NEAR(5.99e6, p_n, error);
	
}

TEST(Task_2, QP_EULER) {

	const zadacha_1 iniz1;

	double V = 4 * iniz1.q / (M_PI * pow(iniz1.d_m, 2));    // скорость в трубопровде
	double Re = V * iniz1.d_m / (iniz1.nu * pow(10, -6));   // число Рейнольдса  
	double e = iniz1.nu * pow(10, -2) / iniz1.d;            // шероховатость 


	double resistance = hydraulic_resistance_isaev(Re, e);  // ф-ция расчета значения лямбды


	double tw = resistance / 8 * iniz1.density * pow(V, 2);

	double h = iniz1.L / iniz1.n;

	double pressure_current;
	double pressure_previous = iniz1.p_k;

	for (int i = 1; i <= iniz1.n; i++) {

		pressure_current = pressure_previous - iniz1.h * (-4 * tw / iniz1.d_m - iniz1.density * M_G * ((iniz1.z_l - iniz1.z_0) / ((iniz1.n  - 1) * iniz1.h)));
		
		pressure_current = pressure_previous;

	}

	std::cout << "Задача QP по Эйлеру: " << pressure_current << " Па" << std::endl;

	double error = 0.01e6;
	EXPECT_NEAR(0.6e6, pressure_current, error);
}



/// @brief Задача 2, PP
/// @param iniz1, ссылка на данные по задаче 1
/// @param iniz2, ссылка на данные по задаче 2
/// @param V, скорость [м/с]
/// @param Re, число Рейнольдса
/// @param index, указатель счетчика 
/// @param p_n, давление в начале трубопровода [МПа]

TEST(Task_3, calculateFlowAndIterations) {
		zadacha_1 iniz1;
		zadacha_2 iniz2(iniz1);
		double V, Re;
		int index = 0;

		// Нахождение lamd_2 методом последовательных приближений 
		do {
			iniz2.lamd = iniz2.lamd_2; //предположили, что начльное значение lamd_2 = 0.02

			V = pow(((iniz1.d_m * 2 * M_G / iniz1.L * ((iniz2.pp_n - iniz2.pp_k) / (iniz1.density * M_G)
				+ iniz1.z_0 - iniz1.z_l)) / iniz2.lamd), 0.5);

			Re = V * iniz1.d_m / iniz1.nu / pow(10, -6);

			iniz2.lamd_2 = 0.11 * pow(((iniz2.sher + 68 / Re)), 0.25);	//уточняем наше предположение 

			index++;

		} while (abs(iniz2.lamd_2 - iniz2.lamd) > iniz2.epsilon);

		std::cout << "Скорость: " << V << std::endl;
		std::cout << "Количество итераций: " << index << std::endl;

		double Q_s = (M_PI * pow(iniz1.d_m, 2)) * V / 4;		//	объемный расход, м^3/c
		double Q_h = Q_s * 3600;								//  объемный расход, м^3/ч		

		std::cout << "Расход: " << Q_h << std::endl;

	
}


TEST (Task_4, Newton){

	setlocale(LC_ALL, "Russian");

	zadacha_1 iniz1;
	zadacha_2 iniz2(iniz1);

/// @brief задача РР Ньютон. Класс, для системы размерности <1>
// <1> - Размерность системы уравнений
	class Newton : public fixed_system_t<1> {

		/// @brief ссылка на переменные в структурах 

		const zadacha_1 iniz1;
		const zadacha_2 iniz2;

		using fixed_system_t<1>::var_type;

	public:

		/// @brief конструктор
		/// @param iniz1 данные по здаче 1
		/// @param iniz2 данные по задаче 2
		Newton(const zadacha_1& iniz1, const zadacha_2& iniz2) : iniz1(iniz1), iniz2(iniz2) {}

		/// @brief Задание функции невязок
		/// @param V неизвестная величина
		/// @return ищем неизвестную V из уравнения Бернулли 
		var_type residuals(const var_type& V)

		{
			double Re = V * iniz1.d_m / (iniz1.nu * pow(10, -6));
			double resistance = hydraulic_resistance_altshul(Re, iniz2.sher);

			return
			{
			 (iniz2.pp_n / (iniz1.density * M_G) + iniz1.z_0) - (iniz2.pp_k / (iniz1.density * M_G) + iniz1.z_l)
				- resistance * iniz1.L / iniz1.d_m * pow(V,2) / 2 / M_G
			};
		};
	};
		// Создание экземпляра класса, который и будет решаемой системой

		Newton test(iniz1, iniz2);

		// Задание настроек решателя по умолчанию

		fixed_solver_parameters_t<1, 0> parameters;

		// Создание структуры для записи результатов расчета

		fixed_solver_result_t<1> result;

		// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
		// { 0.4 } - Начальное приближение скорости 1

		fixed_newton_raphson<1>::solve_dense(test, { 0.1 }, parameters, & result);

		std::cout << "Скорость по Ньютону = " << result.argument << std::endl;

		double Q_s = (M_PI * pow(iniz1.d_m, 2)) * result.argument / 4;		//	объемный расход, м^3/c
		double Q_h = Q_s * 3600;											//  объемный расход, м^3/ч	

		std::cout << "Расход по Ньютону = " << Q_h << std::endl;	

		
		double error = 0.01;

		EXPECT_NEAR(1.99, result.argument, error);
};


TEST(zadacha_5, PP_ABOVE_EULER_ON_NEWTON) {

	setlocale(LC_ALL, "Russian");

	zadacha_1 iniz1;
	zadacha_2 iniz2(iniz1);

	iniz1.h = iniz1.L / iniz1.n;
	
	massiv pressure;
	
	pressure.parameter = vector<double>(iniz1.n); 

	class PP_ABOVE_EULER_ON_NEWTON : public fixed_system_t<1>
	{
		/// @brief Ссылка на структуру с параметрами трубы 
		const zadacha_1 iniz1;
		const zadacha_2 iniz2;
		massiv& massiv_dannye;
		
	using fixed_system_t<1>::var_type; public:

		
		PP_ABOVE_EULER_ON_NEWTON(const zadacha_1& iniz1, const zadacha_2& iniz2, massiv& massiv_dannye)
			: iniz1(iniz1), iniz2(iniz2), massiv_dannye(massiv_dannye) {}

		/// @brief Функция невязок - все члены уравнения Бернулли в правой части 
		/// @param v - скорость течения нефти, [м/с] 
		/// @return Значение функции невязок при заданной скорости 
		var_type residuals(const var_type& V)
		{ //V - искомая скорость
			double Re = V * iniz1.d_m / (iniz1.nu * pow(10, -6));
			double resistance = hydraulic_resistance_altshul(Re, iniz2.sher);
			double t_w = resistance / 8 * iniz1.density * pow(V, 2);
			double p_n;
			double p_k = iniz2.pp_k;
			// расчет p_n по методу Эйлера 

			for (int i = 0; i < iniz1.n; ++i) {
				p_n = p_k - iniz1.h * (-4 / iniz1.d_m * t_w - iniz1.density * M_G * (iniz1.z_l - iniz1.z_0) / ((iniz1.n - 1) * iniz1.h));
				p_k = p_n;
				massiv_dannye.parameter[i] = p_n;

			};

			double delta_p_0;

			return
			{

			   delta_p_0 = iniz2.pp_n - p_n

			};

		};
	};
	
	// Создание экземпляра класса, который и будет решаемой системой
	PP_ABOVE_EULER_ON_NEWTON test(iniz1, iniz2, pressure);
	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
	// { 1} - Начальное приближение скорости 1 
	fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);
	std::cout << "Классическая задача PP поверх Эйлера на методе Ньютона " << " V = " << result.argument << "\n";

	double abs_error = 1;
	EXPECT_NEAR(1.99, result.argument, abs_error);

}











