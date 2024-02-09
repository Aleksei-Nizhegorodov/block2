#include <iomanip>					   // Необходимо для устранения ошибок 'setprecision': is not a member of 'std' ; 'setprecision': identifier not found (из-за бага в pde_solvers ) 
#include <iostream> 
#include <fstream> 
#include <string.h> 
#include <stdio.h> 
#include <vector> 
#include <cmath>                       // для использования функции pow
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>  // подключение библиотеки pde_solvers (В.В. Южанин)

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

	zadacha_1() {

		L = 80 * pow(10, 3);	// км
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
	}
};


/// @brief Данные по задаче 2
struct zadacha_2 {

	double Z; 
	double lamd_2;
	double epsilon;		// допустимая погрешность 
	double lamd; 
	double sher;		// относительная шероховатость 
	double pp_n;
	double pp_k; 

	zadacha_2(zadacha_1& iniz1) {

		Z = 0.15; // мм
		lamd_2 = 0.02; 
		epsilon = 1 * pow(10, (- 3.21)); //1e-5;
		pp_n = 5 * pow(10, 6);
		pp_k = 0.8 * pow(10, 6); 

		sher = Z / iniz1.d;
	}
};


/// @brief Задача 1, PP
/// @param iniz1 ссылка на данные по задаче 1
/// @return 
double calculatePressure(zadacha_1& iniz1, zadacha_2& iniz2) {

	double V = 4 * iniz1.q / (M_PI * pow(iniz1.d_m, 2));    // скорость в трубопровде
	double Re = V * iniz1.d_m / (iniz1.nu * pow(10, -6));   // число Рейнольдса  
	double e = iniz1.nu * pow(10, -2) / iniz1.d;            // шероховатость 
	

	double resistance = hydraulic_resistance_isaev(Re, e);  // ф-ция расчета значения лямбды

	double p_n = (iniz1.p_k / (iniz1.density * M_G) + iniz1.z_0 - iniz1.z_l + resistance
		* (iniz1.L /iniz1.d_m * pow(V, 2) / 2 / M_G)) * (iniz1.density * M_G);	// значение распределения начльного давления, МПа

	double pogr = (6 * pow(10, 6) - p_n) /(6 * pow(10,6)) ; // Вычисление расхождения в расчетах программы и М.В.Лурье
	
	cout << "Результат: " << p_n << "Па" << endl;

	return p_n;
}


/// @brief Задача 2, QP
/// @param iniz1, ссылка на данные по задаче 1
/// @param iniz2, ссылка на данные по задаче 2
/// @param V, скорость [м/с]
/// @param Re, число Рейнольдса
/// @param index, указатель счетчика 
/// @param p_n, давление в начале трубопровода [МПа]
void calculateFlowAndIterations(zadacha_1& iniz1, zadacha_2& iniz2, 
	double& V, double& Re, int& index, double& p_n) {
	
	//double v = iniz1.d_m * 2 * M_G / iniz1.L * ((p_n - iniz1.p_k) / (iniz1.density * M_G) + iniz1.z_0 - iniz1.z_l); для проверки по Лурье 

	// Нахождение lamd_2 методом последовательных приближений 
	do {
		iniz2.lamd = iniz2.lamd_2; //предположили, что начльное значение lamd_2 = 0.02

		V = pow(((iniz1.d_m * 2 * M_G / iniz1.L * ((iniz2.pp_n - iniz2.pp_k) / (iniz1.density * M_G) 
			+ iniz1.z_0 - iniz1.z_l)) / iniz2.lamd), 0.5);
		Re = V * iniz1.d_m / iniz1.nu / pow(10, -6);
		iniz2.lamd_2 = 0.11 * pow(((iniz2.sher + 68 / Re)), 0.25);	//уточняем наше предположение (можно ли здесь использовать формулу исаева??????)

		index++;

	} while (abs(iniz2.lamd_2 - iniz2.lamd) > iniz2.epsilon);

	std::cout << "Скорость: " << V << std::endl;
	std::cout << "Количество итераций: " << index << std::endl;

	double Q_s = (M_PI * pow(iniz1.d_m, 2)) * V / 4;		//	объемный расход, м^3/c
	double Q_h = Q_s * 3600;								// объемный расход, м^3/ч		

	std::cout << "Расход: " << Q_h << std::endl; 

}
/// @brief 
/// @return 


int main() {

	setlocale(LC_ALL, "Russian"); // Корректный вывод руского текста

	zadacha_1 iniz1;			  //объявление переменной iniz1 (хранение данных о трубопроводе) типа zadacha_1
	zadacha_2 iniz2(iniz1);		  //объявление переменной iniz2(iniz1) (хранение данных о трубопроводе, ) типа zadacha_2

	double V, Re, p_n;

	int index = 0;

	p_n = calculatePressure(iniz1, iniz2);	// Находим знчение начального давления в задаче 1

	calculateFlowAndIterations(iniz1, iniz2, V, Re, index, p_n); // решаем задачу 2, обращаясь к значению начального давления в задаче 1

	return 0;
}
