#include <iostream> 
#include <fstream> 
#include <string.h> 
#include <stdio.h> 
#include <vector> 
#include <cmath> // для использования функции pow
using namespace std;

struct zadacha_1 {
	double L; //длина трубопровода
	double D; //внешний диаметр трубопровода
	double b; // толщина стенки трубопровода
	double z_0; // начальная высотная отметка профиля трубопровода 
	double z_l; // конечная высотная отметка профиля трубопровода 
	double density; // плоность нефти
	double nu; //кинематическая вязкость нефти
	double Q; // объемный расход 
	double p_k; // давление в конце трубопровода
	double g;
	double d;
	double d_m;
	double q;

	zadacha_1() {

		L = 80 * pow(10, 3); // км
		D = 720; // мм
		b = 10; // мм
		z_0 = 50; //м 
		z_l = 100; // м
		density = 870; // кг/м3
		nu = 15; // сСт
		Q = 3500; // м3/ч
		p_k = 0.6 * pow(10, 6); // Па
		g = 9.81;

		d = D - 2 * b; // мм
		d_m = d / 1000; // м
		q = Q / 3600; // м3/с
	}
};


struct zadacha_2 {

	double Z;
	double lamd_2;
	double epsilon;
	double lamd;
	double sher;

	zadacha_2(zadacha_1& iniz1) {

		Z = 0.15;
		lamd_2 = 0.02;
		epsilon = 1e-5;

		sher = Z / iniz1.d;
	}
};




int main() {

	setlocale(LC_ALL, "Russian"); // Корректный вывод руского текста

	zadacha_1 iniz1; //объявление переменной iniz1 (хранение данных о трубопроводе) типа zadacha_1
	zadacha_2 iniz2(iniz1);


	double d_1 = (iniz1.D - 2 * iniz1.b) * pow(10, -3);
	double V = 4 * iniz1.q / (3.14 * pow(d_1, 2));
	double Re = V * d_1 / (iniz1.nu * pow(10, -6));
	double e = iniz1.nu * pow(10, -2) / iniz1.d;
	double Lyam;
	if (Re < 2300)
		Lyam = 64 / Re;
	else if (Re > 500 / e)
		Lyam = 0.11 * pow(e, 0.25);
	else if (Re > 2300)
		Lyam = 0.11 * pow((e + 68 / Re), 0.25);

	double p_n = (iniz1.p_k / (iniz1.density * 9.81) + iniz1.z_0 - iniz1.z_l + Lyam * (iniz1.L / d_1 * pow(V, 2) / 2 / 9.91)) * (iniz1.density * 9.81);

	cout << "Результат: " << p_n << "Па";


	int index = 0;



	do {
		iniz2.lamd = iniz2.lamd_2;
		V = pow(((((p_n - iniz1.p_k) / (iniz1.density * iniz1.g) + iniz1.z_0 - iniz1.z_l) * 2 * iniz1.g * iniz1.d_m) / (iniz2.lamd * iniz1.L)), 0.5);
		Re = V * iniz1.d_m / iniz1.nu / pow(10, -6);
		iniz2.lamd_2 = 0.11 * pow(((iniz2.sher + 68 / Re)), 0.25);

		index++;


	} while (abs(iniz2.lamd_2 - iniz2.lamd) > iniz2.epsilon);


	std::cout << "Скорость: " << V << std::endl;
	std::cout << "Количество итераций: " << index << std::endl;

	double Q = 3.14 * pow(iniz1.d_m, 2) / V;

	std::cout << "Расход: " << Q << std::endl; //123
}