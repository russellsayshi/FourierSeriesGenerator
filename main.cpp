#include <lexer.h>
#include <interpreter.h>
#include <iostream>
#include <unordered_map>
#include <iomanip>
#include <limits>

/* Make sure to compile me with the MathParseC- repo given here: https://github.com/russellsayshi/MathParseC- */

double numeric_integral(interpreter& inter, std::unordered_map<std::string, double>& map, std::string variable, int num_sums, double left_bound, double right_bound) {
	double bound_difference = right_bound - left_bound;
	double delta = bound_difference / num_sums;
	double position = left_bound;
	double running_left_sum = 0;
	double running_right_sum = 0;
	map[variable] = position;
	running_left_sum += inter.interpret(map) * delta;
	position += delta;
	for(int i = 1; i < num_sums - 1; i++) {
		map[variable] = position;
		double rectArea = inter.interpret(map) * delta;
		running_left_sum += rectArea;
		running_right_sum += rectArea;
		position += delta;
	}
	map[variable] = position;
	running_right_sum += inter.interpret(map) * delta;
	return (running_left_sum + running_right_sum)/2;
}

int main(int argc, char** argv) {
	//Fetch user input
	std::cout << "Enter the function to create the fourier series (e.g. 2*sin(x)): ";
	std::string input;
	getline(std::cin, input);
	std::cout << "What bound do you want on the function (Fourier series will be generated from -b to b): ";
	double bound;
	std::cin >> bound;
	std::cout << "How many terms do you want in your fourier series? ";
	int num_fourier_terms;
	std::cin >> num_fourier_terms;
	std::cout << "How many terms do you want in your Riemann sum? ";
	int num_riemann_terms;
	std::cin >> num_riemann_terms;

	//Get geared up to calculate
	lexer lex(input);
	interpreter inter;
	inter.fetch_tokens(lex);

	std::string cos_term_str = "(" + input + ")*cos(n*pi*x/b)";
	std::string sin_term_str = "(" + input + ")*sin(n*pi*x/b)";
	std::cout << "Cosine func: " << cos_term_str << std::endl;
	std::cout << "Sine func: " << sin_term_str << std::endl;
	lexer cos_term_lex(cos_term_str);
	lexer sin_term_lex(sin_term_str);
	interpreter inter_cos;
	interpreter inter_sin;
	inter_cos.fetch_tokens(cos_term_lex);
	inter_sin.fetch_tokens(sin_term_lex);

	std::unordered_map<std::string, double> map;
	map["e"] = 2.718281828459045235360287471352662497757247093699959574;
	double pi = 3.141592653589793238462643383279502884197169399375105821;
	map["pi"] = pi;
	map["b"] = bound;

	//Calculate offset of series
	double offset = numeric_integral(inter, map, "x", num_riemann_terms, -bound, bound)/(2 * bound);
	std::cout << std::fixed;
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << offset;

	for(int i = 0; i < num_fourier_terms; i++) {
		map["n"] = i/2+1;
		double cos_coeff = numeric_integral(inter_cos, map, "x", num_riemann_terms, -bound, bound)/bound;
		double np_over_l = (i/2+1)*pi/bound;
		std::cout << " + (" << cos_coeff << ")*cos(x*" << np_over_l << ")";
		std::cout.flush();
		if(++i < num_fourier_terms) {
			double sin_coeff = numeric_integral(inter_sin, map, "x", num_riemann_terms, -bound, bound)/bound;
			std::cout << " + (" << sin_coeff << ")*sin(x*" << np_over_l << ")";
			std::cout.flush();
		}
	}
	std::cout << std::endl;
}
