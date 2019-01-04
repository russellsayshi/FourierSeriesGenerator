#include <lexer.h>
#include <interpreter.h>
#include <iostream>
#include <unordered_map>
#include <iomanip>
#include <limits>
#include <thread>

/* Make sure to compile me with the MathParseC- repo given here: https://github.com/russellsayshi/MathParseC- */

void numeric_integral(const double* const riemann_values, const double* const sin_values, const double* const cos_values, double* const sin_out, double* const cos_out, const long num_riemann_terms, const long n) {
	double sin_sum = 0;
	double cos_sum = 0;

	for(long i = 0; i < num_riemann_terms; i++) {
		long pos = i;
		if(i >= num_riemann_terms/2) {
			pos -= num_riemann_terms/2;
		} else {
			pos = num_riemann_terms/2 - pos;
		}
		//std::cout << std::endl << pos;
		long nth_trig_term = (pos * n) % (num_riemann_terms);
		//std::cout << " " << nth_trig_term;
		double sin_val, cos_val;
		if(nth_trig_term < num_riemann_terms/2) {
			sin_val = sin_values[nth_trig_term];
			cos_val = cos_values[nth_trig_term];
		} else {
			//std::cout << " fetching sin: -(" << (nth_trig_term - num_riemann_terms/2);
			sin_val = -sin_values[nth_trig_term - num_riemann_terms/2];
			//std::cout << ") fetching cos: (" << (num_riemann_terms - nth_trig_term) << ")";
			cos_val = cos_values[num_riemann_terms - nth_trig_term];
		}
		if(i < num_riemann_terms/2) {
			sin_val *= -1;
		}
		//std::cout << " sin: " << sin_val << ", cos: " << cos_val;
		sin_sum += sin_val * riemann_values[i];
		cos_sum += cos_val * riemann_values[i];
	}
	*sin_out = sin_sum;
	*cos_out = cos_sum;
}

void integrate_multithreaded(const double* const riemann_values, const double* const sin_values, const double* const cos_values, const long num_riemann_terms, double* const fourier_coeffs_sin, double* const fourier_coeffs_cos, int start_inc, int end_exc) {
	for(int i = start_inc; i < end_exc; i++) {
		numeric_integral(riemann_values, sin_values, cos_values, &(fourier_coeffs_sin[i]), &(fourier_coeffs_cos[i]), num_riemann_terms, i);
	}
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

	std::unordered_map<std::string, double> map;
	map["e"] = 2.718281828459045235360287471352662497757247093699959574;
	double pi = 3.141592653589793238462643383279502884197169399375105821;
	map["pi"] = pi;
	map["b"] = bound;

	double* riemann_values = (double*)malloc(sizeof(double) * num_riemann_terms);
	if(riemann_values == nullptr) {
		std::cerr << "Cannot allocate enough space to perform computations (#1)." << std::endl;
		return 1;
	}
	double* cos_values = (double*)malloc(sizeof(double) * (num_riemann_terms/2+1));
	if(cos_values == nullptr) {
		std::cerr << "Cannot allocate enough space to perform computations (#2)." << std::endl;
		return 2;
	}
	double* sin_values = (double*)malloc(sizeof(double) * (num_riemann_terms/2+1));
	if(sin_values == nullptr) {
		std::cerr << "Cannot allocate enough space to perform computations (#3)." << std::endl;
		return 3;
	}
	double x_pos = -bound;
	double delta = 2*bound/num_riemann_terms;
	double offset = 0;
	bool skip_failed = false;
	for(int i = 0; i < num_riemann_terms; i++) {
		map["x"] = x_pos;
		double tval = inter.interpret(map) * delta;
		if(isinf(tval) || isnan(tval)) {
			if(skip_failed) {
				tval = 0;
			} else {
				std::cerr << "Failed to evaluate function 'f(x) = " << input << "' at x = " << x_pos << ". Value: " << tval << "." << std::endl;
				std::cerr << "Do you want to ignore this (and future) calculation fails? [y/N] ";
				std::cin.clear();
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::string ignore_eval;
				getline(std::cin, ignore_eval);
				if(ignore_eval.size() > 0 && (ignore_eval[0] == 'y' || ignore_eval[0] == 'Y')) {
					std::cout << "Continuing." << std::endl;
					tval = 0;
					skip_failed = true;
				} else {
					std::cerr << "Aborting." << std::endl;
					return 4;
				}
			}
		}
		offset += tval;
		riemann_values[i] = tval;
		if(i >= num_riemann_terms/2) {
			double npxol = pi * x_pos / bound;
			cos_values[i-num_riemann_terms/2] = cos(npxol);
			sin_values[i-num_riemann_terms/2] = sin(npxol);
		}
		x_pos += delta;
	}
	sin_values[num_riemann_terms/2] = sin_values[0];
	cos_values[num_riemann_terms/2] = -cos_values[0];
	offset /= (2 * bound);
	std::cout << "Evaluated at all x." << std::endl;

	//Calculate offset of series
	std::cout << std::fixed;
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << offset;

	double* fourier_coeffs_cos = (double*)malloc(sizeof(double) * ((num_fourier_terms+1)/2));
	if(fourier_coeffs_cos == nullptr) {
		std::cerr << "Unable to allocate #4." << std::endl;
		return 6;
	}
	double* fourier_coeffs_sin = (double*)malloc(sizeof(double) * ((num_fourier_terms+1)/2));
	if(fourier_coeffs_sin == nullptr) {
		std::cerr << "Unable to allocate #5." << std::endl;
		return 7;
	}

	unsigned int max_threads = std::thread::hardware_concurrency();
	std::thread threads[max_threads];

	int num_per_thread = ((num_fourier_terms+1)/2)/max_threads;
	for(int i = 0; i < max_threads; i++) {
		int num_this_thread = (i == max_threads - 1) ? ((num_fourier_terms+1)/2) - (max_threads - 1) * num_per_thread : num_per_thread;
		int tt_lower_bound = num_per_thread * i;
		int tt_upper_bound = tt_lower_bound + num_this_thread;
		threads[i] = std::thread(integrate_multithreaded, (const double* const)riemann_values, (const double* const)sin_values, (const double* const)cos_values, num_riemann_terms, fourier_coeffs_sin, fourier_coeffs_cos, tt_lower_bound, tt_upper_bound);
	}
	for(int i = 0; i < max_threads; i++) {
		threads[i].join();
	}

	for(int i = 1; i <= (num_fourier_terms+1)/2; i++) {
		double x_coeff = i * pi / bound;
		std::cout << " + (" << (fourier_coeffs_cos[i-1]/bound) << ")*cos(x*" << x_coeff << ")";
		if(!(i == (num_fourier_terms+1)/2 && (num_fourier_terms & 1) == 1)) {
			std::cout << " + (" << (fourier_coeffs_sin[i-1]/bound) << ")*sin(x*" << x_coeff << ")";
		}
	}
	std::cout << std::endl;
	free(riemann_values);
	free(sin_values);
	free(cos_values);
	free(fourier_coeffs_cos);
	free(fourier_coeffs_sin);
	return 0;
}
