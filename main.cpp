#include "includes/dtoa.h"

#include <iostream>
#include <memory>
#include <array>
#include <cstring>
#include <bitset>
#include <cmath>
#include <float.h>
#include <random>
#include <charconv>
#include <chrono>

// https://baseconvert.com/ieee-754-floating-point

std::string d_to_s_v(double a_value) {
	std::string lr = d_to_s(a_value);
	double read = s_to_d(lr);
	std::string rl = d_to_s(read);
	std::string er = d_to_s(a_value - read);

	std::cout.precision(18);
	std::cout << std::scientific
		<< "dtos [" << lr << "]" << std::endl
		<< "cout [" << a_value << "]" << std::endl
		<< "stod [" << rl << "]" << std::endl
		<< "errr [" << er << "]" << std::endl << std::endl;

	return lr;
}

[[noreturn]] void print_mul_map()
{
	uint64_t mul = 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5;
	int mul_idx = 13;
	uint64_t v_low = 1ull << 31;
	uint64_t v_hig = (1ull << 32) - 1ull;
	for(int i = 1; i < 32 && mul > 1; i++) {
		// get highest bit
		int v_low_idx = 0;
		for(v_low_idx = 0; (((v_low * mul) << v_low_idx) >> 63) != 1; v_low_idx++);
		int v_hig_idx = 0;
		for(v_hig_idx = 0; (((v_hig * mul) << v_hig_idx) >> 63) != 1; v_hig_idx++);

		if((v_low * mul) / mul != v_low) break;
		if((v_hig * mul) / mul != v_hig) break;

		// std::cout << "------------------" << std::endl;
		// std::cout << std::bitset<32>(v_low >> 32) << ":" << std::bitset<32>(v_low) << std::endl;
		// std::cout << std::bitset<32>(v_hig >> 32) << ":" << std::bitset<32>(v_hig) << std::endl;
		// std::cout << std::bitset<32>((v_low * mul) >> 32) << ":" << std::bitset<32>(v_low * mul) << std::endl;
		// std::cout << std::bitset<32>((v_hig * mul) >> 32) << ":" << std::bitset<32>(v_hig * mul) << std::endl;
		// std::cout << mul << ", l=" << v_low_idx << ", h=" << v_hig_idx << std::endl;

		std::printf("% 10lu, h=%2d, l=%2d, i=%2d\n", mul, 32 - v_hig_idx + mul_idx, 32 - v_low_idx + mul_idx, mul_idx);
		// std::cout << mul << ", h=" << v_hig_idx << ", l=" << v_low_idx << std::endl;

		mul = mul / 5;
		mul_idx -= 1;
	}

	std::exit(0);
}

[[noreturn]] void print_div_map()
{
	uint64_t div = 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5;
	int div_idx = 13;
	uint64_t v_low = 1ull << 63;
	uint64_t v_hig = ~0ull;
	for(int i = 1; i < 32 && div > 1; i++) {
		if((v_low / div) == 0) break;
		if((v_hig / div) == 0) break;

		// get highest bit
		int v_low_idx = 0;
		for(v_low_idx = 0; (((v_low / div) << v_low_idx) >> 63) != 1; v_low_idx++);
		int v_hig_idx = 0;
		for(v_hig_idx = 0; (((v_hig / div) << v_hig_idx) >> 63) != 1; v_hig_idx++);

		// std::cout << "------------------" << std::endl;
		// std::cout << std::bitset<32>(v_low >> 32) << ":" << std::bitset<32>(v_low) << std::endl;
		// std::cout << std::bitset<32>(v_hig >> 32) << ":" << std::bitset<32>(v_hig) << std::endl;
		// std::cout << std::bitset<32>((v_low / div) >> 32) << ":" << std::bitset<32>(v_low / div) << std::endl;
		// std::cout << std::bitset<32>((v_hig / div) >> 32) << ":" << std::bitset<32>(v_hig / div) << std::endl;
		// std::cout << div << ", l=" << v_low_idx << ", h=" << v_hig_idx << std::endl;

		std::printf("% 10lu, h=%2d, l=%2d, i=%2d\n", div, v_hig_idx + div_idx, v_low_idx + div_idx, div_idx);
		// std::cout << div << ", h=" << v_hig_idx << ", l=" << v_low_idx << std::endl;

		div = div / 5;
		div_idx -= 1;
	}

	std::exit(0);
}

int main(int argc, char** argv)
{
	// print_mul_map();
	// print_div_map();
/*
	d_to_s_v(DBL_MAX);
	d_to_s_v(DBL_MIN);
	d_to_s_v(9007199254740990.0 * (1ull << 12));
	d_to_s_v(9007199254740990.0 * (1ull << 13));
	d_to_s_v(9007199254740990.0 * (1ull << 14));
	d_to_s_v(9007199254740990.0 * (1ull << 15));
	d_to_s_v(9007199254740990.0 * (1ull << 16));
*/
	s_to_d(d_to_s_v(9007199254740990.0 * (1ull << 11)));

	d_to_s_v(9007199254740990.0 * (1ull << 11));
	d_to_s_v(9007199254740990.0 * (1ull << 12));
	d_to_s_v(9007199254740990.0 * (1ull << 13));
	d_to_s_v(9007199254740990.0 * (1ull << 14));
	d_to_s_v(9007199254740990.0 * (1ull << 15));

	d_to_s_v(9007199254740990.0 / (1ull << 11));
	d_to_s_v(9007199254740990.0 / (1ull << 10));
	d_to_s_v(9007199254740990.0 / (1ull << 9));
	d_to_s_v(9007199254740990.0 / (1ull << 63));
	d_to_s_v(DBL_MIN);
	d_to_s_v(4503599627370495.0 / 8);
	d_to_s_v(4503599627370495.0 / 4);
	d_to_s_v(4503599627370495.0 / 2);
	d_to_s_v(4503599627370495.0);
	d_to_s_v(9007199254740990.0);
	d_to_s_v(9007199254740990.5);
	d_to_s_v(9007199254740991.0);
	d_to_s_v(100);
	d_to_s_v(DBL_MIN);
	d_to_s_v(DBL_MAX);
	d_to_s_v(12431234237423482734982734982734982734287381.128342783482);
	d_to_s_v(100.314159265359);
	d_to_s_v(0.314159265359);
	d_to_s_v(1243123423781.128342783482);
	d_to_s_v(100);
	d_to_s_v(1234);
	d_to_s_v(1.5);
	d_to_s_v(1.99999);
	d_to_s_v(1.00001);
	d_to_s_v(0.99999);
	d_to_s_v(1.0);
	d_to_s_v(0.0);
	d_to_s_v(0.5);
	d_to_s_v(0.1);
	d_to_s_v(0.01);
	d_to_s_v(0.001);
	d_to_s_v(0.0001);
	d_to_s_v(0.00001);
	d_to_s_v(0.000001);
	d_to_s_v(18'446'744'073'709'551'615.0);
	d_to_s_v(20'446'744'073'709'551'615.0);
	d_to_s_v(4.5865e+199);
	/**/

	//std::uniform_real_distribution<double> unif(1e-100,1e-26);
	std::uniform_int_distribution<int> uexp(-307, 308);
	std::uniform_real_distribution<double> upar(0.5,1.0);
	std::default_random_engine re;


	std::cout.precision(18);
	std::cout << std::scientific;
	constexpr size_t buff = 2'000'000ull;
	constexpr size_t iter = 1'000'000'000'000ull;
	constexpr bool check_std = false;
	std::vector<double> data(buff);
	for(size_t i = 0; i < buff; i++) {
		double v{};
		v = upar(re) * std::pow(10.0, uexp(re));
		data[i] = v;
	}

	auto start = std::chrono::high_resolution_clock::now();
	if(check_std) {
		char buffer[4096]{};
		// 2460900 ops / sec
		for(size_t i = 0; i < iter; i++) {
			double v = data[i % buff];
			double c{};
			auto res = std::to_chars(buffer, buffer + 4095, v);
			std::from_chars(buffer, res.ptr, c);
			*res.ptr = '\0';

			if((i % 100000) == 0) {
				auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(
					std::chrono::high_resolution_clock::now() - start
				);
				auto ratio = static_cast<double>(i) / (static_cast<double>(time.count()) / 1000000000.0);
				std::cout << i << " : " << (int) ratio << " ops / sec" << std::endl;
			}
			if((c - v) != 0) {
				std::string rv = d_to_s(c);
				std::string dv = d_to_s(c - v);
				std::cout << v << " : " << buffer << " : " << rv << ": " << dv << std::endl;
				break;
			}
		}
	} else {
		std::array<char, 32> buffer{};
		// 1335000 ops / sec
		// 2930000 ops / sec -O3  20% faster than charconv
		// 198'564'200'000 : 2923844 ops / sec
		// 12167976
		//  3641439

		// 2x slower than charconv
		for(size_t i = 0; i < iter; i++) {
			double v = data[i % buff];

			int len = d_to_s_a(buffer, v);
			double c = s_to_d_a(buffer.data(), buffer.data() + len);
			buffer[len] = '\0';

			if((i % 100000) == 0) {
				auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(
					std::chrono::high_resolution_clock::now() - start
				);
				auto ratio = static_cast<double>(i) / (static_cast<double>(time.count()) / 1000000000.0);
				std::cout << i << " : " << (int) ratio << " ops / sec" << std::endl;
			}
			if((c - v) != 0) {
				std::string rv = d_to_s(c);
				std::string dv = d_to_s(c - v);
				std::cout << i << ":" << v << " : " << buffer.data() << " : " << rv << ": " << dv << std::endl;
				break;
			}
		}
	}

	// d_to_s(0.0/0.0);
	// d_to_s(NAN);
	// d_to_s(INFINITY);
	// d_to_s(-INFINITY);
	// d_to_s(HUGE_VAL);

	return 0;
}

// g++ -g main.cpp && gdb -ex 'set pagination off' -ex run ./a.out
