#ifndef DTOA_H
#define DTOA_H

#include <iostream>
#include <memory>
#include <array>
#include <cstring>
#include <bitset>
#include <cmath>
#include <float.h>

// https://baseconvert.com/ieee-754-floating-point

struct d_norm {
	uint64_t value;
	int exp;
};

constexpr std::array<uint32_t, 14 * 2> NEG_LUT{
	1, 0,
	5, 3,
	25, 5,
	125, 7,
	625, 10,
	3125, 12,
	15625, 14,
	78125, 17,
	390625, 19,
	1953125, 21,
	9765625, 24,
	48828125, 26,
	244140625, 28,
	1220703125, 31};
inline d_norm d_to_s_norm_64_opt(uint64_t raw, int pow_2) {
/*
	1220703125, h=44, l=43, i=13
	 244140625, h=40, l=39, i=12
	  48828125, h=37, l=36, i=11
		9765625, h=34, l=33, i=10
		1953125, h=30, l=29, i= 9
		 390625, h=27, l=26, i= 8
		  78125, h=24, l=23, i= 7
		  15625, h=20, l=19, i= 6
			3125, h=17, l=16, i= 5
			 625, h=14, l=13, i= 4
			 125, h=10, l= 9, i= 3
			  25, h= 7, l= 6, i= 2
				5, h= 4, l= 3, i= 1
*/
	int pow_10 = 0;
	uint64_t c = 0;

	// i = (-pow_2 - 3) / 10 + (-pow_2 - 3) % 10 >= 6 + (-pow_2 - 3) % 10 >= 10
	// i = ((-pow_2 - 4 + 3) - ((-pow_2 - 3) / 10)) / 3
	// i = (2 - pow_2 * 9) / 30

	if (pow_2 < 0) {
		while (pow_2 <= -44) { // 1220703125, low=2, hig=1
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;
			c = (c * 1220703125ull);
			b = (b * 1220703125ull) + (c >> 32);
			a = (a * 1220703125ull) + (b >> 32);
			b = b & 0xffffffff;
			c = c & 0xffffffff;
			pow_10 -= 13;
			if (a >> 62) { // hig
				raw = (a << 1) | (b >> 31);
				c = ((b << 1) | (c >> 31)) & 0xffffffff;
				pow_2 += 31 + 13;
			} else { // low
				raw = (a << 2) | (b >> 30);
				c = ((b << 2) | (c >> 30)) & 0xffffffff;
				pow_2 += 30 + 13;
			}
		}
		// i = (2 - pow_2 * 9) / 30
		// (1 1 1) (2 2 2) (3 3 3 3) (4 4 4) (5 5 5) (6 6 6 6) ...
		// if (pow_2 <= -4) {
			const int idx = (2 - pow_2 * 9) / 30;
			// const int idx = (2 - (pow_2 << 3) - pow_2) / 30;
			const uint64_t mul = NEG_LUT[(idx << 1) | 0];
			const uint64_t shf = NEG_LUT[(idx << 1) | 1];
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;
			c = (c * mul);
			b = (b * mul) + (c >> 32);
			a = (a * mul) + (b >> 32);
			b = b & 0xffffffff;
			pow_10 -= idx;
			if (a >> (30 + shf)) {
				raw = (a << (32 - shf)) | (b >> shf);
				pow_2 += shf + idx;
			} else {
				raw = (a << (31 - shf)) | (b >> (shf - 1));
				pow_2 += shf + idx - 1;
			}
		// }
		raw >>= -pow_2;
		// pow_2 = 0;
	} else {
		while (pow_2 >= 44) {
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;

			uint64_t av = (a / 1220703125ull);
			uint64_t ar = (a % 1220703125ull);
			uint64_t bv = ((b + (ar << 32)) / 1220703125ull);
			uint64_t br = ((b + (ar << 32)) % 1220703125ull);
			uint64_t cv = ((c + (br << 32)) / 1220703125ull);

			if (av >> 1) {
				// 30 bit shift
				av <<= 30;
				bv <<= 30;
				cv <<= 30;
				pow_2 -= 30 + 13;
			} else {
				// 31 bit shift
				av <<= 31;
				bv <<= 31;
				cv <<= 31;
				pow_2 -= 31 + 13;
			}

			// enable c to make b reuse old bits
			c = (cv & 0xffffffff);
			raw = (av << 32) + (bv << 0) + (cv >> 32);
			pow_10 += 13;
			// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
		}
/*
		1220703125, h=30, l=31, i=13
		 244140625, h=27, l=28, i=12
		  48828125, h=25, l=26, i=11
			9765625, h=23, l=24, i=10
			1953125, h=20, l=21, i= 9
			 390625, h=18, l=19, i= 8
			  78125, h=16, l=17, i= 7
			  15625, h=13, l=14, i= 6
				3125, h=11, l=12, i= 5
				 625, h= 9, l=10, i= 4
				 125, h= 6, l= 7, i= 3
				  25, h= 4, l= 5, i= 2
					5, h= 2, l= 3, i= 1
*/
		// if (pow_2 > 0) {
			const int idx = ((2 + pow_2 * 9) / 30) + 1;
			// const int idx = (2 - (pow_2 << 3) - pow_2) / 30;
			const uint64_t div = NEG_LUT[(idx << 1) | 0];
			const uint64_t shf = NEG_LUT[(idx << 1) | 1];

			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;

			uint64_t av = (a / div);
			uint64_t ar = (a % div);
			uint64_t bv = ((b + (ar << 32)) / div);
			uint64_t br = ((b + (ar << 32)) % div);
			uint64_t cv = ((c + (br << 32)) / div);

			if (av >> (32 - shf)) {
				// 2 bit shift
				av <<= shf - 1;
				bv <<= shf - 1;
				cv <<= shf - 1;
				pow_2 -= shf + idx - 1;
			} else {
				// 3 bit shift
				av <<= shf;
				bv <<= shf;
				cv <<= shf;
				pow_2 -= shf + idx;
			}

			// enable c to make b reuse old bits
			c = (cv & 0xffffffff);
			raw = (av << 32) + (bv << 0) + (cv >> 32);
			pow_10 += idx;
			// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
		// }
		// if (pow_2 < 0) {
			raw >>= -pow_2;
			pow_2 = 0;
		// }
	}

	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;

	// Make sure raw is 19 digits
	if (raw >= 10'000'000'000'000'000'000ull) {
		raw /= 10;
		pow_10 += 1;
	}

	return { raw, pow_10 };
}

inline d_norm d_to_s_norm_64(uint64_t raw, int pow_2) {
	constexpr bool enable_c = true;
	uint64_t c = 0;
	int pow_10 = 0;

	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << std::endl;

	while (pow_2 <= -44) { // < 0 works
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;
		// enable c to make b reuse old bits
		c = (c * 1220703125ull) * (enable_c);
		b = (b * 1220703125ull) + (c >> 32);
		a = (a * 1220703125ull) + (b >> 32);
		b = b & 0xffffffff;
		c = c & 0xffffffff;
		pow_10 -= 13;

		// std::cout
		//	 << std::bitset<32>(a >> 32) << ":"
		//	 << std::bitset<32>(a) << ":"
		//	 << std::bitset<32>(b) << std::endl;


		// (a >> 34) // checks if 2 bit needs to be shifted
		// (a >> 62) // checks if 31 bits need to be shifted
		// (a >> 34) -> (32 - 2 - 1) : (0 - 2 - 1)
		// (a >> 62) -> (32 - 30 - 1) : (0 - 30 - 1)
		// (a >> 62) -> (32 - 29 - 1) : (0 - 29 - 1)

		if (a >> 62) {
			raw = (a << 1) | (b >> 31);
			c = ((b << 1) | (c >> 31)) & 0xffffffff;
			pow_2 += 31 + 13;
		} else {
			raw = (a << 2) | (b >> 30);
			c = ((b << 2) | (c >> 30)) & 0xffffffff;
			pow_2 += 30 + 13;
		}

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	while (pow_2 <= -4) { // < 0 works
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;
		// enable c to make b reuse old bits
		c = (c * 5ull) * (enable_c);
		b = (b * 5ull) + (c >> 32);
		a = (a * 5ull) + (b >> 32);
		b = b & 0xffffffff;
		pow_10 -= 1;

		// std::cout
		//	<< std::bitset<3>(a >> 32) << ":"
		//	<< std::bitset<32>(a) << ":"
		//	<< std::bitset<32>(b) << std::endl;

		if (a >> 34) {
			raw = (a << 29) | (b >> 3);
			c = ((b << 29) | (c >> 3)) & 0xffffffff;
			pow_2 += 3 + 1;
		} else {
			raw = (a << 30) | (b >> 2);
			c = ((b << 30) | (c >> 2)) & 0xffffffff;
			pow_2 += 2 + 1;
		}

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	while (pow_2 >= 44) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;

		uint64_t av = (a / 1220703125ull);
		uint64_t ar = (a % 1220703125ull);
		uint64_t bv = ((b + (ar << 32)) / 1220703125ull);
		uint64_t br = ((b + (ar << 32)) % 1220703125ull);
		uint64_t cv = ((c + (br << 32)) / 1220703125ull);

		if (av >> 1) {
			// 30 bit shift
			av <<= 30;
			bv <<= 30;
			cv <<= 30;
			pow_2 -= 30 + 13;
		} else {
			// 31 bit shift
			av <<= 31;
			bv <<= 31;
			cv <<= 31;
			pow_2 -= 31 + 13;
		}

		// enable c to make b reuse old bits
		c = (cv & 0xffffffff) * (enable_c);
		raw = (av << 32) + (bv << 0) + (cv >> 32);
		pow_10 += 13;

		/*
		std::cout
			<< std::bitset<32>(a >> 32) << ":"
			<< std::bitset<32>(a) << ":"
			<< std::bitset<32>(b) << ":"
			<< std::bitset<32>(c) << std::endl;
		std::cout
			<< std::bitset<32>(av >> 32) << ":"
			<< std::bitset<32>(av) << ":"
			<< std::bitset<32>(bv) << ":"
			<< std::bitset<32>(cv) << std::endl;
		 */

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	while (pow_2 > 0) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;

		uint64_t av = (a / 5ull);
		uint64_t ar = (a % 5ull);
		uint64_t bv = ((b + (ar << 32)) / 5ull);
		uint64_t br = ((b + (ar << 32)) % 5ull);
		uint64_t cv = ((c + (br << 32)) / 5ull);

		if (av >> 29) { // 31 - shf
			// 2 bit shift
			av <<= 2; // shf
			bv <<= 2; // shf
			cv <<= 2; // shf
			pow_2 -= 2 + 1; // idx + shf
		} else {
			// 3 bit shift
			av <<= 3; // shf + 1
			bv <<= 3; // shf + 1
			cv <<= 3; // shf + 1
			pow_2 -= 3 + 1; // idx + shf + 1
		}

		// enable c to make b reuse old bits
		c = (cv & 0xffffffff) * (enable_c);
		raw = (av << 32) + (bv << 0) + (cv >> 32);
		pow_10 += 1;

		/*
		std::cout
			<< std::bitset<3>(a >> 32) << ":"
			<< std::bitset<32>(a) << ":"
			<< std::bitset<32>(b) << std::endl;
		*/

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	if (pow_2 < 0) {
		raw >>= -pow_2;
		pow_2 = 0;
	}
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;

	// Make sure raw is 19 digits
	if (raw >= 10'000'000'000'000'000'000ull) {
		raw /= 10;
		pow_10 += 1;
	}

	return { raw, pow_10 };
}

inline d_norm s_to_d_norm_64_opt(uint64_t raw, int pow_10) {
	constexpr bool enable_c = true;
	uint64_t c = 0;
	int pow_2 = 0;
	
	while(raw > 0 && (raw >> 63) != 1) {
		raw <<= 1;
		pow_2 -= 1;
	}
 
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << std::endl;
	if (pow_10 > 0) {
		while (pow_10 >= 13) {
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;
			// enable c to make b reuse old bits
			c = (c * 1220703125ull) * (enable_c);
			b = (b * 1220703125ull) + (c >> 32);
			a = (a * 1220703125ull) + (b >> 32);
			b = b & 0xffffffff;
			c = c & 0xffffffff;
			pow_10 -= 13;
		
			// std::cout << std::bitset<32>(a >> 32) << ":" << std::bitset<32>(a) << ":" << std::bitset<32>(b) << std::endl;
		
		
			// (a >> 34) // checks if 2 bit needs to be shifted
			// (a >> 62) // checks if 31 bits need to be shifted
			// (a >> 34) -> (32 - 2 - 1) : (0 - 2 - 1)
			// (a >> 62) -> (32 - 30 - 1) : (0 - 30 - 1)
			// (a >> 62) -> (32 - 29 - 1) : (0 - 29 - 1)
		
			if (a >> 62) {
				raw = (a << 1) | (b >> 31);
				c = ((b << 1) | (c >> 31)) & 0xffffffff;
				pow_2 += 31 + 13;
			} else {
				raw = (a << 2) | (b >> 30);
				c = ((b << 2) | (c >> 30)) & 0xffffffff;
				pow_2 += 30 + 13;
			}
		
			// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
		}
		// if (pow_10 > 0) {
			const int idx = pow_10;
			// const int idx = (2 - (pow_2 << 3) - pow_2) / 30;
			const uint64_t mul = NEG_LUT[(idx << 1) | 0];
			const uint64_t shf = NEG_LUT[(idx << 1) | 1];
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;
			c = (c * mul);
			b = (b * mul) + (c >> 32);
			a = (a * mul) + (b >> 32);
			b = b & 0xffffffff;
			pow_10 -= idx;
			// std::cout << mul << ", " << shf << std::endl;
			// std::cout << std::bitset<64>(raw) << std::endl;
			// std::cout << std::bitset<32>(a >> 32) << ":" << std::bitset<32>(a) << ":" << std::bitset<32>(b) << std::endl;

			if (a >> (30 + shf)) { // TODO - should be 31 + shf some math is wrong here
				raw = (a << (32 - shf)) | (b >> shf);
				pow_2 += shf + idx;
			} else {
				raw = (a << (31 - shf)) | (b >> (shf - 1));
				pow_2 += shf + idx - 1;
			}
			// std::cout << std::bitset<64>(raw) << std::endl;
		// }
	} else {
		while (pow_10 <= -13) {
			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;
	
			uint64_t av = (a / 1220703125ull);
			uint64_t ar = (a % 1220703125ull);
			uint64_t bv = ((b + (ar << 32)) / 1220703125ull);
			uint64_t br = ((b + (ar << 32)) % 1220703125ull);
			uint64_t cv = ((c + (br << 32)) / 1220703125ull);
	
			if (av >> 1) {
				// 30 bit shift
				av <<= 30;
				bv <<= 30;
				cv <<= 30;
				pow_2 -= 30 + 13;
			} else {
				// 31 bit shift
				av <<= 31;
				bv <<= 31;
				cv <<= 31;
				pow_2 -= 31 + 13;
			}
	
			// enable c to make b reuse old bits
			c = (cv & 0xffffffff) * (enable_c);
			raw = (av << 32) + (bv << 0) + (cv >> 32);
			pow_10 += 13;
		}
		// if (pow_10 < 0) {
			const int idx = -pow_10; // ((2 + pow_2 * 9) / 30) + 1;
			// const int idx = (2 - (pow_2 << 3) - pow_2) / 30;
			const uint64_t div = NEG_LUT[(idx << 1) | 0];
			const uint64_t shf = NEG_LUT[(idx << 1) | 1];

			uint64_t a = raw >> 32;
			uint64_t b = raw & 0xffffffff;

			uint64_t av = (a / div);
			uint64_t ar = (a % div);
			uint64_t bv = ((b + (ar << 32)) / div);
			uint64_t br = ((b + (ar << 32)) % div);
			uint64_t cv = ((c + (br << 32)) / div);

			if (av >> (32 - shf)) {
				// 2 bit shift
				av <<= shf - 1;
				bv <<= shf - 1;
				cv <<= shf - 1;
				pow_2 -= shf + idx - 1;
			} else {
				// 3 bit shift
				av <<= shf;
				bv <<= shf;
				cv <<= shf;
				pow_2 -= shf + idx;
			}

			// enable c to make b reuse old bits
			c = (cv & 0xffffffff);
			raw = (av << 32) + (bv << 0) + (cv >> 32);
			pow_10 += idx;	
		// }
	}
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
 
	while (raw > 0 && (raw >> 63) != 1) {
		raw <<= 1;
		raw |= (c >> 31) & 1;
		c <<= 1;
		pow_2 -= 1;
	}
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << " E" << std::endl;
 
	return { raw, pow_2 };
}

inline d_norm s_to_d_norm_64(uint64_t raw, int pow_10) {
	constexpr bool enable_c = true;
	uint64_t c = 0;
	int pow_2 = 0;

	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << std::endl;
	while (pow_10 >= 13) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;
		// enable c to make b reuse old bits
		c = (c * 1220703125ull) * (enable_c);
		b = (b * 1220703125ull) + (c >> 32);
		a = (a * 1220703125ull) + (b >> 32);
		b = b & 0xffffffff;
		c = c & 0xffffffff;
		pow_10 -= 13;

		// std::cout
		//	 << std::bitset<32>(a >> 32) << ":"
		//	 << std::bitset<32>(a) << ":"
		//	 << std::bitset<32>(b) << std::endl;


		// (a >> 34) // checks if 2 bit needs to be shifted
		// (a >> 62) // checks if 31 bits need to be shifted
		// (a >> 34) -> (32 - 2 - 1) : (0 - 2 - 1)
		// (a >> 62) -> (32 - 30 - 1) : (0 - 30 - 1)
		// (a >> 62) -> (32 - 29 - 1) : (0 - 29 - 1)

		if (a >> 62) {
			raw = (a << 1) | (b >> 31);
			c = ((b << 1) | (c >> 31)) & 0xffffffff;
			pow_2 += 31 + 13;
		} else {
			raw = (a << 2) | (b >> 30);
			c = ((b << 2) | (c >> 30)) & 0xffffffff;
			pow_2 += 30 + 13;
		}

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	while (pow_10 > 0) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;
		// enable c to make b reuse old bits
		c = (c * 5ull) * (enable_c);
		b = (b * 5ull) + (c >> 32);
		a = (a * 5ull) + (b >> 32);
		b = b & 0xffffffff;
		pow_10 -= 1;

		/*std::cout
			<< std::bitset<3>(a >> 32) << ":"
			<< std::bitset<32>(a) << ":"
			<< std::bitset<32>(b) << std::endl;
		*/

		if (a >> 34) {
			raw = (a << 29) | (b >> 3);
			c = ((b << 29) | (c >> 3)) & 0xffffffff;
			pow_2 += 3 + 1;
		} else {
			raw = (a << 30) | (b >> 2);
			c = ((b << 30) | (c >> 2)) & 0xffffffff;
			pow_2 += 2 + 1;
		}

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	while (pow_10 <= -13) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;

		uint64_t av = (a / 1220703125ull);
		uint64_t ar = (a % 1220703125ull);
		uint64_t bv = ((b + (ar << 32)) / 1220703125ull);
		uint64_t br = ((b + (ar << 32)) % 1220703125ull);
		uint64_t cv = ((c + (br << 32)) / 1220703125ull);

		if (av >> 1) {
			// 30 bit shift
			av <<= 30;
			bv <<= 30;
			cv <<= 30;
			pow_2 -= 30 + 13;
		} else {
			// 31 bit shift
			av <<= 31;
			bv <<= 31;
			cv <<= 31;
			pow_2 -= 31 + 13;
		}

		// enable c to make b reuse old bits
		c = (cv & 0xffffffff) * (enable_c);
		raw = (av << 32) + (bv << 0) + (cv >> 32);
		pow_10 += 13;
	}
	while (pow_10 < 0) {
		uint64_t a = raw >> 32;
		uint64_t b = raw & 0xffffffff;

		uint64_t av = (a / 5ull);
		uint64_t ar = (a % 5ull);
		uint64_t bv = ((b + (ar << 32)) / 5ull);
		uint64_t br = ((b + (ar << 32)) % 5ull);
		uint64_t cv = ((c + (br << 32)) / 5ull);

		if (av >> 29) {
			// 2 bit shift
			av <<= 2;
			bv <<= 2;
			cv <<= 2;
			pow_2 -= 2 + 1;
		} else {
			// 3 bit shift
			av <<= 3;
			bv <<= 3;
			cv <<= 3;
			pow_2 -= 3 + 1;
		}

		// enable c to make b reuse old bits
		c = (cv & 0xffffffff) * (enable_c);
		raw = (av << 32) + (bv << 0) + (cv >> 32);
		pow_10 += 1;

		/*
		std::cout
			<< std::bitset<3>(a >> 32) << ":"
			<< std::bitset<32>(a) << ":"
			<< std::bitset<32>(b) << std::endl;
		*/

		// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;
	}
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << std::endl;

	while(raw > 0 && (raw >> 63) != 1) {
		raw <<= 1;
		raw |= (c >> 31) & 1;
		c <<= 1;
		pow_2 -= 1;
	}
	// std::cout << raw << ":" << std::bitset<64>(raw) << ":" << pow_2 << " / " << pow_10 << " E" << std::endl;

	return { raw, pow_2 };
}

int d_to_s_a(std::array<char, 32>& buffer, double a_value) {
	uint64_t repr{};
	std::memcpy(&repr, &a_value, 8);
	uint64_t raw_mantissa = (repr & ((1ull << 52) - 1ull));
	uint64_t raw_exponent = (repr << 1) >> 53;
	uint64_t sign = (repr >> 63);

	int index = 0;
	if (raw_exponent == 0 && raw_mantissa == 0) {
		if (sign != 0) {
			buffer[index++] = '-';
		}
		buffer[index++] = '0';
		buffer[index++] = '.';
		buffer[index++] = '0';
		buffer[index++] = 'e';
		buffer[index++] = '+';
		buffer[index++] = '0';
		return index;
	}
	if (raw_exponent == 2047) {
		if (raw_mantissa == 0) {
			if (sign != 0) {
				buffer[index++] = '-';
			}

			buffer[index++] = 'i';
			buffer[index++] = 'n';
			buffer[index++] = 'f';
		} else {
			buffer[index++] = 'N';
			buffer[index++] = 'a';
			buffer[index++] = 'N';
		}
		return index;
	}

	uint64_t mantissa = (1ull << 52) | raw_mantissa;
	int32_t exponent = static_cast<int32_t>(raw_exponent) - (1ull << 10) - 1 - 50;

	// (2 ** exponent) * mantissa = value
	// we want the 19 digits in base 10 that describe value such that
	// A * (10 ** B) ~= (2 ** exponent) * mantissa
	/*
	if (exp >= 0) {
		// We just have 'mantissa' followed by 'exp' zero bits
		div = (mantissa << 11) / 10
		rem = (mantissa << 11) % 10

		// we could save precision on rem
	}
	*/

	if (sign != 0) {
		buffer[index++] = '-';
	}

	auto [raw, pow10] = d_to_s_norm_64_opt(mantissa << 11, exponent - 11);

	// Calculate how many digits are in raw
	int digits = 1;
	uint64_t compare = 10'000'000'000'000'000'000ull;
	//                  9'223'372'036'854'775'807ull
	//                 18'446'744'073'709'551'615ull
	for (int i = 20; i > 1; --i) {
		if (raw >= compare) {
			digits = i;
			break;
		}
		compare /= 10;
	}

	for(int i = digits - 1; i >= 0; --i) {
		if (i == 0) {
			buffer[index] = (raw % 10) + '0';
			buffer[index + 1] = '.';
			break;
		}
		char digit = (raw % 10) + '0';
		raw /= 10;
		buffer[index + i + 1] = digit;
	}
	pow10 += digits - 1;
	index += digits + 1;

	// Add scientific
	buffer[index++] = 'e';
	buffer[index++] = pow10 < 0 ? '-' : '+';
	pow10 = pow10 < 0 ? -pow10 : pow10;

	if(pow10 < 100) {
		buffer[index++] = '0' + (pow10 / 10);
		buffer[index++] = '0' + (pow10 % 10);
	} else {
		buffer[index++] = '0' + (pow10 / 100);
		buffer[index++] = '0' + ((pow10 / 10) % 10);
		buffer[index++] = '0' + (pow10 % 10);
	}
	return index;
}

std::string d_to_s(double a_value) {
	std::array<char, 32> buffer{};
	int end = d_to_s_a(buffer, a_value);
	buffer[end] = '\0';
	return std::string(buffer.data());
}

double s_to_d_a(const char* start, const char* end) {
	bool has_frac = false;
	bool skip_dig = false;
	int pow_10 = 0;
	int sign = 0;

	if(start == end) {
		return 0.0 / 0.0; // NaN
	}

	if(*start == '-') {
		sign = 1;
		start++;
	} else if(*start == '+') {
		start++;
	}

	// value can hold up to 19 digits
	uint64_t value = 0;
	while(start != end) {
		char c = *(start++);
		if(c >= '0' && c <= '9') {
			if(value < 1'000'000'000'000'000'000ull) {
				value = value * 10 + (c - '0');
			} else {
				pow_10 += 1;
			}

			if(has_frac) {
				pow_10 -= 1;
			} else if (value > 9) {
				pow_10 += 1;
			}
		} else if(c == '.') {
			// invalid
			if(has_frac) {
				break;
			}

			has_frac = true;
		} else if(c == 'e') {
			int exp_val = 0;
			int exp_sign = 1;
			if(*start == '-') {
				start++;
				exp_sign = -1;
			} else if(*start == '+') {
				start++;
			}

			while(start != end) {
				c = *(start++);
				if(c < '0' || c > '9') {
					break;
				}
				exp_val = exp_val * 10 + (c - '0');
			}

			pow_10 += exp_val * exp_sign;
			break;
		} else {
			break;
		}
	}

	// Make value 19 digits
	while (value > 0 && value < 1'000'000'000'000'000'000ull) {
		value *= 10;
		pow_10 -= 1;
	}

	// Represent float exactly
	auto [v, e] = s_to_d_norm_64_opt(value, pow_10);

	// Round value up
	if((v & 0b111'1111'1111) > 0b100'0000'0000) {
		if(v + 0b1000'0000'0000 > v) { // not wrap
			v &= ~0b0111'1111'1111;
			v +=  0b1000'0000'0000;
		} else { // wrap
			e += 1;
			v &= ~0b111'1111'1111;
			v >>= 1;
			v +=  0b100'0000'0000;
		}
	}
	// auto [va, ea] = d_to_s_norm_64((v >> 11), e + 11);
	// std::cout << value << ":" << pow_10 << std::endl;
	// std::cout << va << ":" << ea << std::endl;
	double result = (v >> 11); // Exact representation
	result *= std::pow(2.0, e + 11);
	// e += 11;
	// for(int i = 0; i < e; ++i) result *= 2.0;
	// for(int i = 0; i > e; --i) result *= 0.5;

	return sign ? -result : result;
}

double s_to_d(std::string a_value) {
	return s_to_d_a(a_value.c_str(), a_value.c_str() + a_value.size());
}

#endif  // DTOA_H
