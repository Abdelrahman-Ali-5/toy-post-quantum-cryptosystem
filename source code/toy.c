/*
 * toy.c
 *
 *  Created on: Dec 18 , 2023
 *      Author: Abdelrahman Ali
 *
 */



#include "toy.h"
#include <stdlib.h>
#include "stdio.h"


// polynomial multiplication in Z97[X]/(X^4+1)
//if  ( add_to_res = true : res += a*b; if false: res = a*b;

static void toy_polmul_naive(short* res, const short* a, const short* b, int add_to_res){
	short temp[TK_N]; // Temporary array to hold the result before deciding whether to add or replace res

	// Perform polynomial multiplication using coefficients of a and b
	temp[0] = (a[0] * b[0] - a[3] * b[1] - a[2] * b[2] - a[1] * b[3]) % TK_Q;
	temp[1] = (a[1] * b[0] + a[0] * b[1] - a[3] * b[2] - a[2] * b[3]) % TK_Q;
	temp[2] = (a[2] * b[0] + a[1] * b[1] + a[0] * b[2] - a[3] * b[3]) % TK_Q;
	temp[3] = (a[3] * b[0] + a[2] * b[1] + a[1] * b[2] + a[0] * b[3]) % TK_Q;

	// Ensure results are within the modulo range
	for (int i = 0; i < TK_N; ++i) {
		if (temp[i] < 0)
			temp[i] += TK_Q; 
		// If the result is negative, adjust it to fit within the range [0, 96]
		// this is an extra check 
	}

	// Decide whether to add or replace res based on the add_to_res flag
	if (add_to_res) {
		for (int i = 0; i < TK_N; ++i) {
			res[i] += temp[i];  // Add the result to the existing contents of res
			res[i] %= TK_Q;     // Ensure the result remains within the modulo range
		}
	}
	else {
		for (int i = 0; i < TK_N; ++i)
			res[i] = temp[i];   // Replace the contents of res with the result from the temporary array
	}
}


void multiplyAndAccumulate(const short* A, const short* s, short A_RowCount, short* result) {
	int count = 0;

	// Perform pairwise multiplication and accumulation
	for (int i = 0; i < A_RowCount; i++)
	{
		short temporaryResult[TK_N] = { 0 }; // Temporary array to store the multiplication result
		for (int j = 0; j < 12; j += 4)
		{
			toy_polmul_naive(temporaryResult, &A[j + 12 * i], &s[j], 1); // Perform multiplication
			count++;

			// Accumulate every 3 multiplications into result array
			if (count == 3)
			{
				for (int k = 0; k < TK_N; k++)
				{
					result[k + i * TK_N] = 0;
					result[k + i * TK_N] = temporaryResult[k];
				}
				count = 0;
			}
		}
	}
}


void toy_gen(short* A, short* t, short* s)
{
	// Fill K*K-matrix A with uniformly random numbers mod q
	for (int i = 0; i < (TK_K * TK_K * TK_N); ++i) {
		A[i] = rand() % TK_Q;
	}

	// Generate small random numbers mod q for e
	short e[TK_K * TK_N];
	for (int i = 0; i < (TK_K * TK_N); ++i) {
		short val = rand() & 3;
		e[i] = (val & 1) - ((val >> 1) & 1);
	}

	// Generate small random numbers mod q for s
	for (int i = 0; i < (TK_K * TK_N); ++i) {
		short val = rand() & 3;
		s[i] = (val & 1) - ((val >> 1) & 1);
	}

	short result[TK_K * TK_N] = { 0 };

	// Perform matrix-vector multiplication: result = A.s
	multiplyAndAccumulate(A, s, 3, result);

	// Perform vector addition: t = result + e
	for (int i = 0; i < (TK_K * TK_N); ++i)
	{
		t[i] = 0;
		t[i] += result[i] + e[i]; // Add the corresponding e value
		t[i] %= TK_Q;             // Ensure the result remains within the modulo range
	}
}

// function to get the transpose of a matrix given number of rows and coloumns
void transpose(short* A, short* A_transposed, int k, int n) {
	// Iterate over the original array
	for (int i = 0; i < k * k; ++i) {
		// Calculate indices for the transposed array
		int row = i / k;
		int col = i % k;

		// Calculate the index in the transposed array
		int transposed_index = col * k * n + row * n;

		// Copy the polynomial from A to its transposed position in A_transposed
		for (int j = 0; j < n; ++j) {
			A_transposed[transposed_index + j] = A[i * n + j];
		}
	}
}


void toy_enc(short* A, short* t, int plain, short* u, short* v) {

	// Generate small random numbers mod q for e1
	short e1[TK_K * TK_N];
	for (int i = 0; i < (TK_K * TK_N); ++i) {
		short val = rand() & 3;
		e1[i] = (val & 1) - ((val >> 1) & 1);
	}

	// Generate small random numbers mod q for e2
	short e2[TK_N];
	for (int i = 0; i < TK_N; ++i) {
		short val = rand() & 3;
		e2[i] = (val & 1) - ((val >> 1) & 1);
	}

	// Generate small random numbers mod q for r
	short r[TK_K * TK_N];
	for (int i = 0; i < (TK_K * TK_N); ++i) {
		short val = rand() & 3;
		r[i] = (val & 1) - ((val >> 1) & 1);
	}

	// Transpose matrix A
	short A_transpose[TK_K * TK_K * TK_N];
	transpose(A, A_transpose, TK_K, TK_N);

	short u_result[TK_K * TK_N] = { 0 };

	// Perform matrix-vector multiplication: u_result = A_transpose.r
	multiplyAndAccumulate(A_transpose, r, 3, u_result);

	// Perform vector addition: u = result + e1
	for (int i = 0; i < (TK_K * TK_N); ++i)
	{
		u[i] = u_result[i] + e1[i]; // Add the corresponding e value
		u[i] %= TK_Q;             // Ensure the result remains within the modulo range
	}


	short v_result[TK_N] = { 0 };

	// Perform matrix-vector multiplication: v_result = t_transpose.r
	// Transpose matrix t is the same as t so we use t directly
	multiplyAndAccumulate(t, r, 1, v_result);

	short msg_bits[TK_N] = { 0 };

	for (int i = 0; i < TK_N; ++i)
	{
		msg_bits[i] = (plain >> i) & 1;
	}

	// Perform matrix-vector multiplication: v = tT.r + e2 + Round(plain * q/2)
	for (int i = 0; i < TK_N; ++i)
	{
		v[i] = v_result[i] + e2[i] + (msg_bits[i] * (TK_Q >> 1)); // Add the corresponding e value
		v[i] %= TK_Q; // Modulo operation
	}
}



int toy_dec(const short* s, const short* u, const short* v) {

	int plain = 0;

	short p[TK_N];

	short s_u_result[TK_N] = { 0 };

	// Perform matrix-vector multiplication: s_u_result = s_transpose.u
	// transpose s matrix is the same as s so we use s directly
	multiplyAndAccumulate(s, u, 1, s_u_result);

	// Perform matrix-vector multiplication: p =  v - s_u_result
	for (int i = 0; i < TK_N; ++i)
	{
		p[i] = v[i] - s_u_result[i];
		p[i] %= TK_Q; // Modulo operation
	}

	int bit = 0;

	// extract the plaintext bits
	for (int i = 0; i < TK_N; ++i) {
		if (p[i] > (TK_Q / 2)) {
			p[i] -= TK_Q;
		}
		bit = abs(p[i]) > (TK_Q / 4);
		plain |= bit << i;
	}

	return plain;
}






