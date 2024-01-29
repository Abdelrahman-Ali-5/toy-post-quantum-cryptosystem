/*
 * toy.h
 *
 *  Created on: Dec 18, 2023
 *      Author: Abdelrahman Ali
 */

#ifndef _TOY_H

/**============ toy Post-Quantum Public-Key Cryptosystem defines ================**/
#define TK_K 3     // level of security ( matrix size )
#define TK_N 4    // Coefficient modulus
#define TK_Q 97  // Field modulus


void toy_gen(short* A, short* t, short* s);
void toy_enc(short* A, short* t, int plain, short* u, short* v) ;
int toy_dec(const short* s, const short* u, const short* v);

#endif // !_TOY_H
