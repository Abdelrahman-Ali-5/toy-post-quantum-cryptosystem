/*
 * main.c
 *
 *  Created on: Dec 18, 2023
 *      Author: Abdelrahman Ali
 */

#include "toy.h"

int main ( ){

    short A[TK_K * TK_K * TK_N];
    short t[TK_K * TK_N];
    short s[TK_K * TK_N];

    short u[TK_K * TK_N];
    short v[TK_N];

    // from 0 to 15 it's the same
    // above 15 it will be reduced by modulo into the range of 0 -> 15
    for( int i = 0 ; i < 32 ; i++){

        int plain = i;
        int plain_after_decrypt = 0;

        printf("%d \t",plain);

        // Call toy_gen to initialize A, t, and s
        toy_gen(A, t, s);

        toy_enc(A, t,  plain,  u,  v);

        plain_after_decrypt = toy_dec(s,u,v);

        printf("%d \n",plain_after_decrypt);
    }
}
