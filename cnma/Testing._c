#include "matrix.h"
void test();

void test(){
	char num1[100], num2[100], op[4];
	mpz_t n1, n2, res, mod_min, mod_plus;
	int n;
	clock_t start, end;
	float sec;
	mpz_init(n1); 
	mpz_init(n2); 
	mpz_init(mod_min); 
	mpz_init(mod_plus);	
	mpz_init(res);

	printf("Enter numbers as following:\n");
	printf("num\nnum\nn\n");
	scanf(" %s", num1);
	scanf(" %s", num2);
	scanf(" %d", &n); 

	mpz_set_str (n1, num1, 10);
	mpz_set_str (n2, num2, 10);
	get_mod(mod_min, n);
	get_mod_plus(mod_plus, n);

	printf("%9s%40s%12s%2s%40s%12s\n", "|", "2^n-1", "Time (s)", "|", "2^n+1", "Time (s)");	
	printf("--------|-----------------------------------------------------|-----------------------------------------------------\n");

//MODADD
	start = clock();
	modadd(res, n1, n2, mod_min);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%9s%40Zd%12f%2s", "modadd |", res, sec, "|");
	start = clock();
	modadd(res, n1, n2, mod_plus);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);

//MODSUB
	start = clock();
	modsub(res, n1, n2, mod_min);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%9s%40Zd%12f%2s", "modsub |", res, sec, "|");
	start = clock();
	modsub(res, n1, n2, mod_plus);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);

//MODMUL
	start = clock();
	modmul(res, n1, n2, mod_min);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%9s%40Zd%12f%2s", "modmul |", res, sec, "|");
	start = clock();
	modmul(res, n1, n2, mod_plus);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);
		printf("--------|-----------------------------------------------------|-----------------------------------------------------\n");

//MyADD
	start = clock();
	minadd(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	//gmp_printf("%Zd, %Zd\n", n1, n2);
//printf("\n");
	gmp_printf("%9s%40Zd%12f%2s", "2n add |", res, sec, "|");
	start = clock();
	plusadd(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);
//MySUB
	start = clock();
	minsub(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%9s%40Zd%12f%2s", "2n sub |", res, sec, "|");
	start = clock();
	plussub(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);
//MyMUL
	start = clock();
	minmul(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%9s%40Zd%12f%2s", "2n mul |", res, sec, "|");
	start = clock();
	plusmul(res, n1, n2, n);
	end = clock();
	sec = (float)(end - start) / CLOCKS_PER_SEC;
	gmp_printf("%40Zd%12f\n", res, sec);		

	
	printf("\n");
	
	
}

