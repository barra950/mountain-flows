#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(){
    
int K = 200;
long double alpha = 0.1;
long double visc = 5;  
long double diff = 5;  
long double N = 0.01;    
long double L = 5000; 
long double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

long double tick = 0.1;

long double H(long double y) {
    return ( 300 * (1 + cos(2 * pi * y/L)) );
    //return 0;
}

long double Bsfc(long double y) {
    //return ( 0.1 * (1 + cos(2 * pi * y/L)) );
    return 0.1;
}

// printf("%lf",H(0.5));
// printf("%lf",Bsfc(0.5));

long double R;
long double Q;
long double S1;
long double S2;
long double phi;
long double Lk;
long double m1;
long double complex m2;
long double complex m3;


// int Aki[10*K];
// int Cki[10*K];
// int Dki[10*K];
long double final_system[10*K];
long double b[10*K];  
for (int q = -K; q <= K; q++){
    long double complex equation1[6*K+2];
    long double complex equation2[6*K+2];
    long double complex equation3[6*K+2];
    int ctr = 0;
    int ctr2 = 0;
    int ctr3 = 0;
    for (int k = -K; k <= K; k++){
        R = 2 * pow(N,2) * pow(cos(alpha),2) / (visc * diff) * pow ((k * pi / L),2);

        Q = pow(N,2) * pow(sin(alpha),2) / (3 * visc * diff);

        S1 = pow(fabs(R + sqrt(pow(Q,3) + pow(R,2)) ),(1.0/3));
        S2 = - pow(fabs( sqrt(pow(Q,3) + pow(R,2)) - R ),(1.0/3));

        phi = sqrt(pow(S1,2) + pow(S2,2) - S1*S2);

        Lk = acos(- (S1 + S2)/ (2 * phi) );

        m1 = - sqrt(S1 + S2);
        m2 = - sqrt(phi) * cexp(I * Lk/2);
        m3 = conjf(m2);

        long double f1(long double y){
            return (exp(m1 * H(y)) * cos(2 * (q - k) * pi * y / L) );
        }
        long double complex f2(long double y){
            return (cexp(m2 * H(y)) * ccos(2 * (q - k) * pi * y / L) );
        }
        //printf("%.30lf\n",f1(2.45));
        //printf("%.30lf + %.30lfi\n",creal(f2(2.45)),cimag(f2(2.45)));

        long double complex gamma1 = 0.0;
        long double complex gamma2 = 0.0;
        for (long double i = 0; i <= L/2; i = i + tick){
            gamma1 = 2/L * f1(i)*tick + gamma1;
            gamma2 = 2/L * f2(i)*tick + gamma2;
        }
        //printf("%.30lf + %.30lfi\n",creal(gamma2),cimag(gamma2));

        if (k == 0){
            equation1[ctr] = (2 * creal(gamma2));
            //printf("%.30lf + %.30lfi\n",creal(equation1[0]),cimag(equation1[0]));
            //printf("%d\n",ctr);
            ctr++;
            //Cki[k] = k;
            equation1[ctr] = (-2 * cimag(gamma2));
            ctr++;
            //Dki[k] = k;
        }
        else{
            equation1[ctr] = gamma1;
            ctr++;
            //Aki[k] = k;
            equation1[ctr] = (2 * creal(gamma2));
            ctr++;
            //Cki[k]= k;
            equation1[ctr] = (-2 * cimag(gamma2));
            ctr++;
            //Dki[k] = k;
        }

        if (q != 0){
            
            if (k == 0){
                equation2[ctr2] = 0;
                ctr2++;
                equation2[ctr2] = 0;
                ctr2++;
            }
            else {
                equation2[ctr2] = k * gamma1 / pow(m1,3);
                ctr2++;
                equation2[ctr2] = (2 * k * creal(gamma2 / cpow(m2,3) ));
                ctr2++;
                equation2[ctr2] = (-2 * k * cimag(gamma2 / cpow(m2,3) ));
                ctr2++;
            }

        }

        if (k == 0){
            equation3[ctr3] = (2 * creal(cpow(m2,2) * gamma2));
            ctr3++;
            equation3[ctr3] = (-2 * cimag(cpow(m2,2) * gamma2));
            ctr3++;
        }
        else{
            equation3[ctr3] = cpow(m1,2) * gamma1;
            ctr3++;
            equation3[ctr3] = (2 * creal(cpow(m2,2) * gamma2));
            ctr3++;
            equation3[ctr3] = (-2 * cimag(cpow(m2,2) * gamma2));
            ctr3++;
        }
    
    }
    //printf("%.30lf + %.30lfi",creal(equation2[0]),cimag(equation2[0]));
    for (int k = 0; k <= 6*K + 1; k++){
        //printf("%.30lf  %.30lf oo ",creal(equation1[k]),cimag(equation1[k]));
        printf("%.100Lfoo",creall(equation1[k]));
    }
    printf("\n");
    
    long double f4(long double y){
        return (Bsfc(y) * cos(2 * q * pi * y / L) );
    }
    long double gamma4 = 0.0;
    for (long double i = 0; i <= L/2; i = i + tick){
        gamma4 = 2/L * f4(i)*tick + gamma4;
    }
    printf("%.60Lf\n",gamma4);

    if (q != 0){
        for (int k = 0; k <= 6*K + 1; k++){
        // printf("%.30lf  %.30lf oo ",creal(equation2[k]),cimag(equation2[k]));
        printf("%.100Lfoo",creall(equation2[k]));
        }
    printf("\n");

    printf("%d\n",0);
    }

    for (int k = 0; k <= 6*K + 1; k++){
        //printf("%.30lf  %.30lf oo ",creal(equation3[k]),cimag(equation3[k]));
        printf("%.100Lfoo",creall(equation3[k]));
    }
    printf("\n");
    printf("%d\n",0);

}




//printf("%.30lf\n",Q);
//printf("%.30lf + %.30lfi\n",creal(m2),cimag(m2));
//printf("%.30lf + %.30lfi\n",creal(hello(R)),cimag(hello(R)));


// double real = 1.3,imag = 4.9;
// double complex z = CMPLX(real, imag);

// //double complex rato = 2*z;
// double complex rato = ccos(z);  

// printf(
//         "z = %f% + fi\n",
//         creal(z), cimag(z));

//     printf(
//         "rato = %f% + fi\n",
//         creal(rato), cimag(rato));


return 0;

}
