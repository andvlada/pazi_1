#include "jacobi.h"
#include <stdlib.h>

void init_point(struct Point *P, mpz_t X, mpz_t Y, mpz_t T, mpz_t Z){
    mpz_init_set_ui(P->x, X);
    mpz_init_set_ui(P->y, Y);
    mpz_init_set_ui(P->t, T);
    mpz_init_set_ui(P->z, Z);
}

void init_curve(struct Curve *C){
    mpz_init_set_str(C->p, p_s, 10);
    mpz_init_set_str(C->q, q_s, 10);
    mpz_init_set_str(C->b, b_s, 10);
}

struct Point add_points(struct Point p_1, struct Point p_2){
    struct Point new_p;
    struct Curve C;
    init_curve(&C);
    mpz_t help;
    mpz_init(help);
    init_point(&new_p, help, help, help, help);

    //x3 = x1*z1*y2*t2 + y1*t1*x2*z2
    mpz_mul(new_p.x, p_1.x, p_1.z);
    mpz_mul(new_p.x, new_p.x, p_2.y);
    mpz_mul(new_p.x, new_p.x, p_2.t);
    mpz_mul(help, p_1.y, p_1.t);
    mpz_mul(help, help, p_2.x);
    mpz_mul(help, help, p_2.z);
    mpz_add(new_p.x, new_p.x, help);
    mpz_mod(new_p.x, new_p.x, C.p);

    //y3 = y1*z1*y2*z2 - x1*t1*x2*t2
    mpz_set_str(help, "0", 10);
    mpz_mul(new_p.y, p_1.y, p_1.z);
    mpz_mul(new_p.y, new_p.y, p_2.y);
    mpz_mul(new_p.y, new_p.y, p_2.z);
    mpz_mul(help, p_1.x, p_1.t);
    mpz_mul(help, help, p_2.x);
    mpz_mul(help, help, p_2.t);
    mpz_sub(new_p.y, new_p.y, help);
    mpz_mod(new_p.y, new_p.y, C.p);

    //t3 = t1*z1*t2*z2 - b*x1*y1*x2*y2
    mpz_set_str(help, "0", 10);
    mpz_mul(new_p.t, p_1.t, p_1.z);
    mpz_mul(new_p.t, new_p.t, p_2.t);
    mpz_mul(new_p.t, new_p.t, p_2.z);
    mpz_mul(help, C.b, p_1.x);
    mpz_mul(help, help, p_1.y);
    mpz_mul(help, help, p_2.x);
    mpz_mul(help, help, p_2.y);
    mpz_sub(new_p.t, new_p.t, help);
    mpz_mod(new_p.t, new_p.t, C.p);

    //z3 = z1^2*y2^2 + x2^2*t1^2
    mpz_set_str(help, "0", 10);
    mpz_mul(new_p.z, p_1.z, p_1.z);
    mpz_mul(new_p.z, new_p.z, p_2.y);
    mpz_mul(new_p.z, new_p.z, p_2.y);
    mpz_mul(help, p_2.x, p_2.x);
    mpz_mul(help, help, p_1.t);
    mpz_mul(help, help, p_1.t);
    mpz_add(new_p.z, new_p.z, help);
    mpz_mod(new_p.z, new_p.z, C.p);

    mpz_clear(help);

    return new_p;
}

struct Point double_point(struct Point P){
    mpz_t yz, tz, yt;
    mpz_init(yz);
    mpz_init(tz);
    mpz_init(yt);
    struct Point new_p;
    struct Curve C;
    init_curve(&C);
    init_point(&new_p, yz, yz, yz, yz);

    //x3 = 2*x*y*t*z
    mpz_set_str(yz, "2", 10);
    mpz_mul(new_p.x, P.x, P.y);
    mpz_mul(new_p.x, new_p.x, P.t);
    mpz_mul(new_p.x, new_p.x, P.z);
    mpz_mul(new_p.x, new_p.x, yz);
    mpz_mod(new_p.x, new_p.x, C.p);

    //yz = (y*z)^2
    mpz_mul(yz, P.y, P.y);
    mpz_mul(yz, yz, P.z);
    mpz_mul(yz, yz, P.z);

    //tz = (t*z)^2
    mpz_mul(tz, P.t, P.t);
    mpz_mul(tz, tz, P.z);
    mpz_mul(tz, tz, P.z);

    //yt = (y*t)^2
    mpz_mul(yt, P.y, P.y);
    mpz_mul(yt, yt, P.t);
    mpz_mul(yt, yt, P.t);

    //y3 = (y*z)^2 - (t*z)^2 + (y*t)^2
    mpz_sub(new_p.y, yz, tz);
    mpz_add(new_p.y, new_p.y, yt);
    mpz_mod(new_p.y, new_p.y, C.p);

    //t3 = (t*z)^2 - (y*z)^2 + (y*t)^2
    mpz_sub(new_p.t, tz, yz);
    mpz_add(new_p.t, new_p.t, yt);
    mpz_mod(new_p.t, new_p.t, C.p);

    //z3 = (t*z)^2 + (y*z)^2 - (y*t)^2
    mpz_add(new_p.z, tz, yz);
    mpz_sub(new_p.z, new_p.z, yt);
    mpz_mod(new_p.z, new_p.z, C.p);

    mpz_clear(yz);
    mpz_clear(tz);
    mpz_clear(yt);

    return new_p;
}

struct Point multiple_point(struct Point P, mpz_t k){
    struct Point Q, R;
    struct Curve C;
    mpz_t help;
    init_curve(&C);
    mpz_init_set_str(help, "1", 10);
    init_point(&Q, 0, help, help, help);
    init_point(&R, P.x, P.y, P.t, P.z);
    int n = mpz_sizeinbase(k, 2);
    for (int i = n - 1; i >= 0; i--){
        if (mpz_tstbit(k, i) == 0){
            R = add_points(R, Q);
            Q = double_point(Q);
        }
        else {
            Q = add_points(Q, R);
            R = double_point(R);
        }
    }

    mpz_clear(help);

    return Q;
}

int if_contains(struct Point P){
    mpz_t left1, right, left2, help;
    mpz_init(left1);
    mpz_init(right);
    mpz_init(left2);
    mpz_init(help);
    struct Curve C;
    init_curve(&C);

    //left1
    mpz_mul(left1, P.x, P.x);
    mpz_set(left2, left1);
    mpz_mul(help, P.y, P.y);
    mpz_add(left1, left1, help);
    mpz_mod(left1, left1, C.p);

    //left2
    mpz_mul(left2, left2, C.b);
    mpz_mul(help, P.t, P.t);
    mpz_add(left2, left2, help);
    mpz_mod(left2, left2, C.p);

    //right
    mpz_mul(right, P.z, P.z);
    mpz_mod(right, right, C.p);

    if(mpz_cmp(left1, right) == 0 && mpz_cmp(left2, right) == 0)
        return 1;
    else
        return 0;

    mpz_clear(left1);
    mpz_clear(right);
    mpz_clear(left2);
    mpz_clear(help);
}
