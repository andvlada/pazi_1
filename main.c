#include <stdio.h>
#include "jacobi.h"

int main()
{
    struct Curve C;
    struct Point P;
    init_curve(&C);

    mpz_t x1, y1, t1, z1;
    mpz_init_set_str(x1, x_s, 10);
    mpz_init_set_str(y1, y_s, 10);
    mpz_init_set_str(z1, t_s, 10);
    mpz_init_set_str(t1, z_s, 10);
    init_point(&P, x1, y1, t1, z1);

    //тест 1
    mpz_t k, n;
    mpz_init(k);
    mpz_init(n);
    mpz_set_str(n, "100000000000000000000000000", 10);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_urandomm(k, state, n);

    printf("1 тест. Является ли точка [k]P частью кривой: ");
    if (if_contains(multiple_point(P, k)) == 1)
        printf("да\n");
    else
        printf("нет\n");

    //тест 2
    struct Point N;
    mpz_set_str(x1, "1", 10);
    mpz_set_str(y1, "0", 10);
    init_point(&N, y1, x1, x1, x1);

    struct Point Test_2 = multiple_point(P, C.q);

    printf("2 тест. Равны ли нейтральная точка и [q]P: ");
    if (N.x == Test_2.x && N.y == Test_2.y && N.t == Test_2.t && N.z == Test_2.z)
        printf("да\n");
    else
        printf("нет\n");

    //тест 3
    printf("3.1 тест. Верно ли, что [q+1]P = P: ");
    mpz_t qq;
    mpz_init(qq);
    mpz_set_str(qq, "1", 10);
    mpz_add(qq, qq, C.q);
    struct Point Test_3 = multiple_point(P, qq);
    if (P.x == Test_3.x && P.y == Test_3.y && P.t == Test_3.t && P.z == Test_3.z)
        printf("да\n");
    else
        printf("нет\n");

    printf("3.2 тест. Верно ли, что [q-1]P = -P: ");
    mpz_set_str(qq, "1", 10);
    mpz_sub(qq, C.q, qq);
    Test_3 = multiple_point(P, qq);
    mpz_set_str(qq, "-1", 10);
    mpz_mul(qq, qq, P.x);
    if (qq == Test_3.x && P.y == Test_3.y && P.t == Test_3.t && P.z == Test_3.z)
        printf("да\n");
    else
        printf("нет\n");

    //тест 4
    mpz_t k1, k2;
    mpz_init(k1);
    mpz_init(k2);
    mpz_urandomm(k1, state, n);
    mpz_urandomm(k2, state, n);
    struct Point Right, Left;
    Right = multiple_point(P, k2);
    Left = multiple_point(P, k1);
    Left = add_points(Left, Right);
    mpz_add(k1, k1, k2);
    Right = multiple_point(P, k1);
    printf("4 тест. Верно ли, что [k1]P + [k2]P = [k1+k2]P: ");
    if (Left.x == Right.x && Left.y == Right.y && Left.t == Right.t && Left.z == Right.z)
        printf("да\n");
    else
        printf("нет\n");

    mpz_clear(x1);
    mpz_clear(y1);
    mpz_clear(t1);
    mpz_clear(z1);
    mpz_clear(k);
    mpz_clear(n);
    mpz_clear(qq);
    mpz_clear(k1);
    mpz_clear(k2);
}
