# Программная реализация алгоритма поиска кратной точки эллиптической кривой в форме пересечения Якоби (с использованием gmp)
Необходимо
1. Построить/выбрать точку P на кривой
2. Выбрать случайное значение k.
3. Реализовать операцию вычисления кратной точки Q = [k]P.
4. Провести тестирование программы.

Для проведения тестирования необходимо:
1. Проверить, что результирующая точка Q лежит на кривой.
2. Проверить, что [q]P = O, где q – порядок группы точек.
3. Проверить, что [q + 1]P = P и [q − 1]P = −P.
4. Для двух случайных k1, k2 проверить, что [k1]P + [k2]P = [k1 + k2]P.

 
 OS name: Ubuntu 20.04.1 LTS
 
 OS Type: 64 bit
