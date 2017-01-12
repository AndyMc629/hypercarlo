/*
 * mathematical functions not supplied by C++ libraries (or improvements of 
 * them).
 */
//returns correct modulus whether value is negative or not
//library supplied % does not work for negative numerator.
int mod(int x, int N) {
    return ((x %= N) < 0) ? x+N : x;
}