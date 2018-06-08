#ifndef __CONVERSION_H__

#define __CONVERSION_H__

/**
 * @brief Permet la convertion d'un tableau d'unsigned short en float
 *
 * @Param in: tableau a convertir
 *        cols, rows: dimension du tableau
 **/
float*
convertUnsignedShort2Float(unsigned short* in, int cols, int rows);

/**
 * @brief Permet la convertion d'un tableau de float en unsigned short
 *
 * @Param in: tableau a convertir
 *        out: tableau converti
 *        cols, rows: dimension du tableau in
 **/
void
convertFloat2UnsignedShort(float* in, unsigned short* out, int cols, int rows);

/**
 * @brief Calcule la mutiplication entre une matrice3X3 de convertion et un vecteur des trois cannaux
 *
 * @Param cols, rows: dimension du tableau
 *        A: matrice de convertion
 *        B: tableau des trois cannaux
 *        output: résultat de la multiplication
 **/
void
multiple(int cols, int rows, float A[3][3], float* B, float* output);

/**
 * @brief Converti un tableau RGB en YIQ
 *
 * @Param cols, rows: dimension du tableau
 *        in: tableau a convertir
 **/
float*
convertRGB2YIQ(int cols, int rows, float* in);

/**
 * @brief Converti un tableau YIQ en RGB
 *
 * @Param cols, rows: dimension du tableau
 *        in: tableau a convertir
 **/
float*
convertYIQ2RGB(int cols, int rows, float* in);

/**
 * @brief Sort le channel Y, I, Q choisit
 *
 * @Param cols, rows: dimension du tableau
 *        yiq: tableau contenant les channel yiq
 *        channel: 0,1,2 pour récupérer  le channel y,i,q
 **/
float*
channel(int cols, int rows, float* yiq, int channel);

/**
 * @brief Met dans un tableau les channels y,i,q
 *
 * @Param cols, rows: dimension du tableau
 *        y, i, q: Cannaux à mettre dans le tableau
 **/
float*
putChannel(int cols, int rows, float* y, float* I, float* q);

#endif
