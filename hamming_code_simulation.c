/* hamming_code_simulation.c */
#include <math.h>
#include <stdio.h>

/* Function prototypes */
double uniform_random_number(void);
void normal_gaussian_noise(double n_random[]);
void receive(double sigma, double y[]);
void calculate_log_likelihood(double sigma, double y[],
                              double log_likelihood[]);
int select_max_likelihood_codeword(double log_likelihood[]);  // which_the_1st

int main(void) {
  double sigma = 1.0;
  double y[8];
  double log_likelihood[16];

  /* 受信系列を3回発生し、最尤復号法を3回行う */
  for (int event = 0; event < 3; event++) {
    // 受信系列を受信
    receive(sigma, y);
    printf("受信系列: \n");
    for (int i = 0; i < 8; i++) {
      printf("%f ", y[i]);
    }
    printf("\n");

    // すべての復号語に対して対数尤度を計算
    calculate_log_likelihood(sigma, y, log_likelihood);
    printf("符号語の対数尤度: \n");
    for (int i = 0; i < 16; i++) {
      printf(" %d %f \n", i, log_likelihood[i]);
    }

    printf("最大尤度の符号語は %d 番目 \n",
           select_max_likelihood_codeword(log_likelihood));
  }

  return 0;
}

/**
 * @brief 一様乱数を発生させる関数
 *
 * @return double
 */
double uniform_random_number(void) {
  static int seed = 17;
  double yy;
  seed = (seed * 23) % 10000001;
  yy = seed * 1.0e-07;
  return yy;
}

/**
 * @brief 平均0, 分散1 のガウス乱数 (正規乱数)
 *
 * @param n_random
 */
void normal_gaussian_noise(double n_random[]) {
  double eta1, eta2, v1, v2, k;
  do {
    eta1 = uniform_random_number();
    eta2 = uniform_random_number();

    v1 = -1.0 + eta1 * 2.0;
    v2 = -1.0 + eta2 * 2.0;
    k = v1 * v1 + v2 * v2;
  } while (k >= 1.0);
  n_random[0] = v1 * sqrt(-2.0 * log(k) / k);
  n_random[1] = n_random[0] * (v2 / v1);
}

/**
 * @brief 7ビットの受信系列を生成して y[] に格納
 *
 * @param sigma
 * @param y
 */
void receive(double sigma, double y[]) {
  double n_random[2];
  for (int i = 0; i < 4; i++) {
    normal_gaussian_noise(n_random);
    y[2 * i] = 1.0 + n_random[0] * sigma;
    if (2 * i + 1 < 8) {
      y[2 * i + 1] = 1.0 + n_random[1] * sigma;
    }
  }
}

/**
 * @brief sn比σ, 平均1, 分散σ^2 のガウス分布に従う受信系列 y[]
 * に対して、対数尤度を計算
 *
 * @param sigma
 * @param y
 * @param log_likelihood
 */
void calculate_log_likelihood(double sigma, double y[],
                              double log_likelihood[]) {
  double bb;
  int s[16][7];
  int c[16][7] = {
      {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 1, 1}, {0, 0, 1, 0, 1, 1, 0},
      {0, 0, 1, 0, 1, 1, 0}, {0, 1, 0, 0, 1, 1, 1}, {0, 1, 1, 1, 1, 0, 0},
      {0, 1, 1, 0, 0, 0, 1}, {0, 1, 1, 1, 0, 1, 0}, {1, 0, 0, 0, 1, 0, 1},
      {1, 0, 0, 1, 1, 1, 1}, {1, 0, 1, 0, 0, 1, 1}, {1, 0, 1, 1, 0, 0, 1},
      {1, 1, 0, 0, 0, 1, 0}, {1, 1, 0, 1, 0, 0, 1}, {1, 1, 1, 0, 1, 0, 0},
      {1, 1, 1, 1, 1, 1, 0},
  };

  bb = log(1 / sqrt(2 * M_PI * sigma * sigma));

  /* バイナリ・バイポーラ変換 */
  for (int codeword = 0; codeword < 16; codeword++) {
    for (int i = 0; i < 7; i++) {
      if (c[codeword][i] == 0) {
        s[codeword][i] = 1;
      } else {
        s[codeword][i] = -1;
      }
    }
  }

  /* 16個の符号語に対して, 対数尤度の計算 */
  for (int codeword = 0; codeword < 16; codeword++) {
    log_likelihood[codeword] = 0.0;

    for (int i = 0; i < 7; i++) {
      log_likelihood[codeword] =
          log_likelihood[codeword] + bb -
          pow(y[i] - s[codeword][i], 2) / (2 * sigma * sigma);
    }
  }
}

/**
 * @brief 最大尤度の符号語のインデックスを返す
 * 
 * @param log_likelihood 
 * @return int imax 
 */
int select_max_likelihood_codeword(double log_likelihood[]) {
  int imax = 0;

  for (int i = 1; i < 16; i++) {
    if (log_likelihood[i] > log_likelihood[imax]) {
      imax = i;
    }
  }

  return imax;
}
