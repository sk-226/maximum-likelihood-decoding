/**
 * @file 1.c
 * @author Suguru Kurita (g2221034@tcu.ac.jp)
 * @brief 1. SN比 0.0(dB) から 4.0(dB) まで、0.5(dB) おきに $\sigma$ を求めよ。
 * @version 0.1
 * @date 2024-12-18
 */
#include <math.h>
#include <stdio.h>

/**
 * @brief SN比から $\sigma$ を計算する関数
 *
 * @param snr
 * @return double
 */
double calculate_sigma(double snr) {
  double r = 1.0;
  return sqrt(pow(10, -snr / 10) / 2 * r);
}

int main(void) {
  double snr;      // SN比
  double sigma;    // ノイズの標準偏差
  double r = 1.0;  // 符号化率

  for (snr = 0.0; snr <= 4.0; snr += 0.5) {
    sigma = sqrt(pow(10, -snr / 10) / 2 * r);
    printf("snr = %lf, sigma = %lf\n", snr, sigma);
  }

  return 0;
}
