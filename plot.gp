set terminal svg size 800, 600  # svg形式で出力
set output "error_rate.svg" # 出力ファイル名
set title "Decoding Error vs. SNR"  # タイトル
set xlabel "SNR (dB)" # x軸ラベル
set ylabel "Error Rate (log_{10})" # y軸ラベル
set logscale y 10 # y軸をlogスケールにする
set mytics
set grid xtics ytics mytics  # グリッド表示
plot "error_rate.csv" using 1:2 with linespoints title "Error Rate"
