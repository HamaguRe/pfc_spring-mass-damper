# グラフ作成

import csv
import matplotlib.pyplot as plt

t = []  # 時刻
r = []  # 目標値
y = []  # 制御対象出力
u = []  # 制御入力
limit_under = []  # 制御入力上限
limit_upper = []  # 制御入力下限

# CSVからデータを読み出して配列に追加
with open('./result.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        nums = [float(v) for v in row]

        t.append(nums[0])
        r.append(nums[1])
        y.append(nums[2])
        u.append(nums[3])
        limit_under.append(nums[4])
        limit_upper.append(nums[5])


# --- 描画 --- #
fig1 = plt.figure(figsize = (10, 7))
plt.suptitle('Spring-Mass-Damper System',fontsize=20)

ax1 = fig1.add_subplot(211)
ax1.step(t, r, color="royalblue", label="Set-point", linestyle = "--")
ax1.step(t, y, color="limegreen", label="Plant output")
ax1.set_xlim(t[0], t[-1])
ax1.set_ylabel("Displacement [m]", fontsize=15)
ax1.tick_params(labelsize=13)  # 軸目盛の大きさ
ax1.legend(fontsize=15, loc="upper left")

ax2 = fig1.add_subplot(212)
ax2.step(t, limit_upper, color="red", linestyle="--", label="Limit")
ax2.step(t, limit_under, color="red", linestyle="--")
ax2.step(t, u, color="orange", label="Control input")
ax2.set_xlim(t[0], t[-1])
ax2.set_xlabel("Time [s]", fontsize=15)
ax2.set_ylabel("Force [N]", fontsize=15)
ax2.tick_params(labelsize=13)
ax2.legend(fontsize=15, loc="upper left")

plt.show()