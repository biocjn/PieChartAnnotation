import matplotlib.pyplot as plt
import numpy as np

# 1. 准备数据
labels = ['苹果', '香蕉', '橙子', '葡萄', '芒果']
sizes = [15, 30, 25, 10, 20]  # 各部分大小
colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99', '#c2c2f0']  # 颜色

# 2. 创建图形和坐标轴
fig, ax = plt.subplots(figsize=(10, 8))
fig.set_facecolor('white')

# 3. 绘制饼图（不显示默认标签）
wedges, texts, autotexts = ax.pie(
    sizes,
    colors=colors,
    startangle=90,
    wedgeprops={'edgecolor': 'black', 'linewidth': 1},
    autopct='%1.1f%%',
    pctdistance=0.85  # 百分比标签离中心的距离
)

# 4. 设置百分比标签的样式
plt.setp(autotexts, size=10, weight='bold', color='black')

# 5. 添加带引导线的外部标签
for i, (wedge, label) in enumerate(zip(wedges, labels)):
    # 计算引导线起点（扇形边缘）
    ang = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1  # 扇形中心角度
    x = np.cos(np.deg2rad(ang))  # 角度转弧度
    y = np.sin(np.deg2rad(ang))
    
    # 确定引导线终点位置（根据方向调整位置）
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = f"angle,angleA=0,angleB={ang}"
    
    # 计算标签位置（饼图外部）
    label_distance = 1.3  # 标签距离中心的倍数
    label_x = label_distance * x
    label_y = label_distance * y
    
    # 添加带引导线的注释
    ax.annotate(
        f"{label}: {sizes[i]} ({sizes[i]/sum(sizes):.1%})",  # 标签文本
        xy=(x * 0.8, y * 0.8),  # 引导线起点（扇形边缘内）
        xytext=(label_x, label_y),  # 标签位置
        horizontalalignment=horizontalalignment,
        verticalalignment='center',
        arrowprops=dict(
            arrowstyle="-",  # 直线箭头
            color="black",
            connectionstyle=connectionstyle,
            linewidth=1
        ),
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor=colors[i],
            alpha=0.7,
            edgecolor='black'
        )
    )

# 6. 设置标题并显示
plt.title('水果分布比例', fontsize=14, pad=20)
plt.tight_layout()
plt.show()