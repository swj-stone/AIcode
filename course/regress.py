# Feature 数据
X = [[10, 3], [20, 3], [25, 3], [28, 2.5], [30, 2], [35, 2.5], [40, 2.5]]
y = [60, 85, 100, 120, 140, 145, 163]  # Label 数据
# 初始化参数
w = [0.0, 0.0, 0.0]  # w0, w1, w2
lr = 0.0001  # 学习率
num_iterations = 10000  # 迭代次数
# 梯度下降
for i in range(num_iterations):
    # 预测值
    y_pred = [w[0] + w[1] * x[0] + w[2] * x[1] for x in X]
    # 计算损失
    loss = sum((y_pred[j] - y[j]) ** 2 for j in range(len(y))) / len(y)
    # 计算梯度
    grad_w0 = 2 * sum(y_pred[j] - y[j] for j in range(len(y))) / len(y)
    grad_w1 = 2 * sum((y_pred[j] - y[j]) * X[j][0] for j in range(len(y))) / len(y)
    grad_w2 = 2 * sum((y_pred[j] - y[j]) * X[j][1] for j in range(len(y))) / len(y)
    # 更新参数
    w[0] -= lr * grad_w0
    w[1] -= lr * grad_w1
    w[2] -= lr * grad_w2
    # 打印损失
    if i % 100 == 0:
        print(f"Iteration {i}: Loss = {loss}")
# 输出最终参数
print(f"Final parameters: w0 = {w[0]}, w1 = {w[1]}, w2 = {w[2]}")


'''如果模型有多个特征，还可以构造出多个特征之间相乘的高次项。比如你要预测房价，
你收集的feature里有房子的长度，有房子的宽度。那么你可以构造一个新的feature，它的值等于长度乘以宽度。也就是房子的面积。
通过本节的学习，我们知道对于特征和Label之间的非线性问题，我们可以通过构造高次特征来解决。一般的做法是先从二次项开始，
逐步增加，直到达到我们满意的效果。假如你在构造特征的二次项，需要注意的是，构造的特征不光可以是一个特征的平方，也可以是任意两个特征之间的乘积。'''
#特征的排列组合，或者是任意次的多项式组合（Taylar展开和有限区间多项式拟合）