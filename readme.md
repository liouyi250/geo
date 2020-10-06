## 初等变换法求矩阵逆矩阵
1. 写出行增广矩阵 MP
2. 使用初等变换对第一列进行消元，使得除去主元位置外，其余位置都为零（如果还没有消元前，发现主元位置处的数值为零，需要进行换行操作）
3. 循环对2到n列进行消元，如果在某一列主元位置元素为零，且其下方的列元素全为零，则该矩阵没有逆矩阵。
4. 将M的对角矩阵变为单位矩阵
5. 返回其逆矩阵

<hr>
![Image](https://github.com/liouyi250/geo/blob/master/pic/1.PNG)
![Image](https://github.com/liouyi250/geo/blob/master/pic/2.PNG)

