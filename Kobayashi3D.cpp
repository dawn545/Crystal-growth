#include "Kobayashi3D.h"
#include <cfloat> // 包含 FLT_EPSILON，用于浮点数比较，防止除以零

// ==========================================
// 构造函数与初始化
// ==========================================

Kobayashi::Kobayashi(int x, int y, int z, float timeStep) {
    _objectCount = { x, y, z }; // 3D 网格大小，例如 100x100x100
    _dx = 0.03f; // x 方向空间步长
    _dy = 0.03f; // y 方向空间步长
    _dz = 0.03f; // z 方向空间步长
    _dt = timeStep; // 时间步长：每次模拟迭代推进的时间量
    _renderSlice = z / 2; // 默认渲染中间切片

    _initParams();  // 初始化物理参数
    _vectorInit();  // 分配内存并设置初始条件
}

Kobayashi::~Kobayashi() {
    // 析构函数：程序退出时清理显存中的纹理资源
    if (_textureID) glDeleteTextures(1, &_textureID);
}

// 初始化 Kobayashi 晶体生长模型的物理常数
// 这些参数决定了晶体的形状（六角形、树枝状等）和生长速度
// 参考论文：Kobayashi, A. (1993) "Modeling and numerical simulations of dendritic crystal growth"
void Kobayashi::_initParams() {
    // 弛豫时间 τ，控制相变速度（单位：无量纲）
    // 物理意义：越小相变越快，晶体生长越迅速
    _tau = 0.0003f;

    // 平均各向异性强度 ε̄，决定晶界厚度（单位：无量纲）
    // 用于各向异性系数的基础值：ε(θ) = ε̄(1 + δcos(nθ))
    _epsilonBar = 0.010f;

    // 恒相位移动系数 M_η，控制相场演化速度（单位：无量纲）
    // 物理意义：M_η = 1/τ，τ 为弛豫时间，M_η 越大相变越快
    M_eta = 1/_tau;

    // 潜热系数 K，控制温度场对相变的影响（单位：无量纲）
    // 物理意义：相变释放的潜热量，影响周围温度场
    _K = 1.6f;

    // 各向异性强度 δ，值越大晶体越有棱角（范围：0.0-1.0）
    // 公式中的 δ 参数，控制各向异性的幅度
    _delta = 0.05f;

    // 折叠对称性参数 j_fold，决定晶体的对称性阶数
    // 物理意义：晶体在旋转 2π/j_fold 角度后会重复自身的形状
    // 常用值：
    //   j_fold = 4.0 → 四重对称（正方形晶体，如某些金属）
    //   j_fold = 6.0 → 六重对称（六角形晶体，如雪花、石墨）
    //   j_fold = 8.0 → 八重对称（八角形晶体）
    // 在各向异性系数公式中：ε(θ) = ε̄(1 + δ·cos(j_fold·(θ₀ - θ)))
    j_fold = 6.0f;

    // 过冷度系数 α，用于计算驱动力（单位：无量纲）
    // 在 m = (α/π)arctan(γ(T_eq - T)) 中使用
    _alpha = 0.9f;

    // 过冷度灵敏度 γ，控制温度对驱动力的影响程度
    // 值越大，温度变化对生长速度的影响越敏感
    _gamma = 10.0f;

    // 平衡温度 T_eq，晶体与液体的相平衡温度（单位：无量纲）
    _tEq = 1.0f;

    // 热扩散系数 a²，控制温度场的扩散速度（单位：无量纲）
    // 物理意义：热量在材料中传播的快慢
    // 典型值：0.5 - 2.0，值越大温度扩散越快
    _alpha_T = 1.0f;

    // 取向场驱动力系数 H，控制取向场对晶体生长的影响强度
    // 物理意义：取向场梯度越大，对相场演化的驱动力越强
    _H = 0.0f; // 默认为 0，可根据需要调整

    // 取向场迁移率 M_ori，控制取向场演化的速度
    // 物理意义：值越大，取向场调整越快
    M_ori = 1.0f;
}

// 分配内存并重置模拟状态
void Kobayashi::_vectorInit() {
    size_t vSize = _objectCount.x * _objectCount.y * _objectCount.z;

    // _phi: 相场变量 (0=液, 1=固)
    // _t: 温度场
    _phi.assign(vSize, 0.0f);
    _t.assign(vSize, 0.0f);

    // 相场梯度和拉普拉斯算子
    _gradPhiX.assign(vSize, 0.0f);
    _gradPhiY.assign(vSize, 0.0f);
    _gradPhiZ.assign(vSize, 0.0f);
    _lapPhi.assign(vSize, 0.0f);
    _lapT.assign(vSize, 0.0f);

    // 相场梯度模和局部方向角
    _gradPhiMag.assign(vSize, 0.0f);
    _tau_field.assign(vSize, 0.0f);
    _theta.assign(vSize, 0.0f);
    _phi_angle.assign(vSize, 0.0f);

    // 各向异性系数及其导数
    _epsilon.assign(vSize, 0.0f);
    _epsilonDerivTheta.assign(vSize, 0.0f);
    _epsilonDerivPhi.assign(vSize, 0.0f);

    // 取向场：Ω_ori 用单位球上的点 (x, y, z) 表示
    // 初始化为指向 z 轴正方向 (0, 0, 1)
    _omega_ori_x.assign(vSize, 0.0f);
    _omega_ori_y.assign(vSize, 0.0f);
    _omega_ori_z.assign(vSize, 1.0f);

    // 取向场的局部极坐标表示
    _rho_x_plus.assign(vSize, 0.0f);
    _rho_x_minus.assign(vSize, 0.0f);
    _rho_y_plus.assign(vSize, 0.0f);
    _rho_y_minus.assign(vSize, 0.0f);
    _rho_z_plus.assign(vSize, 0.0f);
    _rho_z_minus.assign(vSize, 0.0f);

    _lambda_x_plus.assign(vSize, 0.0f);
    _lambda_x_minus.assign(vSize, 0.0f);
    _lambda_y_plus.assign(vSize, 0.0f);
    _lambda_y_minus.assign(vSize, 0.0f);
    _lambda_z_plus.assign(vSize, 0.0f);
    _lambda_z_minus.assign(vSize, 0.0f);

    // 取向场梯度模
    _gradOmegaOriMag.assign(vSize, 0.0f);

    // 像素缓冲区：用于渲染切片，每个像素4个字节 (R,G,B,A)
    _pixelBuffer.assign(_objectCount.x * _objectCount.y * 4, 0);

    // 在中心创建一个初始晶核
    _createNucleus(_objectCount.x / 2, _objectCount.y / 2, _objectCount.z / 2);

    // 立即更新一次纹理，确保初始画面不是黑的
    _updateTexture();
}

// ==========================================
// 辅助函数：计算取向场相关的几何量
// ==========================================

// 计算两个单位向量之间的大圆距离（中心角）
// 输入：两个单位向量 ω_p = (x_p, y_p, z_p) 和 ω_q = (x_q, y_q, z_q)
// 输出：中心角 ρ ∈ [0, π]
inline float centralAngle(float x_p, float y_p, float z_p, float x_q, float y_q, float z_q) {
    // 使用点积计算：cos(ρ) = ω_p · ω_q
    float dot = x_p * x_q + y_p * y_q + z_p * z_q;
    // 限制在 [-1, 1] 范围内以避免数值误差
    dot = fmax(-1.0f, fmin(1.0f, dot));
    return acos(dot);
}

// 立体投影：将单位球上的点 ω_q 投影到以 ω_p 为极点的平面上
// 然后计算投影点在局部极坐标中的角度 λ
// 输入：ω_p (极点), ω_q (要投影的点)
// 输出：角度 λ ∈ [0, 2π)
inline float stereographicAngle(float x_p, float y_p, float z_p, float x_q, float y_q, float z_q) {
    // 建立局部坐标系：
    // 极点为 ω_p，需要找到两个正交的切向量作为基

    // 第一个切向量：选择一个与 ω_p 不平行的向量，然后叉乘
    float ref_x = 0.0f, ref_y = 0.0f, ref_z = 1.0f;
    if (fabs(z_p) > 0.9f) {
        // 如果 ω_p 接近 z 轴，使用 x 轴作为参考
        ref_x = 1.0f; ref_y = 0.0f; ref_z = 0.0f;
    }

    // 第一个切向量：e1 = ref × ω_p（归一化）
    float e1_x = ref_y * z_p - ref_z * y_p;
    float e1_y = ref_z * x_p - ref_x * z_p;
    float e1_z = ref_x * y_p - ref_y * x_p;
    float e1_len = sqrt(e1_x * e1_x + e1_y * e1_y + e1_z * e1_z);
    if (e1_len < FLT_EPSILON) return 0.0f;
    e1_x /= e1_len; e1_y /= e1_len; e1_z /= e1_len;

    // 第二个切向量：e2 = ω_p × e1（已经归一化）
    float e2_x = y_p * e1_z - z_p * e1_y;
    float e2_y = z_p * e1_x - x_p * e1_z;
    float e2_z = x_p * e1_y - y_p * e1_x;

    // 计算从 ω_p 到 ω_q 的向量在切平面上的投影
    // 投影向量 = ω_q - (ω_q · ω_p) * ω_p
    float dot_pq = x_p * x_q + y_p * y_q + z_p * z_q;
    float proj_x = x_q - dot_pq * x_p;
    float proj_y = y_q - dot_pq * y_p;
    float proj_z = z_q - dot_pq * z_p;

    // 将投影向量表示在 (e1, e2) 基下
    float coord_e1 = proj_x * e1_x + proj_y * e1_y + proj_z * e1_z;
    float coord_e2 = proj_x * e2_x + proj_y * e2_y + proj_z * e2_z;

    // 计算极坐标角度
    float lambda = atan2(coord_e2, coord_e1);
    if (lambda < 0.0f) lambda += 2.0f * PI_F;

    return lambda;
}

// 在网格中心放置一个微小的”种子”，让晶体开始生长
void Kobayashi::_createNucleus(int x, int y, int z)
{
    // 在3D空间中创建一个小球形晶核
    // 将中心及周围的点设为 1.0 (固体)
    _phi[_INDEX(x, y, z)] = 1.0f;
    _phi[_INDEX(x - 1, y, z)] = 1.0f;
    _phi[_INDEX(x + 1, y, z)] = 1.0f;
    _phi[_INDEX(x, y - 1, z)] = 1.0f;
    _phi[_INDEX(x, y + 1, z)] = 1.0f;
    _phi[_INDEX(x, y, z - 1)] = 1.0f;
    _phi[_INDEX(x, y, z + 1)] = 1.0f;
}

// ==========================================
// 物理模拟核心：计算导数
// ==========================================

// 计算空间导数（梯度和拉普拉斯算子）- 3D版本
// 这是有限差分法的核心：通过邻居格子的值来推算当前的斜率和曲率
// 参考：有限差分法 (Finite Difference Method, FDM)
void Kobayashi::_computeGradientLaplacian()
{
    for (int k = 0; k < _objectCount.z; k++)
    {
        for (int j = 0; j < _objectCount.y; j++)
        {
            for (int i = 0; i < _objectCount.x; i++)
            {
                // 周期性边界条件 (Periodic Boundary Condition)
                int i_plus = (i + 1) % _objectCount.x;
                int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
                int j_plus = (j + 1) % _objectCount.y;
                int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;
                int k_plus = (k + 1) % _objectCount.z;
                int k_minus = ((k - 1) + _objectCount.z) % _objectCount.z;

                int idx = _INDEX(i, j, k);

                // ========== 1. 计算相场梯度 (Gradient / 梯度) ==========
                // 物理意义：相场变化的”坡度”，用于确定界面法线方向
                // 公式：∂η/∂x ≈ (η(i+1,j,k) - η(i-1,j,k)) / (2·Δx)  [中心差分法]
                _gradPhiX[idx] = (_phi[_INDEX(i_plus, j, k)] - _phi[_INDEX(i_minus, j, k)]) / (2.0f * _dx);
                _gradPhiY[idx] = (_phi[_INDEX(i, j_plus, k)] - _phi[_INDEX(i, j_minus, k)]) / (2.0f * _dy);
                _gradPhiZ[idx] = (_phi[_INDEX(i, j, k_plus)] - _phi[_INDEX(i, j, k_minus)]) / (2.0f * _dz);

                // 计算梯度模：|∇η| = sqrt((∂η/∂x)² + (∂η/∂y)² + (∂η/∂z)²)
                _gradPhiMag[idx] = sqrt(_gradPhiX[idx] * _gradPhiX[idx]
                                      + _gradPhiY[idx] * _gradPhiY[idx]
                                      + _gradPhiZ[idx] * _gradPhiZ[idx]);

                // 计算 τ = sqrt((∂η/∂x)² + (∂η/∂y)²)
                // 注意：τ 只包含 x 和 y 方向的梯度，不包含 z 方向
                _tau_field[idx] = sqrt(_gradPhiX[idx] * _gradPhiX[idx]
                                     + _gradPhiY[idx] * _gradPhiY[idx]);

                // ========== 2. 计算局部相位前沿方向 Ω = (θ, φ) ==========
                // 物理意义：界面法线方向的球坐标表示
                // 公式：θ = -arccos((∂η/∂z)/|∇η|)
                // 公式：φ = -arctan((∂η/∂y)/(∂η/∂x))
                if (_gradPhiMag[idx] > FLT_EPSILON) {
                    _theta[idx] = -acos(_gradPhiZ[idx] / _gradPhiMag[idx]);
                } else {
                    _theta[idx] = 0.0f;
                }

                if (fabs(_gradPhiX[idx]) > FLT_EPSILON) {
                    _phi_angle[idx] = -atan2(_gradPhiY[idx], _gradPhiX[idx]);
                } else {
                    _phi_angle[idx] = 0.0f;
                }

                // ========== 3. 计算拉普拉斯算子 (Laplacian) - 3D 7点模板 ==========
                // 物理意义：场的”扩散”趋势
                // 公式：∇²η = ∂²η/∂x² + ∂²η/∂y² + ∂²η/∂z²
                // 使用7点中心差分：∇²η ≈ (η_E + η_W + η_N + η_S + η_U + η_D - 6η_C) / Δx²
                _lapPhi[idx] = (_phi[_INDEX(i_plus, j, k)] + _phi[_INDEX(i_minus, j, k)]
                              + _phi[_INDEX(i, j_plus, k)] + _phi[_INDEX(i, j_minus, k)]
                              + _phi[_INDEX(i, j, k_plus)] + _phi[_INDEX(i, j, k_minus)]
                              - 6.0f * _phi[idx]) / (_dx * _dx);

                // 同样计算温度场的拉普拉斯算子
                _lapT[idx] = (_t[_INDEX(i_plus, j, k)] + _t[_INDEX(i_minus, j, k)]
                            + _t[_INDEX(i, j_plus, k)] + _t[_INDEX(i, j_minus, k)]
                            + _t[_INDEX(i, j, k_plus)] + _t[_INDEX(i, j, k_minus)]
                            - 6.0f * _t[idx]) / (_dx * _dx);

                // ========== 4. 算法1：计算取向场梯度 ∇Ω_ori ==========
                // 使用立体投影和局部极坐标方法
                // 参考：Algorithm 1 - Calculation of ∇Ω_ori

                // 当前点的取向 ω_p
                float omega_p_x = _omega_ori_x[idx];
                float omega_p_y = _omega_ori_y[idx];
                float omega_p_z = _omega_ori_z[idx];

                // 对每个邻居计算 (ρ, λ)
                // x+ 方向
                {
                    int idx_q = _INDEX(i_plus, j, k);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_x_plus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_x_plus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // x- 方向
                {
                    int idx_q = _INDEX(i_minus, j, k);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_x_minus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_x_minus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // y+ 方向
                {
                    int idx_q = _INDEX(i, j_plus, k);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_y_plus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_y_plus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // y- 方向
                {
                    int idx_q = _INDEX(i, j_minus, k);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_y_minus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_y_minus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // z+ 方向
                {
                    int idx_q = _INDEX(i, j, k_plus);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_z_plus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_z_plus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // z- 方向
                {
                    int idx_q = _INDEX(i, j, k_minus);
                    float omega_q_x = _omega_ori_x[idx_q];
                    float omega_q_y = _omega_ori_y[idx_q];
                    float omega_q_z = _omega_ori_z[idx_q];

                    _rho_z_minus[idx] = centralAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                    _lambda_z_minus[idx] = stereographicAngle(omega_p_x, omega_p_y, omega_p_z, omega_q_x, omega_q_y, omega_q_z);
                }

                // 计算 ∇Ω_ori 的模
                // 使用 (ρ, λ) 场的梯度来近似
                // ||∇Ω_ori|| ≈ sqrt((∂ρ/∂x)² + (∂ρ/∂y)² + (∂ρ/∂z)² + (∂λ/∂x)² + (∂λ/∂y)² + (∂λ/∂z)²)
                float grad_rho_x = (_rho_x_plus[idx] - _rho_x_minus[idx]) / (2.0f * _dx);
                float grad_rho_y = (_rho_y_plus[idx] - _rho_y_minus[idx]) / (2.0f * _dy);
                float grad_rho_z = (_rho_z_plus[idx] - _rho_z_minus[idx]) / (2.0f * _dz);

                float grad_lambda_x = (_lambda_x_plus[idx] - _lambda_x_minus[idx]) / (2.0f * _dx);
                float grad_lambda_y = (_lambda_y_plus[idx] - _lambda_y_minus[idx]) / (2.0f * _dy);
                float grad_lambda_z = (_lambda_z_plus[idx] - _lambda_z_minus[idx]) / (2.0f * _dz);

                _gradOmegaOriMag[idx] = sqrt(grad_rho_x * grad_rho_x + grad_rho_y * grad_rho_y + grad_rho_z * grad_rho_z
                                           + grad_lambda_x * grad_lambda_x + grad_lambda_y * grad_lambda_y + grad_lambda_z * grad_lambda_z);

                // ========== 5. 计算各向异性系数 ε(Ω, Ω_ori) 及其导数 ==========
                // 物理意义：晶体在不同方向生长速度不同
                // 计算局部相位前沿方向 Ω 和取向场 Ω_ori 之间的夹角

                // Ω = -∇η (归一化)
                float omega_x = 0.0f, omega_y = 0.0f, omega_z = 0.0f;
                if (_gradPhiMag[idx] > FLT_EPSILON) {
                    omega_x = -_gradPhiX[idx] / _gradPhiMag[idx];
                    omega_y = -_gradPhiY[idx] / _gradPhiMag[idx];
                    omega_z = -_gradPhiZ[idx] / _gradPhiMag[idx];
                }

                // 计算 Ω 和 Ω_ori 之间的夹角
                float dot_omega = omega_x * omega_p_x + omega_y * omega_p_y + omega_z * omega_p_z;
                dot_omega = fmax(-1.0f, fmin(1.0f, dot_omega));
                float angle_diff = acos(dot_omega);

                // 各向异性系数：ε = ε̄(1 + δ·cos(j_fold·angle_diff))
                _epsilon[idx] = _epsilonBar * (1.0f + _delta * cos(j_fold * angle_diff));

                // 各向异性系数对角度的导数：∂ε/∂angle = -ε̄·j_fold·δ·sin(j_fold·angle_diff)
                _epsilonDerivTheta[idx] = -_epsilonBar * j_fold * _delta * sin(j_fold * angle_diff);

                // 对 φ 的导数（简化为0）
                _epsilonDerivPhi[idx] = 0.0f;
            }
        }
    }
}

// ==========================================
// 物理模拟核心：时间演化 - 3D版本
// ==========================================

// 根据微分方程更新 _phi, _t, _theta_ori, _phi_ori 的值
// 实现公式(17)和(18)
// 参考：Kobayashi, A. (1993) "Modeling and numerical simulations of dendritic crystal growth"
// ==========================================
// 物理模拟核心：时间演化 - 3D版本
// ==========================================

// 根据微分方程更新 _phi, _t, _theta_ori, _phi_ori 的值
// 实现公式(17)和(18)
void Kobayashi::_evolution()
{
    for (int k = 0; k < _objectCount.z; k++)
    {
        for (int j = 0; j < _objectCount.y; j++)
        {
            for (int i = 0; i < _objectCount.x; i++)
            {
                // 周期性边界条件
                int i_plus = (i + 1) % _objectCount.x;
                int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
                int j_plus = (j + 1) % _objectCount.y;
                int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;
                int k_plus = (k + 1) % _objectCount.z;
                int k_minus = ((k - 1) + _objectCount.z) % _objectCount.z;

                int idx = _INDEX(i, j, k);

                // ========== 计算驱动力 (Driving Force) ==========
                // 驱动力由过冷度（T_eq - T）决定
                // 公式：m = (α/π)·arctan(γ·(T_eq - T))
                float m = _alpha / PI_F * atan(_gamma * (_tEq - _t[idx]));

                // 保存旧值
                float oldPhi = _phi[idx];
                float oldT = _t[idx];
                float oldOmegaX = _omega_ori_x[idx];
                float oldOmegaY = _omega_ori_y[idx];
                float oldOmegaZ = _omega_ori_z[idx];

                // 获取当前点的各向异性参数
                float eps = _epsilon[idx];
                float eps_theta = _epsilonDerivTheta[idx];
                float eps_phi = _epsilonDerivPhi[idx];
                float tau = _tau_field[idx];
                float gradPhiMag = _gradPhiMag[idx];

                // ========== 公式(17)：相场演化方程 ==========
                // ∂η/∂t = M_η[∇·(ε²∇η) + ∂/∂z(...) + ∂/∂y(...) - ∂/∂z(ε·∂ε/∂θ·τ) - g'(η) - p'(η)(f_s - f_t + f_ori)]

                // 第一项：∇·(ε²∇η) = ε²∇²η + ∇(ε²)·∇η
                float term_diffusion = eps * eps * _lapPhi[idx];

                // ∇(ε²)·∇η
                float gradEps2_x = (_epsilon[_INDEX(i_plus, j, k)] * _epsilon[_INDEX(i_plus, j, k)]
                                  - _epsilon[_INDEX(i_minus, j, k)] * _epsilon[_INDEX(i_minus, j, k)]) / (2.0f * _dx);
                float gradEps2_y = (_epsilon[_INDEX(i, j_plus, k)] * _epsilon[_INDEX(i, j_plus, k)]
                                  - _epsilon[_INDEX(i, j_minus, k)] * _epsilon[_INDEX(i, j_minus, k)]) / (2.0f * _dy);
                float gradEps2_z = (_epsilon[_INDEX(i, j, k_plus)] * _epsilon[_INDEX(i, j, k_plus)]
                                  - _epsilon[_INDEX(i, j, k_minus)] * _epsilon[_INDEX(i, j, k_minus)]) / (2.0f * _dz);

                float term_grad_eps2 = gradEps2_x * _gradPhiX[idx]
                                     + gradEps2_y * _gradPhiY[idx]
                                     + gradEps2_z * _gradPhiZ[idx];

                // 第二项：∂/∂z[ε/τ·∂ε/∂θ·∂η/∂x - ε/τ²·∂ε/∂φ·|∇η|²·∂η/∂y]
                // 计算 z+ 和 z- 处的值
                float val_z_plus, val_z_minus;
                {
                    int idx_zp = _INDEX(i, j, k_plus);
                    float eps_zp = _epsilon[idx_zp];
                    float tau_zp = _tau_field[idx_zp];
                    float eps_theta_zp = _epsilonDerivTheta[idx_zp];
                    float eps_phi_zp = _epsilonDerivPhi[idx_zp];
                    float gradPhiMag2_zp = _gradPhiMag[idx_zp] * _gradPhiMag[idx_zp];

                    val_z_plus = (eps_zp / (tau_zp + FLT_EPSILON)) * eps_theta_zp * _gradPhiX[idx_zp]
                               - (eps_zp / ((tau_zp * tau_zp) + FLT_EPSILON)) * eps_phi_zp * gradPhiMag2_zp * _gradPhiY[idx_zp];
                }
                {
                    int idx_zm = _INDEX(i, j, k_minus);
                    float eps_zm = _epsilon[idx_zm];
                    float tau_zm = _tau_field[idx_zm];
                    float eps_theta_zm = _epsilonDerivTheta[idx_zm];
                    float eps_phi_zm = _epsilonDerivPhi[idx_zm];
                    float gradPhiMag2_zm = _gradPhiMag[idx_zm] * _gradPhiMag[idx_zm];

                    val_z_minus = (eps_zm / (tau_zm + FLT_EPSILON)) * eps_theta_zm * _gradPhiX[idx_zm]
                                - (eps_zm / ((tau_zm * tau_zm) + FLT_EPSILON)) * eps_phi_zm * gradPhiMag2_zm * _gradPhiY[idx_zm];
                }
                float term_z = (val_z_plus - val_z_minus) / (2.0f * _dz);

                // 第三项：∂/∂y[ε/τ·∂ε/∂θ·∂η/∂z + ε/τ²·∂ε/∂φ·|∇η|²·∂η/∂x]
                float val_y_plus, val_y_minus;
                {
                    int idx_yp = _INDEX(i, j_plus, k);
                    float eps_yp = _epsilon[idx_yp];
                    float tau_yp = _tau_field[idx_yp];
                    float eps_theta_yp = _epsilonDerivTheta[idx_yp];
                    float eps_phi_yp = _epsilonDerivPhi[idx_yp];
                    float gradPhiMag2_yp = _gradPhiMag[idx_yp] * _gradPhiMag[idx_yp];

                    val_y_plus = (eps_yp / (tau_yp + FLT_EPSILON)) * eps_theta_yp * _gradPhiZ[idx_yp]
                               + (eps_yp / ((tau_yp * tau_yp) + FLT_EPSILON)) * eps_phi_yp * gradPhiMag2_yp * _gradPhiX[idx_yp];
                }
                {
                    int idx_ym = _INDEX(i, j_minus, k);
                    float eps_ym = _epsilon[idx_ym];
                    float tau_ym = _tau_field[idx_ym];
                    float eps_theta_ym = _epsilonDerivTheta[idx_ym];
                    float eps_phi_ym = _epsilonDerivPhi[idx_ym];
                    float gradPhiMag2_ym = _gradPhiMag[idx_ym] * _gradPhiMag[idx_ym];

                    val_y_minus = (eps_ym / (tau_ym + FLT_EPSILON)) * eps_theta_ym * _gradPhiZ[idx_ym]
                                + (eps_ym / ((tau_ym * tau_ym) + FLT_EPSILON)) * eps_phi_ym * gradPhiMag2_ym * _gradPhiX[idx_ym];
                }
                float term_y = (val_y_plus - val_y_minus) / (2.0f * _dy);

                // 第四项：-∂/∂z(ε·∂ε/∂θ·τ)
                float val_eps_tau_z_plus = _epsilon[_INDEX(i, j, k_plus)] * _epsilonDerivTheta[_INDEX(i, j, k_plus)] * _tau_field[_INDEX(i, j, k_plus)];
                float val_eps_tau_z_minus = _epsilon[_INDEX(i, j, k_minus)] * _epsilonDerivTheta[_INDEX(i, j, k_minus)] * _tau_field[_INDEX(i, j, k_minus)];
                float term_eps_tau = -(val_eps_tau_z_plus - val_eps_tau_z_minus) / (2.0f * _dz);

                // 第五项：-g'(η)，其中 g(η) = η²(η-1)²/4
                // g'(η) = η(η-1)(η-0.5)
                float g_prime = oldPhi * (oldPhi - 1.0f) * (oldPhi - 0.5f);

                // 第六项：-p'(η)(f_s - f_t + f_ori)
                // p(η) = η²(3-2η), p'(η) = 6η(1-η)
                float p_prime = 6.0f * oldPhi * (1.0f - oldPhi);

                // 使用已经计算好的取向场梯度模
                float gradOmegaOriMag = _gradOmegaOriMag[idx];

                // f_s - f_t + f_ori
                // f_t = 0 (液相自由能), f_s = -m/6 (固相自由能), f_ori = H·||∇Ω_ori||
                float f_diff = -m / 6.0f + _H * gradOmegaOriMag;

                // 计算 ∂η/∂t
                float dPhiDt = M_eta * (term_diffusion + term_grad_eps2 + term_z + term_y + term_eps_tau
                                      - g_prime - p_prime * f_diff);

                // 更新相场
                _phi[idx] = oldPhi + dPhiDt * _dt;

                // ========== 热传导方程 + 潜热释放 ==========
                // ∂T/∂t = a²·∇²T + K·∂η/∂t
                _t[idx] = oldT + (_alpha_T * _lapT[idx] + _K * dPhiDt) * _dt;

                // ========== 公式(18)：取向场演化方程 ==========
                // ∂Ω_ori/∂t = -M_ori·H·(1-p(η))·∇·[p(η)·∇Ω_ori/||∇Ω_ori||]
                // Ω_ori 用单位向量 (x, y, z) 表示
                // 这个方程需要在切空间中计算

                float p_eta = oldPhi * oldPhi * (3.0f - 2.0f * oldPhi);

                // 计算取向场的拉普拉斯算子（对每个分量）
                float lapOmegaX = (_omega_ori_x[_INDEX(i_plus, j, k)] + _omega_ori_x[_INDEX(i_minus, j, k)]
                                 + _omega_ori_x[_INDEX(i, j_plus, k)] + _omega_ori_x[_INDEX(i, j_minus, k)]
                                 + _omega_ori_x[_INDEX(i, j, k_plus)] + _omega_ori_x[_INDEX(i, j, k_minus)]
                                 - 6.0f * oldOmegaX) / (_dx * _dx);

                float lapOmegaY = (_omega_ori_y[_INDEX(i_plus, j, k)] + _omega_ori_y[_INDEX(i_minus, j, k)]
                                 + _omega_ori_y[_INDEX(i, j_plus, k)] + _omega_ori_y[_INDEX(i, j_minus, k)]
                                 + _omega_ori_y[_INDEX(i, j, k_plus)] + _omega_ori_y[_INDEX(i, j, k_minus)]
                                 - 6.0f * oldOmegaY) / (_dx * _dx);

                float lapOmegaZ = (_omega_ori_z[_INDEX(i_plus, j, k)] + _omega_ori_z[_INDEX(i_minus, j, k)]
                                 + _omega_ori_z[_INDEX(i, j_plus, k)] + _omega_ori_z[_INDEX(i, j_minus, k)]
                                 + _omega_ori_z[_INDEX(i, j, k_plus)] + _omega_ori_z[_INDEX(i, j, k_minus)]
                                 - 6.0f * oldOmegaZ) / (_dx * _dx);

                // 投影到切空间：去除法向分量
                // 切空间拉普拉斯 = 拉普拉斯 - (拉普拉斯·ω)ω
                float lap_dot_omega = lapOmegaX * oldOmegaX + lapOmegaY * oldOmegaY + lapOmegaZ * oldOmegaZ;
                float lapOmegaX_tangent = lapOmegaX - lap_dot_omega * oldOmegaX;
                float lapOmegaY_tangent = lapOmegaY - lap_dot_omega * oldOmegaY;
                float lapOmegaZ_tangent = lapOmegaZ - lap_dot_omega * oldOmegaZ;

                // 计算演化速率
                float coeff = -M_ori * _H * (1.0f - p_eta) * p_eta;
                if (gradOmegaOriMag > FLT_EPSILON) {
                    coeff /= gradOmegaOriMag;
                } else {
                    coeff = 0.0f;
                }

                float dOmegaXDt = coeff * lapOmegaX_tangent;
                float dOmegaYDt = coeff * lapOmegaY_tangent;
                float dOmegaZDt = coeff * lapOmegaZ_tangent;

                // 更新取向场
                float newOmegaX = oldOmegaX + dOmegaXDt * _dt;
                float newOmegaY = oldOmegaY + dOmegaYDt * _dt;
                float newOmegaZ = oldOmegaZ + dOmegaZDt * _dt;

                // 重新归一化到单位球面
                float omegaNorm = sqrt(newOmegaX * newOmegaX + newOmegaY * newOmegaY + newOmegaZ * newOmegaZ);
                if (omegaNorm > FLT_EPSILON) {
                    _omega_ori_x[idx] = newOmegaX / omegaNorm;
                    _omega_ori_y[idx] = newOmegaY / omegaNorm;
                    _omega_ori_z[idx] = newOmegaZ / omegaNorm;
                } else {
                    // 如果归一化失败，保持原值
                    _omega_ori_x[idx] = oldOmegaX;
                    _omega_ori_y[idx] = oldOmegaY;
                    _omega_ori_z[idx] = oldOmegaZ;
                }
            }
        }
    }
}

// 主更新循环
void Kobayashi::update() {
    if (!_updateFlag) return; // 如果暂停则不计算

    // 为了加快视觉效果，每一帧渲染前，我们计算 10 次物理步骤
    for (int i = 0; i < 10; i++) {
        _computeGradientLaplacian();
        _evolution();
    }
    _updateTexture(); // 计算完后，准备将数据传给显卡
}

void Kobayashi::reset() {
    _vectorInit();
}

// ==========================================
// 图形渲染部分 (OpenGL)
// ==========================================

// 初始化 OpenGL 设置
void Kobayashi::glInit() {
    // 启用深度测试（3D渲染必需）
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // 启用混合（用于半透明效果）
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // 设置点的大小
    glPointSize(2.0f);

    // 设置背景颜色为深灰色
    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
}

// 更新纹理（3D版本 - 实际上不需要纹理，这里保留接口）
void Kobayashi::_updateTexture()
{
    // 3D渲染不需要预先生成纹理
    // 直接在 glRender() 中绘制体素
}

// 3D渲染函数 - 使用点云渲染晶体
void Kobayashi::glRender() {
    // 清除颜色和深度缓冲区
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // 设置投影矩阵
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = 1.0f; // 假设窗口是正方形
    gluPerspective(45.0, aspect, 0.1, 100.0);

    // 设置模型视图矩阵
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // 相机位置：从远处观察晶体
    static float angle = 0.0f;
    angle += 0.5f; // 自动旋转
    float camDist = 3.0f;
    gluLookAt(camDist * cos(angle * PI_F / 180.0f), camDist * 0.5f, camDist * sin(angle * PI_F / 180.0f),  // 相机位置
              0.0f, 0.0f, 0.0f,   // 看向原点
              0.0f, 1.0f, 0.0f);  // 上方向

    // 将晶体居中并缩放到合适大小
    glTranslatef(-0.5f, -0.5f, -0.5f); // 移动到原点
    float scale = 1.0f / (float)_objectCount.x;
    glScalef(scale, scale, scale);

    // 绘制晶体体素
    glBegin(GL_POINTS);

    // 定义颜色映射
    auto getColor = [](float phi) -> void {
        if (phi < 0.1f) {
            // 液相：不绘制或半透明蓝色
            return;
        } else if (phi < 0.5f) {
            // 界面区域：蓝色到青色渐变
            float t = (phi - 0.1f) / 0.4f;
            glColor4f(0.2f, 0.5f + 0.5f * t, 1.0f, 0.3f + 0.4f * t);
        } else if (phi < 0.9f) {
            // 过渡区域：青色到白色
            float t = (phi - 0.5f) / 0.4f;
            glColor4f(0.5f + 0.5f * t, 0.8f + 0.2f * t, 1.0f, 0.7f + 0.3f * t);
        } else {
            // 固相中心：白色
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        }
    };

    // 遍历所有体素，只绘制相场值大于阈值的点
    int skip = 1; // 采样间隔，可以调整以提高性能
    for (int k = 0; k < _objectCount.z; k += skip) {
        for (int j = 0; j < _objectCount.y; j += skip) {
            for (int i = 0; i < _objectCount.x; i += skip) {
                int idx = _INDEX(i, j, k);
                float phi = _phi[idx];

                if (phi > 0.1f) { // 只绘制相场值大于0.1的点
                    getColor(phi);
                    glVertex3f((float)i, (float)j, (float)k);
                }
            }
        }
    }

    glEnd();

    // 交换缓冲区
    glutSwapBuffers();
}
