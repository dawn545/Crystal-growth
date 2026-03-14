#include "Kobayashi3D.h"
#include <cfloat> // 包含 FLT_EPSILON，用于浮点数比较，防止除以零

// ==========================================
// 构造函数与初始化
// ==========================================

Kobayashi::Kobayashi(int x, int y, float timeStep) {
    _objectCount = { x, y }; // 网格大小，例如 300x300
    _dx = 0.03f; // 空间步长：模拟世界中每个格子代表的物理距离
    _dy = 0.03f; 
    _dt = timeStep; // 时间步长：每次模拟迭代推进的时间量
    
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

    M_ori = 1/_tau;           // 界面迁移率系数

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

    // 各向异性参考角度 θ_ori，决定晶体的初始朝向（单位：弧度）
    // 物理意义：控制晶体主轴的旋转角度
    // θ_ori = 0.0 表示晶体的一个主轴沿 x 轴正方向
    // θ_ori = π/6 会让六角形晶体旋转 30 度
    theta_ori = 0.0f;

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
}

// 分配内存并重置模拟状态
void Kobayashi::_vectorInit() {
    size_t vSize = _objectCount.x * _objectCount.y;
    
    // _phi: 相场变量 (0=液, 1=固)
    // _t: 温度场
    _phi.assign(vSize, 0.0f); 
    _t.assign(vSize, 0.0f);
    
    // 梯度(Gradient)和拉普拉斯(Laplacian)缓存，用于计算变化率
    _gradPhiX.assign(vSize, 0.0f); 
    _gradPhiY.assign(vSize, 0.0f);
    _lapPhi.assign(vSize, 0.0f); 
    _lapT.assign(vSize, 0.0f);
    
    // 辅助变量：角度、各向异性系数及其导数
    _angl.assign(vSize, 0.0f); 
    _epsilon.assign(vSize, 0.0f); 
    _epsilonDeriv.assign(vSize, 0.0f);
    
    // 像素缓冲区：这是我们要传给显卡的数据，每个像素4个字节 (R,G,B,A)
    _pixelBuffer.assign(vSize * 4, 0);
    
    // 在中心创建一个初始晶核
    _createNucleus(_objectCount.x / 2, _objectCount.y / 2);
    
    // 立即更新一次纹理，确保初始画面不是黑的
    _updateTexture();
}

// 在网格中心放置一个微小的“种子”，让晶体开始生长
void Kobayashi::_createNucleus(int x, int y)
{
    // 将中心及周围几个点设为 1.0 (固体)
    _phi[_INDEX(x, y)] = 1.0f;
    _phi[_INDEX(x - 1, y)] = 1.0f;
    _phi[_INDEX(x + 1, y)] = 1.0f;
    _phi[_INDEX(x, y - 1)] = 1.0f;
    _phi[_INDEX(x, y + 1)] = 1.0f;
}

// ==========================================
// 物理模拟核心：计算导数
// ==========================================

// 计算空间导数（梯度和拉普拉斯算子）
// 这是有限差分法的核心：通过邻居格子的值来推算当前的斜率和曲率
// 参考：有限差分法 (Finite Difference Method, FDM)
void Kobayashi::_computeGradientLaplacian()
{
    for (int j = 0; j < _objectCount.y; j++)
    {
        for (int i = 0; i < _objectCount.x; i++)
        {
            // 周期性边界条件 (Periodic Boundary Condition)
            // 意味着如果你走出右边界，会从左边界回来（像《贪吃蛇》或地球仪）
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            // ========== 1. 计算一阶导数 (Gradient / 梯度) ==========
            // 物理意义：相场变化的”坡度”，用于确定界面法线方向
            // 公式：∂φ/∂x ≈ (φ(i+1,j) - φ(i-1,j)) / (2·Δx)  [中心差分法]
            // 公式：∂φ/∂y ≈ (φ(i,j+1) - φ(i,j-1)) / (2·Δy)
            _gradPhiX[_INDEX(i, j)] = (_phi[_INDEX(i_plus, j)] - _phi[_INDEX(i_minus, j)]) / _dx;
            _gradPhiY[_INDEX(i, j)] = (_phi[_INDEX(i, j_plus)] - _phi[_INDEX(i, j_minus)]) / _dy;

            // ========== 2. 计算二阶导数 (Laplacian / 拉普拉斯算子) ==========
            // 物理意义：场的”扩散”趋势。如果是凸的，数值会扩散出去；如果是凹的，会汇聚进来。
            // 这里使用的是 9点差分法 (9-point stencil)
            // 公式：∇²φ = ∂²φ/∂x² + ∂²φ/∂y²
            // 展开式：∇²φ ≈ [2(φ_E + φ_W + φ_N + φ_S) + φ_NE + φ_NW + φ_SE + φ_SW - 12φ_C] / (3·Δx²)
            // 其中 E=East(右), W=West(左), N=North(上), S=South(下), C=Center(中心)
            _lapPhi[_INDEX(i, j)] = (2.0f * (_phi[_INDEX(i_plus, j)] + _phi[_INDEX(i_minus, j)] + _phi[_INDEX(i, j_plus)] + _phi[_INDEX(i, j_minus)])
                + _phi[_INDEX(i_plus, j_plus)] + _phi[_INDEX(i_minus, j_minus)] + _phi[_INDEX(i_minus, j_plus)] + _phi[_INDEX(i_plus, j_minus)]
                - 12.0f * _phi[_INDEX(i, j)]) / (3.0f * _dx * _dx);

            // 同样计算温度场的拉普拉斯算子
            _lapT[_INDEX(i, j)] = (2.0f * (_t[_INDEX(i_plus, j)] + _t[_INDEX(i_minus, j)] + _t[_INDEX(i, j_plus)] + _t[_INDEX(i, j_minus)])
                + _t[_INDEX(i_plus, j_plus)] + _t[_INDEX(i_minus, j_minus)] + _t[_INDEX(i_minus, j_plus)] + _t[_INDEX(i_plus, j_minus)]
                - 12.0f * _t[_INDEX(i, j)]) / (3.0f * _dx * _dx);

            // ========== 3. 计算界面法线角度 (Angle / 方向角) ==========
            // 物理意义：晶体界面的法线方向，用于确定各向异性系数
            // 公式：θ = atan2(∂φ/∂y, ∂φ/∂x)
            // 注意：需要处理梯度为零的特殊情况（即平坦区域，没有明确的角度）
            // FLT_EPSILON 是浮点数最小精度，用于判断是否接近零
            if (_gradPhiX[_INDEX(i, j)] <= +FLT_EPSILON && _gradPhiX[_INDEX(i, j)] >= -FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = -0.5f * PI_F;
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = 0.5f * PI_F;

            if (_gradPhiX[_INDEX(i, j)] > +FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = 2.0f * PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            if (_gradPhiX[_INDEX(i, j)] < -FLT_EPSILON)
                _angl[_INDEX(i, j)] = PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            // ========== 4. 计算各向异性系数 (Anisotropy Coefficient) ==========
            // 物理意义：晶体不是圆球，它在不同方向生长速度不同（这就形成了雪花形状）。
            // 这个公式让 epsilon 随角度变化，从而产生各向异性。
            // 公式：ε(θ) = ε̄(1 + δ·cos(j_fold·(θ_ori - θ)))
            // 其中：ε̄ = 平均各向异性强度，δ = 各向异性强度，j_fold = 折叠对称性（6=六角形），θ_ori = 参考角度
            _epsilon[_INDEX(i, j)] = _epsilonBar * (1.0f + _delta * cos(j_fold * (theta_ori - _angl[_INDEX(i, j)])));

            // 各向异性系数的导数，用于计算界面能的梯度
            // 公式：dε/dθ = ε̄·j_fold·δ·sin(j_fold·(θ_ori - θ))
            _epsilonDeriv[_INDEX(i, j)] = _epsilonBar * j_fold * _delta * sin(j_fold * (theta_ori - _angl[_INDEX(i, j)]));
        }
    }
}

// ==========================================
// 物理模拟核心：时间演化
// ==========================================

// 根据微分方程更新 _phi 和 _t 的值
// 这是 Kobayashi 模型的核心演化方程
// 参考：Kobayashi, A. (1993) "Modeling and numerical simulations of dendritic crystal growth"
void Kobayashi::_evolution()
{
    for (int j = 0; j < _objectCount.y; j++)
    {
        for (int i = 0; i < _objectCount.x; i++)
        {
            // 获取邻居索引（周期性边界条件）
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            // ========== 计算各向异性界面能的贡献项 ==========
            // 这些项来自于变分法推导，用于处理各向异性界面能
            // 公式：∇·(ε²∇φ) 的展开式

            // term1 和 term2：各向异性系数及其导数的贡献
            // 这些项确保晶体在不同方向有不同的生长速度
            float gradEpsPowX = (_epsilon[_INDEX(i_plus, j)] * _epsilon[_INDEX(i_plus, j)] - _epsilon[_INDEX(i_minus, j)] * _epsilon[_INDEX(i_minus, j)]) / _dx;
            float gradEpsPowY = (_epsilon[_INDEX(i, j_plus)] * _epsilon[_INDEX(i, j_plus)] - _epsilon[_INDEX(i, j_minus)] * _epsilon[_INDEX(i, j_minus)]) / _dy;

            // 公式：term1 = ∂/∂y[ε·dε/dθ·∂φ/∂x]
            float term1 = (_epsilon[_INDEX(i, j_plus)] * _epsilonDeriv[_INDEX(i, j_plus)] * _gradPhiX[_INDEX(i, j_plus)]
                - _epsilon[_INDEX(i, j_minus)] * _epsilonDeriv[_INDEX(i, j_minus)] * _gradPhiX[_INDEX(i, j_minus)]) / _dy;

            // 公式：term2 = -∂/∂x[ε·dε/dθ·∂φ/∂y]
            float term2 = -(_epsilon[_INDEX(i_plus, j)] * _epsilonDeriv[_INDEX(i_plus, j)] * _gradPhiY[_INDEX(i_plus, j)]
                - _epsilon[_INDEX(i_minus, j)] * _epsilonDeriv[_INDEX(i_minus, j)] * _gradPhiY[_INDEX(i_minus, j)]) / _dx;

            // 公式：term3 = ∇(ε²)·∇φ
            float term3 = gradEpsPowX * _gradPhiX[_INDEX(i, j)] + gradEpsPowY * _gradPhiY[_INDEX(i, j)];

            // ========== 计算驱动力 (Driving Force) ==========
            // 驱动力由过冷度（T_eq - T）决定
            // 公式：m = (α/π)·arctan(γ·(T_eq - T))
            // 物理意义：温度越低（过冷度越大），驱动力越大，晶体生长越快
            float m = _alpha / PI_F * atan(_gamma * (_tEq - _t[_INDEX(i, j)]));

            // 保存旧值，用于计算潜热释放
            float oldPhi = _phi[_INDEX(i, j)];
            float oldT = _t[_INDEX(i, j)];

            // ========== 核心更新公式 1：相场演化方程 (Allen-Cahn Equation) ==========
            // 这是一个非线性偏微分方程，描述相场随时间的变化
            // 公式：∂φ/∂t = M_ori*[ε²∇²φ + term1 + term2 + term3 + φ(1-φ)(φ-0.5+m)]
            // 其中：
            //   - ε²∇²φ：各向异性界面能的扩散项
            //   - term1, term2, term3：各向异性系数的贡献
            //   - φ(1-φ)(φ-0.5+m)：双势阱势能项 + 驱动力项
            // 时间积分：使用显式欧拉法 φ_new = φ_old + (∂φ/∂t)·Δt

            // 先计算 ∂φ/∂t
            float dPhiDt = (term1 + term2 + _epsilon[_INDEX(i, j)] * _epsilon[_INDEX(i, j)] * _lapPhi[_INDEX(i, j)] + term3
                    + oldPhi * (1.0f - oldPhi) * (oldPhi - 0.5f + m)) * M_ori;

            // 更新相场
            _phi[_INDEX(i, j)] = oldPhi + dPhiDt * _dt;

            // ========== 核心更新公式 2：热传导方程 + 潜热释放 ==========
            // 这个方程描述温度场如何随时间变化
            // 公式：∂T/∂t = a²·∇²T + K·∂φ/∂t
            // 其中：
            //   - a²·∇²T：热扩散项（温度向周围扩散，a² 是热扩散系数）
            //   - K·∂φ/∂t：潜热释放项（相变速率，注意这里用的是 ∂φ/∂t 而不是 Δφ）
            // 物理意义：当液体凝固成固体时（φ 从 0 变为 1），会释放潜热，加热周围
            //          这个热量会减缓进一步的晶体生长（自我调节机制）
            // 时间积分：使用显式欧拉法 T_new = T_old + (∂T/∂t)·Δt
            _t[_INDEX(i, j)] = oldT + (_alpha_T * _lapT[_INDEX(i, j)] + _K * dPhiDt) * _dt;
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

// 初始化 OpenGL 纹理设置
void Kobayashi::glInit() {
    glEnable(GL_TEXTURE_2D); // 开启 2D 纹理功能
    glGenTextures(1, &_textureID); // 申请一个纹理 ID
    glBindTexture(GL_TEXTURE_2D, _textureID); // 绑定这个 ID 为当前操作对象

    // 设置纹理过滤方式 (Linear = 线性插值，图像放大时会模糊平滑，而不是马赛克)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // 设置纹理包裹方式 (Clamp = 边缘拉伸，防止纹理重复平铺)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
}

// 将模拟数据 (_phi) 转换为颜色数据 (_pixelBuffer) 并上传到显卡
// 这个函数实现了从物理场到可视化的映射
void Kobayashi::_updateTexture()
{
    // 定义颜色 (RGB格式, 范围 0.0-1.0)
    // 这些颜色用于表示晶体的不同状态
    struct float3 { float x, y, z; };
    float3 c0 = { 0.0f, 0.0f, 0.0f };             // 黑色 (背景/液体，φ ≈ 0)
    float3 c1 = { 0.25f, 0.50f, 0.98f };          // 蓝色 (晶体边缘，φ ≈ 0.9)
    float3 c2 = { 0.36f, 1.00f, 0.98f };          // 青色 (过渡区域，φ ≈ 0.99)
    float3 c3 = { 0.90f, 1.00f, 0.98f };          // 白色 (晶体中心，φ ≈ 1.0)

    // 定义颜色分界线（相场值的阈值）
    // 这些值决定了不同颜色的过渡点
    float c1Boundary = 0.9f;   // φ < 0.9：从黑色过渡到蓝色
    float c2Boundary = 0.99f;  // 0.9 < φ < 0.99：从蓝色过渡到青色
                               // φ > 0.99：从青色过渡到白色

    size_t size = _phi.size();

    // 遍历每一个格子
    for (size_t k = 0; k < size; k++)
    {
        float phi = _phi[k]; // 获取当前格子的相场值 (0.0 = 液体, 1.0 = 固体)
        float3 color;
        float ratio;

        // ========== 颜色映射逻辑 (Color Mapping) ==========
        // 根据 phi 的值，在 c0, c1, c2, c3 之间进行线性插值 (Linear Interpolation / Lerp)
        // 公式：color = c_old + (c_new - c_old) * ratio，其中 ratio ∈ [0, 1]
        if (phi <= c1Boundary)
        {
            // 区间 1：φ ∈ [0, 0.9]，从黑色 (c0) 过渡到蓝色 (c1)
            ratio = phi * (1.0f / c1Boundary);  // 归一化到 [0, 1]
            color.x = c0.x * (1.0f - ratio) + c1.x * ratio;
            color.y = c0.y * (1.0f - ratio) + c1.y * ratio;
            color.z = c0.z * (1.0f - ratio) + c1.z * ratio;
        }
        else if (phi > c1Boundary && phi <= c2Boundary)
        {
            // 区间 2：φ ∈ (0.9, 0.99]，从蓝色 (c1) 过渡到青色 (c2)
            ratio = (phi - c1Boundary) * (1.0f / (c2Boundary - c1Boundary));  // 归一化到 [0, 1]
            color.x = c1.x * (1.0f - ratio) + c2.x * ratio;
            color.y = c1.y * (1.0f - ratio) + c2.y * ratio;
            color.z = c1.z * (1.0f - ratio) + c2.z * ratio;
        }
        else
        {
            // 区间 3：φ ∈ (0.99, 1.0]，从青色 (c2) 过渡到白色 (c3)
            float c3Boundary = 1.0f;
            ratio = (phi - c2Boundary) * (1.0f / (c3Boundary - c2Boundary));  // 归一化到 [0, 1]
            color.x = c2.x * (1.0f - ratio) + c3.x * ratio;
            color.y = c2.y * (1.0f - ratio) + c3.y * ratio;
            color.z = c2.z * (1.0f - ratio) + c3.z * ratio;
        }

        // ========== 数据打包 (Data Packing) ==========
        // 将浮点颜色 (0.0-1.0) 转换为字节 (0-255)，以便 GPU 处理
        auto toByte = [](float v) -> unsigned char {
            if (v < 0.0f) return 0;
            if (v > 1.0f) return 255;
            return static_cast<unsigned char>(v * 255.0f);
        };

        // 填充 RGBA 像素数据 (Red, Green, Blue, Alpha)
        // 每个像素占 4 个字节，按 RGBA 顺序排列
        _pixelBuffer[k * 4 + 0] = toByte(color.x);  // Red 通道
        _pixelBuffer[k * 4 + 1] = toByte(color.y);  // Green 通道
        _pixelBuffer[k * 4 + 2] = toByte(color.z);  // Blue 通道
        _pixelBuffer[k * 4 + 3] = 255;              // Alpha 通道 = 255 (完全不透明)
    }

    // ========== 上传纹理到 GPU ==========
    glBindTexture(GL_TEXTURE_2D, _textureID);

    // glTexImage2D 函数：将 CPU 内存中的像素数据上传到 GPU 显存
    // 参数解释：
    //   GL_TEXTURE_2D：目标纹理类型（2D 纹理）
    //   0：Mipmap 层级（0 表示原始分辨率，不使用 Mipmap）
    //   GL_RGBA：显卡内部存储格式（4 通道，每通道 8 位）
    //   _objectCount.x, _objectCount.y：纹理尺寸（宽度和高度，单位：像素）
    //   0：边框宽度（必须是 0，OpenGL 规范要求）
    //   GL_RGBA：我们提供的数据格式（与内部格式相同）
    //   GL_UNSIGNED_BYTE：我们提供的数据类型（无符号字节，范围 0-255）
    //   _pixelBuffer.data()：指向像素数据的指针
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _objectCount.x, _objectCount.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, _pixelBuffer.data());
}

// 实际的绘制函数，每帧调用一次
// 这个函数负责将计算好的纹理显示到屏幕上
void Kobayashi::glRender() {
    // 1. 清除屏幕颜色缓冲区
    // 这会用当前的清除颜色（通常是黑色）填充整个屏幕
    glClear(GL_COLOR_BUFFER_BIT);

    // 2. 启用 2D 纹理功能并绑定我们要使用的纹理
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _textureID);

    // 设置绘制颜色为白色 (1.0, 1.0, 1.0)
    // 这样纹理会保持原色，不会被额外的颜色调制（Modulation）
    glColor3f(1.0f, 1.0f, 1.0f);

    // 3. 绘制一个四边形 (Quad)，充满整个视口
    // OpenGL 的坐标系：
    //   - 中心是 (0, 0)
    //   - 左下是 (-1, -1)
    //   - 右上是 (1, 1)
    // 纹理坐标系：
    //   - 左下是 (0, 0)
    //   - 右上是 (1, 1)
    // 我们的任务是将纹理的四个角对应到四边形的四个角上
    glBegin(GL_QUADS);
        // 左下角：纹理坐标 (0, 0) → 屏幕坐标 (-1, -1)
        glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);

        // 右下角：纹理坐标 (1, 0) → 屏幕坐标 (1, -1)
        glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);

        // 右上角：纹理坐标 (1, 1) → 屏幕坐标 (1, 1)
        glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);

        // 左上角：纹理坐标 (0, 1) → 屏幕坐标 (-1, 1)
        glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();

    // 4. 交换缓冲区 (Double Buffering)
    // OpenGL 使用双缓冲技术：
    //   - 我们在一个隐藏的后台缓冲区 (Back Buffer) 画画
    //   - 画好后瞬间交换到前台缓冲区 (Front Buffer) 显示
    //   - 这样可以防止画面闪烁和撕裂 (Tearing)
    glutSwapBuffers();
}
