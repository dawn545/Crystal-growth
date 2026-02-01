#include "Kobayashi.h"
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
void Kobayashi::_initParams() {
    _tau = 0.0003f;       // 弛豫时间，控制相变速度
    _epsilonBar = 0.010f; // 平均各向异性强度，决定晶界厚度
    _mu = 1.0f;           
    _K = 1.6f;            // 潜热系数，控制温度场对相变的影响
    _delta = 0.05f;       // 各向异性强度，值越大晶体越有棱角
    _anisotropy = 6.0f;   // 各向异性模数：6.0 代表六角形对称（雪花状）
    _alpha = 0.9f; 
    _gamma = 10.0f; 
    _tEq = 1.0f;          // 平衡温度
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
void Kobayashi::_computeGradientLaplacian()
{
    for (int j = 0; j < _objectCount.y; j++)
    {
        for (int i = 0; i < _objectCount.x; i++)
        {
            // 周期性边界条件 (Periodic Boundary)
            // 意味着如果你走出右边界，会从左边界回来（像《贪吃蛇》或地球仪）
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            // 1. 计算一阶导数 (Gradient / 梯度)
            // 物理意义：相场变化的“坡度”，用于确定界面法线方向
            _gradPhiX[_INDEX(i, j)] = (_phi[_INDEX(i_plus, j)] - _phi[_INDEX(i_minus, j)]) / _dx;
            _gradPhiY[_INDEX(i, j)] = (_phi[_INDEX(i, j_plus)] - _phi[_INDEX(i, j_minus)]) / _dy;

            // 2. 计算二阶导数 (Laplacian / 拉普拉斯算子)
            // 物理意义：场的“扩散”趋势。如果是凸的，数值会扩散出去；如果是凹的，会汇聚进来。
            // 这里使用的是 5点差分法 或者 9点差分法 的变种
            _lapPhi[_INDEX(i, j)] = (2.0f * (_phi[_INDEX(i_plus, j)] + _phi[_INDEX(i_minus, j)] + _phi[_INDEX(i, j_plus)] + _phi[_INDEX(i, j_minus)])
                + _phi[_INDEX(i_plus, j_plus)] + _phi[_INDEX(i_minus, j_minus)] + _phi[_INDEX(i_minus, j_plus)] + _phi[_INDEX(i_plus, j_minus)]
                - 12.0f * _phi[_INDEX(i, j)]) / (3.0f * _dx * _dx);
            
            _lapT[_INDEX(i, j)] = (2.0f * (_t[_INDEX(i_plus, j)] + _t[_INDEX(i_minus, j)] + _t[_INDEX(i, j_plus)] + _t[_INDEX(i, j_minus)])
                + _t[_INDEX(i_plus, j_plus)] + _t[_INDEX(i_minus, j_minus)] + _t[_INDEX(i_minus, j_plus)] + _t[_INDEX(i_plus, j_minus)]
                - 12.0f * _t[_INDEX(i, j)]) / (3.0f * _dx * _dx);

            // 3. 计算界面法线角度 (Angle)
            // 通过梯度 (dx, dy) 使用 atan2 计算角度。
            // FLT_EPSILON 用于处理梯度接近 0 的情况（即平坦区域，没有角度）
            if (_gradPhiX[_INDEX(i, j)] <= +FLT_EPSILON && _gradPhiX[_INDEX(i, j)] >= -FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = -0.5f * PI_F;
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = 0.5f * PI_F;

            if (_gradPhiX[_INDEX(i, j)] > +FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = 2.0f * PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            if (_gradPhiX[_INDEX(i, j)] < -FLT_EPSILON)
                _angl[_INDEX(i, j)] = PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            // 4. 计算各向异性系数 (Epsilon)
            // 晶体不是圆球，它在不同方向生长速度不同（这就形成了雪花形状）。
            // 这个公式让 epsilon 随角度变化。
            _epsilon[_INDEX(i, j)] = _epsilonBar * (1.0f + _delta * cos(_anisotropy * _angl[_INDEX(i, j)]));
            _epsilonDeriv[_INDEX(i, j)] = -_epsilonBar * _anisotropy * _delta * sin(_anisotropy * _angl[_INDEX(i, j)]);
        }
    }
}

// ==========================================
// 物理模拟核心：时间演化
// ==========================================

// 根据微分方程更新 _phi 和 _t 的值
void Kobayashi::_evolution()
{
    for (int j = 0; j < _objectCount.y; j++)
    {
        for (int i = 0; i < _objectCount.x; i++)
        {
            // 获取邻居索引
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            // 计算 epsilon 的梯度的平方项（复杂的物理推导结果，用于处理各向异性界面能）
            float gradEpsPowX = (_epsilon[_INDEX(i_plus, j)] * _epsilon[_INDEX(i_plus, j)] - _epsilon[_INDEX(i_minus, j)] * _epsilon[_INDEX(i_minus, j)]) / _dx;
            float gradEpsPowY = (_epsilon[_INDEX(i, j_plus)] * _epsilon[_INDEX(i, j_plus)] - _epsilon[_INDEX(i, j_minus)] * _epsilon[_INDEX(i, j_minus)]) / _dy;

            float term1 = (_epsilon[_INDEX(i, j_plus)] * _epsilonDeriv[_INDEX(i, j_plus)] * _gradPhiX[_INDEX(i, j_plus)]
                - _epsilon[_INDEX(i, j_minus)] * _epsilonDeriv[_INDEX(i, j_minus)] * _gradPhiX[_INDEX(i, j_minus)]) / _dy;

            float term2 = -(_epsilon[_INDEX(i_plus, j)] * _epsilonDeriv[_INDEX(i_plus, j)] * _gradPhiY[_INDEX(i_plus, j)]
                - _epsilon[_INDEX(i_minus, j)] * _epsilonDeriv[_INDEX(i_minus, j)] * _gradPhiY[_INDEX(i_minus, j)]) / _dx;
            
            float term3 = gradEpsPowX * _gradPhiX[_INDEX(i, j)] + gradEpsPowY * _gradPhiY[_INDEX(i, j)];

            // m项：驱动力项，受过冷度（T_eq - T）控制
            float m = _alpha / PI_F * atan(_gamma * (_tEq - _t[_INDEX(i, j)]));

            float oldPhi = _phi[_INDEX(i, j)];
            float oldT = _t[_INDEX(i, j)];

            // === 核心更新公式 ===
            // 1. 更新相场 phi (Allen-Cahn equation 变体)
            _phi[_INDEX(i, j)] = _phi[_INDEX(i, j)] +
                (term1 + term2 + _epsilon[_INDEX(i, j)] * _epsilon[_INDEX(i, j)] * _lapPhi[_INDEX(i, j)] + term3
                    + oldPhi * (1.0f - oldPhi) * (oldPhi - 0.5f + m)) * _dt / _tau;
            
            // 2. 更新温度场 T (热传导方程 + 潜热释放)
            // _K * (phi_new - phi_old) 代表相变时释放的潜热，这会加热周围，减缓进一步生长
            _t[_INDEX(i, j)] = oldT + _lapT[_INDEX(i, j)] * _dt + _K * (_phi[_INDEX(i, j)] - oldPhi);
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
void Kobayashi::_updateTexture()
{
    // 定义颜色 (RGB格式, 范围 0.0-1.0)
    struct float3 { float x, y, z; };
    float3 c0 = { 0.0f, 0.0f, 0.0f };             // 黑色 (背景/液体)
    float3 c1 = { 0.25f, 0.50f, 0.98f };          // 蓝色 (边缘)
    float3 c2 = { 0.36f, 1.00f, 0.98f };          // 青色 (过渡)
    float3 c3 = { 0.90f, 1.00f, 0.98f };          // 白色 (晶体中心)

    // 定义颜色分界线
    float c1Boundary = 0.9f;
    float c2Boundary = 0.99f;

    size_t size = _phi.size();
    
    // 遍历每一个格子
    for (size_t k = 0; k < size; k++)
    {
        float phi = _phi[k]; // 获取当前格子的状态 (0.0 - 1.0)
        float3 color;
        float ratio;

        // --- 颜色映射逻辑 (Color Mapping) ---
        // 根据 phi 的值，在 c0, c1, c2, c3 之间进行线性插值 (Lerp)
        if (phi <= c1Boundary)
        {
            ratio = phi * (1.0f / c1Boundary);
            color.x = c0.x * (1.0f - ratio) + c1.x * ratio;
            color.y = c0.y * (1.0f - ratio) + c1.y * ratio;
            color.z = c0.z * (1.0f - ratio) + c1.z * ratio;
        }
        else if (phi > c1Boundary && phi <= c2Boundary)
        {
            ratio = (phi - c1Boundary) * (1.0f / (c2Boundary - c1Boundary));
            color.x = c1.x * (1.0f - ratio) + c2.x * ratio;
            color.y = c1.y * (1.0f - ratio) + c2.y * ratio;
            color.z = c1.z * (1.0f - ratio) + c2.z * ratio;
        }
        else
        {
            float c3Boundary = 1.0f;
            ratio = (phi - c2Boundary) * (1.0f / (c3Boundary - c2Boundary));
            color.x = c2.x * (1.0f - ratio) + c3.x * ratio;
            color.y = c2.y * (1.0f - ratio) + c3.y * ratio;
            color.z = c2.z * (1.0f - ratio) + c3.z * ratio;
        }

        // --- 数据打包 ---
        // 将浮点颜色 (0.0-1.0) 转换为字节 (0-255)
        auto toByte = [](float v) -> unsigned char {
            if (v < 0.0f) return 0;
            if (v > 1.0f) return 255;
            return static_cast<unsigned char>(v * 255.0f);
        };

        // 填充 RGBA (Red, Green, Blue, Alpha)
        _pixelBuffer[k * 4 + 0] = toByte(color.x);
        _pixelBuffer[k * 4 + 1] = toByte(color.y);
        _pixelBuffer[k * 4 + 2] = toByte(color.z);
        _pixelBuffer[k * 4 + 3] = 255; // Alpha = 255 (不透明)
    }

    // --- 上传纹理到 GPU ---
    glBindTexture(GL_TEXTURE_2D, _textureID);
    // glTexImage2D 参数解释：
    // GL_TEXTURE_2D: 目标类型
    // 0: Mipmap 层级 (0是原图)
    // GL_RGBA: 显卡内部存储格式
    // width, height: 纹理尺寸
    // 0: 边框 (必须是0)
    // GL_RGBA: 我们提供的数据格式
    // GL_UNSIGNED_BYTE: 我们提供的数据类型 (uchar)
    // data: 数据指针
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _objectCount.x, _objectCount.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, _pixelBuffer.data());
}

// 实际的绘制函数，每帧调用一次
void Kobayashi::glRender() {
    // 1. 清除屏幕颜色
    glClear(GL_COLOR_BUFFER_BIT);
    
    // 2. 绑定我们要使用的纹理
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _textureID);
    
    // 设置绘制颜色为白色 (这样纹理会保持原色，不会被染色)
    glColor3f(1.0f, 1.0f, 1.0f);

    // 3. 绘制一个四边形 (Quad)，充满整个视口
    // OpenGL 的坐标系：中心是 (0,0)，左下是 (-1,-1)，右上是 (1,1)
    // 纹理坐标 (TexCoord)：左下是 (0,0)，右上是 (1,1)
    // 我们的任务是将纹理的四个角对应到四边形的四个角上
    glBegin(GL_QUADS);
        // 左下角
        glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
        // 右下角
        glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);
        // 右上角
        glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);
        // 左上角
        glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();
    
    // 4. 交换缓冲区 (Double Buffering)
    // 我们在一个隐藏的缓冲区画画，画好后瞬间交换到前台显示，防止闪烁
    glutSwapBuffers(); 
}