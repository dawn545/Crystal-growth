#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>
#include <GL/freeglut.h>

const float PI_F = 3.14159265358979f;

class Kobayashi
{
public:
    Kobayashi(int x, int y, int z, float timeStep);
    ~Kobayashi();

    // 核心模拟逻辑
    void update();
    void reset();

    // 渲染逻辑
    void glInit();
    void glRender();

    // 简单的控制接口
    void togglePause() { _updateFlag = !_updateFlag; }
    bool isPaused() const { return !_updateFlag; }

private:
    // 3D 网格参数
    struct int3 { int x; int y; int z; };
    int3 _objectCount = { 0, 0, 0 };
    inline int _INDEX(int i, int j, int k) { return (i + _objectCount.x * (j + _objectCount.y * k)); };

    float _dx, _dy, _dz, _dt;
    float _tau, M_eta, _K, _alpha, _gamma, _tEq;
    float _alpha_T; // 热扩散系数 a² (thermal diffusivity)
    float _H; // 取向场驱动力系数
    float M_ori; // 取向场迁移率
    float _c1, _c2; // 公式(19)中的各向异性系数：ε_o(n) = c1 + c2*(sin⁴θ̃(sin⁴φ̃ + cos⁴φ̃) + cos⁴θ̃)

    // 相场、温度场
    std::vector<float> _phi, _t;

    // 相场梯度和拉普拉斯算子
    std::vector<float> _gradPhiX, _gradPhiY, _gradPhiZ;
    std::vector<float> _lapPhi, _lapT;

    // 相场梯度模和局部方向角
    std::vector<float> _gradPhiMag; // |∇η|
    std::vector<float> _tau_field;  // τ = sqrt((∂η/∂x)² + (∂η/∂y)²)
    std::vector<float> _theta, _phi_angle; // 局部相位前沿方向 Ω = (θ, φ)

    // 各向异性系数及其导数
    std::vector<float> _epsilon, _epsilonDerivTheta, _epsilonDerivPhi;

    // 取向场：Ω_ori 用单位球上的点表示
    // 使用笛卡尔坐标 (x, y, z) 存储单位向量
    std::vector<float> _omega_ori_x, _omega_ori_y, _omega_ori_z;

    // 取向场的局部极坐标表示 (ρ, λ)
    // ρ: 中心角（大圆距离）
    // λ: 立体投影后的极坐标角度
    std::vector<float> _rho_x_plus, _rho_x_minus, _rho_y_plus, _rho_y_minus, _rho_z_plus, _rho_z_minus;
    std::vector<float> _lambda_x_plus, _lambda_x_minus, _lambda_y_plus, _lambda_y_minus, _lambda_z_plus, _lambda_z_minus;

    // 取向场梯度：∇Ω_ori（使用 (ρ, λ) 计算）
    std::vector<float> _gradOmegaOriMag; // ||∇Ω_ori||

    // OpenGL 相关
    bool _updateFlag = true;

    void _initParams();
    void _vectorInit();
    void _createNucleus(int x, int y, int z);
    void _computeGradientLaplacian();
    void _evolution();
};
