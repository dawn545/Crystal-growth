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
    float _tau, _epsilonBar, M_eta, _K, _delta, j_fold, _alpha, _gamma, _tEq;
    float _alpha_T; // 热扩散系数 a² (thermal diffusivity)
    float _H; // 取向场驱动力系数
    float M_ori; // 取向场迁移率

    // 相场、温度场
    std::vector<float> _phi, _t;

    // 取向场：Ω_ori = (θ_ori, φ_ori)
    std::vector<float> _theta_ori, _phi_ori;

    // 相场梯度和拉普拉斯算子
    std::vector<float> _gradPhiX, _gradPhiY, _gradPhiZ;
    std::vector<float> _lapPhi, _lapT;

    // 相场梯度模和局部方向角
    std::vector<float> _gradPhiMag; // |∇η|
    std::vector<float> _tau_field;  // τ = sqrt((∂η/∂x)² + (∂η/∂y)²)
    std::vector<float> _theta, _phi_angle; // 局部相位前沿方向 Ω = (θ, φ)

    // 各向异性系数及其导数
    std::vector<float> _epsilon, _epsilonDerivTheta, _epsilonDerivPhi;

    // 取向场梯度：∇Ω_ori
    std::vector<float> _gradThetaOriX, _gradThetaOriY, _gradThetaOriZ;
    std::vector<float> _gradPhiOriX, _gradPhiOriY, _gradPhiOriZ;

    // OpenGL 纹理（显示中间切片）
    std::vector<unsigned char> _pixelBuffer;
    GLuint _textureID = 0;
    bool _updateFlag = true;
    int _renderSlice; // 渲染哪个z切片

    void _initParams();
    void _vectorInit();
    void _createNucleus(int x, int y, int z);
    void _computeGradientLaplacian();
    void _evolution();
    void _updateTexture();
};
