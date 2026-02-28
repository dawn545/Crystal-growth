#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <GL/freeglut.h>

const float PI_F = 3.14159265358979f;


class Kobayashi3D
{
public:
    Kobayashi3D(int x, int y, int z, float timeStep);
    ~Kobayashi3D();

    // 核心模拟逻辑
    void update();
    void reset();

    // 3D 渲染逻辑
    void glInit();
    void glRender();

    // 控制接口
    void togglePause() { _updateFlag = !_updateFlag; }
    bool isPaused() const { return !_updateFlag; }

    // 取向场控制接口
    void resetOrientationField() { _resetOrientationField(); }
    void paintOrientation(int screenX, int screenY, float angle, float radius = 20.0f, float blendFactor = 0.6f);

    // 相机控制
    void rotateCamera(float deltaX, float deltaY);
    void zoomCamera(float delta);

private:
    struct int3 { int x; int y; int z; };
    int3 _objectCount = { 0, 0, 0 };

    // 3D 索引宏：将 (x, y, z) 转换为一维数组索引
    inline int IDX(int x, int y, int z) const {
        return z * _objectCount.y * _objectCount.x + y * _objectCount.x + x;
    }

    // 3D 边界处理：钳位到有效范围（Neumann 零通量边界条件）
    inline int clampX(int x) const { return (x < 0) ? 0 : ((x >= _objectCount.x) ? _objectCount.x - 1 : x); }
    inline int clampY(int y) const { return (y < 0) ? 0 : ((y >= _objectCount.y) ? _objectCount.y - 1 : y); }
    inline int clampZ(int z) const { return (z < 0) ? 0 : ((z >= _objectCount.z) ? _objectCount.z - 1 : z); }

    float _dx, _dy, _dz, _dt;
    float _tau, _epsilonBar, _mu, _K, _delta, _anisotropy, _alpha, _gamma, _tEq;

    // 3D 场数据：使用一维平铺数组（Flat Array）确保内存连续性
    std::vector<double> _phi;              // 相场 η(x,y,z)
    std::vector<double> _t;                // 温度场 T(x,y,z)
    std::vector<double> _orientationField; // 取向场 Ω(x,y,z) - 单位球上的向量 [θ, φ]

    // 辅助计算场
    std::vector<double> _epsilon, _epsilonDeriv;
    std::vector<double> _gradPhiX, _gradPhiY, _gradPhiZ;  // 3D 梯度分量
    std::vector<double> _lapPhi, _lapT;
    std::vector<double> _gradPhiMag;   // |∇η| - 梯度模长
    std::vector<double> _theta;        // θ = -arccos(∂η/∂z / |∇η|) - 极角
    std::vector<double> _phi_angle;    // φ = -arctan2(∂η/∂y, ∂η/∂x) - 方位角
    std::vector<double> _tau_field;    // τ = sqrt((∂η/∂x)² + (∂η/∂y)²) - 辅助算子
    std::vector<double> _angl;         // 保留用于调试


    // OpenGL 纹理和渲染（用于2D切片可视化）
    std::vector<unsigned char> _pixelBuffer;
    GLuint _textureID = 0;
    bool _updateFlag = true;


    // 3D 相机参数
    float _cameraRotX = 30.0f;
    float _cameraRotY = 45.0f;
    float _cameraDistance = 3.0f;
    
    float _crystalThickness = 0.03f;  // 冰晶最大厚度

    // 内部方法
    void _initParams();
    void _vectorInit();
    void _createNucleus(int x, int y, int z, float radius);
    void _computeGradientLaplacian();
    void _evolution();
    void _updateTexture();
    void _render3DCrystal();

    // 取向场管理 (Ren et al. 2018)
    void _resetOrientationField();      // 重置为全局统一方向
};
