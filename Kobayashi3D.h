#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <GL/freeglut.h>

const float PI_F = 3.14159265358979f;

class Kobayashi3D
{
public:
    Kobayashi3D(int x, int y, float timeStep);
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

    // 相机控制
    void rotateCamera(float deltaX, float deltaY);
    void zoomCamera(float delta);

private:
    // 模拟参数
    struct int2 { int x; int y; };
    int2 _objectCount = { 0, 0 };
    inline int _INDEX(int i, int j) { return (i + _objectCount.x * j); };

    float _dx, _dy, _dt;
    float _tau, _epsilonBar, _mu, _K, _delta, _anisotropy, _alpha, _gamma, _tEq;

    // 相场数据
    std::vector<float> _phi, _t, _epsilon, _epsilonDeriv, _gradPhiX, _gradPhiY, _lapPhi, _lapT, _angl;

    // OpenGL 纹理和渲染
    std::vector<unsigned char> _pixelBuffer;
    GLuint _textureID = 0;
    bool _updateFlag = true;

    // 3D 相机参数
    float _cameraRotX = 30.0f;
    float _cameraRotY = 45.0f;
    float _cameraDistance = 3.0f;

    // 厚度参数
    float _crystalThickness = 0.03f;  // 冰晶最大厚度

    // 内部方法
    void _initParams();
    void _vectorInit();
    void _createNucleus(int x, int y);
    void _computeGradientLaplacian();
    void _evolution();
    void _updateTexture();
    void _render3DCrystal();
};
