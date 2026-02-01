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
    Kobayashi(int x, int y, float timeStep);
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
    // 模拟参数保持不变
    struct int2 { int x; int y; };
    int2 _objectCount = { 0, 0 };
    inline int _INDEX(int i, int j) { return (i + _objectCount.x * j); };

    float _dx, _dy, _dt;
    float _tau, _epsilonBar, _mu, _K, _delta, _anisotropy, _alpha, _gamma, _tEq;

    std::vector<float> _phi, _t, _epsilon, _epsilonDeriv, _gradPhiX, _gradPhiY, _lapPhi, _lapT, _angl;
    
    // OpenGL 纹理
    std::vector<unsigned char> _pixelBuffer;
    GLuint _textureID = 0;
    bool _updateFlag = true;

    void _initParams();
    void _vectorInit();
    void _createNucleus(int x, int y);
    void _computeGradientLaplacian();
    void _evolution();
    void _updateTexture();
};