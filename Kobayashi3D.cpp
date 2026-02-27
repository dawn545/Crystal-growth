#include "Kobayashi3D.h"
#include <cfloat>

// ==========================================
// 构造函数与初始化
// ==========================================

Kobayashi3D::Kobayashi3D(int x, int y, float timeStep) {
    _objectCount = { x, y };
    _dx = 0.03f;
    _dy = 0.03f;
    _dt = timeStep;

    _initParams();
    _vectorInit();
}

Kobayashi3D::~Kobayashi3D() {
    if (_textureID) glDeleteTextures(1, &_textureID);
}

void Kobayashi3D::_initParams() {
    _tau = 0.0003f;
    _epsilonBar = 0.010f;
    _mu = 1.0f;
    _K = 1.6f;
    _delta = 0.05f;
    _anisotropy = 6.0f;
    _alpha = 0.9f;
    _gamma = 10.0f;
    _tEq = 1.0f;
}

void Kobayashi3D::_vectorInit() {
    size_t vSize = _objectCount.x * _objectCount.y;

    _phi.assign(vSize, 0.0f);
    _t.assign(vSize, 0.0f);
    _gradPhiX.assign(vSize, 0.0f);
    _gradPhiY.assign(vSize, 0.0f);
    _lapPhi.assign(vSize, 0.0f);
    _lapT.assign(vSize, 0.0f);
    _angl.assign(vSize, 0.0f);
    _epsilon.assign(vSize, 0.0f);
    _epsilonDeriv.assign(vSize, 0.0f);
    _pixelBuffer.assign(vSize * 4, 0);

    _createNucleus(_objectCount.x / 2, _objectCount.y / 2);
    _updateTexture();
}

void Kobayashi3D::_createNucleus(int x, int y) {
    _phi[_INDEX(x, y)] = 1.0f;
    _phi[_INDEX(x - 1, y)] = 1.0f;
    _phi[_INDEX(x + 1, y)] = 1.0f;
    _phi[_INDEX(x, y - 1)] = 1.0f;
    _phi[_INDEX(x, y + 1)] = 1.0f;
}

// ==========================================
// 物理模拟核心（与原版相同）
// ==========================================

void Kobayashi3D::_computeGradientLaplacian() {
    for (int j = 0; j < _objectCount.y; j++) {
        for (int i = 0; i < _objectCount.x; i++) {
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            _gradPhiX[_INDEX(i, j)] = (_phi[_INDEX(i_plus, j)] - _phi[_INDEX(i_minus, j)]) / _dx;
            _gradPhiY[_INDEX(i, j)] = (_phi[_INDEX(i, j_plus)] - _phi[_INDEX(i, j_minus)]) / _dy;

            _lapPhi[_INDEX(i, j)] = (2.0f * (_phi[_INDEX(i_plus, j)] + _phi[_INDEX(i_minus, j)] + _phi[_INDEX(i, j_plus)] + _phi[_INDEX(i, j_minus)])
                + _phi[_INDEX(i_plus, j_plus)] + _phi[_INDEX(i_minus, j_minus)] + _phi[_INDEX(i_minus, j_plus)] + _phi[_INDEX(i_plus, j_minus)]
                - 12.0f * _phi[_INDEX(i, j)]) / (3.0f * _dx * _dx);

            _lapT[_INDEX(i, j)] = (2.0f * (_t[_INDEX(i_plus, j)] + _t[_INDEX(i_minus, j)] + _t[_INDEX(i, j_plus)] + _t[_INDEX(i, j_minus)])
                + _t[_INDEX(i_plus, j_plus)] + _t[_INDEX(i_minus, j_minus)] + _t[_INDEX(i_minus, j_plus)] + _t[_INDEX(i_plus, j_minus)]
                - 12.0f * _t[_INDEX(i, j)]) / (3.0f * _dx * _dx);

            if (_gradPhiX[_INDEX(i, j)] <= +FLT_EPSILON && _gradPhiX[_INDEX(i, j)] >= -FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = -0.5f * PI_F;
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = 0.5f * PI_F;

            if (_gradPhiX[_INDEX(i, j)] > +FLT_EPSILON)
                if (_gradPhiY[_INDEX(i, j)] < -FLT_EPSILON) _angl[_INDEX(i, j)] = 2.0f * PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);
                else if (_gradPhiY[_INDEX(i, j)] > +FLT_EPSILON) _angl[_INDEX(i, j)] = atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            if (_gradPhiX[_INDEX(i, j)] < -FLT_EPSILON)
                _angl[_INDEX(i, j)] = PI_F + atan(_gradPhiY[_INDEX(i, j)] / _gradPhiX[_INDEX(i, j)]);

            // 计算各向异性系数
            _epsilon[_INDEX(i, j)] = _epsilonBar * (1.0f + _delta * cos(_anisotropy * _angl[_INDEX(i, j)]));
            _epsilonDeriv[_INDEX(i, j)] = -_epsilonBar * _anisotropy * _delta * sin(_anisotropy * _angl[_INDEX(i, j)]);
        }
    }
}

void Kobayashi3D::_evolution() {
    for (int j = 0; j < _objectCount.y; j++) {
        for (int i = 0; i < _objectCount.x; i++) {
            int i_plus = (i + 1) % _objectCount.x;
            int i_minus = ((i - 1) + _objectCount.x) % _objectCount.x;
            int j_plus = (j + 1) % _objectCount.y;
            int j_minus = ((j - 1) + _objectCount.y) % _objectCount.y;

            float gradEpsPowX = (_epsilon[_INDEX(i_plus, j)] * _epsilon[_INDEX(i_plus, j)] - _epsilon[_INDEX(i_minus, j)] * _epsilon[_INDEX(i_minus, j)]) / _dx;
            float gradEpsPowY = (_epsilon[_INDEX(i, j_plus)] * _epsilon[_INDEX(i, j_plus)] - _epsilon[_INDEX(i, j_minus)] * _epsilon[_INDEX(i, j_minus)]) / _dy;

            float term1 = (_epsilon[_INDEX(i, j_plus)] * _epsilonDeriv[_INDEX(i, j_plus)] * _gradPhiX[_INDEX(i, j_plus)]
                - _epsilon[_INDEX(i, j_minus)] * _epsilonDeriv[_INDEX(i, j_minus)] * _gradPhiX[_INDEX(i, j_minus)]) / _dy;

            float term2 = -(_epsilon[_INDEX(i_plus, j)] * _epsilonDeriv[_INDEX(i_plus, j)] * _gradPhiY[_INDEX(i_plus, j)]
                - _epsilon[_INDEX(i_minus, j)] * _epsilonDeriv[_INDEX(i_minus, j)] * _gradPhiY[_INDEX(i_minus, j)]) / _dx;

            float term3 = gradEpsPowX * _gradPhiX[_INDEX(i, j)] + gradEpsPowY * _gradPhiY[_INDEX(i, j)];
            float m = _alpha / PI_F * atan(_gamma * (_tEq - _t[_INDEX(i, j)]));

            float oldPhi = _phi[_INDEX(i, j)];
            float oldT = _t[_INDEX(i, j)];

            _phi[_INDEX(i, j)] = _phi[_INDEX(i, j)] +
                (term1 + term2 + _epsilon[_INDEX(i, j)] * _epsilon[_INDEX(i, j)] * _lapPhi[_INDEX(i, j)] + term3
                    + oldPhi * (1.0f - oldPhi) * (oldPhi - 0.5f + m)) * _dt / _tau;

            _t[_INDEX(i, j)] = oldT + _lapT[_INDEX(i, j)] * _dt + _K * (_phi[_INDEX(i, j)] - oldPhi);
        }
    }
}

void Kobayashi3D::update() {
    if (!_updateFlag) return;

    for (int i = 0; i < 10; i++) {
        _computeGradientLaplacian();
        _evolution();
    }
    _updateTexture();
}

void Kobayashi3D::reset() {
    _vectorInit();
}

// ==========================================
// 3D 渲染部分
// ==========================================

void Kobayashi3D::glInit() {
    // 启用深度测试和光照
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    // 设置主光源（从上方照射）
    GLfloat light0_pos[] = { 0.0f, 0.0f, 5.0f, 1.0f };
    GLfloat light0_ambient[] = { 0.2f, 0.2f, 0.3f, 1.0f };
    GLfloat light0_diffuse[] = { 0.8f, 0.8f, 1.0f, 1.0f };
    GLfloat light0_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);

    // 添加第二个光源（侧面补光）
    glEnable(GL_LIGHT1);
    GLfloat light1_pos[] = { 3.0f, 2.0f, 2.0f, 1.0f };
    GLfloat light1_ambient[] = { 0.1f, 0.1f, 0.15f, 1.0f };
    GLfloat light1_diffuse[] = { 0.4f, 0.5f, 0.6f, 1.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);

    // 设置材质属性
    GLfloat mat_specular[] = { 0.8f, 0.8f, 1.0f, 1.0f };
    GLfloat mat_shininess[] = { 50.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    // 启用平滑着色
    glShadeModel(GL_SMOOTH);

    // 初始化纹理
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &_textureID);
    glBindTexture(GL_TEXTURE_2D, _textureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
}

void Kobayashi3D::_updateTexture() {
    struct float3 { float x, y, z; };
    float3 c0 = { 0.0f, 0.0f, 0.0f };
    float3 c1 = { 0.25f, 0.50f, 0.98f };
    float3 c2 = { 0.36f, 1.00f, 0.98f };
    float3 c3 = { 0.90f, 1.00f, 0.98f };

    float c1Boundary = 0.9f;
    float c2Boundary = 0.99f;

    size_t size = _phi.size();

    for (size_t k = 0; k < size; k++) {
        float phi = _phi[k];
        float3 color;
        float ratio;

        if (phi <= c1Boundary) {
            ratio = phi * (1.0f / c1Boundary);
            color.x = c0.x * (1.0f - ratio) + c1.x * ratio;
            color.y = c0.y * (1.0f - ratio) + c1.y * ratio;
            color.z = c0.z * (1.0f - ratio) + c1.z * ratio;
        }
        else if (phi > c1Boundary && phi <= c2Boundary) {
            ratio = (phi - c1Boundary) * (1.0f / (c2Boundary - c1Boundary));
            color.x = c1.x * (1.0f - ratio) + c2.x * ratio;
            color.y = c1.y * (1.0f - ratio) + c2.y * ratio;
            color.z = c1.z * (1.0f - ratio) + c2.z * ratio;
        }
        else {
            float c3Boundary = 1.0f;
            ratio = (phi - c2Boundary) * (1.0f / (c3Boundary - c2Boundary));
            color.x = c2.x * (1.0f - ratio) + c3.x * ratio;
            color.y = c2.y * (1.0f - ratio) + c3.y * ratio;
            color.z = c2.z * (1.0f - ratio) + c3.z * ratio;
        }

        auto toByte = [](float v) -> unsigned char {
            if (v < 0.0f) return 0;
            if (v > 1.0f) return 255;
            return static_cast<unsigned char>(v * 255.0f);
        };

        _pixelBuffer[k * 4 + 0] = toByte(color.x);
        _pixelBuffer[k * 4 + 1] = toByte(color.y);
        _pixelBuffer[k * 4 + 2] = toByte(color.z);
        _pixelBuffer[k * 4 + 3] = 255;
    }

    glBindTexture(GL_TEXTURE_2D, _textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _objectCount.x, _objectCount.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, _pixelBuffer.data());
}

void Kobayashi3D::_render3DCrystal() {
    // 计算网格单元大小
    float cellWidth = 2.0f / _objectCount.x;
    float cellHeight = 2.0f / _objectCount.y;

    glEnable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    // 绘制每个网格单元为立方体
    for (int j = 0; j < _objectCount.y - 1; j++) {
        for (int i = 0; i < _objectCount.x - 1; i++) {
            int idx = _INDEX(i, j);
            float phi = _phi[idx];

            // 只绘制有晶体的区域
            if (phi < 0.01f) continue;

            // 计算位置
            float x = -1.0f + i * cellWidth;
            float y = -1.0f + j * cellHeight;
            float z = phi * _crystalThickness;  // 根据相场值确定高度

            // 获取颜色
            int pixelIdx = idx * 4;
            float r = _pixelBuffer[pixelIdx + 0] / 255.0f;
            float g = _pixelBuffer[pixelIdx + 1] / 255.0f;
            float b = _pixelBuffer[pixelIdx + 2] / 255.0f;

            // 设置材质颜色
            glColor3f(r, g, b);

            // 绘制立方体（顶面）
            glBegin(GL_QUADS);
            // 顶面
            glNormal3f(0.0f, 0.0f, 1.0f);
            glVertex3f(x, y, z);
            glVertex3f(x + cellWidth, y, z);
            glVertex3f(x + cellWidth, y + cellHeight, z);
            glVertex3f(x, y + cellHeight, z);
            glEnd();

            // 如果相场值较高，绘制侧面增强立体感
            if (phi > 0.3f) {
                // 稍微变暗的颜色用于侧面
                glColor3f(r * 0.7f, g * 0.7f, b * 0.7f);

                glBegin(GL_QUADS);
                // 前面
                glNormal3f(0.0f, -1.0f, 0.0f);
                glVertex3f(x, y, 0.0f);
                glVertex3f(x + cellWidth, y, 0.0f);
                glVertex3f(x + cellWidth, y, z);
                glVertex3f(x, y, z);

                // 右面
                glNormal3f(1.0f, 0.0f, 0.0f);
                glVertex3f(x + cellWidth, y, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, z);
                glVertex3f(x + cellWidth, y, z);

                // 后面
                glNormal3f(0.0f, 1.0f, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, 0.0f);
                glVertex3f(x, y + cellHeight, 0.0f);
                glVertex3f(x, y + cellHeight, z);
                glVertex3f(x + cellWidth, y + cellHeight, z);

                // 左面
                glNormal3f(-1.0f, 0.0f, 0.0f);
                glVertex3f(x, y + cellHeight, 0.0f);
                glVertex3f(x, y, 0.0f);
                glVertex3f(x, y, z);
                glVertex3f(x, y + cellHeight, z);
                glEnd();
            }
        }
    }

    // 绘制底部平面作为参考
    glDisable(GL_LIGHTING);
    glColor4f(0.1f, 0.1f, 0.15f, 0.5f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_QUADS);
    glVertex3f(-1.0f, -1.0f, 0.0f);
    glVertex3f( 1.0f, -1.0f, 0.0f);
    glVertex3f( 1.0f,  1.0f, 0.0f);
    glVertex3f(-1.0f,  1.0f, 0.0f);
    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void Kobayashi3D::glRender() {
    // 设置深蓝色背景，模拟冰雪环境
    glClearColor(0.05f, 0.08f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // 设置投影矩阵
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 0.1, 100.0);

    // 设置模型视图矩阵
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // 相机位置
    gluLookAt(
        0.0, 0.0, _cameraDistance,  // 相机位置
        0.0, 0.0, 0.0,               // 看向原点
        0.0, 1.0, 0.0                // 上方向
    );

    // 应用旋转
    glRotatef(_cameraRotX, 1.0f, 0.0f, 0.0f);
    glRotatef(_cameraRotY, 0.0f, 1.0f, 0.0f);

    _render3DCrystal();

    glutSwapBuffers();
}

// ==========================================
// 相机控制
// ==========================================

void Kobayashi3D::rotateCamera(float deltaX, float deltaY) {
    _cameraRotY += deltaX;
    _cameraRotX += deltaY;

    // 限制 X 轴旋转角度
    if (_cameraRotX > 89.0f) _cameraRotX = 89.0f;
    if (_cameraRotX < -89.0f) _cameraRotX = -89.0f;
}

void Kobayashi3D::zoomCamera(float delta) {
    _cameraDistance += delta;
    if (_cameraDistance < 1.0f) _cameraDistance = 1.0f;
    if (_cameraDistance > 10.0f) _cameraDistance = 10.0f;
}
