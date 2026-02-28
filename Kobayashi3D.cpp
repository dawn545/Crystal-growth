#include "Kobayashi3D.h"
#include <cfloat>

// ==========================================
// 构造函数与初始化
// ==========================================

Kobayashi3D::Kobayashi3D(int x, int y, int z, float timeStep) {
    _objectCount = { x, y, z };
    // 恢复原始网格间距
    _dx = 0.03f;
    _dy = 0.03f;
    _dz = 0.03f;
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
    size_t vSize = _objectCount.x * _objectCount.y * _objectCount.z;

    // 主要场
    _phi.assign(vSize, 0.0);
    _t.assign(vSize, 0.0);
    _orientationField.assign(vSize, 0.0);

    // 辅助场
    _gradPhiX.assign(vSize, 0.0);
    _gradPhiY.assign(vSize, 0.0);
    _gradPhiZ.assign(vSize, 0.0);
    _lapPhi.assign(vSize, 0.0);
    _lapT.assign(vSize, 0.0);
    _gradPhiMag.assign(vSize, 0.0);
    _epsilon.assign(vSize, 0.0);
    _epsilonDeriv.assign(vSize, 0.0);

    _pixelBuffer.assign(_objectCount.x * _objectCount.y * 4, 0);

    // 创建初始晶核
    int centerX = _objectCount.x / 2;
    int centerY = _objectCount.y / 2;
    int centerZ = _objectCount.z / 2;
    float seedRadius = 3.0f;

    _createNucleus(centerX, centerY, centerZ, seedRadius);
    _updateTexture();
}

void Kobayashi3D::_createNucleus(int x, int y, int z, float radius) {
    // 在 3D 空间中创建球形晶核
    int iRadius = static_cast<int>(radius) + 1;

    for (int k = z - iRadius; k <= z + iRadius; k++) {
        for (int j = y - iRadius; j <= y + iRadius; j++) {
            for (int i = x - iRadius; i <= x + iRadius; i++) {
                if (i < 0 || i >= _objectCount.x ||
                    j < 0 || j >= _objectCount.y ||
                    k < 0 || k >= _objectCount.z) {
                    continue;
                }

                float dx = static_cast<float>(i - x);
                float dy = static_cast<float>(j - y);
                float dz = static_cast<float>(k - z);
                float dist = sqrtf(dx * dx + dy * dy + dz * dz);

                if (dist < radius) {
                    _phi[IDX(i, j, k)] = 1.0;
                }
            }
        }
    }
}

// ==========================================
// 物理模拟核心（与原版相同）
// ==========================================

void Kobayashi3D::_computeGradientLaplacian() {
    // 在中间 z 切片计算（2D 模拟）
    const double EPSILON_GRAD = 1e-10;
    int sliceZ = _objectCount.z / 2;

    for (int j = 0; j < _objectCount.y; j++) {
        for (int i = 0; i < _objectCount.x; i++) {
            int idx = IDX(i, j, sliceZ);

            int i_plus = clampX(i + 1);
            int i_minus = clampX(i - 1);
            int j_plus = clampY(j + 1);
            int j_minus = clampY(j - 1);

            // 计算梯度
            _gradPhiX[idx] = (_phi[IDX(i_plus, j, sliceZ)] - _phi[IDX(i_minus, j, sliceZ)]) / (2.0 * _dx);
            _gradPhiY[idx] = (_phi[IDX(i, j_plus, sliceZ)] - _phi[IDX(i, j_minus, sliceZ)]) / (2.0 * _dy);

            double gx = _gradPhiX[idx];
            double gy = _gradPhiY[idx];
            _gradPhiMag[idx] = sqrt(gx * gx + gy * gy);

            // 计算拉普拉斯（9点模板）
            _lapPhi[idx] = (2.0 * (_phi[IDX(i_plus, j, sliceZ)] + _phi[IDX(i_minus, j, sliceZ)] +
                                   _phi[IDX(i, j_plus, sliceZ)] + _phi[IDX(i, j_minus, sliceZ)])
                + _phi[IDX(i_plus, j_plus, sliceZ)] + _phi[IDX(i_minus, j_minus, sliceZ)] +
                  _phi[IDX(i_minus, j_plus, sliceZ)] + _phi[IDX(i_plus, j_minus, sliceZ)]
                - 12.0 * _phi[idx]) / (3.0 * _dx * _dx);

            _lapT[idx] = (2.0 * (_t[IDX(i_plus, j, sliceZ)] + _t[IDX(i_minus, j, sliceZ)] +
                                 _t[IDX(i, j_plus, sliceZ)] + _t[IDX(i, j_minus, sliceZ)])
                + _t[IDX(i_plus, j_plus, sliceZ)] + _t[IDX(i_minus, j_minus, sliceZ)] +
                  _t[IDX(i_minus, j_plus, sliceZ)] + _t[IDX(i_plus, j_minus, sliceZ)]
                - 12.0 * _t[idx]) / (3.0 * _dx * _dx);

            // 各向异性
            double theta = atan2(-gy, -gx);
            double omega = _orientationField[idx];
            double relativeAngle = theta - omega;

            double sigma = 1.0 + _delta * cos(_anisotropy * relativeAngle);
            _epsilon[idx] = _epsilonBar * sigma;
            _epsilonDeriv[idx] = -_epsilonBar * _anisotropy * _delta * sin(_anisotropy * relativeAngle);
        }
    }
}

void Kobayashi3D::_evolution() {
    // 2D 相场演化（在中间 z 切片）
    int sliceZ = _objectCount.z / 2;

    for (int j = 0; j < _objectCount.y; j++) {
        for (int i = 0; i < _objectCount.x; i++) {
            int i_plus = clampX(i + 1);
            int i_minus = clampX(i - 1);
            int j_plus = clampY(j + 1);
            int j_minus = clampY(j - 1);

            double gradEpsPowX = (_epsilon[IDX(i_plus, j, sliceZ)] * _epsilon[IDX(i_plus, j, sliceZ)] -
                                  _epsilon[IDX(i_minus, j, sliceZ)] * _epsilon[IDX(i_minus, j, sliceZ)]) / (2.0 * _dx);
            double gradEpsPowY = (_epsilon[IDX(i, j_plus, sliceZ)] * _epsilon[IDX(i, j_plus, sliceZ)] -
                                  _epsilon[IDX(i, j_minus, sliceZ)] * _epsilon[IDX(i, j_minus, sliceZ)]) / (2.0 * _dy);

            double term1 = (_epsilon[IDX(i, j_plus, sliceZ)] * _epsilonDeriv[IDX(i, j_plus, sliceZ)] * _gradPhiX[IDX(i, j_plus, sliceZ)]
                - _epsilon[IDX(i, j_minus, sliceZ)] * _epsilonDeriv[IDX(i, j_minus, sliceZ)] * _gradPhiX[IDX(i, j_minus, sliceZ)]) / (2.0 * _dy);

            double term2 = -(_epsilon[IDX(i_plus, j, sliceZ)] * _epsilonDeriv[IDX(i_plus, j, sliceZ)] * _gradPhiY[IDX(i_plus, j, sliceZ)]
                - _epsilon[IDX(i_minus, j, sliceZ)] * _epsilonDeriv[IDX(i_minus, j, sliceZ)] * _gradPhiY[IDX(i_minus, j, sliceZ)]) / (2.0 * _dx);

            int idx = IDX(i, j, sliceZ);
            double term3 = gradEpsPowX * _gradPhiX[idx] + gradEpsPowY * _gradPhiY[idx];
            double m = _alpha / PI_F * atan(_gamma * (_tEq - _t[idx]));

            double oldPhi = _phi[idx];
            double oldT = _t[idx];

            _phi[idx] = _phi[idx] +
                (term1 + term2 + _epsilon[idx] * _epsilon[idx] * _lapPhi[idx] + term3
                    + oldPhi * (1.0 - oldPhi) * (oldPhi - 0.5 + m)) * _dt / _tau;

            _t[idx] = oldT + _lapT[idx] * _dt + _K * (_phi[idx] - oldPhi);
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

    // 提取中间 z 切片用于 2D 纹理可视化
    int sliceZ = _objectCount.z / 2;

    for (int j = 0; j < _objectCount.y; j++) {
        for (int i = 0; i < _objectCount.x; i++) {
            int idx3D = IDX(i, j, sliceZ);
            int idx2D = j * _objectCount.x + i;

            float phi = static_cast<float>(_phi[idx3D]);
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

            _pixelBuffer[idx2D * 4 + 0] = toByte(color.x);
            _pixelBuffer[idx2D * 4 + 1] = toByte(color.y);
            _pixelBuffer[idx2D * 4 + 2] = toByte(color.z);
            _pixelBuffer[idx2D * 4 + 3] = 255;
        }
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

    // 渲染中间 z 切片（2D 可视化兼容）
    int sliceZ = _objectCount.z / 2;

    for (int j = 0; j < _objectCount.y - 1; j++) {
        for (int i = 0; i < _objectCount.x - 1; i++) {
            int idx = IDX(i, j, sliceZ);
            float phi = static_cast<float>(_phi[idx]);

            if (phi < 0.01f) continue;

            float x = -1.0f + i * cellWidth;
            float y = -1.0f + j * cellHeight;
            float z = phi * _crystalThickness;

            int idx2D = j * _objectCount.x + i;
            int pixelIdx = idx2D * 4;
            float r = _pixelBuffer[pixelIdx + 0] / 255.0f;
            float g = _pixelBuffer[pixelIdx + 1] / 255.0f;
            float b = _pixelBuffer[pixelIdx + 2] / 255.0f;

            glColor3f(r, g, b);

            glBegin(GL_QUADS);
            glNormal3f(0.0f, 0.0f, 1.0f);
            glVertex3f(x, y, z);
            glVertex3f(x + cellWidth, y, z);
            glVertex3f(x + cellWidth, y + cellHeight, z);
            glVertex3f(x, y + cellHeight, z);
            glEnd();

            if (phi > 0.3f) {
                glColor3f(r * 0.7f, g * 0.7f, b * 0.7f);

                glBegin(GL_QUADS);
                glNormal3f(0.0f, -1.0f, 0.0f);
                glVertex3f(x, y, 0.0f);
                glVertex3f(x + cellWidth, y, 0.0f);
                glVertex3f(x + cellWidth, y, z);
                glVertex3f(x, y, z);

                glNormal3f(1.0f, 0.0f, 0.0f);
                glVertex3f(x + cellWidth, y, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, z);
                glVertex3f(x + cellWidth, y, z);

                glNormal3f(0.0f, 1.0f, 0.0f);
                glVertex3f(x + cellWidth, y + cellHeight, 0.0f);
                glVertex3f(x, y + cellHeight, 0.0f);
                glVertex3f(x, y + cellHeight, z);
                glVertex3f(x + cellWidth, y + cellHeight, z);

                glNormal3f(-1.0f, 0.0f, 0.0f);
                glVertex3f(x, y + cellHeight, 0.0f);
                glVertex3f(x, y, 0.0f);
                glVertex3f(x, y, z);
                glVertex3f(x, y + cellHeight, z);
                glEnd();
            }
        }
    }

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

    // 放宽 X 轴旋转角度限制，允许更大范围的垂直旋转
    if (_cameraRotX > 85.0f) _cameraRotX = 85.0f;
    if (_cameraRotX < -85.0f) _cameraRotX = -85.0f;
}

void Kobayashi3D::zoomCamera(float delta) {
    _cameraDistance += delta;
    if (_cameraDistance < 0.5f) _cameraDistance = 0.5f;  // 允许更近距离观察
    if (_cameraDistance > 15.0f) _cameraDistance = 15.0f;  // 允许更远距离观察
}

// ==========================================
// 取向场管理 (Ren et al. 2018)
// ==========================================

void Kobayashi3D::_resetOrientationField() {
    // 重置取向场为 0
    size_t vSize = _objectCount.x * _objectCount.y * _objectCount.z;
    _orientationField.assign(vSize, 0.0);
}

// ==========================================
// 交互式笔刷工具 (Ren et al. 2018 Section 4.1 - Guiding)
// ==========================================

void Kobayashi3D::paintOrientation(int screenX, int screenY, float angle, float radius, float blendFactor) {
    // 在中间 z 切片绘制取向场（兼容 2D 交互）
    int gridX = static_cast<int>((float)screenX / 800.0f * _objectCount.x);
    int gridY = static_cast<int>((float)screenY / 800.0f * _objectCount.y);
    int sliceZ = _objectCount.z / 2;

    gridY = _objectCount.y - 1 - gridY;

    if (gridX < 0 || gridX >= _objectCount.x || gridY < 0 || gridY >= _objectCount.y) {
        return;
    }

    float gridRadius = radius * _objectCount.x / 800.0f;
    int iRadius = static_cast<int>(gridRadius) + 1;

    for (int dy = -iRadius; dy <= iRadius; dy++) {
        for (int dx = -iRadius; dx <= iRadius; dx++) {
            int targetX = gridX + dx;
            int targetY = gridY + dy;

            if (targetX < 0 || targetX >= _objectCount.x ||
                targetY < 0 || targetY >= _objectCount.y) {
                continue;
            }

            float dist = sqrtf(static_cast<float>(dx * dx + dy * dy));

            if (dist <= gridRadius) {
                int idx = IDX(targetX, targetY, sliceZ);

                float falloff = 1.0f - (dist / gridRadius);
                float effectiveBlend = blendFactor * falloff;

                double oldAngle = _orientationField[idx];
                double newAngle = angle;

                const double PI = 3.14159265358979323846;
                double angleDiff = newAngle - oldAngle;

                while (angleDiff > PI) angleDiff -= 2.0 * PI;
                while (angleDiff < -PI) angleDiff += 2.0 * PI;

                double blendedAngle = oldAngle + angleDiff * effectiveBlend;

                while (blendedAngle < 0.0) blendedAngle += 2.0 * PI;
                while (blendedAngle >= 2.0 * PI) blendedAngle -= 2.0 * PI;

                _orientationField[idx] = blendedAngle;
            }
        }
    }
}
