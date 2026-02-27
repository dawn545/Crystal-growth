#include <GL/freeglut.h>
#include "Kobayashi3D.h"

Kobayashi3D* g_sim = nullptr;

// 鼠标控制变量
int g_lastMouseX = 0;
int g_lastMouseY = 0;
bool g_mouseLeftDown = false;
bool g_mouseRightDown = false;  // 右键用于笔刷工具
float g_brushAngle = 0.0f;      // 笔刷方向角度
float g_brushRadius = 20.0f;    // 笔刷半径（像素）- 增大默认值
float g_brushBlend = 0.6f;      // 笔刷混合强度 - 提高默认强度

// 渲染回调
void display() {
    if (g_sim) g_sim->glRender();
}

// 闲置回调
void idle() {
    if (g_sim) {
        g_sim->update();
        glutPostRedisplay();
    }
}

// 键盘回调
void keyboard(unsigned char key, int x, int y) {
    if (!g_sim) return;

    switch (key) {
    case 27: // ESC 键
        glutLeaveMainLoop();
        break;
    case ' ': // 空格键暂停/播放
        g_sim->togglePause();
        std::cout << (g_sim->isPaused() ? "Paused" : "Running") << std::endl;
        break;
    case 'r': // R 键重置
    case 'R':
        g_sim->reset();
        std::cout << "Reset" << std::endl;
        break;
    case '0': // 数字 0 - 重置为全局统一取向场
        g_sim->resetOrientationField();
        std::cout << "Orientation Field: Uniform (Omega = 0)" << std::endl;
        break;
    case '1': // 数字 1 - 应用逆时针涡旋取向场
        g_sim->applyVortexField(false);
        std::cout << "Orientation Field: Vortex Counter-Clockwise" << std::endl;
        break;
    case '2': // 数字 2 - 应用顺时针涡旋取向场
        g_sim->applyVortexField(true);
        std::cout << "Orientation Field: Vortex Clockwise" << std::endl;
        break;
    case '+': // 增加笔刷大小
    case '=':
        g_brushRadius += 5.0f;
        if (g_brushRadius > 50.0f) g_brushRadius = 50.0f;
        std::cout << "Brush radius: " << g_brushRadius << " pixels" << std::endl;
        break;
    case '-': // 减小笔刷大小
    case '_':
        g_brushRadius -= 5.0f;
        if (g_brushRadius < 5.0f) g_brushRadius = 5.0f;
        std::cout << "Brush radius: " << g_brushRadius << " pixels" << std::endl;
        break;
    case '[': // 减小笔刷强度
        g_brushBlend -= 0.1f;
        if (g_brushBlend < 0.1f) g_brushBlend = 0.1f;
        std::cout << "Brush strength: " << g_brushBlend << std::endl;
        break;
    case ']': // 增加笔刷强度
        g_brushBlend += 0.1f;
        if (g_brushBlend > 1.0f) g_brushBlend = 1.0f;
        std::cout << "Brush strength: " << g_brushBlend << std::endl;
        break;
    }
}

// 鼠标按钮回调
void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            g_mouseLeftDown = true;
            g_lastMouseX = x;
            g_lastMouseY = y;
        }
        else {
            g_mouseLeftDown = false;
        }
    }
    else if (button == GLUT_RIGHT_BUTTON) {
        if (state == GLUT_DOWN) {
            g_mouseRightDown = true;
            g_lastMouseX = x;
            g_lastMouseY = y;
            std::cout << "Brush tool activated - drag to paint orientation" << std::endl;
        }
        else {
            g_mouseRightDown = false;
            std::cout << "Brush tool deactivated" << std::endl;
        }
    }
    // 鼠标滚轮缩放
    else if (button == 3) { // 滚轮向上
        if (g_sim) g_sim->zoomCamera(-0.3f);  // 增加缩放步长
        glutPostRedisplay();
    }
    else if (button == 4) { // 滚轮向下
        if (g_sim) g_sim->zoomCamera(0.3f);
        glutPostRedisplay();
    }
}

// 鼠标移动回调
void motion(int x, int y) {
    if (!g_sim) return;

    // 左键拖拽旋转视角
    if (g_mouseLeftDown) {
        float deltaX = (x - g_lastMouseX) * 0.3f;  // 降低灵敏度，更平滑
        float deltaY = (y - g_lastMouseY) * 0.3f;

        g_sim->rotateCamera(deltaY, deltaX);

        g_lastMouseX = x;
        g_lastMouseY = y;

        glutPostRedisplay();
    }
    // 右键拖拽涂抹取向场
    else if (g_mouseRightDown) {
        // 计算鼠标移动向量
        int dx = x - g_lastMouseX;
        int dy = y - g_lastMouseY;

        // 只有当鼠标移动时才更新角度和涂抹
        if (dx != 0 || dy != 0) {
            // 计算移动方向的角度
            g_brushAngle = atan2f(static_cast<float>(dy), static_cast<float>(dx));

            // 应用笔刷（使用全局笔刷参数）
            g_sim->paintOrientation(x, y, g_brushAngle, g_brushRadius, g_brushBlend);

            g_lastMouseX = x;
            g_lastMouseY = y;

            glutPostRedisplay();
        }
    }
}

int main(int argc, char** argv) {
    // 1. 初始化 GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Kobayashi Crystal 3D - Ice Crystal with Thickness");

    // 启用抗锯齿
    glEnable(GLUT_MULTISAMPLE);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    // 2. 初始化模拟器
    g_sim = new Kobayashi3D(250, 250, 0.0001f);
    g_sim->glInit();

    std::cout << "=== 3D Ice Crystal Simulation ===" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << " [Left Mouse + Drag]: Rotate view" << std::endl;
    std::cout << " [Right Mouse + Drag]: Paint orientation field (Guiding)" << std::endl;
    std::cout << " [Mouse Wheel]: Zoom in/out" << std::endl;
    std::cout << " [Space]: Pause/Play simulation" << std::endl;
    std::cout << " [R]: Reset simulation" << std::endl;
    std::cout << " [0]: Uniform orientation field (default)" << std::endl;
    std::cout << " [1]: Vortex field - Counter-Clockwise" << std::endl;
    std::cout << " [2]: Vortex field - Clockwise" << std::endl;
    std::cout << " [+/-]: Increase/Decrease brush size" << std::endl;
    std::cout << " [[/]]: Decrease/Increase brush strength" << std::endl;
    std::cout << " [ESC]: Quit" << std::endl;

    // 3. 注册回调函数
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);

    // 4. 进入主循环
    glutMainLoop();

    delete g_sim;
    return 0;
}
