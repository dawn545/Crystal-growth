#include <GL/freeglut.h>
#include "Kobayashi3D.h"

Kobayashi3D* g_sim = nullptr;

// 鼠标控制变量
int g_lastMouseX = 0;
int g_lastMouseY = 0;
bool g_mouseLeftDown = false;

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
    // 鼠标滚轮缩放
    else if (button == 3) { // 滚轮向上
        if (g_sim) g_sim->zoomCamera(-0.2f);
        glutPostRedisplay();
    }
    else if (button == 4) { // 滚轮向下
        if (g_sim) g_sim->zoomCamera(0.2f);
        glutPostRedisplay();
    }
}

// 鼠标移动回调
void motion(int x, int y) {
    // 左键拖拽旋转视角
    if (g_mouseLeftDown && g_sim) {
        float deltaX = (x - g_lastMouseX) * 0.5f;
        float deltaY = (y - g_lastMouseY) * 0.5f;

        g_sim->rotateCamera(deltaY, deltaX);

        g_lastMouseX = x;
        g_lastMouseY = y;

        glutPostRedisplay();
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
    std::cout << " [Mouse Wheel]: Zoom in/out" << std::endl;
    std::cout << " [Space]: Pause/Play simulation" << std::endl;
    std::cout << " [R]: Reset simulation" << std::endl;
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
