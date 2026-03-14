#include <GL/freeglut.h>
#include "Kobayashi3D.h"

Kobayashi* g_sim = nullptr;

// 渲染回调
void display() {
    if (g_sim) g_sim->glRender();
}

// 闲置回调（相当于 Update 循环）
void idle() {
    if (g_sim) {
        g_sim->update();
        glutPostRedisplay(); // 请求重绘
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

int main(int argc, char** argv) {
    // 1. 初始化 GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Kobayashi 3D Crystal Growth");

    // 启用深度测试（3D渲染必需）
    glEnable(GL_DEPTH_TEST);

    // 2. 初始化模拟器（3D网格：100x100x100）
    g_sim = new Kobayashi(100, 100, 100, 0.0001f);
    g_sim->glInit();

    std::cout << "Controls:\n [Space]: Pause/Play\n [R]: Reset\n [ESC]: Quit" << std::endl;

    // 3. 注册回调函数
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);

    // 4. 进入主循环
    glutMainLoop();

    delete g_sim;
    return 0;
}
