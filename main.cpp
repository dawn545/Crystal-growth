#include <GL/freeglut.h>
#include "Kobayashi.h"

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
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Kobayashi Crystal (FreeGLUT)");

    // 2. 初始化模拟器
    g_sim = new Kobayashi(250, 250, 0.0001f);
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