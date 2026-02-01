# CrystalGrowth

This project simulates dendritic crystal growth using the Kobayashi model. It was initially based on the [jklae/CrystalGrowth](https://github.com/jklae/CrystalGrowth) project, which used DXViewer for rendering. I modified it to use FreeGLUT for cross-platform compatibility and improved accessibility.

## Features

- **Dendritic Crystal Growth Simulation**: Models crystal growth in a 2D grid using the Kobayashi model, incorporating anisotropic material properties and thermal dynamics.
- **FreeGLUT Integration**: Replaces DXViewer with FreeGLUT for rendering and interaction, ensuring compatibility with a wider range of platforms.
- **Physical Modeling**: Includes calculations for gradient, Laplacian, and anisotropy effects on crystal growth, temperature field evolution, and phase transitions.

## Reference
Kim, Y., & Lin, S. (2003). Visual Simulation of Ice Crystal Growth. Eurographics/SIGGRAPH Symposium on Computer Animation.

## Compilation Command

To compile the project, use the following command:

```bash
g++ main.cpp Kobayashi.cpp -I. -lopengl32 -lfreeglut -o main.exe
